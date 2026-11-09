/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_GEOMECH_CENTRE_STRESS_INTERPOLATION_HPP
#define OPM_GEOMECH_CENTRE_STRESS_INTERPOLATION_HPP

#include <dune/common/fvector.hh>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Fault.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaultCollection.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaultFace.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace Opm::centre_stress
{

using Vec3 = Dune::FieldVector<double, 3>;
using Sym6 = Dune::FieldVector<double, 6>;

enum class Mode { None, Lateral, InPlane };

inline Mode parseMode(const std::string& name)
{
    if (name == "none") {
        return Mode::None;
    }
    if (name == "centre_lateral") {
        return Mode::Lateral;
    }
    if (name == "centre_in_plane") {
        return Mode::InPlane;
    }
    throw std::runtime_error("Unknown reservoir.stress_interpolation: " + name);
}

// Why a cell side may or may not be interpolated across.
enum class SideKind {
    Clean,              // one neighbour on each side, adjacent in logical indices
    Boundary,           // grid or process boundary, or inactive cell on part of the side
    NonConforming,      // several neighbours: throw, or a coarse/fine interface
    LogicalOffset,      // single neighbour in another layer or column: throw, pinch-out
    DeclaredFault,      // listed in FAULTS, whatever the throw
    NoTransmissibility, // barrier to flow, so the temperature step is real
    NoMechanics         // flow cell outside a mechanics sector, its stress is zero
};
constexpr std::size_t numSideKinds = 7;

inline bool isFaultSide(SideKind kind)
{
    return kind == SideKind::NonConforming || kind == SideKind::LogicalOffset
        || kind == SideKind::DeclaredFault || kind == SideKind::NoTransmissibility;
}

struct Neighbour {
    int cell = -1;
    double area = 0.0;
};

struct CellSide {
    SideKind kind = SideKind::Boundary;
    std::vector<Neighbour> neighbours;
};

// Sides are indexed as Intersection::indexInInside(): I-, I+, J-, J+, K-, K+.
using CellSides = std::array<CellSide, 6>;

inline int faceTag(FaceDir::DirEnum dir)
{
    switch (dir) {
    case FaceDir::XMinus: return 0;
    case FaceDir::XPlus: return 1;
    case FaceDir::YMinus: return 2;
    case FaceDir::YPlus: return 3;
    case FaceDir::ZMinus: return 4;
    case FaceDir::ZPlus: return 5;
    default: return -1;
    }
}

class DeclaredFaults
{
public:
    DeclaredFaults() = default;

    explicit DeclaredFaults(const FaultCollection& faults)
    {
        for (std::size_t f = 0; f < faults.size(); ++f) {
            for (const auto& face : faults.getFault(f)) {
                const int tag = faceTag(face.getDir());
                if (tag < 0) {
                    continue;
                }
                for (const auto cartesianIdx : face) {
                    add(cartesianIdx, tag);
                }
            }
        }
    }

    void add(std::size_t cartesianIdx, int tag)
    {
        faces_.insert(key(cartesianIdx, tag));
    }

    // A deck may list a fault on either of the two cells sharing the face.
    bool contains(std::size_t cartesianIdx, int tag, std::size_t otherCartesianIdx) const
    {
        return faces_.count(key(cartesianIdx, tag)) > 0
            || faces_.count(key(otherCartesianIdx, tag ^ 1)) > 0;
    }

private:
    static std::uint64_t key(std::size_t cartesianIdx, int tag)
    {
        return (static_cast<std::uint64_t>(cartesianIdx) << 3) | static_cast<std::uint64_t>(tag);
    }

    std::unordered_set<std::uint64_t> faces_;
};

// Single neighbour test: adjacent along the side's axis and level on the other two.
// A difference of zero along the axis occurs between cells refined from one parent.
inline bool logicallyAdjacent(const std::array<int, 3>& cell, const std::array<int, 3>& other, int tag)
{
    const int axis = tag / 2;
    const int expected = (tag % 2 == 0) ? -1 : 1;
    for (int d = 0; d < 3; ++d) {
        const int diff = other[d] - cell[d];
        if (d == axis) {
            if (diff != expected && diff != 0) {
                return false;
            }
        } else if (diff != 0) {
            return false;
        }
    }
    return true;
}

struct RawSide {
    bool hasBoundary = false;
    std::vector<Neighbour> neighbours;
};

template <class GridView, class Element>
std::array<RawSide, 6> rawSides(const GridView& gridView, const Element& element)
{
    std::array<RawSide, 6> sides;
    for (const auto& is : intersections(gridView, element)) {
        const int tag = is.indexInInside();
        if (tag < 0 || tag > 5) {
            continue;
        }
        if (is.neighbor()) {
            const int other = gridView.indexSet().index(is.outside());
            auto& list = sides[tag].neighbours;
            const auto hit = std::find_if(list.begin(), list.end(),
                                          [other](const Neighbour& n) { return n.cell == other; });
            if (hit == list.end()) {
                list.push_back({other, is.geometry().volume()});
            } else {
                hit->area += is.geometry().volume();
            }
        } else {
            sides[tag].hasBoundary = true;
        }
    }
    return sides;
}

// Decide per side of one cell.  `logical(cell)` gives (i,j,k), `cartesian(cell)` the
// index used by FAULTS, `flows(a,b)` false marks a flow barrier, `hasMech(cell)` a
// solved stress.  With `acrossFaults` only boundaries and cells without mechanics stop
// the interpolation.
template <class GridView, class Element, class Logical, class Cartesian, class Flows, class HasMech>
CellSides classifySides(const GridView& gridView, const Element& element,
                        const Logical& logical, const Cartesian& cartesian,
                        const DeclaredFaults& faults, const Flows& flows,
                        const HasMech& hasMech, bool acrossFaults)
{
    CellSides result;
    const int cell = gridView.indexSet().index(element);
    const auto raw = rawSides(gridView, element);
    for (int tag = 0; tag < 6; ++tag) {
        auto& side = result[tag];
        side.neighbours = raw[tag].neighbours;
        if (side.neighbours.empty() || raw[tag].hasBoundary) {
            side.kind = SideKind::Boundary;
            continue;
        }
        if (std::any_of(side.neighbours.begin(), side.neighbours.end(),
                        [&hasMech](const Neighbour& n) { return !hasMech(n.cell); })) {
            side.kind = SideKind::NoMechanics;
            continue;
        }
        if (acrossFaults) {
            side.kind = SideKind::Clean;
            continue;
        }
        if (side.neighbours.size() > 1) {
            side.kind = SideKind::NonConforming;
            continue;
        }
        const int other = side.neighbours.front().cell;
        if (!logicallyAdjacent(logical(cell), logical(other), tag)) {
            side.kind = SideKind::LogicalOffset;
            continue;
        }
        // One neighbour here can still be one of several seen from the other cell.
        bool mirrored = true;
        for (const auto& is : intersections(gridView, element)) {
            if (is.indexInInside() == tag && is.neighbor()) {
                const auto mirror = rawSides(gridView, is.outside())[tag ^ 1];
                mirrored = !mirror.hasBoundary && mirror.neighbours.size() == 1
                    && mirror.neighbours.front().cell == cell;
                break;
            }
        }
        if (!mirrored) {
            side.kind = SideKind::NonConforming;
            continue;
        }
        if (faults.contains(cartesian(cell), tag, cartesian(other))) {
            side.kind = SideKind::DeclaredFault;
            continue;
        }
        if (!flows(cell, other)) {
            side.kind = SideKind::NoTransmissibility;
            continue;
        }
        side.kind = SideKind::Clean;
    }
    return result;
}

// One usable neighbour value on one logical axis.
struct AxisSample {
    int axis = 0; // 0 = I, 1 = J, 2 = K
    Vec3 centre;
    Sym6 stress;
};

// Placeholder for a fault-specific value on a fault side; empty keeps the side flat.
template <class CentreOf, class StressOf>
std::optional<AxisSample> faultSideSample(const CellSide&, int /*tag*/, const CentreOf&, const StressOf&)
{
    return std::nullopt;
}

// Area-weighted neighbour value per usable side.
template <class CentreOf, class StressOf>
std::vector<AxisSample> samplesFromSides(const CellSides& sides, const CentreOf& centreOf, const StressOf& stressOf)
{
    std::vector<AxisSample> samples;
    for (int tag = 0; tag < 6; ++tag) {
        const auto& side = sides[tag];
        if (isFaultSide(side.kind)) {
            if (auto sample = faultSideSample(side, tag, centreOf, stressOf)) {
                samples.push_back(*sample);
            }
            continue;
        }
        if (side.kind != SideKind::Clean) {
            continue;
        }
        AxisSample sample;
        sample.axis = tag / 2;
        sample.centre = 0.0;
        sample.stress = 0.0;
        double area = 0.0;
        for (const auto& n : side.neighbours) {
            sample.centre.axpy(n.area, Vec3(centreOf(n.cell)));
            sample.stress.axpy(n.area, Sym6(stressOf(n.cell)));
            area += n.area;
        }
        if (area <= 0.0) {
            continue;
        }
        sample.centre /= area;
        sample.stress /= area;
        samples.push_back(sample);
    }
    return samples;
}

// Stress at x in a cell.  Keeps the cell's own value where x projects on the cell centre,
// is linear towards each neighbour centre, and keeps the normal traction within the values
// of the cells used.  Offsets are measured in the fracture plane, so neighbours mainly
// across the plane are not used.  `traction(stress)` is the normal traction.
template <class Traction>
Sym6 interpolate(const Vec3& x, const Vec3& normal, const Vec3& centre, const Sym6& stress,
                 const std::vector<AxisSample>& samples, Mode mode, const Traction& traction)
{
    if (mode == Mode::None) {
        return stress;
    }
    const auto inPlane = [&normal](Vec3 v) {
        v.axpy(-(v * normal), normal);
        return v;
    };
    Vec3 delta = x;
    delta -= centre;
    delta = inPlane(delta);

    struct Candidate {
        Vec3 d;
        const Sym6* value;
        double alignment;
    };
    // Per axis the neighbour on the side of x, then the two most in-plane axes.
    std::array<std::optional<Candidate>, 3> perAxis;
    for (const auto& s : samples) {
        if (mode == Mode::Lateral && s.axis == 2) {
            continue;
        }
        Vec3 full = s.centre;
        full -= centre;
        const Vec3 d = inPlane(full);
        const double inPlaneLength = d.two_norm();
        if (inPlaneLength <= 0.0 || std::abs(full * normal) > inPlaneLength) {
            continue;
        }
        if (delta * d <= 0.0) {
            continue;
        }
        perAxis[s.axis] = Candidate{d, &s.stress, inPlaneLength / full.two_norm()};
    }
    std::vector<Candidate> used;
    for (const auto& c : perAxis) {
        if (c) {
            used.push_back(*c);
        }
    }
    std::sort(used.begin(), used.end(),
              [](const Candidate& a, const Candidate& b) { return a.alignment > b.alignment; });
    if (used.size() > 2) {
        used.resize(2);
    }

    std::array<double, 2> t {0.0, 0.0};
    if (used.size() == 2) {
        const double a11 = used[0].d * used[0].d;
        const double a12 = used[0].d * used[1].d;
        const double a22 = used[1].d * used[1].d;
        const double det = a11 * a22 - a12 * a12;
        if (det > 1e-6 * a11 * a22) {
            const double b1 = delta * used[0].d;
            const double b2 = delta * used[1].d;
            t[0] = (a22 * b1 - a12 * b2) / det;
            t[1] = (a11 * b2 - a12 * b1) / det;
        } else {
            used.resize(1);
        }
    }
    if (used.size() == 1) {
        t[0] = (delta * used[0].d) / (used[0].d * used[0].d);
    }

    Sym6 result = stress;
    const double own = traction(stress);
    double low = own;
    double high = own;
    for (std::size_t a = 0; a < used.size(); ++a) {
        const double w = std::clamp(t[a], 0.0, 1.0);
        Sym6 diff = *used[a].value;
        diff -= stress;
        result.axpy(w, diff);
        const double other = traction(*used[a].value);
        low = std::min(low, other);
        high = std::max(high, other);
    }
    const double value = traction(result);
    if (value < low || value > high) {
        const double bound = (value < low) ? low : high;
        const double theta = (bound - own) / (value - own);
        Sym6 diff = result;
        diff -= stress;
        result = stress;
        result.axpy(theta, diff);
    }
    return result;
}

} // namespace Opm::centre_stress

#endif
