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

#include <config.h>
#include <opm/geomech/CentreStressInterpolation.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm::centre_stress;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool close(double a, double b, double tol = 1e-9)
{
    return std::abs(a - b) <= tol * std::max(1.0, std::abs(a));
}

// Isotropic tensors, so the normal traction is the value itself.
Sym6 iso(double v)
{
    Sym6 s(0.0);
    s[0] = s[1] = s[2] = v;
    return s;
}

double traction(const Sym6& s)
{
    return s[0];
}

const Vec3 normalX {1.0, 0.0, 0.0};

double at(const Vec3& x, const Vec3& normal, const Vec3& centre, double own,
          const std::vector<AxisSample>& samples, Mode mode)
{
    return traction(interpolate(x, normal, centre, iso(own), samples, mode, traction));
}

void testOperator()
{
    const Vec3 c {0.0, 0.0, 0.0};
    // Stretched spacing along J, as on a 9 m row next to an 88 m row.
    const std::vector<AxisSample> samples {
        {1, Vec3 {0.0, 48.5, 0.0}, iso(10.0)},
        {1, Vec3 {0.0, -88.0, 0.0}, iso(30.0)},
        {2, Vec3 {0.0, 0.0, 5.0}, iso(25.0)},
        {2, Vec3 {0.0, 0.0, -5.0}, iso(15.0)},
    };

    require(close(at(c, normalX, c, 20.0, samples, Mode::InPlane), 20.0), "centre value kept");
    require(close(at(Vec3 {7.0, 0.0, 0.0}, normalX, c, 20.0, samples, Mode::InPlane), 20.0),
                  "offset along the normal has no effect");
    require(close(at(Vec3 {0.0, 24.25, 0.0}, normalX, c, 20.0, samples, Mode::InPlane), 15.0), "half way to the neighbour");
    require(close(at(Vec3 {0.0, -44.0, 0.0}, normalX, c, 20.0, samples, Mode::InPlane), 25.0), "other side");
    require(close(at(Vec3 {0.0, 500.0, 0.0}, normalX, c, 20.0, samples, Mode::InPlane), 10.0), "stays at the neighbour value");
    require(close(at(Vec3 {0.0, 0.0, 2.5}, normalX, c, 20.0, samples, Mode::Lateral), 20.0), "lateral ignores K");
    require(close(at(Vec3 {0.0, 0.0, 2.5}, normalX, c, 20.0, samples, Mode::InPlane), 22.5), "in-plane uses K");
    require(close(at(Vec3 {0.0, 24.25, 0.0}, normalX, c, 20.0, samples, Mode::None), 20.0), "none is the cell value");
    require(close(at(Vec3 {0.0, 24.25, 0.0}, normalX, c, 20.0, {}, Mode::InPlane), 20.0), "no usable side is flat");

    // No overshoot when both neighbours pull the same way.
    const std::vector<AxisSample> both {
        {1, Vec3 {0.0, 10.0, 0.0}, iso(10.0)},
        {2, Vec3 {0.0, 0.0, 10.0}, iso(10.0)},
    };
    require(at(Vec3 {0.0, 9.0, 9.0}, normalX, c, 20.0, both, Mode::InPlane) >= 10.0 - 1e-12, "bounded");

    // A neighbour mainly across the plane is not used.
    const std::vector<AxisSample> across {{0, Vec3 {10.0, 1.0, 0.0}, iso(99.0)}};
    require(close(at(Vec3 {0.0, 1.0, 0.0}, normalX, c, 20.0, across, Mode::InPlane), 20.0), "across the plane ignored");

    // Linear fields are reproduced on a tilted plane with non-orthogonal axes.
    Vec3 tilted {1.0, 1.0, 0.3};
    tilted /= tilted.two_norm();
    const auto field = [](const Vec3& x) { return 100.0 + 2.0 * x[0] - 3.0 * x[1] + 0.5 * x[2]; };
    const Vec3 dj {-30.0, 40.0, 3.0};
    const Vec3 dk {-1.0, 0.5, 8.0};
    const std::vector<AxisSample> skew {{1, dj, iso(field(dj))}, {2, dk, iso(field(dk))}};
    const auto project = [&tilted](Vec3 v) {
        v.axpy(-(v * tilted), tilted);
        return v;
    };
    Vec3 x = project(dj);
    x *= 0.2;
    Vec3 x2 = project(dk);
    x2 *= 0.3;
    x += x2;
    // Value along the plane: neighbour offsets are projected, so compare with the same projection.
    const double slope1 = (field(dj) - field(c)) * 0.2;
    const double slope2 = (field(dk) - field(c)) * 0.3;
    require(close(at(x, tilted, c, field(c), skew, Mode::InPlane), field(c) + slope1 + slope2), "skew axes");

    // Continuous across a clean face between unequal cells.
    const Vec3 ca {0.0, 96.5, 0.0};
    const Vec3 cb {0.0, 0.0, 0.0};
    const Vec3 face {0.0, 44.0, 0.0};
    const double fromB = at(face, normalX, cb, 20.0, {{1, ca, iso(50.0)}}, Mode::InPlane);
    const double fromA = at(face, normalX, ca, 50.0, {{1, cb, iso(20.0)}}, Mode::InPlane);
    require(close(fromB, fromA), "continuous at a clean face");
}

constexpr int nx = 4;
constexpr int ny = 1;
constexpr int nz = 3;

// Layers 10 m thick; the columns from i = 2 on are moved down by `slip`.
Dune::CpGrid makeGrid(double slip, const std::vector<int>& actnum, std::vector<double>& coord, std::vector<double>& zcorn)
{
    coord.clear();
    zcorn.assign(8 * nx * ny * nz, 0.0);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = 100.0 * i;
            const double y = 100.0 * j;
            coord.insert(coord.end(), {x, y, 0.0, x, y, 100.0});
        }
    }
    for (int k = 0; k < 2 * nz; ++k) {
        for (int j = 0; j < 2 * ny; ++j) {
            for (int i = 0; i < 2 * nx; ++i) {
                const int cellI = i / 2;
                const double depth = 10.0 * ((k + 1) / 2);
                zcorn[k * 4 * nx * ny + j * 2 * nx + i] = depth + (cellI >= 2 ? slip : 0.0);
            }
        }
    }
    grdecl raw;
    raw.dims[0] = nx;
    raw.dims[1] = ny;
    raw.dims[2] = nz;
    raw.coord = coord.data();
    raw.zcorn = zcorn.data();
    raw.actnum = actnum.data();
    Dune::CpGrid grid;
    grid.processEclipseFormat(raw, false);
    return grid;
}

// Kinds of the I+ side of (1,0,1), the cell next to the slip surface.
CellSide sideAtSlip(const Dune::CpGrid& grid, const DeclaredFaults& faults, bool flowing, bool acrossFaults)
{
    const auto gv = grid.leafGridView();
    const auto& global = grid.globalCell();
    const int wanted = 1 + nx * ny * 1;
    for (const auto& element : elements(gv)) {
        const int c = gv.indexSet().index(element);
        if (global[c] != wanted) {
            continue;
        }
        const auto logical = [&grid](int cc) {
            std::array<int, 3> ijk;
            grid.getIJK(cc, ijk);
            return ijk;
        };
        const auto cartesian = [&global](int cc) { return std::size_t(global[cc]); };
        const auto flows = [flowing](int, int) { return flowing; };
        const auto hasMech = [](int) { return true; };
        return classifySides(gv, element, logical, cartesian, faults, flows, hasMech, acrossFaults)[1];
    }
    throw std::runtime_error("cell not found");
}

void testFaults()
{
    std::vector<double> coord;
    std::vector<double> zcorn;
    const std::vector<int> all(nx * ny * nz, 1);
    const DeclaredFaults none;
    const int cellLeft = 1 + nx * ny * 1;
    const int cellRight = 2 + nx * ny * 1;

    {
        auto grid = makeGrid(0.0, all, coord, zcorn);
        require(sideAtSlip(grid, none, true, false).kind == SideKind::Clean, "no fault is clean");

        DeclaredFaults left;
        left.add(cellLeft, faceTag(Opm::FaceDir::XPlus));
        require(sideAtSlip(grid, left, true, false).kind == SideKind::DeclaredFault, "FAULTS, zero throw, this cell");

        DeclaredFaults right;
        right.add(cellRight, faceTag(Opm::FaceDir::XMinus));
        require(sideAtSlip(grid, right, true, false).kind == SideKind::DeclaredFault, "FAULTS, zero throw, other cell");

        DeclaredFaults elsewhere;
        elsewhere.add(cellLeft, faceTag(Opm::FaceDir::XMinus));
        require(sideAtSlip(grid, elsewhere, true, false).kind == SideKind::Clean, "fault on the other side of the cell");

        require(sideAtSlip(grid, none, false, false).kind == SideKind::NoTransmissibility, "flow barrier");
        require(sideAtSlip(grid, left, true, true).kind == SideKind::Clean, "option: across a declared fault");
    }
    {
        // Half a layer of throw: two neighbours on the side.
        auto grid = makeGrid(5.0, all, coord, zcorn);
        const auto stopped = sideAtSlip(grid, none, true, false);
        require(stopped.kind == SideKind::NonConforming, "half-layer throw");
        require(stopped.neighbours.size() == 2, "half-layer throw has two neighbours");

        const auto across = sideAtSlip(grid, none, true, true);
        require(across.kind == SideKind::Clean, "option: across a throw");
        require(close(across.neighbours[0].area, across.neighbours[1].area), "equal overlap areas");

        // Weighted by overlap, so a skewed neighbour pair gives their mean here.
        CellSides sides;
        sides[1] = across;
        const auto samples = samplesFromSides(
            sides, [](int c) { return Vec3 {double(c), 0.0, 0.0}; }, [&across](int c) { return iso(c == across.neighbours[0].cell ? 10.0 : 20.0); });
        require(samples.size() == 1 && samples[0].axis == 0, "one sample on the I axis");
        require(close(traction(samples[0].stress), 15.0), "area-weighted neighbour stress");
    }
    {
        // A whole layer of throw: one neighbour, but in another layer.
        auto grid = makeGrid(10.0, all, coord, zcorn);
        const auto side = sideAtSlip(grid, none, true, false);
        require(side.kind == SideKind::LogicalOffset, "full-layer throw");
        require(side.neighbours.size() == 1, "full-layer throw has one neighbour");
    }
    {
        // Only the lower half of the side has a neighbour left.
        auto actnum = all;
        actnum[2 + nx * ny * 1] = 0;
        auto grid = makeGrid(5.0, actnum, coord, zcorn);
        const auto side = sideAtSlip(grid, none, true, false);
        require(side.kind != SideKind::Clean, "partly inactive side is not clean");
    }
    {
        // Neighbour inactive, no throw.
        auto actnum = all;
        actnum[2 + nx * ny * 1] = 0;
        auto grid = makeGrid(0.0, actnum, coord, zcorn);
        require(sideAtSlip(grid, none, true, false).kind == SideKind::Boundary, "inactive neighbour");
    }
}

void testFaultSidesStayFlat()
{
    const auto centreOf = [](int) { return Vec3(1.0); };
    const auto stressOf = [](int) { return Sym6(1.0); };
    for (const auto kind : {SideKind::NonConforming, SideKind::LogicalOffset, SideKind::DeclaredFault,
                            SideKind::NoTransmissibility}) {
        CellSides sides;
        sides[0] = CellSide{kind, {Neighbour{1, 1.0}}};
        require(isFaultSide(kind), "fault kind");
        require(samplesFromSides(sides, centreOf, stressOf).empty(), "fault side stays flat");
    }
    CellSides sides;
    sides[0] = CellSide{SideKind::Clean, {Neighbour{1, 1.0}}};
    require(samplesFromSides(sides, centreOf, stressOf).size() == 1, "clean side is used");
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    try {
        testOperator();
        testFaults();
        testFaultSidesStayFlat();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "ok" << std::endl;
    return 0;
}
