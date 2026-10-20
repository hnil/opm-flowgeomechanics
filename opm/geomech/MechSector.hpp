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
#ifndef OPM_MECH_SECTOR_HPP
#define OPM_MECH_SECTOR_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <utility>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

/// A mechanics grid cut out of the deck's corner-point grid and padded, built at
/// run time so that the deck itself carries no padding.
namespace Opm::MechSector
{

struct Grdecl
{
    std::array<int,3> dims{};
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<double> bend;    // per pillar: vertical above bend[2p], below bend[2p+1]; empty: straight
};

/// `cells` cells over `size` metres, each `growth` times the one before, the
/// smallest next to the box. Vertically `depth` (if set) replaces `size`: the
/// padding then reaches that depth in every column.
struct Padding
{
    int cells = 0;
    double size = 0.0;
    double growth = 1.0;
    double depth = std::numeric_limits<double>::quiet_NaN();
    bool verticalPillars = false;    // lateral only
};

struct Spec
{
    std::array<int,6> box{};    // I1 I2 J1 J2 K1 K2, zero-based, end exclusive
    bool boxFromFlow = false;   // box = the flow cells' bounding box plus `margin`
    std::array<int,3> margin{};
    Padding lateral, top, bottom;
    bool straighten = false;
    bool bendPillars = false;   // box pillars vertical above and below the box
    std::string dump;           // GRDECL file of the generated grid, if set
};

/// Vertical pillars through each pillar's point at the box's mean depth; the
/// depths stay. Returns the largest horizontal move of a corner.
inline double straighten(Grdecl& g)
{
    double zref = 0.0;
    for (const double z : g.zcorn) {
        zref += z;
    }
    zref /= g.zcorn.size();
    const auto [zmin, zmax] = std::minmax_element(g.zcorn.begin(), g.zcorn.end());
    double moved = 0.0;
    for (std::size_t p = 0; p < g.coord.size(); p += 6) {
        double* c = &g.coord[p];
        const double dz = c[5] - c[2];
        if (std::abs(dz) < 1e-12) {
            continue;
        }
        const auto at = [c, dz](double z, int a) { return c[a] + (z - c[2])/dz*(c[3 + a] - c[a]); };
        const double x = at(zref, 0), y = at(zref, 1);
        for (const double z : {*zmin, *zmax}) {
            moved = std::max(moved, std::hypot(at(z, 0) - x, at(z, 1) - y));
        }
        c[0] = c[3] = x;
        c[1] = c[4] = y;
    }
    return moved;
}

/// The box cell (I1,J1,K1) sits at `origin` in the result.
struct Result
{
    Grdecl grid;
    std::array<int,3> origin{};
    double moved = 0.0;    // largest corner move by straightening
    int folded = 0;        // columns whose pillars cross or fold within the grid's depths
};

inline std::size_t zIndex(const std::array<int,3>& d, int i, int j, int k,
                          int ci, int cj, int ck)
{
    const std::size_t nx = d[0], ny = d[1];
    return (2*i + ci) + 2*nx*(2*j + cj) + 4*nx*ny*(2*std::size_t(k) + ck);
}

inline std::vector<double> widths(const Padding& p)
{
    std::vector<double> w(p.cells);
    if (p.cells == 0) {
        return w;
    }
    double sum = 0.0;
    for (int m = 0; m < p.cells; ++m) {
        w[m] = std::pow(p.growth, m);
        sum += w[m];
    }
    for (auto& x : w) {
        x /= sum;    // fractions; scaled by the size, or by each column's distance to `depth`
    }
    return w;
}

inline Grdecl extractBox(const Grdecl& in, const std::array<int,6>& b)
{
    Grdecl out;
    out.dims = {b[1] - b[0], b[3] - b[2], b[5] - b[4]};
    const int px = in.dims[0] + 1;
    for (int j = b[2]; j <= b[3]; ++j) {
        for (int i = b[0]; i <= b[1]; ++i) {
            const auto* p = &in.coord[6*(i + px*j)];
            out.coord.insert(out.coord.end(), p, p + 6);
        }
    }
    out.zcorn.resize(8*std::size_t(out.dims[0])*out.dims[1]*out.dims[2]);
    for (int k = 0; k < out.dims[2]; ++k)
    for (int j = 0; j < out.dims[1]; ++j)
    for (int i = 0; i < out.dims[0]; ++i)
    for (int c = 0; c < 8; ++c) {
        const int ci = c & 1, cj = (c >> 1) & 1, ck = c >> 2;
        out.zcorn[zIndex(out.dims, i, j, k, ci, cj, ck)] =
            in.zcorn[zIndex(in.dims, i + b[0], j + b[2], k + b[4], ci, cj, ck)];
    }
    return out;
}

/// The grid as GRDECL, with one extra integer array per entry of `ints`.
inline void writeGrdecl(const std::string& file, const Grdecl& g, const std::vector<int>& actnum,
                        const std::vector<std::pair<std::string, std::vector<int>>>& ints)
{
    std::ofstream os(file);
    if (!os) {
        OPM_THROW(std::runtime_error, "Cannot write " + file);
    }
    os.precision(10);
    const auto write = [&os](const std::string& name, const auto& v) {
        os << name << "\n";
        for (std::size_t i = 0; i < v.size(); ++i) {
            os << v[i] << ((i % 6 == 5) ? "\n" : " ");
        }
        os << "/\n\n";
    };
    os << "SPECGRID\n" << g.dims[0] << " " << g.dims[1] << " " << g.dims[2] << " 1 F /\n\n";
    write("COORD", g.coord);
    write("ZCORN", g.zcorn);
    write("ACTNUM", actnum);
    for (const auto& [name, v] : ints) {
        write(name, v);
    }
}

inline std::array<double,2> pillarAt(const double* c, double z, const double* bend = nullptr)
{
    if (bend != nullptr) {
        if (bend[0] == bend[0]) {
            z = std::max(z, bend[0]);
        }
        if (bend[1] == bend[1]) {
            z = std::min(z, bend[1]);
        }
    }
    const double dz = c[5] - c[2];
    const double t = std::abs(dz) < 1e-12 ? 0.0 : (z - c[2])/dz;
    return {c[0] + t*(c[3] - c[0]), c[1] + t*(c[4] - c[1])};
}

/// Pad along axis 0 (I) or 1 (J). New pillars are the boundary pillars moved
/// outward horizontally; the boundary layers continue flat, so the faces
/// against the box keep their corners. With vertical pillars, the first one
/// clears the boundary pillar by its full offset anywhere in [zlo, zhi].
inline Grdecl padLateral(const Grdecl& in, int axis, const Padding& pad,
                         double zlo, double zhi)
{
    const auto frac = widths(pad);
    const int n = pad.cells;
    const auto& d = in.dims;
    Grdecl out;
    out.dims = d;
    out.dims[axis] += 2*n;

    // Distance of new pillar q (outside the box) from the boundary pillar.
    std::vector<double> dist(n + 1, 0.0);
    for (int m = 0; m < n; ++m) {
        dist[m + 1] = dist[m] + frac[m]*pad.size;
    }

    double zref = 0.0;
    for (const double z : in.zcorn) {
        zref += z;
    }
    zref /= in.zcorn.size();

    const int pxIn = d[0] + 1, pxOut = out.dims[0] + 1;
    const auto bendOf = [&in](int p) { return in.bend.empty() ? nullptr : &in.bend[2*p]; };
    out.coord.resize(6*std::size_t(pxOut)*(out.dims[1] + 1));
    if (!in.bend.empty()) {
        out.bend.assign(2*std::size_t(pxOut)*(out.dims[1] + 1),
                        std::numeric_limits<double>::quiet_NaN());
    }
    for (int pj = 0; pj <= out.dims[1]; ++pj)
    for (int pi = 0; pi <= out.dims[0]; ++pi) {
        std::array<int,2> q{pi, pj};
        std::array<int,2> s = q;
        s[axis] = q[axis] - n;
        double shift = 0.0;
        std::array<int,2> inner = s;
        if (s[axis] < 0) {
            shift = dist[-s[axis]];
            s[axis] = 0;
            inner = s; inner[axis] = 1;
        } else if (s[axis] > d[axis]) {
            shift = dist[s[axis] - d[axis]];
            s[axis] = d[axis];
            inner = s; inner[axis] = d[axis] - 1;
        }
        const int ps_ = s[0] + pxIn*s[1];
        const double* src = &in.coord[6*ps_];
        double* dst = &out.coord[6*(q[0] + pxOut*q[1])];
        std::copy(src, src + 6, dst);
        if (!in.bend.empty() && !(shift > 0.0 && pad.verticalPillars)) {
            std::copy_n(&in.bend[2*ps_], 2, &out.bend[2*(q[0] + pxOut*q[1])]);
        }
        if (shift > 0.0) {
            const double* nb = &in.coord[6*(inner[0] + pxIn*inner[1])];
            const int pn_ = inner[0] + pxIn*inner[1];
            const auto ps = pillarAt(src, zref, bendOf(ps_)), pn = pillarAt(nb, zref, bendOf(pn_));
            double ux = pad.verticalPillars ? ps[0] - pn[0] : src[0] - nb[0];
            double uy = pad.verticalPillars ? ps[1] - pn[1] : src[1] - nb[1];
            const double len = std::hypot(ux, uy);
            if (!(len > 0.0)) {
                OPM_THROW(std::runtime_error, "Mechanics padding: zero pillar spacing at the box edge");
            }
            if (pad.verticalPillars) {
                ux /= len;
                uy /= len;
                double drift = 0.0;
                for (const double z : {zlo, zhi}) {
                    const auto p = pillarAt(src, z, bendOf(ps_));
                    drift = std::max(drift, (p[0] - ps[0])*ux + (p[1] - ps[1])*uy);
                }
                dst[0] = dst[3] = ps[0] + (drift + shift)*ux;
                dst[1] = dst[4] = ps[1] + (drift + shift)*uy;
            } else {
                ux *= shift/len;
                uy *= shift/len;
                dst[0] += ux; dst[1] += uy;
                dst[3] += ux; dst[4] += uy;
            }
        }
    }

    out.zcorn.resize(8*std::size_t(out.dims[0])*out.dims[1]*out.dims[2]);
    for (int k = 0; k < out.dims[2]; ++k)
    for (int j = 0; j < out.dims[1]; ++j)
    for (int i = 0; i < out.dims[0]; ++i) {
        std::array<int,2> c{i, j};
        std::array<int,2> s = c;
        s[axis] = c[axis] - n;
        int corner = -1;    // the boundary corner every padding corner copies
        if (s[axis] < 0) {
            s[axis] = 0; corner = 0;
        } else if (s[axis] >= d[axis]) {
            s[axis] = d[axis] - 1; corner = 1;
        }
        for (int cc = 0; cc < 8; ++cc) {
            std::array<int,3> cn{cc & 1, (cc >> 1) & 1, cc >> 2};
            std::array<int,3> sc = cn;
            if (corner >= 0) {
                sc[axis] = corner;
            }
            out.zcorn[zIndex(out.dims, i, j, k, cn[0], cn[1], cn[2])] =
                in.zcorn[zIndex(d, s[0], s[1], k, sc[0], sc[1], sc[2])];
        }
    }
    return out;
}

/// Layers above and below that follow the box's top and bottom surfaces.
inline Grdecl padVertical(const Grdecl& in, const Padding& top, const Padding& bottom)
{
    const auto ft = widths(top), fb = widths(bottom);
    const int nt = top.cells, nb = bottom.cells;
    const auto& d = in.dims;
    Grdecl out;
    out.dims = {d[0], d[1], d[2] + nt + nb};
    out.zcorn.resize(8*std::size_t(out.dims[0])*out.dims[1]*out.dims[2]);

    const auto thickness = [](const Padding& p, double z, bool up) {
        if (std::isnan(p.depth)) {
            return p.size;
        }
        const double t = up ? z - p.depth : p.depth - z;
        if (!(t > 0.0)) {
            OPM_THROW(std::runtime_error, "Mechanics padding: depth " + std::to_string(p.depth)
                      + " is not " + (up ? "above" : "below") + " the box at depth "
                      + std::to_string(z));
        }
        return t;
    };

    // Per pillar, the shallowest box top and deepest box bottom of the columns
    // around it. Only the first padding layer follows each column; the layers
    // beyond it are shared by all columns at a pillar, so faults stop there.
    const int px = d[0] + 1, py = d[1] + 1;
    std::vector<double> topP(px*py, std::numeric_limits<double>::max());
    std::vector<double> botP(px*py, std::numeric_limits<double>::lowest());
    for (int j = 0; j < d[1]; ++j)
    for (int i = 0; i < d[0]; ++i)
    for (int ci = 0; ci < 2; ++ci)
    for (int cj = 0; cj < 2; ++cj) {
        const int p = (i + ci) + px*(j + cj);
        topP[p] = std::min(topP[p], in.zcorn[zIndex(d, i, j, 0, ci, cj, 0)]);
        botP[p] = std::max(botP[p], in.zcorn[zIndex(d, i, j, d[2] - 1, ci, cj, 1)]);
    }
    // Boundaries of the shared layers at each pillar: z[m] is the far side of layer m.
    const auto shared = [&](const std::vector<double>& surf, const Padding& pad,
                            const std::vector<double>& f, bool up) {
        std::vector<std::vector<double>> z(surf.size());
        for (std::size_t p = 0; p < surf.size(); ++p) {
            const double t = thickness(pad, surf[p], up);
            double zz = surf[p];
            for (int m = 0; m < pad.cells; ++m) {
                zz += (up ? -1.0 : 1.0)*f[m]*t;
                z[p].push_back(zz);
            }
        }
        return z;
    };
    const auto zt = nt > 0 ? shared(topP, top, ft, true) : std::vector<std::vector<double>>{};
    const auto zb = nb > 0 ? shared(botP, bottom, fb, false) : std::vector<std::vector<double>>{};

    for (int j = 0; j < d[1]; ++j)
    for (int i = 0; i < d[0]; ++i)
    for (int ci = 0; ci < 2; ++ci)
    for (int cj = 0; cj < 2; ++cj) {
        for (int k = 0; k < d[2]; ++k) {
            for (int ck = 0; ck < 2; ++ck) {
                out.zcorn[zIndex(out.dims, i, j, k + nt, ci, cj, ck)] =
                    in.zcorn[zIndex(d, i, j, k, ci, cj, ck)];
            }
        }
        const int p = (i + ci) + px*(j + cj);
        double z = in.zcorn[zIndex(d, i, j, 0, ci, cj, 0)];
        for (int m = 0; m < nt; ++m) {    // m = 0 touches the box
            const int k = nt - 1 - m;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 1)] = z;
            z = zt[p][m];
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 0)] = z;
        }
        z = in.zcorn[zIndex(d, i, j, d[2] - 1, ci, cj, 1)];
        for (int m = 0; m < nb; ++m) {
            const int k = nt + d[2] + m;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 0)] = z;
            z = zb[p][m];
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 1)] = z;
        }
    }

    // Pillars must span the new depths.
    const auto [zmin, zmax] = std::minmax_element(out.zcorn.begin(), out.zcorn.end());
    out.coord = in.coord;
    out.bend = in.bend;
    for (std::size_t p = 0; p < out.coord.size(); p += 6) {
        double* c = &out.coord[p];
        const double dz = c[5] - c[2];
        if (std::abs(dz) < 1e-12) {
            continue;
        }
        const double z1 = *zmin - 1.0, z2 = *zmax + 1.0;
        const double t1 = (z1 - c[2])/dz, t2 = (z2 - c[2])/dz;
        const std::array<double,6> old{c[0], c[1], c[2], c[3], c[4], c[5]};
        for (int a = 0; a < 3; ++a) {
            c[a] = old[a] + t1*(old[3 + a] - old[a]);
            c[3 + a] = old[a] + t2*(old[3 + a] - old[a]);
        }
    }
    return out;
}

/// Per pillar, the shallowest top and deepest bottom of the columns around it:
/// where padVertical's first layer ends, so the pillars bend on a fault-free surface.
inline std::vector<double> boxBends(const Grdecl& g)
{
    const auto& d = g.dims;
    const int px = d[0] + 1, py = d[1] + 1;
    std::vector<double> b(2*std::size_t(px)*py);
    for (int p = 0; p < px*py; ++p) {
        b[2*p] = std::numeric_limits<double>::max();
        b[2*p + 1] = std::numeric_limits<double>::lowest();
    }
    for (int j = 0; j < d[1]; ++j)
    for (int i = 0; i < d[0]; ++i)
    for (int ci = 0; ci < 2; ++ci)
    for (int cj = 0; cj < 2; ++cj) {
        const int p = (i + ci) + px*(j + cj);
        b[2*p] = std::min(b[2*p], g.zcorn[zIndex(d, i, j, 0, ci, cj, 0)]);
        b[2*p + 1] = std::max(b[2*p + 1], g.zcorn[zIndex(d, i, j, d[2] - 1, ci, cj, 1)]);
    }
    return b;
}

/// Columns whose four pillars, taken at the column's own depths, do not keep
/// the orientation they have at the column's mid-depth.
inline int foldedColumns(const Grdecl& g)
{
    const auto& d = g.dims;
    const int px = d[0] + 1;
    int folded = 0;
    for (int j = 0; j < d[1]; ++j)
    for (int i = 0; i < d[0]; ++i) {
        double zmin = std::numeric_limits<double>::max(), zmax = std::numeric_limits<double>::lowest();
        for (int k = 0; k < d[2]; ++k)
        for (int c = 0; c < 8; ++c) {
            const double z = g.zcorn[zIndex(d, i, j, k, c & 1, (c >> 1) & 1, c >> 2)];
            zmin = std::min(zmin, z);
            zmax = std::max(zmax, z);
        }
        const std::array<int,4> ring{i + px*j, i + 1 + px*j, i + 1 + px*(j + 1), i + px*(j + 1)};
        const auto corners = [&](double z, std::array<double,4>& cr) {
            std::array<std::array<double,2>,4> p;
            for (int a = 0; a < 4; ++a) {
                p[a] = pillarAt(&g.coord[6*ring[a]], z,
                                g.bend.empty() ? nullptr : &g.bend[2*ring[a]]);
            }
            for (int a = 0; a < 4; ++a) {
                const auto& o = p[a];
                const auto& n = p[(a + 1) % 4];
                const auto& b = p[(a + 3) % 4];
                cr[a] = (n[0] - o[0])*(b[1] - o[1]) - (n[1] - o[1])*(b[0] - o[0]);
            }
        };
        std::array<double,4> mid, cr;
        corners(0.5*(zmin + zmax), mid);
        bool bad = false;
        for (int s = 0; s <= 4 && !bad; ++s) {
            corners(zmin + 0.25*s*(zmax - zmin), cr);
            for (int a = 0; a < 4; ++a) {
                bad = bad || cr[a]*mid[0] <= 0.0;
            }
        }
        folded += bad;
    }
    return folded;
}

inline Result build(const Grdecl& deck, const Spec& spec)
{
    Result r;
    r.grid = extractBox(deck, spec.box);
    if (spec.straighten) {
        r.moved = straighten(r.grid);
    }
    if (spec.bendPillars) {
        r.grid.bend = boxBends(r.grid);
    }
    if (spec.lateral.cells > 0) {
        // Depth range of the finished grid, which the vertical pillars must clear.
        const auto [bmin, bmax] = std::minmax_element(r.grid.zcorn.begin(), r.grid.zcorn.end());
        const auto reach = [](const Padding& p, double z, double sign) {
            return p.cells == 0 ? z : (std::isnan(p.depth) ? z + sign*p.size : p.depth);
        };
        const double zlo = reach(spec.top, *bmin, -1.0), zhi = reach(spec.bottom, *bmax, 1.0);
        r.grid = padLateral(r.grid, 0, spec.lateral, zlo, zhi);
        r.grid = padLateral(r.grid, 1, spec.lateral, zlo, zhi);
    }
    if (spec.top.cells > 0 || spec.bottom.cells > 0) {
        r.grid = padVertical(r.grid, spec.top, spec.bottom);
    }
    r.origin = {spec.lateral.cells, spec.lateral.cells, spec.top.cells};
    r.folded = foldedColumns(r.grid);
    return r;
}

/// The JSON spec. "box" is I1 I2 J1 J2 K1 K2, one-based and inclusive as in a
/// deck; without it the box is the whole deck grid.
inline Spec readSpec(const std::string& file, const std::array<int,3>& deckDims)
{
    const PropertyTree tree(file);
    Spec spec;
    spec.box = {0, deckDims[0], 0, deckDims[1], 0, deckDims[2]};
    spec.boxFromFlow = tree.get<bool>("box_from_flow", false);
    if (const auto margin = tree.get_child_items_as_vector<int>("margin")) {
        if (margin->size() != 3) {
            OPM_THROW(std::runtime_error, "Mechanics sector: \"margin\" needs three numbers");
        }
        spec.margin = {(*margin)[0], (*margin)[1], (*margin)[2]};
    }
    spec.straighten = tree.get<bool>("straighten", false);
    spec.bendPillars = tree.get<bool>("bend_pillars", false);
    spec.dump = tree.get<std::string>("dump", "");
    if (const auto box = tree.get_child_items_as_vector<int>("box")) {
        if (box->size() != 6) {
            OPM_THROW(std::runtime_error, "Mechanics sector: \"box\" needs six numbers");
        }
        for (int a = 0; a < 3; ++a) {
            const int lo = (*box)[2*a] - 1, hi = (*box)[2*a + 1];
            if (lo < 0 || hi > deckDims[a] || lo >= hi) {
                OPM_THROW(std::runtime_error, "Mechanics sector: \"box\" is outside the grid");
            }
            spec.box[2*a] = lo;
            spec.box[2*a + 1] = hi;
        }
    }
    const auto padding = [&tree](const std::string& name, bool vertical) {
        Padding p;
        if (!tree.get_child_optional(name)) {
            return p;
        }
        p.cells = tree.get<int>(name + ".cells", 0);
        p.size = tree.get<double>(name + (vertical ? ".thickness" : ".width"), 0.0);
        p.growth = tree.get<double>(name + ".growth", 1.0);
        if (vertical && tree.get_child_optional(name + ".depth")) {
            p.depth = tree.get<double>(name + ".depth");
        }
        if (!vertical) {
            p.verticalPillars = tree.get<bool>(name + ".vertical_pillars", false);
        }
        if (p.cells < 0 || !(p.growth > 0.0)
            || (p.cells > 0 && !(p.size > 0.0) && std::isnan(p.depth))) {
            OPM_THROW(std::runtime_error, "Mechanics sector: bad \"" + name + "\" padding");
        }
        return p;
    };
    spec.lateral = padding("lateral", false);
    spec.top = padding("top", true);
    spec.bottom = padding("bottom", true);
    return spec;
}

} // namespace Opm::MechSector

#endif // OPM_MECH_SECTOR_HPP
