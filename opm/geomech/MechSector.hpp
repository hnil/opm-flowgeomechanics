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
};

struct Spec
{
    std::array<int,6> box{};    // I1 I2 J1 J2 K1 K2, zero-based, end exclusive
    Padding lateral, top, bottom;
};

/// The box cell (I1,J1,K1) sits at `origin` in the result.
struct Result
{
    Grdecl grid;
    std::array<int,3> origin{};
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

/// Pad along axis 0 (I) or 1 (J). New pillars are the boundary pillars moved
/// outward horizontally; the boundary layers continue flat, so the faces
/// against the box keep their corners.
inline Grdecl padLateral(const Grdecl& in, int axis, const Padding& pad)
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

    const int pxIn = d[0] + 1, pxOut = out.dims[0] + 1;
    out.coord.resize(6*std::size_t(pxOut)*(out.dims[1] + 1));
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
        const double* src = &in.coord[6*(s[0] + pxIn*s[1])];
        double* dst = &out.coord[6*(q[0] + pxOut*q[1])];
        std::copy(src, src + 6, dst);
        if (shift > 0.0) {
            const double* nb = &in.coord[6*(inner[0] + pxIn*inner[1])];
            double ux = src[0] - nb[0], uy = src[1] - nb[1];
            const double len = std::hypot(ux, uy);
            if (!(len > 0.0)) {
                OPM_THROW(std::runtime_error, "Mechanics padding: zero pillar spacing at the box edge");
            }
            ux *= shift/len;
            uy *= shift/len;
            dst[0] += ux; dst[1] += uy;
            dst[3] += ux; dst[4] += uy;
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
        const double ztop = in.zcorn[zIndex(d, i, j, 0, ci, cj, 0)];
        const double zbot = in.zcorn[zIndex(d, i, j, d[2] - 1, ci, cj, 1)];
        const double tt = nt > 0 ? thickness(top, ztop, true) : 0.0;
        const double tb = nb > 0 ? thickness(bottom, zbot, false) : 0.0;
        double z = ztop;
        for (int m = 0; m < nt; ++m) {    // m = 0 touches the box
            const int k = nt - 1 - m;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 1)] = z;
            z -= ft[m]*tt;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 0)] = z;
        }
        z = zbot;
        for (int m = 0; m < nb; ++m) {
            const int k = nt + d[2] + m;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 0)] = z;
            z += fb[m]*tb;
            out.zcorn[zIndex(out.dims, i, j, k, ci, cj, 1)] = z;
        }
    }

    // Pillars must span the new depths.
    const auto [zmin, zmax] = std::minmax_element(out.zcorn.begin(), out.zcorn.end());
    out.coord = in.coord;
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

inline Result build(const Grdecl& deck, const Spec& spec)
{
    Result r;
    r.grid = extractBox(deck, spec.box);
    if (spec.lateral.cells > 0) {
        r.grid = padLateral(r.grid, 0, spec.lateral);
        r.grid = padLateral(r.grid, 1, spec.lateral);
    }
    if (spec.top.cells > 0 || spec.bottom.cells > 0) {
        r.grid = padVertical(r.grid, spec.top, spec.bottom);
    }
    r.origin = {spec.lateral.cells, spec.lateral.cells, spec.top.cells};
    return r;
}

/// The JSON spec. "box" is I1 I2 J1 J2 K1 K2, one-based and inclusive as in a
/// deck; without it the box is the whole deck grid.
inline Spec readSpec(const std::string& file, const std::array<int,3>& deckDims)
{
    const PropertyTree tree(file);
    Spec spec;
    spec.box = {0, deckDims[0], 0, deckDims[1], 0, deckDims[2]};
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
