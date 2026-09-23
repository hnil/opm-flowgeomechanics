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

// Check the face orientation of a grid with an LGR: node order, stored normal
// and face-to-cell order must all agree with the geometry, on every level.
#include "config.h"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace
{

struct Tensor
{
    std::array<int, 3> dims;
    std::vector<double> coord, zcorn;
    std::vector<int> actnum;
};

Tensor tensorGrid(int nx, int ny, int nz, double dx, double dz)
{
    Tensor g;
    g.dims = {nx, ny, nz};
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            for (double z : {0.0, nz*dz}) {
                g.coord.insert(g.coord.end(), {i*dx, j*dx, z});
            }
        }
    }
    g.zcorn.assign(8ull*nx*ny*nz, 0.0);
    for (int k = 0; k < nz; ++k) {
        for (int dk = 0; dk < 2; ++dk) {
            for (int j = 0; j < 2*ny; ++j) {
                for (int i = 0; i < 2*nx; ++i) {
                    g.zcorn[i + 2ull*nx*j + 4ull*nx*ny*(2*k + dk)] = (k + dk)*dz;
                }
            }
        }
    }
    g.actnum.assign(1ull*nx*ny*nz, 1);
    return g;
}

double dot(const Dune::FieldVector<double, 3>& a, const Dune::FieldVector<double, 3>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Counts on the leaf, through the CpGrid index interface the VEM uses.
void checkLeaf(const Dune::CpGrid& grid)
{
    const auto& gv = grid.leafGridView();
    std::vector<int> level(gv.size(0));
    for (const auto& e : elements(gv)) {
        level[gv.indexSet().index(e)] = e.level();
    }

    // per face: normal vs node order, normal vs faceCell order
    int nodeVsNormal = 0, normalVsCells = 0, faces = 0;
    int kinds[3][3] = {};    // [level of cell0 side][level of cell1 side]
    int shown = 0;
    for (int f = 0; f < grid.numFaces(); ++f) {
        const int nv = grid.numFaceVertices(f);
        const auto& c = grid.faceCentroid(f);
        Dune::FieldVector<double, 3> poly(0.0);
        for (int v = 0; v < nv; ++v) {
            auto a = grid.vertexPosition(grid.faceVertex(f, v));
            auto b = grid.vertexPosition(grid.faceVertex(f, (v + 1) % nv));
            a -= c;
            b -= c;
            poly[0] += a[1]*b[2] - a[2]*b[1];
            poly[1] += a[2]*b[0] - a[0]*b[2];
            poly[2] += a[0]*b[1] - a[1]*b[0];
        }
        const auto& n = grid.faceNormal(f);
        const int c0 = grid.faceCell(f, 0), c1 = grid.faceCell(f, 1);
        // the normal should point from c0 to c1
        Dune::FieldVector<double, 3> d(0.0);
        if (c0 >= 0) { d = c; d -= grid.cellCentroid(c0); }
        else         { d = grid.cellCentroid(c1); d -= c; }
        const bool badNodes = dot(poly, n) < 0.0;
        const bool badCells = dot(n, d) < 0.0;
        nodeVsNormal += badNodes;
        normalVsCells += badCells;
        ++faces;
        if (badNodes || badCells) {
            const int l0 = c0 >= 0 ? std::min(level[c0], 1) + 1 : 0;
            const int l1 = c1 >= 0 ? std::min(level[c1], 1) + 1 : 0;
            ++kinds[l0][l1];
            if (shown++ < 8) {
                std::cout << "  face " << f << " nodes " << nv << " cells " << c0 << "/" << c1
                          << " (levels " << (c0 >= 0 ? level[c0] : -1) << "/"
                          << (c1 >= 0 ? level[c1] : -1) << ")"
                          << (badNodes ? " node order against normal" : "")
                          << (badCells ? " normal against faceCell order" : "")
                          << "  normal " << n << '\n';
            }
        }
    }
    std::cout << "leaf: " << faces << " faces, node order against normal: " << nodeVsNormal
              << ", normal against faceCell order: " << normalVsCells << '\n';
    const char* name[3] = {"boundary", "coarse", "refined"};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            if (kinds[a][b]) {
                std::cout << "  " << name[a] << " | " << name[b] << ": " << kinds[a][b] << '\n';
            }
        }
    }
}

// Through the Dune interface, on any grid view: the outer normal of each
// intersection has to point away from the cell.
template <class GV>
void checkView(const GV& gv, const std::string& what)
{
    int bad = 0, total = 0;
    for (const auto& e : elements(gv)) {
        const auto cc = e.geometry().center();
        for (const auto& is : intersections(gv, e)) {
            auto d = is.geometry().center();
            d -= cc;
            bad += dot(is.centerUnitOuterNormal(), d) < 0.0;
            ++total;
        }
    }
    std::cout << what << ": " << total << " intersections, outer normal pointing in: " << bad
              << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    // SIMPLE-like refinement: a 3x3x3 box of an 11x11x25 tensor grid, 9x9x3 each.
    const int nx = 11, ny = 11, nz = 25;
    auto g = tensorGrid(nx, ny, nz, 181.0, 10.0);
    grdecl input{};
    input.dims[0] = nx; input.dims[1] = ny; input.dims[2] = nz;
    input.coord = g.coord.data();
    input.zcorn = g.zcorn.data();
    input.actnum = g.actnum.data();

    const std::array<int, 3> refine {9, 9, 3};
    const std::array<int, 3> lo {4, 4, 14}, hi {7, 7, 17};
    const bool edgeConformal = argc > 1 && std::string(argv[1]) == "ec";

    Dune::CpGrid grid;
    grid.processEclipseFormat(input, false, false, edgeConformal);
    checkLeaf(grid);
    grid.addLgrsUpdateLeafView({refine}, {lo}, {hi}, {"LGR1"});

    for (int l = 0; l <= grid.maxLevel(); ++l) {
        checkView(grid.levelGridView(l), "level " + std::to_string(l));
    }
    checkView(grid.leafGridView(), "leaf");
    checkLeaf(grid);
    return 0;
}
