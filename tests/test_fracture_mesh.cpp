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

#include <opm/geomech/FractureMesh.hpp>
#include <opm/geomech/FractureMeshFoamGrid.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <cmath>
#include <cstring>
#include <iostream>
#include <span>
#include <string>
#include <vector>

using namespace Opm;

namespace
{
int failures = 0;

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
bool same(const FractureMesh::Vec3& a, const FractureMesh::Vec3& b)
{
    return same(a[0], b[0]) && same(a[1], b[1]) && same(a[2], b[2]);
}
bool same(int a, int b) { return a == b; }
bool same(const FractureMesh::Triangle& a, const FractureMesh::Triangle& b)
{
    return same(a[0], b[0]) && same(a[1], b[1]) && same(a[2], b[2]);
}
template <class T>
bool same(std::span<const T> a, std::span<const T> b)
{
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same(a[i], b[i]))
            return false;
    return true;
}

void check(bool ok, const std::string& what, int i)
{
    if (!ok && failures++ < 20)
        std::cerr << "differs: " << what << " at " << i << '\n';
}

// Both builders must give the same numbering and the same bits.
void compare(const Triangulation& tri, const std::string& name)
{
    const auto grid = RegularTrimesh::createDuneGrid(tri);
    const auto g = buildFractureMesh(*grid);
    const auto t = buildFractureMesh(tri.nodes, tri.triangles);
    const int before = failures;
    check(g.numCells() == t.numCells() && g.numFaces() == t.numFaces() && g.numNodes() == t.numNodes(),
          name + " sizes", 0);
    if (failures > before)
        return;
    for (int c = 0; c < g.numCells(); ++c) {
        check(same(g.cellArea(c), t.cellArea(c)), name + " cellArea", c);
        check(same(g.cellCentroid(c), t.cellCentroid(c)), name + " cellCentroid", c);
        check(same(g.cellNormal(c), t.cellNormal(c)), name + " cellNormal", c);
        check(g.cellLevel(c) == t.cellLevel(c), name + " cellLevel", c);
        check(same(g.cellTriangles(c), t.cellTriangles(c)), name + " cellTriangles", c);
        check(same(g.cellNodes(c), t.cellNodes(c)), name + " cellNodes", c);
        check(same(g.cellCorners(c), t.cellCorners(c)), name + " cellCorners", c);
        check(same(g.cellFaces(c), t.cellFaces(c)), name + " cellFaces", c);
        check(same(g.cellNeighbours(c), t.cellNeighbours(c)), name + " cellNeighbours", c);
        for (int j = 0; j < static_cast<int>(g.cellFaces(c).size()); ++j) {
            check(same(g.cellFaceCentroid(c, j), t.cellFaceCentroid(c, j)), name + " cellFaceCentroid", c);
            check(same(g.cellFaceLength(c, j), t.cellFaceLength(c, j)), name + " cellFaceLength", c);
            check(same(g.cellFaceNormal(c, j), t.cellFaceNormal(c, j)), name + " cellFaceNormal", c);
            check(same(g.cellFaceCorners(c, j)[0], t.cellFaceCorners(c, j)[0])
                      && same(g.cellFaceCorners(c, j)[1], t.cellFaceCorners(c, j)[1]),
                  name + " cellFaceCorners", c);
        }
    }
    for (int f = 0; f < g.numFaces(); ++f) {
        check(same(g.faceCells(f), t.faceCells(f)), name + " faceCells", f);
        check(g.faceNodes(f) == t.faceNodes(f), name + " faceNodes", f);
        check(same(g.faceCentroid(f), t.faceCentroid(f)), name + " faceCentroid", f);
        check(same(g.faceLength(f), t.faceLength(f)), name + " faceLength", f);
    }
    for (int n = 0; n < g.numNodes(); ++n)
        check(same(g.nodeCoord(n), t.nodeCoord(n)), name + " nodeCoord", n);
    std::cout << name << ": " << g.numCells() << " cells, " << g.numFaces() << " faces, "
              << (failures == before ? "identical" : "DIFFER") << '\n';
}
} // namespace

int main()
{
    // large offsets and an oblique plane, so round-off is not trivially exact
    const std::array<double, 3> origin {456789.123, 7123456.789, 2534.5};
    const double s = 1.0 / std::sqrt(2.0);
    const std::array<double, 3> ax1 {s, s, 0.0};
    const std::array<double, 3> ax2 {0.5 * s, 0.5 * s, std::sqrt(3.0) / 2};
    RegularTrimesh mesh(60.0, origin, ax1, ax2, {7.3, 7.3});
    mesh.removeSawtooths();
    mesh.expandGrid(CellRef {9, 0, 0});
    mesh.expandGrid(CellRef {-3, 8, 1});
    mesh.removeSawtooths();
    const auto ring = RegularTrimesh::inner_ring_cells();

    compare(mesh.createTriangulation(1, ring), "seed level 1");
    compare(mesh.createTriangulation(0, {}, true), "fine, smoothing");
    compare(mesh.createTriangulation(3, ring, false, 10, 1), "coarsen 3");
    compare(mesh.createTriangulation(3, ring, true, 4, 2), "coarsen 3, smoothing, band 2");

    std::cout << (failures == 0 ? "PASSED" : "FAILED") << '\n';
    return failures == 0 ? 0 : 1;
}
