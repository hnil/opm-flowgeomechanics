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

#include "config.h"

#include <opm/geomech/FractureMeshFoamGrid.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/geomech/Math.hpp>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <stdexcept>
#include <utility>

namespace Opm
{

FractureMesh
buildFractureMesh(const Dune::FoamGrid<2, 3>& grid)
{
    using Vec3 = FractureMesh::Vec3;
    using GridView = Dune::FoamGrid<2, 3>::LeafGridView;
    const GridView gv = grid.leafGridView();
    const auto& iset = gv.indexSet();
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(gv, Dune::mcmgElementLayout());
    const int nc = gv.size(0);
    const int nf = gv.size(1);
    const int nn = gv.size(2);

    FractureMesh::Arrays a;
    a.node_coord.resize(nn);
    for (const auto& v : vertices(gv)) {
        a.node_coord[iset.index(v)] = v.geometry().center();
    }

    a.cell_area.reserve(nc);
    a.cell_centroid.reserve(nc);
    a.cell_normal.reserve(nc);
    a.cell_level.assign(nc, 0);
    std::vector<std::vector<int>> face_cells(nf);
    std::vector<int> face_halfface(nf, -1);
    std::vector<std::array<int, 2>> face_nodes(nf);

    int c = 0;
    for (const auto& e : elements(gv)) {
        // physics loops over 0..n-1 in place of grid iteration; both orders must agree
        if (static_cast<int>(mapper.index(e)) != c) {
            OPM_THROW(std::logic_error, "buildFractureMesh: element iteration is not in mapper order");
        }
        const auto geom = e.geometry();
        a.cell_area.push_back(geom.volume());
        a.cell_centroid.push_back(geom.center());

        FractureMesh::Triangle tri;
        const int nv = e.subEntities(2);
        if (nv != 3) {
            OPM_THROW(std::logic_error, "buildFractureMesh: only triangle cells are supported");
        }
        for (int i = 0; i < nv; ++i) {
            const int n = iset.subIndex(e, i, 2);
            a.cell_nodes.push_back(n);
            a.cell_corners.push_back(geom.corner(i));
            tri[i] = a.node_coord[n];
        }
        a.cell_nodepos.push_back(a.cell_nodes.size());
        a.triangles.push_back(tri);
        a.cell_tripos.push_back(a.triangles.size());

        // same operations as ddm::normalOfElement
        const Vec3 e01 = tri[1] - tri[0];
        const Vec3 e02 = tri[2] - tri[0];
        Vec3 normal = crossProduct(e01, e02);
        normal /= normal.two_norm();
        a.cell_normal.push_back(normal);

        const auto& ref = Dune::ReferenceElements<double, 2>::general(e.type());
        for (const auto& is : intersections(gv, e)) {
            if (is.neighbor() == is.boundary()) {
                OPM_THROW(std::logic_error, "buildFractureMesh: intersection is neither interior nor boundary");
            }
            const int lf = is.indexInInside();
            const int f = iset.subIndex(e, lf, 1);
            const int hf = a.cell_faces.size();
            a.cell_faces.push_back(f);
            a.cell_neighbours.push_back(is.neighbor() ? static_cast<int>(mapper.index(is.outside())) : -1);
            const auto ig = is.geometry();
            a.hf_centroid.push_back(ig.center());
            a.hf_length.push_back(ig.volume());
            a.hf_normal.push_back(is.centerUnitOuterNormal());
            a.hf_corners.push_back({ig.corner(0), ig.corner(1)});
            if (face_halfface[f] < 0) {
                face_halfface[f] = hf;
                face_nodes[f] = {static_cast<int>(iset.subIndex(e, ref.subEntity(lf, 1, 0, 2), 2)),
                                 static_cast<int>(iset.subIndex(e, ref.subEntity(lf, 1, 1, 2), 2))};
            }
            face_cells[f].push_back(c);
        }
        a.cell_facepos.push_back(a.cell_faces.size());
        ++c;
    }

    for (int f = 0; f < nf; ++f) {
        if (face_cells[f].empty() || face_cells[f].size() > 2) {
            OPM_THROW(std::logic_error, "buildFractureMesh: a face needs one or two cells (no junctions)");
        }
        a.face_cells.insert(a.face_cells.end(), face_cells[f].begin(), face_cells[f].end());
        a.face_cellpos.push_back(a.face_cells.size());
    }
    a.face_halfface = std::move(face_halfface);
    a.face_nodes = std::move(face_nodes);
    return FractureMesh(std::move(a));
}

} // namespace Opm
