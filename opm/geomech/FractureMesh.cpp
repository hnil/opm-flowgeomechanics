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

#include <opm/geomech/FractureMesh.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/geomech/Math.hpp>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>

namespace Opm
{

FractureMesh::FractureMesh(Arrays a)
    : a_(std::move(a))
{
    const std::size_t nc = a_.cell_area.size();
    const std::size_t nf = a_.face_halfface.size();
    const std::size_t nhf = a_.cell_faces.size();
    const bool ok = a_.cell_centroid.size() == nc && a_.cell_normal.size() == nc && a_.cell_level.size() == nc
        && a_.cell_nodepos.size() == nc + 1 && a_.cell_tripos.size() == nc + 1
        && a_.cell_facepos.size() == nc + 1 && a_.cell_nodes.size() == std::size_t(a_.cell_nodepos.back())
        && a_.cell_corners.size() == a_.cell_nodes.size() && a_.triangles.size() == std::size_t(a_.cell_tripos.back())
        && nhf == std::size_t(a_.cell_facepos.back()) && a_.cell_neighbours.size() == nhf
        && a_.hf_centroid.size() == nhf && a_.hf_length.size() == nhf && a_.hf_normal.size() == nhf
        && a_.hf_corners.size() == nhf && a_.face_cellpos.size() == nf + 1
        && a_.face_cells.size() == std::size_t(a_.face_cellpos.back()) && a_.face_nodes.size() == nf;
    if (!ok) {
        OPM_THROW(std::logic_error, "FractureMesh: inconsistent array sizes");
    }
}

FractureMesh::Vec3
FractureMesh::faceNormal(int f, int c) const
{
    const auto faces = cellFaces(c);
    for (std::size_t j = 0; j < faces.size(); ++j) {
        if (faces[j] == f) {
            return cellFaceNormal(c, static_cast<int>(j));
        }
    }
    OPM_THROW(std::logic_error, "FractureMesh::faceNormal: face is not a face of the cell");
}

namespace
{
// FoamGrid's unit outer normal of edge v0-v1 of triangle p.
FractureMesh::Vec3
outerNormal(const FractureMesh::Triangle& p, int v0, int v1)
{
    using Vec3 = FractureMesh::Vec3;
    int v2 = (v1 + 1) % 3;
    if (v2 == v0)
        v2 = (v0 + 1) % 3;
    Vec3 facet = p[v0] - p[v1];
    Vec3 otherEdge = p[v2] - p[v1];
    Vec3 scaledElementNormal;
    Vec3 n;
    scaledElementNormal[0] = facet[1]*otherEdge[2] - facet[2]*otherEdge[1];
    scaledElementNormal[1] = facet[2]*otherEdge[0] - facet[0]*otherEdge[2];
    scaledElementNormal[2] = facet[0]*otherEdge[1] - facet[1]*otherEdge[0];
    n[0] = facet[1]*scaledElementNormal[2] - facet[2]*scaledElementNormal[1];
    n[1] = facet[2]*scaledElementNormal[0] - facet[0]*scaledElementNormal[2];
    n[2] = facet[0]*scaledElementNormal[1] - facet[1]*scaledElementNormal[0];
    otherEdge = p[v0] - p[v2];
    if (otherEdge * n < 0)
        n *= -1.0;
    n /= n.two_norm();
    return n;
}
} // namespace

FractureMesh
buildFractureMesh(const std::vector<std::array<double, 3>>& nodes,
                  const std::vector<std::array<unsigned int, 3>>& triangles)
{
    using Vec3 = FractureMesh::Vec3;
    const auto type = Dune::GeometryTypes::simplex(2);
    const auto& ref = Dune::ReferenceElements<double, 2>::general(type);
    const int nc = static_cast<int>(triangles.size());

    FractureMesh::Arrays a;
    a.node_coord.reserve(nodes.size());
    for (const auto& x : nodes)
        a.node_coord.push_back(Vec3 {x[0], x[1], x[2]});

    // faces numbered in order of first appearance over cells and local faces, as FoamGrid does
    std::map<std::pair<unsigned int, unsigned int>, int> face_of;
    std::vector<std::vector<int>> face_cells;
    std::vector<int> face_halfface;
    std::vector<std::array<int, 2>> face_nodes;
    for (int c = 0; c < nc; ++c) {
        for (int lf = 0; lf < 3; ++lf) {
            unsigned int n0 = triangles[c][ref.subEntity(lf, 1, 0, 2)];
            unsigned int n1 = triangles[c][ref.subEntity(lf, 1, 1, 2)];
            const auto key = std::minmax(n0, n1);
            auto [it, inserted] = face_of.try_emplace({key.first, key.second}, static_cast<int>(face_cells.size()));
            const int f = it->second;
            if (inserted) {
                face_cells.emplace_back();
                face_halfface.push_back(3 * c + lf);
                face_nodes.push_back({static_cast<int>(n0), static_cast<int>(n1)});
            }
            a.cell_faces.push_back(f);
            face_cells[f].push_back(c);
        }
    }

    a.cell_area.reserve(nc);
    a.cell_centroid.reserve(nc);
    a.cell_normal.reserve(nc);
    a.cell_level.assign(nc, 0);
    for (int c = 0; c < nc; ++c) {
        FractureMesh::Triangle tri;
        for (int i = 0; i < 3; ++i) {
            const int n = static_cast<int>(triangles[c][i]);
            a.cell_nodes.push_back(n);
            tri[i] = a.node_coord[n];
        }
        const Dune::AffineGeometry<double, 2, 3> geom(type, tri);
        a.cell_area.push_back(geom.volume());
        a.cell_centroid.push_back(geom.center());
        for (int i = 0; i < 3; ++i)
            a.cell_corners.push_back(geom.corner(i));
        a.cell_nodepos.push_back(a.cell_nodes.size());
        a.triangles.push_back(tri);
        a.cell_tripos.push_back(a.triangles.size());

        // same operations as ddm::normalOfElement
        const Vec3 e01 = tri[1] - tri[0];
        const Vec3 e02 = tri[2] - tri[0];
        Vec3 normal = crossProduct(e01, e02);
        normal /= normal.two_norm();
        a.cell_normal.push_back(normal);

        for (int lf = 0; lf < 3; ++lf) {
            const int f = a.cell_faces[3 * c + lf];
            const auto& fc = face_cells[f];
            if (fc.size() > 2) {
                OPM_THROW(std::logic_error, "buildFractureMesh: a face needs one or two cells (no junctions)");
            }
            a.cell_neighbours.push_back(fc.size() == 2 ? (fc[0] == c ? fc[1] : fc[0]) : -1);
            const int v0 = ref.subEntity(lf, 1, 0, 2);
            const int v1 = ref.subEntity(lf, 1, 1, 2);
            const std::array<Vec3, 2> ends {tri[v0], tri[v1]};
            const Dune::AffineGeometry<double, 1, 3> ig(Dune::GeometryTypes::simplex(1), ends);
            a.hf_centroid.push_back(ig.center());
            a.hf_length.push_back(ig.volume());
            a.hf_normal.push_back(outerNormal(tri, v0, v1));
            a.hf_corners.push_back({ig.corner(0), ig.corner(1)});
        }
        a.cell_facepos.push_back(3 * (c + 1));
    }

    for (const auto& fc : face_cells) {
        a.face_cells.insert(a.face_cells.end(), fc.begin(), fc.end());
        a.face_cellpos.push_back(a.face_cells.size());
    }
    a.face_halfface = std::move(face_halfface);
    a.face_nodes = std::move(face_nodes);
    return FractureMesh(std::move(a));
}

} // namespace Opm
