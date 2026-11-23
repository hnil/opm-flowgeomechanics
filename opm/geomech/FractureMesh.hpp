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

#pragma once

#include <dune/common/fvector.hh>

#include <array>
#include <span>
#include <vector>

namespace Opm
{

// Read-only fracture surface mesh in CSR arrays: cells, faces (edges) and nodes, all by index.
// A half-face is entry j of cellFaces(c); its geometry is as seen from c.
class FractureMesh
{
public:
    using Vec3 = Dune::FieldVector<double, 3>;
    using Triangle = std::array<Vec3, 3>;

    struct Arrays
    {
        std::vector<double> cell_area;
        std::vector<Vec3> cell_centroid;
        std::vector<Vec3> cell_normal;
        std::vector<int> cell_level;
        std::vector<int> cell_nodepos {0}; // CSR into cell_nodes and cell_corners
        std::vector<int> cell_nodes;
        std::vector<Vec3> cell_corners; // corners of the cell's geometry map, one per cell node
        std::vector<int> cell_tripos {0}; // CSR into triangles
        std::vector<Triangle> triangles;
        std::vector<int> cell_facepos {0}; // CSR into the half-face arrays
        std::vector<int> cell_faces;
        std::vector<int> cell_neighbours; // -1 on the boundary
        std::vector<Vec3> hf_centroid;
        std::vector<double> hf_length;
        std::vector<Vec3> hf_normal; // unit, in-plane, outward from the cell
        std::vector<std::array<Vec3, 2>> hf_corners;
        std::vector<int> face_cellpos {0}; // CSR into face_cells
        std::vector<int> face_cells;
        std::vector<int> face_halfface; // half-face of face_cells[face_cellpos[f]]
        std::vector<std::array<int, 2>> face_nodes;
        std::vector<Vec3> node_coord;
    };

    FractureMesh() = default;
    explicit FractureMesh(Arrays a);

    int numCells() const { return static_cast<int>(a_.cell_area.size()); }
    int numFaces() const { return static_cast<int>(a_.face_halfface.size()); }
    int numNodes() const { return static_cast<int>(a_.node_coord.size()); }

    double cellArea(int c) const { return a_.cell_area[c]; }
    const Vec3& cellCentroid(int c) const { return a_.cell_centroid[c]; }
    const Vec3& cellNormal(int c) const { return a_.cell_normal[c]; }
    int cellLevel(int c) const { return a_.cell_level[c]; }
    std::span<const Triangle> cellTriangles(int c) const { return csr(a_.triangles, a_.cell_tripos, c); }
    std::span<const int> cellNodes(int c) const { return csr(a_.cell_nodes, a_.cell_nodepos, c); }
    std::span<const Vec3> cellCorners(int c) const { return csr(a_.cell_corners, a_.cell_nodepos, c); }

    std::span<const int> cellFaces(int c) const { return csr(a_.cell_faces, a_.cell_facepos, c); }
    std::span<const int> cellNeighbours(int c) const { return csr(a_.cell_neighbours, a_.cell_facepos, c); }
    const Vec3& cellFaceCentroid(int c, int j) const { return a_.hf_centroid[a_.cell_facepos[c] + j]; }
    double cellFaceLength(int c, int j) const { return a_.hf_length[a_.cell_facepos[c] + j]; }
    const Vec3& cellFaceNormal(int c, int j) const { return a_.hf_normal[a_.cell_facepos[c] + j]; }
    const std::array<Vec3, 2>& cellFaceCorners(int c, int j) const
    {
        return a_.hf_corners[a_.cell_facepos[c] + j];
    }

    std::span<const int> faceCells(int f) const { return csr(a_.face_cells, a_.face_cellpos, f); }
    bool isBoundaryFace(int f) const { return faceCells(f).size() == 1; }
    const std::array<int, 2>& faceNodes(int f) const { return a_.face_nodes[f]; }
    const Vec3& faceCentroid(int f) const { return a_.hf_centroid[a_.face_halfface[f]]; }
    double faceLength(int f) const { return a_.hf_length[a_.face_halfface[f]]; }
    Vec3 faceNormal(int f, int c) const;

    const Vec3& nodeCoord(int n) const { return a_.node_coord[n]; }

private:
    template <class T>
    static std::span<const T> csr(const std::vector<T>& v, const std::vector<int>& pos, int i)
    {
        return {v.data() + pos[i], static_cast<std::size_t>(pos[i + 1] - pos[i])};
    }

    Arrays a_;
};

// Same numbering and floating-point operations as buildFractureMesh on a FoamGrid made from these
// vertices and triangles (insertion order), so the two agree bitwise.
FractureMesh buildFractureMesh(const std::vector<std::array<double, 3>>& nodes,
                               const std::vector<std::array<unsigned int, 3>>& triangles);

} // namespace Opm
