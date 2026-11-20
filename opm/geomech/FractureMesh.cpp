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

} // namespace Opm
