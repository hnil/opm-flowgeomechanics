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

#include <opm/geomech/FractureMesh.hpp>

#include <dune/foamgrid/foamgrid.hh>

namespace Opm
{

// Cells in element-mapper order, faces and nodes in leaf index-set order; geometry values are
// the grid's own (element and intersection geometries), so results built on either agree bitwise.
FractureMesh buildFractureMesh(const Dune::FoamGrid<2, 3>& grid);

} // namespace Opm
