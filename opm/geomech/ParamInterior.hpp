/*
  Copyright 2024 Equinor ASA.

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

#include <vector>
#include <array>

namespace Opm{

  typedef std::array<std::array<double, 3>, 3> Axis3D;

  // project a set of points on a 2D plane embedded in 3D down to that plane in
  // a 2D representation
  Axis3D project_to_2D(const std::vector<double>& p3d, std::vector<double>& p2d);

  // lift a set of points in 2D to a plane in 3D
  void lift_to_3D(const std::vector<double>& p2d,
                  const Axis3D& axis,
                  std::vector<double>& p3d);

  // redistribute a set of 2D points on a loop so that they are equally
  // distributed around the loop.  The first point is the "anchor point", which will
  // remain in place.
  void redistribute_2D(const std::vector<double>& points,
                       std::vector<double>& result);
  
  // parametrize a set of points in terms of the points on the boundary of a 2D
  // polygon
  void parametrize_interior_2D(const std::vector<double>& bpoints,
                               const std::vector<double>& ipoints,
                               std::vector<double>& result);
}; // end namespace Opm
