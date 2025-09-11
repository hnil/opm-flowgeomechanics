/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_CUTDE_HPP_INCLUDED
#define OPM_CUTDE_HPP_INCLUDED

#include <array>

namespace ddm
{
using Real = double;

struct Real2
{
    Real x {};
    Real y {};
};

struct Real3
{
    Real x {};
    Real y {};
    Real z {};
};

struct Real6
{
    Real x {};
    Real y {};
    Real z {};
    Real a {};
    Real b {};
    Real c {};
};

Real3 make3(Real x, Real y, Real z);
Real6 make6(Real x, Real y, Real z, Real a, Real b, Real c);
Real3 disp_fs(const Real3& obs, const std::array<Real3, 3>& tri, const Real3& slip, Real nu);
Real6 strain_fs(const Real3& obs, const std::array<Real3, 3>& tri, const Real3& slip, Real nu);

} // namespace ddm

#endif // OPM_CUTDE_HPP_INCLUDED
