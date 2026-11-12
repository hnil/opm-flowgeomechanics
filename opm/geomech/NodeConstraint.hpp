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

#ifndef OPM_GEOMECH_NODE_CONSTRAINT_HPP
#define OPM_GEOMECH_NODE_CONSTRAINT_HPP

#include <opm/input/eclipse/Schedule/BCState.hpp>

#include <array>
#include <cstddef>
#include <tuple>
#include <vector>

namespace Opm::Elasticity {

//! Displacement constraints of one mechanics node: u . frame[i] = value[i]
//! for every i with fixed[i]. The rows of frame are orthonormal.
struct NodeConstraint
{
    std::size_t node = 0;
    std::array<std::array<double, 3>, 3> frame{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    std::array<bool, 3> fixed{};
    std::array<double, 3> value{};

    bool axisAligned() const
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (frame[i][j] != (i == j ? 1.0 : 0.0)) {
                    return false;
                }
            }
        }
        return true;
    }
};

//! BCCON masks constrain the coordinate axes.
inline std::vector<NodeConstraint>
toNodeConstraints(const std::vector<std::tuple<std::size_t, MechBCValue>>& bc_nodes)
{
    std::vector<NodeConstraint> result;
    result.reserve(bc_nodes.size());
    for (const auto& [node, bc] : bc_nodes) {
        NodeConstraint c;
        c.node = node;
        c.fixed = bc.fixeddir;
        c.value = bc.disp;
        result.push_back(c);
    }
    return result;
}

} // namespace Opm::Elasticity

#endif
