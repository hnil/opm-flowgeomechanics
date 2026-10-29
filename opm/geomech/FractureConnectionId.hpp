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

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <algorithm>
#include <cstddef>

namespace Opm {

// Completion indices are unique only within their grid. WSEED parent-cell
// indices are separate identities and must not be used for this lookup.
struct FractureConnectionId
{
    std::size_t global_index{};
    int lgr_level{};

    bool matches(const Connection& connection) const
    {
        return global_index == connection.global_index()
            && lgr_level == connection.get_lgr_level();
    }

    bool operator==(const FractureConnectionId& other) const
    {
        return global_index == other.global_index && lgr_level == other.lgr_level;
    }

    const Connection* find(const WellConnections& connections) const
    {
        const auto it = std::find_if(connections.begin(), connections.end(),
                                    [this](const auto& connection) { return matches(connection); });
        return it == connections.end() ? nullptr : &*it;
    }
};

} // namespace Opm
