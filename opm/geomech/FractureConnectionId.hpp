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
