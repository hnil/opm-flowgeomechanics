#pragma once

#include <type_traits>

namespace Opm {

template <typename T>
inline double diagScalar(const T& d)
{
    if constexpr (std::is_arithmetic_v<std::decay_t<T>>) {
        return d;
    }
    else {
        return d[0][0];
    }
}

} // namespace Opm
