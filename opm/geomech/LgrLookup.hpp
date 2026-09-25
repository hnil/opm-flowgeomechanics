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
#ifndef OPM_GEOMECH_LGR_LOOKUP_HPP
#define OPM_GEOMECH_LGR_LOOKUP_HPP

#include <opm/input/eclipse/EclipseState/Grid/Carfin.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>

#include <array>
#include <cstddef>
#include <string>

namespace Opm
{

/// LGR names and parent cells from the CARFIN records, which every rank
/// holds; the input grid that could answer the same exists on rank 0 only.
class LgrLookup
{
public:
    LgrLookup(const LgrCollection& lgrs, const std::array<int, 3>& levelZeroDims)
        : lgrs_(lgrs), dims_(levelZeroDims)
    {}

    /// Name of LGR grid number n (1-based, as connections carry it).
    std::string name(int n) const
    {
        return (n >= 1 && static_cast<std::size_t>(n) <= lgrs_.size())
            ? lgrs_.getLgr(static_cast<std::size_t>(n) - 1).NAME()
            : std::string{};
    }

    /// Level-zero Cartesian index of the cell that LGR n's cell lgrIndex was
    /// refined from, or -1 where that cannot be read off a uniform box.
    int father(int n, std::size_t lgrIndex) const
    {
        if (n < 1 || static_cast<std::size_t>(n) > lgrs_.size()) {
            return -1;
        }
        const auto& box = lgrs_.getLgr(static_cast<std::size_t>(n) - 1);
        if (box.isGraded() || box.PARENT_NAME() != "GLOBAL") {
            return -1;
        }
        const std::array<int, 3> lo {box.I1(), box.J1(), box.K1()};
        const std::array<int, 3> hi {box.I2(), box.J2(), box.K2()};
        const std::array<int, 3> fine {box.NX(), box.NY(), box.NZ()};
        std::array<int, 3> local {};
        auto rest = lgrIndex;
        for (int d = 0; d < 3; ++d) {
            local[d] = static_cast<int>(rest % fine[d]);
            rest /= fine[d];
        }
        std::array<int, 3> ijk {};
        for (int d = 0; d < 3; ++d) {
            const int perParent = fine[d] / (hi[d] + 1 - lo[d]);
            ijk[d] = lo[d] + local[d] / perParent;
        }
        return ijk[0] + dims_[0] * (ijk[1] + dims_[1] * ijk[2]);
    }

private:
    const LgrCollection& lgrs_;
    std::array<int, 3> dims_;
};

} // namespace Opm

#endif // OPM_GEOMECH_LGR_LOOKUP_HPP
