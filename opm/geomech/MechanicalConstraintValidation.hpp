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

#include <opm/input/eclipse/Schedule/BCState.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace Dune { class CpGrid; }

namespace Opm::Elasticity {

struct ConstrainedMechanicalVertex
{
    std::int64_t id;
    std::array<double, 3> position;
    unsigned mask;
};

struct MechanicalConstraintComponent
{
    std::int64_t id;
    int rank;
};

// Edges describe mechanical connectivity through shared vertex DOFs. Duplicate
// vertices/constraints from overlap partitions are merged by global ID.
std::vector<MechanicalConstraintComponent> rigidBodyConstraintRanks(
    const std::vector<std::pair<std::int64_t, std::int64_t>>& edges,
    const std::vector<ConstrainedMechanicalVertex>& constraints);

// Collective on grid.comm(): every rank receives the same component ranks.
std::vector<MechanicalConstraintComponent> mechanicalConstraintRanks(
    const Dune::CpGrid& grid,
    const std::vector<std::tuple<std::size_t, MechBCValue>>& boundaryNodes);

} // namespace Opm::Elasticity
