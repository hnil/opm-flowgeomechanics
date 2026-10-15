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
