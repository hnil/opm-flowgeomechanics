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

#include <config.h>

#include <opm/geomech/MechanicalPartitionValidation.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;

// Rank 1 owns one octant, so its inner corner node touches a rank 0 cell only at that node.
void balance(Dune::CpGrid& grid, bool corners, bool byVertex)
{
    grid.createCartesian({4, 4, 4}, {1, 1, 1});
    std::vector<int> parts(64);
    for (int k = 0; k < 4; ++k)
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 4; ++i)
                parts[i + 4 * j + 16 * k] = (i >= 2 && j >= 2 && k >= 2) ? 1 : 0;
    grid.setCornerCellsByVertex(byVertex);
    grid.loadBalance(parts, false, corners, 1);
}

void buildComm(const Dune::CpGrid& grid, Comm& comm)
{
    comm.indexSet() = Opm::entityToDofIndexSet(Opm::makeEntityEntityCommunication<3>(grid), 3);
}

template <class Check>
void expectAbort(const char* what, Check&& check)
{
    bool threw = false;
    try {
        check();
    } catch (const std::runtime_error&) {
        threw = true;
    }
    if (!threw) {
        throw std::logic_error(std::string("no abort on ") + what);
    }
}

// Relabels every shared owner dof as copy, which leaves shared nodes without an owner.
void dropOwners(const Dune::CpGrid& grid, Comm& comm)
{
    using Index = Comm::ParallelIndexSet;
    const auto good = Opm::entityToDofIndexSet(Opm::makeEntityEntityCommunication<3>(grid), 3);
    auto& set = comm.indexSet();
    set.beginResize();
    for (const auto& index : good) {
        const auto owner = Dune::OwnerOverlapCopyAttributeSet::owner;
        const auto copy = Dune::OwnerOverlapCopyAttributeSet::copy;
        const bool flip = grid.comm().rank() == 1 && index.global() % 3 == 0
            && index.local().attribute() == owner;
        set.add(index.global(),
                Index::LocalIndex(index.local().local(),
                                  flip ? copy : index.local().attribute(), true));
    }
    set.endResize();
    comm.remoteIndices().rebuild<false>();
}
} // namespace

int main(int argc, char** argv)
{
    const auto& helper = Dune::MPIHelper::instance(argc, argv);
    if (helper.size() != 2) {
        std::cerr << "needs two ranks\n";
        return 1;
    }
    try {
        Dune::CpGrid good;
        balance(good, true, true);
        Opm::Elasticity::validateMechanicalPartition(good);

        // Both overlaps below are short of the cell star at shared nodes.
        Dune::CpGrid faceOnly;
        balance(faceOnly, false, false);
        expectAbort("face-only overlap",
                    [&] { Opm::Elasticity::validateMechanicalPartition(faceOnly); });
        Dune::CpGrid cornersOnly;
        balance(cornersOnly, true, false);
        expectAbort("corners-only overlap",
                    [&] { Opm::Elasticity::validateMechanicalPartition(cornersOnly); });

        Comm wrongOwners(good.comm(), Dune::SolverCategory::overlapping, false);
        buildComm(good, wrongOwners);
        dropOwners(good, wrongOwners);
        expectAbort("nodes without owner",
                    [&] { Opm::Elasticity::validateMechanicalPartition(good, wrongOwners); });

        // Cell labels: an id owned twice, an id not owned, the same id at two positions.
        const auto& cc = good.comm();
        const std::vector<double> origin(3, 0.0);
        const std::vector<double> rankPosition(3, 1.0 * cc.rank());
        if (Opm::Elasticity::detail::countBadPartitionLabels(cc, {7}, {1}, origin) != 1
            || Opm::Elasticity::detail::countBadPartitionLabels(cc, {7}, {0}, origin) != 1
            || Opm::Elasticity::detail::countBadPartitionLabels(
                   cc, {7}, {cc.rank() == 0 ? 1 : 0}, rankPosition) != 1
            || Opm::Elasticity::detail::countBadPartitionLabels(
                   cc, {7}, {cc.rank() == 0 ? 1 : 0}, origin) != 0) {
            throw std::logic_error("label count is wrong");
        }
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
