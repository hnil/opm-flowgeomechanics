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

#include <opm/grid/CpGrid.hpp>

#include <dune/common/fvector.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/istl/bvector.hh>
#if HAVE_MPI
#include <dune/common/parallel/communicator.hh>
#include <dune/istl/owneroverlapcopy.hh>
#endif

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#if HAVE_MPI
#include <opm/geomech/DuneUtilities.hpp>
#include <opm/geomech/DuneCommunicationHelpers.hpp>
#endif
#include <opm/geomech/MechanicalVertexStencilValidation.hpp>

namespace Opm::Elasticity {
#if HAVE_MPI
namespace detail {

//! Number of global ids without exactly one owner or with different positions across ranks.
template <class Communication>
int countBadPartitionLabels(const Communication& cc,
                            const std::vector<int>& ids,
                            const std::vector<int>& owned,
                            const std::vector<double>& xyz)
{
    const int n = static_cast<int>(ids.size());
    std::vector<int> counts(cc.size()), displ(cc.size() + 1, 0);
    cc.gather(&n, counts.data(), 1, 0);
    std::partial_sum(counts.begin(), counts.end(), displ.begin() + 1);
    std::vector<int> xyzCounts(counts), xyzDispl(displ);
    for (std::size_t r = 0; r < counts.size(); ++r) {
        xyzCounts[r] *= 3;
        xyzDispl[r] *= 3;
    }
    const int total = displ.back();
    std::vector<int> allIds(total), allOwned(total);
    std::vector<double> allXyz(3 * total);
    cc.gatherv(ids.data(), n, allIds.data(), counts.data(), displ.data(), 0);
    cc.gatherv(owned.data(), n, allOwned.data(), counts.data(), displ.data(), 0);
    cc.gatherv(xyz.data(), 3 * n, allXyz.data(), xyzCounts.data(), xyzDispl.data(), 0);

    int bad = 0;
    if (cc.rank() == 0) {
        std::vector<int> order(total);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&allIds](int a, int b) { return allIds[a] < allIds[b]; });
        constexpr double positionTolerance = 1e-6; // metres
        for (int begin = 0; begin < total;) {
            int end = begin;
            int owners = 0;
            bool moved = false;
            for (; end < total && allIds[order[end]] == allIds[order[begin]]; ++end) {
                owners += allOwned[order[end]];
                for (int d = 0; d < 3; ++d) {
                    moved = moved
                        || std::abs(allXyz[3 * order[end] + d] - allXyz[3 * order[begin] + d])
                            > positionTolerance;
                }
            }
            bad += (owners != 1 || moved) ? 1 : 0;
            begin = end;
        }
    }
    cc.broadcast(&bad, 1, 0);
    return bad;
}

class CellIdDataHandle : public Dune::CommDataHandleIF<CellIdDataHandle, int>
{
public:
    CellIdDataHandle(std::vector<int>& ids, const Dune::CpGrid::LeafGridView& gv)
        : ids_(ids)
        , gv_(gv)
    {
    }
    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return true; }
    template <class Entity>
    std::size_t size(const Entity&) const { return 1; }
    template <class Buffer, class Entity>
    void gather(Buffer& buffer, const Entity& e) const { buffer.write(ids_[gv_.indexSet().index(e)]); }
    template <class Buffer, class Entity>
    void scatter(Buffer& buffer, const Entity& e, std::size_t) { buffer.read(ids_[gv_.indexSet().index(e)]); }

private:
    std::vector<int>& ids_;
    Dune::CpGrid::LeafGridView gv_;
};

} // namespace detail

//! Throws on every rank unless cell and node labels and the owner-to-copy
//! exchange of both are correct, and every owned node has its full cell star.
inline void validateMechanicalPartition(const Dune::CpGrid& grid,
                                        Dune::OwnerOverlapCopyCommunication<int, int>& comm)
{
    const auto& cc = grid.comm();
    if (cc.size() == 1) {
        return;
    }
    const auto gv = grid.leafGridView();
    const auto& gids = grid.globalIdSet();

    std::vector<int> cellIds(gv.size(0)), cellOwned(gv.size(0));
    std::vector<double> cellXyz(3 * gv.size(0));
    std::vector<int> received(gv.size(0), -1);
    for (const auto& cell : elements(gv)) {
        const int c = gv.indexSet().index(cell);
        const auto x = cell.geometry().center();
        cellIds[c] = gids.id(cell);
        cellOwned[c] = cell.partitionType() == Dune::InteriorEntity;
        std::copy(x.begin(), x.end(), cellXyz.begin() + 3 * c);
        if (cellOwned[c]) {
            received[c] = cellIds[c];
        }
    }
    detail::CellIdDataHandle handle(received, gv);
    gv.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    const int badCellExchange = cc.sum(static_cast<int>(received != cellIds));

    std::vector<int> nodeIds, nodeOwned;
    std::vector<double> nodeXyz;
    std::vector<Dune::FieldVector<double, 3>> position(gv.size(3));
    for (const auto& vertex : vertices(gv)) {
        position[gv.indexSet().index(vertex)] = vertex.geometry().center();
    }
    Dune::BlockVector<Dune::FieldVector<double, 1>> exchanged(3 * gv.size(3));
    exchanged = -1.0;
    using Attribute = Dune::OwnerOverlapCopyAttributeSet;
    for (const auto& index : comm.indexSet()) {
        const int i = index.local().local();
        const bool owner = index.local().attribute() == Attribute::owner;
        if (owner) {
            exchanged[i] = index.global();
        }
        if (i % 3 == 0) {
            nodeIds.push_back(index.global());
            nodeOwned.push_back(owner);
            const auto& x = position[i / 3];
            nodeXyz.insert(nodeXyz.end(), x.begin(), x.end());
        }
    }
    comm.copyOwnerToAll(exchanged, exchanged);
    int badNodeExchange = static_cast<int>(comm.indexSet().size()) != 3 * gv.size(3);
    for (const auto& index : comm.indexSet()) {
        badNodeExchange += exchanged[index.local().local()] != static_cast<double>(index.global());
    }
    badNodeExchange = cc.sum(badNodeExchange);

    const int badCells = detail::countBadPartitionLabels(cc, cellIds, cellOwned, cellXyz);
    const int badNodes = detail::countBadPartitionLabels(cc, nodeIds, nodeOwned, nodeXyz);
    if (badCells + badNodes + badCellExchange + badNodeExchange != 0) {
        throw std::runtime_error("Mechanics partition is inconsistent: " + std::to_string(badCells)
                                 + " cell and " + std::to_string(badNodes)
                                 + " node labels wrong, cell exchange failed on "
                                 + std::to_string(badCellExchange) + " ranks, node exchange at "
                                 + std::to_string(badNodeExchange) + " dofs.");
    }
    validateMechanicalVertexStencils(grid, comm);
}

inline void validateMechanicalPartition(const Dune::CpGrid& grid)
{
    if (grid.comm().size() == 1) {
        return;
    }
    Dune::OwnerOverlapCopyCommunication<int, int> comm(
        grid.comm(), Dune::SolverCategory::overlapping, false);
    comm.indexSet() = Opm::entityToDofIndexSet(Opm::makeEntityEntityCommunication<3>(grid), 3);
    comm.remoteIndices().rebuild<false>();
    validateMechanicalPartition(grid, comm);
}
#endif
} // namespace Opm::Elasticity
