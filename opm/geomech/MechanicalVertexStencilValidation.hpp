#pragma once

#include <opm/grid/CpGrid.hpp>
#include <dune/grid/common/rangegenerators.hh>
#if HAVE_MPI
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/communicator.hh>
#include <opm/geomech/DuneUtilities.hpp>
#include <opm/geomech/DuneCommunicationHelpers.hpp>
#endif
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Elasticity {
#if HAVE_MPI
// The VEM operator is assembled from all local cells. Every owned vertex must
// therefore have its complete incident-cell star, including corner neighbors.
// Count each global cell only on its cell owner and sum over vertex copies;
// compare that global count with the actual cells available to the vertex owner.
// Use face vertices, not only the cell's eight corners: conforming refinements
// may insert additional vertices along faces and edges.
inline void validateMechanicalVertexStencils(
    const Dune::CpGrid& grid,
    Dune::OwnerOverlapCopyCommunication<int, int>& comm)
{
    if (grid.comm().size() == 1) {
        return;
    }
    const auto gv = grid.leafGridView();
    std::vector<int> local(3*gv.size(3), 0), owned(local.size(), 0);
    for (const auto& cell : elements(gv)) {
        const int c = gv.indexSet().index(cell);
        std::vector<int> nodes;
        for (int f = 0; f < grid.numCellFaces(c); ++f) {
            const int face = grid.cellFace(c, f);
            if (grid.numFaceVertices(face) < 3) {
                continue; // NNCs are flow connections, not VEM faces.
            }
            for (int v = 0; v < grid.numFaceVertices(face); ++v) {
                nodes.push_back(grid.faceVertex(face, v));
            }
        }
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        for (const int node : nodes) {
            for (int d = 0; d < 3; ++d) {
                ++local[3*node+d];
                if (cell.partitionType() == Dune::InteriorEntity) {
                    ++owned[3*node+d];
                }
            }
        }
    }
    using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
    Communication::AllSet all;
    Dune::Interface interface;
    interface.build(comm.remoteIndices(), all, all);
    Dune::BufferedCommunicator exchange;
    exchange.build<std::vector<int>>(interface);
    auto global = owned;
    exchange.forward<Dune::FreeAddGatherScatter<std::vector<int>>>(owned, global);
    exchange.free();
    int bad = 0;
    for (const auto& index : comm.indexSet()) {
        const int i = index.local().local();
        if (i % 3 == 0 && index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner
            && local[i] != global[i]) {
            ++bad;
        }
    }
    const int total = grid.comm().sum(bad);
    if (total != 0) {
        throw std::runtime_error("Incomplete VEM vertex stencils at " + std::to_string(total)
            + " owned vertices: the partition overlap omits incident cells. "
              "Use --add-corners=true --add-corners-by-vertex=true.");
    }
}
// Preflight before timestepping: an incomplete overlap cannot be repaired by
// cutting the timestep, so report it outside the adaptive retry loop as well.
inline void validateMechanicalVertexStencils(const Dune::CpGrid& grid)
{
    if (grid.comm().size() == 1) {
        return;
    }
    Dune::OwnerOverlapCopyCommunication<int, int> comm(
        grid.comm(), Dune::SolverCategory::overlapping, false);
    comm.indexSet() = Opm::entityToDofIndexSet(Opm::makeEntityEntityCommunication<3>(grid), 3);
    comm.remoteIndices().rebuild<false>();
    validateMechanicalVertexStencils(grid, comm);
}
#endif
} // namespace Opm::Elasticity
