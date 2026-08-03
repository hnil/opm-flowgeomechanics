/*
  Copyright (C) 2026 SINTEF Digital

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
#ifndef OPM_FRACTURE_AUX_CELLS_IMPL_HPP
#define OPM_FRACTURE_AUX_CELLS_IMPL_HPP

#include <opm/geomech/FractureModel.hpp>

namespace Opm {

template <class TypeTag>
bool
FractureAuxCells<TypeTag>::bind(const FractureModel& fractures)
{
    const auto previousConnections = this->connections_.size();
    const auto previousActive = this->numActive();

    // The fracture's per-cell arrays are rebuilt at different points of its own solve,
    // and part way through a propagation they need not agree with each other.  Binding
    // against such a state attaches connections to the wrong degrees of freedom -- the
    // half-transmissibility list of an older, smaller grid indexed into the new one --
    // which is far worse than being one solve behind.  Keep the previous description
    // and say so.
    for (const auto& wellFractures : fractures.wellFractures()) {
        for (const auto& fracture : wellFractures) {
            const auto n = fracture.numCells();

            std::string bad;
            if (fracture.reservoirCells().size() != n) { bad = "reservoir cells"; }
            else if (fracture.leakOf().size() != n) { bad = "leak-off"; }
            else if (fracture.reservoirMobility().size() != n) { bad = "mobility"; }
            else if (static_cast<std::size_t>(fracture.fractureWidth().size()) != n) { bad = "width"; }

            if (!bad.empty()) {
                OpmLog::info(fmt::format("Embedded fracture flow: fracture state "
                                         "inconsistent ({} vs {} cells); keeping the "
                                         "previous binding", bad, n));
                return false;
            }
        }
    }

    this->connections_.clear();
    this->slotOf_.clear();
    this->wellPerforations_.clear();

    // Every slot starts this round dormant.  A fracture's grid can come back smaller
    // than it was -- the propagation rebuilds it -- and a slot the previous round
    // activated would otherwise keep its flag and its volume with no connections left,
    // an impossible cell that the convergence data and the volume refresh still count.
    std::fill(this->active_.begin(), this->active_.end(), false);
    std::fill(this->bulkVolume_.begin(), this->bulkVolume_.end(), 0.0);

    const auto numGridDof = this->simulator_.model().numGridDof();
    unsigned nextSlot = 0;

    // A slot is claimed once and never moves, so that a cell keeps its unknown from one
    // report step to the next.  Walking the fractures in the order the model holds them
    // is what makes that stable.
    std::size_t fractureIdx = 0;
    for (const auto& wellFractures : fractures.wellFractures()) {
        for (const auto& fracture : wellFractures) {
            const auto numCells = fracture.numCells();

            const auto& reservoirCells = fracture.reservoirCells();
            const auto& leakOf = fracture.leakOf();
            const auto& mobility = fracture.reservoirMobility();
            const auto& width = fracture.fractureWidth();
            const auto areas = fracture.cellAreas();
            const auto depths = fracture.cellDepths();

            if (nextSlot + numCells > this->capacity_) {
                OPM_THROW(std::runtime_error,
                          fmt::format("The fracture cells need more degrees of freedom than "
                                      "were reserved for them: {} in use, {} more wanted, {} "
                                      "reserved. Raise "
                                      "fractureparam.solver.embedded_capacity.",
                                      nextSlot, numCells, this->capacity_));
            }

            const auto firstSlot = nextSlot;
            for (std::size_t cell = 0; cell < numCells; ++cell) {
                const auto slot = firstSlot + cell;
                this->slotOf_.emplace_back(fractureIdx, cell);

                // The aperture is the volume; below the floor a cell is treated as having
                // the floor's aperture, exactly as the fracture's own pressure solve does,
                // so that a closed cell still has a well-defined -- if tiny -- volume.
                const auto aperture = std::max(static_cast<Scalar>(width[cell][0]), this->minWidth_);

                this->bulkVolume_[slot] = static_cast<Scalar>(areas[cell]) * aperture;
                this->depth_[slot] = static_cast<Scalar>(depths[cell]);

                const auto reservoirCell = reservoirCells[cell];
                if (reservoirCell < 0) {
                    // Outside this rank's grid, or not mapped: the cell exists as an
                    // unknown but exchanges nothing.
                    this->active_[slot] = false;
                    continue;
                }

                this->partner_[slot] = static_cast<unsigned>(reservoirCell);
                this->active_[slot] = true;

                // leakOf() carries the reservoir mobility, which the reservoir's own
                // local residual applies again from the upwind cell.  Divide it back out
                // so the connection is a transmissibility and nothing else.
                const auto mob = mobility[cell];
                const auto trans = (mob > 0.0)
                    ? static_cast<Scalar>(leakOf[cell] / mob)
                    : Scalar{0};

                this->connections_.push_back({static_cast<unsigned>(this->localToGlobalDof(slot)),
                                              static_cast<unsigned>(reservoirCell),
                                              trans, 0.0, 0.0});
            }

            // Fracture cell to fracture cell: the cubic law over the two half
            // transmissibilities, the same combination the fracture's own pressure solve
            // forms (FracturePressureAssemblerAD).  Computed fresh from the grid rather
            // than read from the fracture's cache: the cache is rebuilt at particular
            // points of the pressure solve and can lag a re-gridding with indices that
            // are still in range -- a wrong topology no size check catches, and a cell
            // it leaves disconnected has no storage and no path to the reservoir, which
            // is a singular row the moment the ILU eliminates around it.
            const auto freshHalfTrans = fracture.currentHalfTrans();
            for (const auto& [i, j, t1, t2] : freshHalfTrans) {
                const auto slotI = firstSlot + i;
                const auto slotJ = firstSlot + j;

                const auto h1 = std::max(static_cast<Scalar>(width[i][0]), this->minWidth_);
                const auto h2 = std::max(static_cast<Scalar>(width[j][0]), this->minWidth_);

                const auto invTrans = 12.0 / (h1 * h1 * h1 * t1)
                                    + 12.0 / (h2 * h2 * h2 * t2);

                this->connections_.push_back({static_cast<unsigned>(this->localToGlobalDof(slotI)),
                                              static_cast<unsigned>(this->localToGlobalDof(slotJ)),
                                              static_cast<Scalar>(1.0 / invTrans), 0.0, 0.0});
            }

            // Where the well meets the fracture.  perfinj_ is (fracture cell, well index)
            // and is what the fracture's own pressure solve uses to drive itself; the same
            // pair becomes a perforation of the degree of freedom that cell now owns.
            auto& perfs = this->wellPerforations_[fracture.wellInfo().name];
            for (const auto& [cell, wellIndex] : fracture.wellPerforations()) {
                const auto slot = firstSlot + static_cast<std::size_t>(cell);

                RuntimePerforation perf;
                perf.cell = static_cast<int>(this->localToGlobalDof(slot));
                perf.ctf = wellIndex;
                perf.depth = this->depth_[slot];
                perfs.push_back(perf);
            }

            nextSlot = firstSlot + numCells;
            ++fractureIdx;
        }
    }

    // A cell that has just been handed out has been holding a placeholder state; give it
    // the state of the rock it cuts through, at its own depth.
    auto& solution = this->simulator_.model().solution(/*timeIdx=*/0);
    for (unsigned slot = 0; slot < nextSlot; ++slot) {
        if (this->active_[slot]) {
            this->assignStateFromPartner(solution, slot);
        }
    }

    static_cast<void>(numGridDof);

    // Aggregate measures of what this binding feeds the flow, for the outer loop's
    // coupling residual.
    {
        Scalar totalTrans = 0.0;
        const auto gridDofLimit = this->simulator_.model().numGridDof();
        for (const auto& conn : this->connections_) {
            if (conn.dof1 < gridDofLimit || conn.dof2 < gridDofLimit) {
                totalTrans += conn.trans;
            }
        }

        Scalar totalPv = 0.0;
        for (unsigned slot = 0; slot < this->capacity_; ++slot) {
            totalPv += this->bulkVolume_[slot];
        }

        const auto rel = [](const Scalar now, const Scalar before) {
            if (before <= 0.0) {
                return (now > 0.0) ? Scalar{1} : Scalar{0};
            }
            return std::abs(now - before) / before;
        };

        this->lastBindChange_ = std::max(rel(totalTrans, this->lastTotalTrans_),
                                         rel(totalPv, this->lastTotalPv_));
        this->lastTotalTrans_ = totalTrans;
        this->lastTotalPv_ = totalPv;
    }

    const auto active = this->numActive();
    if (this->simulator_.gridView().comm().rank() == 0) {
        OpmLog::info(fmt::format("Embedded fracture flow: {} of {} reserved cells in use, "
                                 "{} connections",
                                 active, this->capacity_, this->connections_.size()));
    }

    // Only a changed connection list needs the sparsity pattern rebuilt; apertures moving
    // is a change of values.
    return (this->connections_.size() != previousConnections) || (active != previousActive);
}

template <class TypeTag>
bool
FractureAuxCells<TypeTag>::layoutMatches(const FractureModel& fractures) const
{
    // Cell counts per fracture, in binding order, reconstructed from the slot registry.
    std::vector<std::size_t> bound;
    for (const auto& [fractureIdx, cell] : this->slotOf_) {
        if (fractureIdx >= bound.size()) {
            bound.resize(fractureIdx + 1, 0);
        }
        bound[fractureIdx] = std::max(bound[fractureIdx], cell + 1);
    }

    std::size_t fractureIdx = 0;
    for (const auto& wellFractures : fractures.wellFractures()) {
        for (const auto& fracture : wellFractures) {
            if (fractureIdx >= bound.size() || fracture.numCells() != bound[fractureIdx]) {
                return false;
            }
            ++fractureIdx;
        }
    }

    return fractureIdx == bound.size();
}

template <class TypeTag>
std::vector<RuntimePerforation>
FractureAuxCells<TypeTag>::wellPerforations(const std::string& wellName) const
{
    auto pos = this->wellPerforations_.find(wellName);

    return (pos == this->wellPerforations_.end())
        ? std::vector<RuntimePerforation> {}
        : pos->second;
}

} // namespace Opm

#endif // OPM_FRACTURE_AUX_CELLS_IMPL_HPP
