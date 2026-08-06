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

#include <algorithm>
#include <cmath>
#include <set>

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

    // Which slots already held a cell.  A slot that did keeps the state it has solved
    // its way to; only a slot that has just been handed out needs one made up for it.
    const auto wasActive = this->active_;

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

            // The floor under the width in the cubic law is the fracture's own
            // (config.min_width), NOT the volume floor above.  The volume floor exists
            // so that a nearly closed cell still has storage -- a well-posedness matter
            // -- but the transmissibility goes with the aperture *cubed*, and flooring
            // it at the (much larger) volume floor short-circuits the fracture: every
            // cell is lifted to the well's pressure, the internal pressure gradient the
            // fracture's own solve computes disappears, and the leak-off is driven far
            // too hard.  The two representations must apply the same law to the same
            // apertures.
            const auto cubicFloor = static_cast<Scalar>(fracture.cubicLawMinWidth());

            for (const auto& [i, j, t1, t2] : freshHalfTrans) {
                const auto slotI = firstSlot + i;
                const auto slotJ = firstSlot + j;

                const auto h1 = std::max(static_cast<Scalar>(width[i][0]), cubicFloor);
                const auto h2 = std::max(static_cast<Scalar>(width[j][0]), cubicFloor);

                const auto invTrans = 12.0 / (h1 * h1 * h1 * t1)
                                    + 12.0 / (h2 * h2 * h2 * t2);

                this->connections_.push_back({static_cast<unsigned>(this->localToGlobalDof(slotI)),
                                              static_cast<unsigned>(this->localToGlobalDof(slotJ)),
                                              static_cast<Scalar>(1.0 / invTrans), 0.0, 0.0});
            }

            // Where the well meets the fracture.  perfinj_ is (fracture cell, well index)
            // and is what the fracture's own pressure solve uses to drive itself; the same
            // cells become perforations of the degrees of freedom they now own.  The well
            // index either carries over as-is -- the current modelling, a constant factor
            // on the cells around the wellbore -- or is estimated as radial flow through
            // the fracture's aperture, with a fixed prescribed width for now.
            auto& perfs = this->wellPerforations_[fracture.wellInfo().name];
            for (const auto& [cell, wellIndex] : fracture.wellPerforations()) {
                const auto slot = firstSlot + static_cast<std::size_t>(cell);

                RuntimePerforation perf;
                perf.cell = static_cast<int>(this->localToGlobalDof(slot));
                perf.depth = this->depth_[slot];

                if (this->perfWiMode_ == PerfWiMode::Estimate) {
                    // Radial flow in the fracture plane towards the wellbore:
                    // WI = 2 pi k_f w / ln(r_e/r_w), with the cubic-law
                    // permeability k_f = w^2/12 and r_e from the cell's area.
                    const auto w = this->perfWidth_;
                    const auto re = std::sqrt(static_cast<Scalar>(areas[cell]) / M_PI);
                    const auto lnTerm = std::log(std::max(re / this->perfRw_, Scalar{1.1}));

                    perf.ctf = 2.0 * M_PI * (w * w * w / 12.0) / lnTerm;
                }
                else {
                    perf.ctf = wellIndex;
                }

                perfs.push_back(perf);
            }

            nextSlot = firstSlot + numCells;
            ++fractureIdx;
        }
    }

    // A cell that has just been handed out has been holding a placeholder state; give it
    // the state of the rock it cuts through, at its own depth.  A slot that was already
    // carrying a cell keeps what it has: the binding is rebuilt at every step boundary,
    // and re-initialising the whole fracture there would overwrite the pressure it has
    // just solved for with the reservoir's, every step, so the fracture could never hold
    // a pressure of its own at all.
    auto& solution = this->simulator_.model().solution(/*timeIdx=*/0);
    for (unsigned slot = 0; slot < nextSlot; ++slot) {
        if (this->active_[slot] && !wasActive[slot]) {
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
void
FractureAuxCells<TypeTag>::leakoffReport(const FractureModel& fractures) const
{
    if constexpr (!Linearizer::assemblesAuxiliaryDofEquations) {
        return;
    }
    else {
        // What the well actually does at the fracture's own degrees of freedom.  This
        // reads only the well state, so it is reported whether or not the binding still
        // describes the fracture -- it is the answer to "is the fracture connected at
        // all", which must not depend on the fracture having stood still.
        this->perforationReport();

        if (!this->layoutMatches(fractures)) {
            OpmLog::info("LEAKOFF-CHECK skipped: binding layout differs from fracture state");
            return;
        }

        const auto& model = this->simulator_.model();
        const auto& problem = this->simulator_.problem();
        const auto& neighborInfo = model.linearizer().getNeighborInfo();
        const auto waterPos = FluidSystem::waterPhaseIdx;

        Scalar qEmb = 0.0;      // reservoir's water rate over the aux connections [sm3/s]
        Scalar qFrac = 0.0;     // fracture solver's own leak-off opinion [m3/s]
        Scalar condEmb = 0.0;   // sum of trans * upwind water mobility
        Scalar condFrac = 0.0;  // sum of the fracture's leakof_ (trans * total mobility)
        Scalar pFracSum = 0.0, pResSum = 0.0;
        Scalar dFracSum = 0.0, dResSum = 0.0, dZgSum = 0.0, dpotSum = 0.0;
        Scalar pMin = 1e30, pMax = -1e30, pPerfSum = 0.0;
        unsigned n = 0, nPerf = 0;

        // The same cell's pressure as the fracture's own solver holds it.  The two
        // representations solve for the same unknown; a difference here is the
        // coupling failing, not a difference of conductance.
        Scalar pOwnSum = 0.0, dpOwnMax = 0.0, pOwnAtMax = 0.0, pEmbAtMax = 0.0;
        std::size_t cellAtMax = 0;

        // Slots the well perforates, by global degree of freedom.
        std::vector<unsigned> perfDofs;
        for (const auto& [wname, perfs] : this->wellPerforations_) {
            static_cast<void>(wname);
            for (const auto& p : perfs) {
                perfDofs.push_back(static_cast<unsigned>(p.cell));
            }
        }

        std::size_t firstSlot = 0;
        for (const auto& wellFractures : fractures.wellFractures()) {
            for (const auto& fracture : wellFractures) {
                const auto numCells = fracture.numCells();
                const auto& leakOf = fracture.leakOf();
                const auto internalRate = fracture.leakOfRate();
                const auto areas = fracture.cellAreas();
                const auto& pOwn = fracture.fracturePressure();
                const auto& widths = fracture.fractureWidth();

                // Aperture statistics: what the cubic law actually works with, and how
                // many cells sit below the volume floor -- each of those is a cell whose
                // conductivity a floor at that level would inflate by (floor/width)^3.
                Scalar wMin = 1e30, wMax = 0.0, wSum = 0.0;
                unsigned nBelowVolumeFloor = 0;
                for (std::size_t cell = 0; cell < numCells && cell < widths.size(); ++cell) {
                    const auto w = static_cast<Scalar>(widths[cell][0]);
                    wMin = std::min(wMin, w);
                    wMax = std::max(wMax, w);
                    wSum += w;
                    if (w < this->minWidth_) {
                        ++nBelowVolumeFloor;
                    }
                }
                OpmLog::info(fmt::format(
                    "LEAKOFF-CHECK apertures: min {:.3g} mean {:.3g} max {:.3g} m  "
                    "cubic-law floor {:.3g} m  volume floor {:.3g} m  "
                    "cells below volume floor {} of {}",
                    wMin, (numCells > 0) ? wSum / numCells : Scalar{0}, wMax,
                    fracture.cubicLawMinWidth(), this->minWidth_,
                    nBelowVolumeFloor, numCells));

                for (std::size_t cell = 0; cell < numCells; ++cell) {
                    const auto slot = firstSlot + cell;
                    if (!this->active_[slot]) {
                        continue;
                    }

                    const auto g = static_cast<unsigned>(this->localToGlobalDof(slot));
                    const auto partner = this->partner_[slot];
                    const auto& iqF = model.intensiveQuantities(g, 0);
                    const auto& iqR = model.intensiveQuantities(partner, 0);

                    pFracSum += getValue(iqF.fluidState().pressure(waterPos));
                    pResSum += getValue(iqR.fluidState().pressure(waterPos));
                    dFracSum += problem.dofCenterDepth(g);
                    dResSum += problem.dofCenterDepth(partner);

                    const auto pF = getValue(iqF.fluidState().pressure(waterPos));
                    pMin = std::min(pMin, pF);
                    pMax = std::max(pMax, pF);
                    if (std::find(perfDofs.begin(), perfDofs.end(), g) != perfDofs.end()) {
                        pPerfSum += pF;
                        ++nPerf;
                    }

                    if (cell < pOwn.size()) {
                        const auto pO = static_cast<Scalar>(pOwn[cell][0]);
                        pOwnSum += pO;
                        if (std::abs(pO - pF) > dpOwnMax) {
                            dpOwnMax = std::abs(pO - pF);
                            pOwnAtMax = pO;
                            pEmbAtMax = pF;
                            cellAtMax = cell;
                        }
                    }

                    for (const auto& nbInfo : neighborInfo[g]) {
                        if (nbInfo.neighbor != partner) {
                            continue;
                        }

                        condEmb += nbInfo.res_nbinfo.trans
                            * getValue(iqF.mobility(waterPos));
                        dZgSum += nbInfo.res_nbinfo.dZg;
                        dpotSum += getValue(iqF.fluidState().pressure(waterPos))
                                 - getValue(iqR.fluidState().pressure(waterPos))
                                 - nbInfo.res_nbinfo.dZg
                                   * getValue(iqF.fluidState().density(waterPos));

                        RateVector flux(0.0);
                        RateVector darcy(0.0);
                        LocalResidual::computeFlux(flux, darcy, g, partner,
                                                   iqF, iqR, nbInfo.res_nbinfo,
                                                   problem.moduleParams());

                        const auto& fsys = iqF.fluidState().fluidSystem();
                        const auto waterEqIdx = Indices::conti0EqIdx
                            + fsys.canonicalToActiveCompIdx(
                                  fsys.solventComponentIndex(waterPos));

                        Scalar rate = getValue(flux[waterEqIdx])
                            * nbInfo.res_nbinfo.faceArea;
                        if constexpr (!getPropValue<TypeTag,
                                      Properties::BlackoilConserveSurfaceVolume>()) {
                            rate /= fsys.referenceDensity(waterPos,
                                                          problem.pvtRegionIndex(g));
                        }
                        qEmb += rate;
                        break;
                    }

                    if (cell < leakOf.size()) {
                        condFrac += static_cast<Scalar>(leakOf[cell]);
                    }
                    if (cell < internalRate.size()) {
                        qFrac += static_cast<Scalar>(internalRate[cell] * areas[cell]);
                    }
                    ++n;
                }

                firstSlot += numCells;
            }
        }

        const Scalar day = 86400.0;
        OpmLog::info(fmt::format(
            "LEAKOFF-CHECK cells {}  qEmb {:.6g} sm3/day  qFracInternal {:.6g} m3/day  "
            "condEmb {:.6g}  condFrac(leakof) {:.6g}  cond ratio emb/frac {:.4g}  "
            "mean pFrac {:.6g} bar  mean pRes {:.6g} bar  "
            "mean depthFrac {:.6g} m  mean depthRes {:.6g} m  "
            "mean dZg {:.6g}  mean dpot {:.6g} bar  "
            "pFrac min {:.6g} max {:.6g} bar  wellPerfs {} meanPatPerf {:.6g} bar",
            n, qEmb * day, qFrac * day, condEmb, condFrac,
            (condFrac > 0.0) ? condEmb / condFrac : Scalar{0},
            (n > 0) ? pFracSum / n / 1e5 : Scalar{0},
            (n > 0) ? pResSum / n / 1e5 : Scalar{0},
            (n > 0) ? dFracSum / n : Scalar{0},
            (n > 0) ? dResSum / n : Scalar{0},
            (n > 0) ? dZgSum / n : Scalar{0},
            (n > 0) ? dpotSum / n / 1e5 : Scalar{0},
            pMin / 1e5, pMax / 1e5, nPerf,
            (nPerf > 0) ? pPerfSum / nPerf / 1e5 : Scalar{0}));

        OpmLog::info(fmt::format(
            "LEAKOFF-CHECK pressure vs fracture solver: mean own {:.6g} bar  "
            "mean embedded {:.6g} bar  max |diff| {:.6g} bar at cell {} "
            "(own {:.6g}, embedded {:.6g})",
            (n > 0) ? pOwnSum / n / 1e5 : Scalar{0},
            (n > 0) ? pFracSum / n / 1e5 : Scalar{0},
            dpOwnMax / 1e5, cellAtMax, pOwnAtMax / 1e5, pEmbAtMax / 1e5));
    }
}

template <class TypeTag>
void
FractureAuxCells<TypeTag>::perforationReport() const
{
    if (this->wellPerforations_.empty()) {
        OpmLog::info("PERF-CHECK no fracture perforations registered");
        return;
    }

    const auto& model = this->simulator_.model();
    const auto& wellModel = this->simulator_.problem().wellModel();
    const auto& wellState = wellModel.wellState();
    const auto waterPos = FluidSystem::waterPhaseIdx;
    const auto waterActive = FluidSystem::canonicalToActivePhaseIdx(waterPos);
    const auto numPhases = wellState.numPhases();
    const Scalar day = 86400.0;

    for (const auto& [wname, perfs] : this->wellPerforations_) {
        if (perfs.empty() || !wellState.has(wname)) {
            continue;
        }

        const auto& ws = wellState[wname];
        const auto& pd = ws.perf_data;

        // The perforations of the fracture are the tail of the well's list, but match
        // on the degree of freedom rather than assuming that: the whole point is to
        // find out whether they are in the list at all.
        std::set<std::size_t> auxDofs;
        for (const auto& p : perfs) {
            auxDofs.insert(static_cast<std::size_t>(p.cell));
        }

        Scalar qAux = 0.0, ctfAux = 0.0, pPerfSum = 0.0;
        Scalar qTotal = 0.0;
        std::size_t nFound = 0;

        for (std::size_t perf = 0; perf < pd.size(); ++perf) {
            const auto rate = (waterActive < numPhases)
                ? pd.phase_rates[perf * numPhases + waterActive] : Scalar{0};
            qTotal += rate;

            if (auxDofs.count(pd.cell_index[perf]) == 0) {
                continue;
            }

            ++nFound;
            qAux += rate;
            ctfAux += pd.connection_transmissibility_factor[perf];
            pPerfSum += pd.pressure[perf];
        }

        // What the flow problem holds at those degrees of freedom, for the drawdown.
        Scalar pCellSum = 0.0;
        for (const auto& p : perfs) {
            const auto& iq = model.intensiveQuantities(static_cast<unsigned>(p.cell), 0);
            pCellSum += getValue(iq.fluidState().pressure(waterPos));
        }

        OpmLog::info(fmt::format(
            "PERF-CHECK {}: registered {}  found in well {} of {} perforations  "
            "bhp {:.6g} bar  mean perf pressure {:.6g} bar  mean cell pressure {:.6g} bar  "
            "sum CTF {:.6g}  water rate through fracture {:.6g} of {:.6g} sm3/day",
            wname, perfs.size(), nFound, pd.size(),
            ws.bhp / 1e5,
            (nFound > 0) ? pPerfSum / nFound / 1e5 : Scalar{0},
            pCellSum / perfs.size() / 1e5,
            ctfAux, qAux * day, qTotal * day));
    }
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
