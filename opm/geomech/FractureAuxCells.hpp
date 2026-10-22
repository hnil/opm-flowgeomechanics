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
#ifndef OPM_FRACTURE_AUX_CELLS_HPP
#define OPM_FRACTURE_AUX_CELLS_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/geomech/FractureMechCoupling.hpp>
#include <opm/simulators/flow/FlowAuxCellModule.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/models/nonlinear/newtonmethodproperties.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class FractureModel;

/*!
 * \brief Fracture cells as degrees of freedom of the reservoir flow problem.
 *
 * The alternative to upscaling the fracture to well indices (\f$WI = q/\Delta p\f$,
 * sign-indefinite): fracture cells carry the reservoir's own conservation equations and
 * connect through ordinary transmissibilities.
 *
 * \section capacity Preallocation
 *
 * A fixed number of DOFs is claimed up front and handed out as cells open. An unused
 * (dormant) cell has no volume or connections and gets an identity row. Running out of
 * capacity stops the run rather than silently capping the fracture.
 */
template <class TypeTag>
class FractureAuxCells : public FlowAuxCellModule<TypeTag>
{
    using ParentType = FlowAuxCellModule<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { dimWorld = GridView::dimensionworld };

public:
    using Connection = typename ParentType::Connection;

    //! How the well index of a well-to-fracture-cell perforation is formed.
    enum class PerfWiMode {
        //! The factor the fracture's own pressure solve uses (perfinj_): the
        //! current modelling -- a constant factor on the cells around the
        //! wellbore.
        Fracture,
        //! A radial-flow estimate through the fracture's aperture,
        //! 2 pi (w^3/12) / ln(r_e/r_w), with a fixed prescribed width.
        Estimate,
    };

    FractureAuxCells(Simulator& simulator,
                     const unsigned capacity,
                     const Scalar minWidth,
                     const PerfWiMode perfWiMode,
                     const Scalar perfWidth,
                     const Scalar perfRw)
        : simulator_(simulator)
        , capacity_(capacity)
        , minWidth_(minWidth)
        , perfWiMode_(perfWiMode)
        , perfWidth_(perfWidth)
        , perfRw_(perfRw)
        , active_(capacity, false)
        , bulkVolume_(capacity, 0.0)
        , depth_(capacity, 0.0)
        , partner_(capacity, 0)
    {}

    unsigned numDofs() const override
    { return this->capacity_; }

    /*!
     * \brief The aperture is the pore volume.
     *
     * A fracture cell is a void, not rock with pores in it, so its pore volume is its
     * whole volume and its porosity is one.  This is also what gives it no rock heat
     * capacity, since the model derives the rock fraction from the porosity.
     */
    Scalar poreVolume(unsigned localIdx) const override
    { return this->bulkVolume_.at(localIdx); }

    Scalar bulkVolume(unsigned localIdx) const override
    { return this->bulkVolume_.at(localIdx); }

    Scalar depth(unsigned localIdx) const override
    { return this->depth_.at(localIdx); }

    //! The fracture holds the same fluid as the rock it cuts through.
    unsigned pvtRegionIndex(unsigned localIdx) const override
    { return this->simulator_.problem().pvtRegionIndex(this->partner_.at(localIdx)); }

    unsigned satRegionIndex(unsigned localIdx) const override
    { return this->simulator_.problem().satnumRegionIndex(this->partner_.at(localIdx)); }

    unsigned initialisationPartner(unsigned localIdx) const override
    { return this->partner_.at(localIdx); }

    bool isActive(unsigned localIdx) const override
    { return this->active_.at(localIdx); }

    /*!
     * \brief Fracture cells stay out of the CNV measure.
     *
     * Their pore volume is an aperture times an area -- minute against any flow worth
     * simulating -- so the volume-scaled residual dwarfs every tolerance while the mass
     * it stands for is negligible.  The material balance still covers them, weighed by
     * that same small mass.
     */
    bool participatesInCnv() const override
    { return false; }

    //! Residual check against the host rock cell's pore volume (opt-in): same units
    //! and tolerance as CNV, and no through-flow in the denominator.
    double cnvReferencePoreVolume(unsigned localIdx) const override
    {
        if (!this->residualCheck_ || localIdx >= this->partner_.size() || !this->active_[localIdx]) {
            return 0.0;
        }
        const auto partner = static_cast<unsigned>(this->partner_[localIdx]);
        return this->simulator_.problem().referencePorosity(partner, /*timeIdx=*/0)
            * this->simulator_.model().dofTotalVolume(partner);
    }
    void setResidualCheck(const bool on) { this->residualCheck_ = on; }

    typename ParentType::LocalStepTolerance localStepTolerance() const override
    { return this->stepTolerance_; }
    void setLocalStepTolerance(const typename ParentType::LocalStepTolerance& tol)
    { this->stepTolerance_ = tol; }

    double maxTemperatureChange() const override
    { return this->maxTempChange_; }
    void setMaxTemperatureChange(const double dT) { this->maxTempChange_ = dT; }

    void connections(std::vector<Connection>& conns) const override
    { conns.insert(conns.end(), this->connections_.begin(), this->connections_.end()); }

    /*!
     * \brief Start every cell from the state of the rock it cuts through.
     *
     * Called once, when the initial solution is applied -- long before any fracture
     * exists.  Every cell is dormant then, and this only gives the rows something
     * well-defined to hold; a cell's real starting state is set when it is handed out,
     * in bind().
     */
    void applyInitial() override
    {
        auto& solution = this->simulator_.model().solution(/*timeIdx=*/0);

        for (unsigned localIdx = 0; localIdx < this->numDofs(); ++localIdx) {
            this->assignStateFromPartner(solution, localIdx);
        }
    }

    /*!
     * \brief Condition the rows of the cells that are not in use.
     *
     * Nothing else writes to them: they have no volume, so the accumulation term
     * vanishes, and no connections, so no flux reaches them.  An identity row leaves the
     * unknown where it started and keeps the matrix invertible.
     */
    void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual) override
    {
        for (unsigned localIdx = 0; localIdx < this->numDofs(); ++localIdx) {
            if (this->active_[localIdx]) {
                continue;
            }

            const auto globalIdx = static_cast<unsigned>(this->localToGlobalDof(localIdx));

            residual[globalIdx] = 0.0;

            auto* block = matrix.blockAddress(globalIdx, globalIdx);
            if (block == nullptr) {
                OPM_THROW(std::logic_error,
                          fmt::format("Fracture cell {} has no diagonal block to condition",
                                      globalIdx));
            }

            *block = 0.0;
            for (unsigned eq = 0; eq < residual[globalIdx].size(); ++eq) {
                (*block)[eq][eq] = 1.0;
            }
        }

    }

    /*!
     * \brief Hand out degrees of freedom to the fracture's cells and describe them.
     *
     * \return Whether the set of connections changed, and so whether the sparsity pattern
     *         has to be rebuilt.  Aperture changes alone do not: they move the value of a
     *         transmissibility, not the shape of the matrix.
     */
    bool bind(const FractureModel& fractures);

    /*!
     * \brief Refresh apertures, pore volumes and transmissibilities over the existing
     *        binding, without changing slots or sparsity.
     *
     * Returns false if a fracture's cell count changed (the trimesh was renumbered);
     * the caller must then rebind.
     */
    bool updateValues(const FractureModel& fractures);

    /*!
     * \brief Global degrees of freedom the last bind() brought into existence.
     *
     * New cells have no start-of-step state; existing ones do, so a mid-step bind
     * must tell them apart.
     */
    const std::vector<unsigned>& newbornDofs() const
    { return this->newbornDofs_; }

    /*!
     * \brief Whether every fracture still has the cell count its slots were handed out for.
     *
     * Used to hold the topology fixed inside a time step.
     */
    bool layoutMatches(const FractureModel& fractures) const;

    /*!
     * \brief Log the leak-off cross-check between the two representations.
     *
     * Per cell: the reservoir's applied conductance and flux against the fracture
     * solver's own leakof_ and leakOfRate(); the aggregate ratio is the calibration number.
     */
    void leakoffReport(const FractureModel& fractures) const;

    /*!
     * \brief Log what the well does at the fracture's degrees of freedom.
     *
     * Per well: carried fracture perforations, their CTFs, BHP and perforation pressures,
     * and the water share through them. Reads the well state alone.
     */
    void perforationReport() const;
    //! Diagnostic: per-cell state of the binding (pressure vs partner and vs the
    //! fracture's own solve, connections, well index); the nWorst cells by
    //! partner-pressure deviation plus every isolated cell. Prints regardless of
    //! whether the binding layout matches the fracture state.
    void cellDump(const FractureModel& fractures, const std::string& where, int nWorst) const;

    /*!
     * \brief Relative change of the binding between the last two binds.
     *
     * Max of the relative changes of total connection transmissibility and total pore
     * volume; aggregates, so per-cell contact chatter cancels. The outer loop's coupling
     * residual in embedded mode.
     */
    Scalar lastBindChange() const
    { return this->lastBindChange_; }

    /*!
     * \brief The well's perforations of this fracture's cells.
     *
     * Handed over directly rather than through the schedule, since an auxiliary cell
     * has no cartesian index.
     */
    std::vector<RuntimePerforation> wellPerforations(const std::string& wellName) const;
    //! Water pressure of the bound cells of one fracture (binding order), per
    //! fracture cell; empty if the fracture is not bound.
    std::vector<Scalar> cellPressures(std::size_t fractureIdx) const
    {
        std::vector<Scalar> p;
        const auto& model = this->simulator_.model();
        // one entry per fracture cell, bound or not, so the caller can tell a
        // cell without a flow DOF (NaN) from a layout that no longer matches
        for (unsigned slot = 0; slot < this->slotOf_.size() && slot < this->capacity_; ++slot) {
            const auto [fidx, cell] = this->slotOf_[slot];
            if (fidx == fractureIdx && cell >= p.size()) {
                p.resize(cell + 1, std::numeric_limits<Scalar>::quiet_NaN());
            }
        }
        for (unsigned slot = 0; slot < this->slotOf_.size() && slot < this->capacity_; ++slot) {
            const auto [fidx, cell] = this->slotOf_[slot];
            if (fidx != fractureIdx || !this->active_[slot]) {
                continue;
            }
            const auto dof = static_cast<unsigned>(this->localToGlobalDof(slot));
            p[cell] = getValue(model.intensiveQuantities(dof, 0).fluidState().pressure(FluidSystem::waterPhaseIdx));
        }
        return p;
    }
    /*!
     * \brief Gather what the fracture <-> mechanics coupling blocks are built from.
     *
     * The flow side of the numbers is this module's: the apertures and the
     * pressures the flow actually holds, the areas and half transmissibilities
     * of the current grid, and the two floors the binding applies.  Pass the
     * result to buildFlowMechCoupling() / buildMechFlowCoupling(), or to
     * checkCouplingMatricesFD() to have them verified.
     */
    MechCouplingInput couplingInput(const Fracture& fracture,
                                    std::size_t fractureIdx,
                                    Scalar dt) const
    {
        MechCouplingInput in;
        const auto nc = fracture.numCells();
        const auto& width = fracture.fractureWidth();
        const auto areas = fracture.cellAreas();
        const auto& mobility = fracture.reservoirMobility();

        in.htrans = fracture.currentHalfTrans();
        in.cubic_law_min_width = fracture.cubicLawMinWidth();
        in.volume_min_width = this->minWidth_;
        in.dt = dt;

        in.aperture.resize(nc, 0.0);
        in.area.resize(nc, 0.0);
        in.face_mobility.resize(nc, 0.0);
        in.density.resize(nc, 1000.0);
        in.open.assign(nc, 1);
        for (std::size_t c = 0; c < nc; ++c) {
            in.aperture[c] = (c < width.size()) ? width[c][0] : 0.0;
            in.area[c] = (c < areas.size()) ? areas[c] : 0.0;
            in.face_mobility[c] = (c < mobility.size()) ? mobility[c] : 0.0;
        }
        for (const auto c : fracture.closedCells()) {
            if ((c >= 0) && (static_cast<std::size_t>(c) < nc)) {
                in.open[c] = 0;
            }
        }

        // the pressures and the volumes are the flow's, not the fracture's
        in.pressure = this->cellPressures(fractureIdx);
        in.pressure.resize(nc, 0.0);
        // cells without a flow DOF come back as NaN; the coupling wants a number
        for (auto& pv : in.pressure) {
            if (!std::isfinite(pv)) {
                pv = 0.0;
            }
        }
        in.volume_prev.resize(nc, 0.0);
        const auto& model = this->simulator_.model();
        for (unsigned slot = 0; slot < this->slotOf_.size() && slot < this->capacity_; ++slot) {
            const auto [fidx, cell] = this->slotOf_[slot];
            if ((fidx != fractureIdx) || !this->active_[slot] || (cell >= nc)) {
                continue;
            }
            const auto dof = static_cast<unsigned>(this->localToGlobalDof(slot));
            in.volume_prev[cell] = this->bulkVolume_[slot];
            in.density[cell] = getValue(model.intensiveQuantities(dof, 0)
                                        .fluidState().density(FluidSystem::waterPhaseIdx));
        }
        return in;
    }

    /*!
     * \brief Verify the coupling blocks of every bound fracture against finite
     *        differences (opt-in; see solver.check_coupling_fd).
     *
     * Costs one residual evaluation per checked column, so a production run
     * should leave \p opt.max_columns small.
     */
    bool checkCoupling(const FractureModel& fractures,
                       const CouplingCheckOptions& opt,
                       Scalar dt) const
    {
        if (!opt.enabled) {
            return true;
        }
        bool ok = true;
        std::size_t fidx = 0;
        for (const auto& wellFractures : fractures.wellFractures()) {
            for (const auto& fracture : wellFractures) {
                const auto in = this->couplingInput(fracture, fidx++, dt);
                if (in.numCells() == 0) {
                    continue;
                }
                std::vector<CouplingCheckReport> reports;
                const bool fine = checkCouplingMatricesFD(in, opt, &reports);
                ok = ok && fine;
                for (const auto& rep : reports) {
                    OpmLog::info(fmt::format("Fracture {}: {}", fracture.name(), rep.summary()));
                }
            }
        }
        return ok;
    }

    //! Global DOF indices of every active cell of the fractures attached to a well.
    std::vector<int> cellsOfWell(const std::string& wellName) const
    {
        std::vector<int> cells;
        if (const auto pos = this->wellCells_.find(wellName); pos != this->wellCells_.end()) {
            cells = pos->second;
        }
        return cells;
    }

    //! Carry the cells' flow state across a regrid geometrically instead of by slot.
    void setRemapState(const bool on) { this->remapState_ = on; }

    //! Cells handed out so far, for the high-water mark in the log.
    unsigned numActive() const
    { return static_cast<unsigned>(std::count(this->active_.begin(), this->active_.end(), true)); }

private:
    //! useCurrentState: the partner's present state (bind-time, needs valid
    //! intensive quantities) rather than its initial one (model init).
    template <class SolutionVector>
    void assignStateFromPartner(SolutionVector& solution, const unsigned localIdx,
                                const bool useCurrentState = false)
    {
        const auto globalIdx = static_cast<unsigned>(this->localToGlobalDof(localIdx));
        const auto partner = this->partner_.at(localIdx);
        const auto& problem = this->simulator_.problem();

        auto fs = problem.initialFluidState(partner);
        // A cell born mid-run starts from the rock's state NOW: at the initial
        // pressure it sits tens of bar below its partner, the well drops to it and
        // cross-flows from its matrix perforations.
        if (useCurrentState) {
            const auto& cur = this->simulator_.model().intensiveQuantities(partner, /*timeIdx=*/0).fluidState();
            for (unsigned phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }
                fs.setPressure(phase, getValue(cur.pressure(phase)));
                // the fracture volume is created by the injected water and the
                // fracture's own solve treats it as water-filled; starting at the
                // rock's Sw would force a 0 -> 1 saturation transient inside the
                // first step, one limiter-capped Newton iteration at a time
                fs.setSaturation(phase, (phase == FluidSystem::waterPhaseIdx) ? 1.0 : 0.0);
            }
            if constexpr (getPropValue<TypeTag, Properties::EnableEnergy>()) {
                fs.setTemperature(getValue(cur.temperature(0)));
            }
            if (FluidSystem::enableDissolvedGas()) {
                fs.setRs(getValue(cur.Rs()));
            }
            if (FluidSystem::enableVaporizedOil()) {
                fs.setRv(getValue(cur.Rv()));
            }
        }

        // Carry the phase pressures to the fracture cell's own depth; the fluid is the
        // rock's, so nothing else about the state changes.
        const auto waterPos = FluidSystem::waterPhaseIdx;
        const auto rho = getValue(fs.density(waterPos));
        const auto gravity = problem.gravity()[dimWorld - 1];
        const auto dz = this->depth_.at(localIdx) - problem.dofCenterDepth(partner);

        for (unsigned phase = 0; phase < FluidSystem::numPhases; ++phase) {
            if (!FluidSystem::phaseIsActive(phase)) {
                continue;
            }

            fs.setPressure(phase, getValue(fs.pressure(phase)) + rho * gravity * dz);
        }

        solution[globalIdx].setPvtRegionIndex(this->pvtRegionIndex(localIdx));
        solution[globalIdx].assignNaive(fs);
    }

    Simulator& simulator_;

    unsigned capacity_{};
    Scalar minWidth_{};
    // Heat exchange between fracture fluid and wall rock (geometric half transmissibilities).
    // Without it a cell with little net flow has an unanchored energy row and the Newton
    // does not converge. solver.embedded_wall_conduction (default true); the rock side uses
    // the leak-off distance unless solver.embedded_thermal_distance (m) is given.
    std::pair<Scalar, Scalar> wallConduction_(const auto& fracture, std::size_t cell,
                                              Scalar area, Scalar width) const
    {
        if (wallConductionOn_ < 0) {
            const auto prm = this->simulator_.problem().getFractureParam();
            wallConductionOn_ = prm.template get<bool>("solver.embedded_wall_conduction", true) ? 1 : 0;
            thermalDistance_ = static_cast<Scalar>(
                prm.template get<double>("solver.embedded_thermal_distance", 0.0));
        }
        const auto& dist = fracture.reservoirDistance();
        if (wallConductionOn_ == 0 || area <= 0 || cell >= dist.size()) {
            return {Scalar{0}, Scalar{0}};
        }
        const auto sides = static_cast<Scalar>(fracture.leakingSides());
        const auto d = (thermalDistance_ > 0) ? thermalDistance_ : static_cast<Scalar>(dist[cell]);
        if (!(d > 0)) {
            return {Scalar{0}, Scalar{0}};
        }
        const auto halfWidth = std::max(width, this->minWidth_) / 2;
        return {sides * area / halfWidth, sides * area / d};
    }
    mutable int wallConductionOn_{-1};
    mutable Scalar thermalDistance_{0};

    PerfWiMode perfWiMode_{PerfWiMode::Fracture};
    Scalar perfWidth_{};
    Scalar perfRw_{};

    std::vector<bool> active_{};

    Scalar lastTotalTrans_ = -1.0;
    Scalar lastTotalPv_ = -1.0;
    Scalar lastBindChange_ = 0.0;
    std::vector<Scalar> bulkVolume_{};
    std::vector<Scalar> depth_{};
    //! Reservoir cell each fracture cell leaks into; also where its initial state
    //! comes from.  Dormant cells keep cell zero, which is only ever read to give the
    //! row a defined state.
    std::vector<unsigned> partner_{};

    std::vector<Connection> connections_{};

    //! (fracture index within the model, cell index within that fracture) -> local slot.
    std::vector<std::pair<std::size_t, std::size_t>> slotOf_{};

    //! Well name -> the perforations of that well's fractures, in degrees of freedom.
    std::map<std::string, std::vector<RuntimePerforation>> wellPerforations_{};
    std::map<std::string, std::vector<int>> wellCells_{}; // active aux DOFs per well
    std::vector<unsigned> newbornDofs_{}; // see newbornDofs()
    std::vector<std::pair<unsigned, std::size_t>> boundLayout_{}; // per fracture: first slot, cells
    bool remapState_{false}; // see setRemapState()
    bool residualCheck_{false}; // see cnvReferencePoreVolume()
    typename ParentType::LocalStepTolerance stepTolerance_{};
    double maxTempChange_{0.0}; // Newton step limit [K]; <= 0: the global one
};

} // namespace Opm

#include <opm/geomech/FractureAuxCells_impl.hpp>

#endif // OPM_FRACTURE_AUX_CELLS_HPP
