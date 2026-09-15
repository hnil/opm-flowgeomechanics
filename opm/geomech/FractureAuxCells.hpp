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
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class FractureModel;

/*!
 * \brief Fracture cells as degrees of freedom of the reservoir flow problem.
 *
 * The alternative to upscaling the fracture to well indices.  In the upscaled form the
 * fracture is not part of the flow problem at all: its conductance is folded into a well
 * index, \f$WI = q/\Delta p\f$, which is sign-indefinite and has to be hard-zeroed when
 * it comes out negative.  Here the fracture cells carry the reservoir's own conservation
 * equations, and what used to be an upscaled index is an ordinary transmissibility on an
 * ordinary connection.
 *
 * \section capacity Preallocation
 *
 * A fracture grows, and a degree of freedom cannot be added to a system that has already
 * been sized.  So a fixed number of them is claimed up front and handed out as cells
 * open.  A cell that has not been handed out is *dormant*: it has no volume and no
 * connections, so nothing in the ordinary assembly writes to its row, and this module
 * puts an identity there to keep the matrix non-singular.  Running out of capacity stops
 * the run rather than silently capping the fracture -- growing the allocation is a
 * follow-up, and a fracture quietly prevented from growing is exactly the failure this
 * whole line of work exists to avoid.
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
     * \brief Refresh what the flow is fed, without touching the shape.
     *
     * Apertures, and with them pore volumes and cubic-law transmissibilities,
     * move at every iteration of the coupled mechanics-pressure solve; the set
     * of cells only changes when the fracture grows.  This recomputes the
     * former over the binding that exists, leaving the slot map, the active
     * set and the sparsity pattern alone, so it is cheap and needs no matrix
     * rebuild.
     *
     * Requires the layout to still match the fracture: once the cell count
     * changes the trimesh has been renumbered and a cell index no longer means
     * the same cell, so there is nothing well defined to update.  Returns false
     * in that case, and the caller must rebind instead.
     */
    bool updateValues(const FractureModel& fractures);

    /*!
     * \brief Global degrees of freedom the last bind() brought into existence.
     *
     * A cell that has just appeared has no start-of-step state; one that was
     * already there does, and it is the reference its accumulation term is
     * measured against.  The two have to be told apart when a bind happens
     * inside a step.
     */
    const std::vector<unsigned>& newbornDofs() const
    { return this->newbornDofs_; }

    /*!
     * \brief Whether the fractures still have the shape of the current binding.
     *
     * True when every fracture has the cell count its slots were handed out for.  Used
     * to hold the topology fixed inside a time step: a fracture solve that grew or
     * reset its grid mid-step changes the shape of the flow system, and re-describing
     * that to a Newton iteration already under way gives it a moving target.  Such
     * changes wait for the step boundary.
     */
    bool layoutMatches(const FractureModel& fractures) const;

    /*!
     * \brief Log the leak-off cross-check between the two representations.
     *
     * For the same fracture, per cell: the conductance the reservoir actually
     * applies (bound transmissibility times the upwind water mobility) against the
     * fracture solver's own leakof_ (transmissibility times the reservoir total
     * mobility), and the two flux opinions -- the reservoir's, evaluated with the
     * same local residual that assembled it, and the fracture's internal
     * leakOfRate().  The aggregate ratio is the calibration number: it says whose
     * conductance the pressure difference between the modes belongs to.
     */
    void leakoffReport(const FractureModel& fractures) const;

    /*!
     * \brief Log what the well does at the fracture's degrees of freedom.
     *
     * Per well: how many of the registered fracture perforations the well actually
     * carries, their transmissibility factors, the bottom-hole and perforation
     * pressures, and the share of the well's water that goes through them.  It reads
     * the well state alone, so it says whether the fracture is connected regardless of
     * whether the binding still describes the fracture's current shape.
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
     * The maximum of the relative changes of total reservoir-connection
     * transmissibility and total pore volume -- aggregates, so that contact chatter of
     * individual cells cancels the way it does in the upscaled well index.  This is the
     * coupling residual the outer loop watches in embedded mode: it measures what the
     * flow is actually fed, where the well-index change list measures a quantity that
     * is never applied.
     */
    Scalar lastBindChange() const
    { return this->lastBindChange_; }

    /*!
     * \brief The well's perforations of this fracture's cells.
     *
     * Each is a degree of freedom of the flow problem and a well index, so the well model
     * treats it exactly as it treats a perforation of a grid cell -- it reads intensive
     * quantities by index, which an auxiliary cell answers as well as any other.  What it
     * cannot do is arrive at one through a cartesian index, which is why these are handed
     * over directly instead of going through the schedule.
     */
    std::vector<RuntimePerforation> wellPerforations(const std::string& wellName) const;
    //! Water pressure of the bound cells of one fracture (binding order), per
    //! fracture cell; empty if the fracture is not bound.
    std::vector<Scalar> cellPressures(std::size_t fractureIdx) const
    {
        std::vector<Scalar> p;
        const auto& model = this->simulator_.model();
        for (unsigned slot = 0; slot < this->slotOf_.size() && slot < this->capacity_; ++slot) {
            const auto [fidx, cell] = this->slotOf_[slot];
            if (fidx != fractureIdx || !this->active_[slot]) {
                continue;
            }
            if (cell >= p.size()) {
                p.resize(cell + 1, std::numeric_limits<Scalar>::quiet_NaN());
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
};

} // namespace Opm

#include <opm/geomech/FractureAuxCells_impl.hpp>

#endif // OPM_FRACTURE_AUX_CELLS_HPP
