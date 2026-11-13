/*
  Copyright 2023 Equinor ASA.

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

#ifndef OPM_GEOMECH_MODEL_HPP
#define OPM_GEOMECH_MODEL_HPP
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/geomech/VemElasticitySolver.hpp>

#include <opm/geomech/FlowGeoMechLinearSolverParameters.hpp>
#include <opm/geomech/FractureMechHost.hpp>
#include <opm/geomech/MechFlowMap.hpp>
#include <opm/geomech/MechGridContext.hpp>
#include <opm/geomech/FractureModel.hpp>
#include <opm/geomech/MechOptions.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>

#include <algorithm>
#include <limits>
#include <array>
#include <cmath>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <sstream>

namespace Opm{
    template<typename TypeTag>
    class GeoMechModel : public BaseAuxiliaryModule<TypeTag>
    {
        using Parent = BaseAuxiliaryModule<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using Stencil = GetPropType<TypeTag, Properties::Stencil>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;

        template <class FluidState>
        static int referencePhaseIdx(const FluidState& fs)
        {
            if (fs.phaseIsActive(FluidSystem::waterPhaseIdx)) {
                return FluidSystem::waterPhaseIdx;
            }
            if (fs.phaseIsActive(FluidSystem::oilPhaseIdx)) {
                return FluidSystem::oilPhaseIdx;
            }
            return FluidSystem::gasPhaseIdx;
        }
    public:
        GeoMechModel(Simulator& simulator):
            simulator_(simulator)
        {

        };
        //ax model things
        void postSolve(GlobalEqVector&) {
        }
        void addNeighbors(std::vector<NeighborSet>&) const {
        };
        void applyInitial(){
         };
        unsigned numDofs() const{return 0;};
        void linearize(SparseMatrixAdapter&, GlobalEqVector&) {
        };

        // model things
        void beginIteration(){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech begin iteration\n");
            }
        }
        void endIteration(){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech end iteration\n");
            }
        }
        void beginTimeStep(){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech begin time step\n");
            }
            solved_in_step_ = false;
            fracHost_.beginTimeStep();
        }
        void endTimeStep(){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech model endstimeStep\n");
            }
            if (solved_in_step_) {
                solved_in_step_ = false;
                return;
            }
            if (this->skipEndOfStepSolve_()) {
                return;
            }
            const double eos_area_fac = simulator_.problem().getGeoMechParam()
                .template get<double>("fractureparam.solver.end_of_step_area_change_fac", -1.0);
            if (eos_area_fac > 0.0 && this->fractureModelActive()) {
                this->fractureModel().setAreaChangeFactorOverride(eos_area_fac);
            }
            try {
                this->solveGeoMechAndFracture();
            } catch (...) {
                if (eos_area_fac > 0.0 && this->fractureModelActive())
                    this->fractureModel().setAreaChangeFactorOverride(-1.0);
                throw;
            }
            if (eos_area_fac > 0.0 && this->fractureModelActive()) {
                this->fractureModel().setAreaChangeFactorOverride(-1.0);
            }
        }

        //! "off" reuses this step's last converged in-loop fracture state; "tolerance" also
        //! needs p and T in the fracture cells to have moved less than the tolerances since.
        bool skipEndOfStepSolve_()
        {
            const auto& prm = simulator_.problem().getGeoMechParam();
            const std::string mode
                = prm.template get<std::string>("fractureparam.solver.end_of_step_fracture", "always");
            if (mode == "always") {
                return false;
            }
            if (mode != "off" && mode != "tolerance") {
                OPM_THROW(std::runtime_error, "Unknown fractureparam.solver.end_of_step_fracture: " + mode);
            }
            const auto& comm = simulator_.gridView().comm();
            int reusable = this->fractureModelActive() && fracHost_.solvesThisStep() > 0
                && !this->fractureModel().anyLastSolveUnconverged();
            reusable = comm.min(reusable);
            double dp = std::numeric_limits<double>::infinity();
            double dT = dp;
            if (reusable) {
                std::tie(dp, dT) = this->fractureModel()
                    .template drivingStateChange<TypeTag, Simulator>(simulator_, /*record*/ false);
                dp = comm.max(dp);
                dT = comm.max(dT);
            }
            bool skip = reusable;
            if (mode == "tolerance") {
                skip = skip
                    && dp <= prm.template get<double>("fractureparam.solver.end_of_step_pressure_tol", 1.0e5)
                    && dT <= prm.template get<double>("fractureparam.solver.end_of_step_temperature_tol", 1.0);
            }
            if (comm.rank() == 0) {
                OpmLog::info(fmt::format("End-of-step fracture solve {} ({}; max dp {:.3g} bar, dT {:.3g} K)",
                                         skip ? "skipped" : "done",
                                         reusable ? "in-loop state converged" : "no converged in-loop state",
                                         dp / 1.0e5, dT));
            }
            return skip;
        }
        struct PostSolveResult {
            bool mech_converged = true;
            bool fracture_converged = true;
        };
        //! PostSolve inside the step's retry scope, so the caller can still reject the step.
        PostSolveResult solveGeoMechAndFractureInStep(){
            solved_in_step_ = true;
            return this->solveGeoMechAndFracture();
        }
        //! solver.postsolve_mech_interval: "step" (default) or "report", which
        //! solves mechanics only on the last substep of a report step.
        bool mechSolveDueAtStepEnd_() const {
            const std::string interval = simulator_.problem().getGeoMechParam()
                .template get<std::string>("solver.postsolve_mech_interval", "step");
            if (interval == "step")
                return true;
            if (interval != "report")
                OPM_THROW(std::invalid_argument, "solver.postsolve_mech_interval must be 'step' or 'report'");
            return simulator_.episodeWillBeOver();
        }
        PostSolveResult solveGeoMechAndFracture(){
            PostSolveResult result;
            const auto& problem = simulator_.problem();
            if (this->mechSolveDueAtStepEnd_()) {
                result.mech_converged = this->solveGeomechanics().converged;
                if (!result.mech_converged && simulator_.gridView().comm().rank() == 0) {
                    OpmLog::warning(fmt::format("End-of-step mechanical solve not converged "
                                                "({} unconverged mechanical solves so far)",
                                                num_mech_unconverged_));
                }
            }
            // the joint solve owns the fractures: growth only on a converged per-well state
            if(problem.hasFractures()
               && !problem.getGeoMechParam().template get<bool>("fractureparam.solver.well_joint_solve", false)){
                this->solveFractures();
                if (this->fractureModelActive()) {
                    result.fracture_converged = !this->fractureModel().anyLastSolveUnconverged();
                }
            }
            return result;
        }

       //! Fracture driving is hosted by FractureMechHost (shared between
       //! mechanics backends); the public surface is forwarded unchanged.
       void solveFractures(){ fracHost_.solveFractures(); }

       void writeFractureSolutionFirst(){ fracHost_.writeFractureSolutionFirst(); }

       void writeFractureSolution(){ fracHost_.writeFractureSolution(); }

       std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname){
           return fracHost_.getExtraWellIndices(wellname);
       }

       std::vector<RuntimePerforation> getFlowWellIndices(const std::string& wellname){
           return fracHost_.getFlowWellIndices(wellname);
       }

       bool fractureModelActive() const{ return fracHost_.fractureModelActive(); }

       void updateFilterCakePropertiesOnFractures(){ fracHost_.updateFilterCakePropertiesOnFractures(); }

       const FractureModel& fractureModel() const{ return fracHost_.fractureModel(); }

       FractureModel& fractureModel() { return fracHost_.fractureModel(); }

       void resetFractureModel(){ fracHost_.resetFractureModel(); }

        //! Mechanics runs on its own grid when one was built, otherwise on
        //! the flow grid.
        const Grid& mechGrid() const {
            return mechCtx_ ? mechCtx_->grid() : simulator_.vanguard().grid();
        }

        //! VEM needs a body without voids, cracks or hanging nodes; checked on
        //! the whole grid on rank 0.
        void checkMechGridIsGeometric_() const
        {
            const std::string mode = Parameters::Get<Parameters::MechGridCheck>();
            if (mode == "off") {
                return;
            }
            if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {
                auto& grid = const_cast<Grid&>(this->mechGrid());
                const auto& comm = simulator_.gridView().comm();
                int fatal = 0;
                int defects = 0;
                std::string message;
                if (comm.rank() == 0) {
                    if (comm.size() > 1) {
                        grid.switchToGlobalView();
                    }
                    const auto check = Dune::cpgrid::checkGeometric(grid);
                    if (comm.size() > 1) {
                        grid.switchToDistributedView();
                    }
                    // VEM needs every cell closed and edge-conformal: no defect is benign.
                    fatal = !check.ok();
                    defects = fatal;
                    message = "Mechanics grid: " + check.summary();
                }
                fatal = comm.max(fatal);
                defects = comm.max(defects);
                if (comm.rank() == 0) {
                    if (fatal && mode == "error") {
                        OpmLog::error(message);
                    } else if (defects) {
                        OpmLog::warning(message);
                    } else {
                        OpmLog::info(message);
                    }
                }
                if (fatal && mode == "error") {
                    OPM_THROW_NOLOG(std::runtime_error,
                                    "The mechanics grid is not a geometric body: voids, cracks, "
                                    "broken cells or hanging nodes (see the log). Process with "
                                    "--edge-conformal=true; --mech-grid-check=warn runs anyway.");
                }
            }
        }

        //! Must be called before init().
        void setMechGrid(const MechGridContext& ctx){
            mechCtx_ = &ctx;
            mechMap_ = ctx.map();
        }

        void updatePotentialForces(bool relative_solve = true){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Update Forces for Geomechanics");
            }
            size_t numDof = simulator_.model().numGridDof();
            const auto& problem = simulator_.problem();
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = simulator_.model().intensiveQuantities(dofIdx,0);
                // pressure part
                const auto& fs = iq.fluidState();
                const int phaseIdx = referencePhaseIdx(fs);
                const auto& press = fs.pressure(phaseIdx);
                //thermal part
                const auto& poelCoef = problem.poelCoef(dofIdx);
                double pressure = Toolbox::value(press);
                double diffpress = pressure;
                if(relative_solve){
                    diffpress -=  problem.initPressure(dofIdx);
                }

		        const auto& pratio = problem.pRatio(dofIdx);
                double fac = (1-pratio)/(1-2*pratio);
                double pcoeff = poelCoef*fac;
		//assert pcoeff == biot
                pressure_[dofIdx] = pressure;
                mechPotentialForce_[dofIdx] = diffpress*pcoeff;
                mechPotentialPressForce_[dofIdx] = diffpress*pcoeff;
                assert(pcoeff<=(1.0+1e-12));
                mechPotentialPressForceFracture_[dofIdx] = diffpress*(1.0-pcoeff);
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();

                if(thermal_expansion){
                    OPM_TIMEBLOCK(addTermalParametersToMech);
                    const auto& temp = fs.temperature(phaseIdx);//NB all phases have equal temperature
                    const auto& thelcoef = problem.thelCoef(dofIdx);
		    // tcoeff = (youngs*tempExp/(1-pratio))*fac;
                    double tcoeff = thelcoef*fac;
                    // assume difftemp = 0 for non termal runs
                    double difftemp = (Toolbox::value(temp));
                    if(relative_solve){
                        difftemp -= problem.initTemperature(dofIdx);
                    }

                    mechPotentialForce_[dofIdx] += difftemp*tcoeff;
                    mechPotentialTempForce_[dofIdx] = difftemp*tcoeff;
                }
                // Flow-only cells (numerical-aquifer cells etc.) carry no
                // mechanical loads: they stay in the stiffness as inert rock
                // but must not push on the mechanics (their pressure/PV/depth
                // are flow bookkeeping, not rock state).
                if (problem.mechFlowOnlyCell(dofIdx)) {
                    mechPotentialForce_[dofIdx] = 0.0;
                    mechPotentialPressForce_[dofIdx] = 0.0;
                    mechPotentialPressForceFracture_[dofIdx] = 0.0;
                    mechPotentialTempForce_[dofIdx] = 0.0;
                }
            }
            const auto& opts = this->mechOptions();
            const bool smooth = opts.smooth_force;
            // Smoothing length in meters: smooths under-resolved fronts over a fixed physical
            // scale (e.g. sqrt(diffusivity*t)) so stress extrema are grid independent.
            // Takes precedence over the grid-dependent one-ring 'smooth_force'.
            const double smooth_length = opts.smooth_force_length;
            // Safeguard: force smoothing averages neighbor values, so the
            // zeroed loads of flow-only cells would bleed into the adjacent
            // rock.  Warn loudly rather than silently distort the loads.
            if ((smooth || smooth_length > 0.0) && problem.hasMechFlowOnlyCells()
                && simulator_.gridView().comm().rank() == 0) {
                OpmLog::warning("smooth_force(_length) is active while flow-only "
                                "(mech-masked) cells exist: smoothing will mix their "
                                "zero loads into neighboring rock cells near the "
                                "aquifer attachment.");
            }
            if(smooth_length > 0.0){
              mechPotentialForce_ = vem::diffuseCellVector(simulator_.vanguard().grid(),mechPotentialForce_,smooth_length);
            } else if(smooth){
              mechPotentialForce_ = vem::smoothCellVector(simulator_.vanguard().grid(),mechPotentialForce_);
            }
            // Smoothing belongs on the flow grid; the mechanics grid sees the
            // volume-weighted mean, which keeps the load integral.
            if(!mechMap_.isIdentity()){
                mechMap_.restrict(mechPotentialForce_, mechPotentialForceMech_);
                // Restriction only reaches cells whose flow children are local.
                mechCtx_->communicate(mechPotentialForceMech_);
            }
        }

        //! The load the mechanics solver assembles: flow cells when the grids
        //! are the same, mechanics cells otherwise.
        const Dune::BlockVector<Dune::FieldVector<double,1>>& mechLoad() const {
            return mechMap_.isIdentity() ? mechPotentialForce_ : mechPotentialForceMech_;
        }

        //! Fields on the mechanics grid's own cells and vertices (output of a separate grid).
        const MechGridContext* mechContext() const { return mechCtx_; }
        const MechFlowMap& mechFlowMap() const { return mechMap_; }
        const Dune::BlockVector<Dune::FieldVector<double,6>>& mechCellLinStress() const { return mechState_.stress; }
        const Dune::BlockVector<Dune::FieldVector<double,6>>& mechCellStrain() const { return mechState_.strain; }
        const Dune::BlockVector<Dune::FieldVector<double,3>>& mechVertexDisplacement() const { return displacement_; }
        const std::vector<double>& mechCellYModule() const { return ymoduleMech_; }
        const std::vector<double>& mechCellPRatio() const { return pratioMech_; }

        //! Mechanics index of a flow cell; the identity while the grids agree.
        std::size_t mechIdx(std::size_t globalIdx) const {
            return static_cast<std::size_t>(mechMap_.mechCell(globalIdx));
        }

        //! False for a flow cell outside a mechanics sector; its mechanics output is zero.
        bool hasMech(std::size_t globalIdx) const {
            return mechMap_.mechCell(globalIdx) >= 0;
        }
        //! Operators once; the matrix (and its preconditioner) again only when the
        //! body force differs from the one it was assembled with.
        void prepareMechSolve_(bool use_body_force){
                const double gravity = use_body_force ? Opm::unit::gravity : 0.0;
                if (assembledGravity_ && *assembledGravity_ == gravity) {
                    return;
                }
                const auto& problem = simulator_.problem();
                const auto& opts = this->mechOptions();
                const std::string& face_quad = opts.face_quadrature;
                if (face_quad == "triangulated") {
                    vem::set_face_quadrature(vem::FaceQuadrature::TRIANGULATED);
                } else if (face_quad == "kite") {
                    vem::set_face_quadrature(vem::FaceQuadrature::KITE);
                } else {
                    OPM_THROW(std::invalid_argument, "vem_face_quadrature must be 'kite' or 'triangulated', got " + face_quad);
                }
                OpmLog::info("VEM face quadrature: " + face_quad);
                mechDisc_->setOptions(opts.vem_stability_choice,
                                      opts.stab_on_stress,
                                      opts.diagonal_scaling);
                OPM_TIMEBLOCK(SetupMechSolver);
                if(use_body_force){
                    OpmLog::info("Using body force in geomechanics to 2000 kg/m^3 \n");
                }
                int nc = this->mechGrid().leafGridView().size(0);
                std::vector<double> density(nc, 2000.0);
                mechDisc_->setBodyForce(gravity, density);
                mechDisc_->fixNodes(problem.bcNodes());
                // a cell that fails to assemble throws on its owner only; the solver setup below is collective
                OPM_BEGIN_PARALLEL_TRY_CATCH();
                mechDisc_->assembleMatrix();
                if (!mechDisc_->hasOperators()) {
                    mechDisc_->buildOperators();
                }
                OPM_END_PARALLEL_TRY_CATCH("Mechanics assembly failed: ", simulator_.vanguard().grid().comm());
                FlowLinearSolverParametersGeoMech p;
                p.init<TypeTag>();
                PropertyTree prm = setupPropertyTree(p, true, true);
                if (p.linear_solver_print_json_definition_) {
                    if(simulator_.gridView().comm().rank() == 0){
                        std::ostringstream os;
                        os << "Property tree for mech linear solver:\n";
                        prm.write_json(os, true);
                        OpmLog::note(os.str());
                    }
                }
                mechLinSolver_->setup(prm, mechDisc_->matrix(), mechDisc_->coords(), mechDisc_->rotationModes());
                mechLinSolver_->comm()->communicator().barrier();
                assembledGravity_ = gravity;
                write_system_ = prm.get<int>("verbosity", 0) > 10;
        }
        void writeMechSystem(){
            OPM_TIMEBLOCK(WriteMechSystem);
                    const auto& problem = simulator_.problem();
                    Opm::Helper::writeMechSystem(simulator_,
                    mechDisc_->matrix(),
                    mechRhs_,
                    mechLinSolver_->comm());
                    {
                        int num_points = this->mechGrid().size(3);
                        Dune::BlockVector<Dune::FieldVector<double,1>> fixed(3*num_points);
                        fixed  = 0.0;
                        const auto& bcnodes = problem.bcNodes();
                        for(const auto& node: bcnodes){
                            for(int i = 0; i < 3; ++i){
                                fixed[3*node.node+i][0] = node.fixed[i];
                            }
                        }
                        Opm::Helper::writeVector(simulator_,
                                                fixed,
                                                "fixed_values_",
                                                mechLinSolver_->comm());
                        // node coordinates in dof order, for offline solver studies
                        Dune::BlockVector<Dune::FieldVector<double,1>> xyz(3*num_points);
                        const auto& gv = this->mechGrid().leafGridView();
                        for (const auto& v : vertices(gv)) {
                            const auto idx = gv.indexSet().index(v);
                            const auto c = v.geometry().center();
                            for (int i = 0; i < 3; ++i) {
                                xyz[3*idx+i][0] = c[i];
                            }
                        }
                        Opm::Helper::writeVector(simulator_, xyz, "coords_", mechLinSolver_->comm());
                    }
        }
      void setOutputPutStress(const Dune::BlockVector<Dune::FieldVector<double,6> >& initstress){
            // used at if initial stress not calculated
            outputstress_ = initstress;
      }
      void calculateOutputQuantitiesMech(){
            OPM_TIMEBLOCK(CalculateOutputQuantitesMech);
            const auto& grid = this->mechGrid();
            const Opm::Elasticity::Vector& field = mechState_.displacement;

            this->makeDisplacement(field);
            // update variables used for output to resinsight
            {
            OPM_TIMEBLOCK(calculateStress);
            mechDisc_->stress(field, mechState_.stress);
            mechDisc_->strain(field, mechState_.strain);
            }
            // Optional patch recovery of the cell-wise stress: reduces grid-dependent
            // jumps at under-resolved fronts, and so the tractions the fracture samples.
            if(this->mechOptions().patch_recovery_stress){
                OPM_TIMEBLOCK(patchRecoveryStress);
                mechState_.stress = vem::patchRecoveryCells6(grid, mechState_.stress);
            }
            /// NB need to be updated after linstress to be correct
            // output stress is saved here in case init stress is changed before
            // output, and is reported per flow cell: it adds that cell's own
            // initial stress and load to the mechanics cell's solution.
            for(size_t dofIdx = 0; dofIdx < simulator_.model().numGridDof(); ++dofIdx){
                outputstress_[dofIdx] = this->stress(dofIdx);
            }
            //NB ssume initial strain is 0

            bool verbose = false;
            if(verbose){
                OPM_TIMEBLOCK(WriteMatrixMarket);
                // debug output to matrixmaket format
                Dune::storeMatrixMarket(mechDisc_->matrix(), "A.mtx");
                Dune::storeMatrixMarket(mechRhs_, "b.mtx");
                Dune::storeMatrixMarket(mechState_.displacement, "u.mtx");
                Dune::storeMatrixMarket(field, "field.mtx");
                Dune::storeMatrixMarket(mechPotentialForce_, "pressforce.mtx");
            }
        }


        //! A failure on any rank is a failure on all: the displacement is shared.
        Opm::Elasticity::MechSolveResult mechSolveOnAllRanks_(bool cold)
        {
            mechDisc_->toLocal(mechState_.displacement);
            auto result = mechLinSolver_->solve(mechRhs_, mechDisc_->scale(),
                                                mechState_.displacement, cold);
            mechDisc_->toGlobal(mechState_.displacement);
            const auto& comm = simulator_.gridView().comm();
            result.converged = comm.min(result.converged ? 1 : 0) == 1;
            result.failed = comm.max(result.failed ? 1 : 0) == 1;
            return result;
        }

        Opm::Elasticity::MechSolveResult solveGeomechanics(bool use_body_force = false, bool relative_solve = true){
            {
                OPM_TIMEBLOCK(endTimeStepMech);
                this->updatePotentialForces(relative_solve);
                this->prepareMechSolve_(use_body_force);
                OPM_TIMEBLOCK(AssembleRhs);
                mechDisc_->rhs(this->mechLoad(), mechRhs_);
            }
            Opm::Elasticity::MechSolveResult result;
            {
                OPM_TIMEBLOCK(SolveMechanicalSystem);
                Dune::Timer solve_timer;
                result = this->mechSolveOnAllRanks_(/*cold*/ false);
                if (!result.converged) {
                    total_mech_linear_iterations_ += result.iterations;
                    if (simulator_.gridView().comm().rank() == 0)
                        OpmLog::warning("Mechanical solve not converged; retrying from a zero start");
                    result = this->mechSolveOnAllRanks_(/*cold*/ true);
                }
                ++num_mech_solves_;
                if (!result.converged)
                    ++num_mech_unconverged_;
                total_mech_linear_iterations_ += result.iterations;
                last_mechanical_solve_time_seconds_ = solve_timer.stop();
                total_mechanical_solve_time_seconds_ += last_mechanical_solve_time_seconds_;
                if(simulator_.gridView().comm().rank() == 0){
                    std::ostringstream os;
                    os << "Mechanical solve stats: solves=" << num_mech_solves_
                       << ", linear_iterations=" << result.iterations
                       << " (total " << total_mech_linear_iterations_ << ")"
                       << ", solve_time_s=" << last_mechanical_solve_time_seconds_
                       << " (total " << total_mechanical_solve_time_seconds_ << ")"
                       << ", converged=" << (result.converged ? "true" : "false");
                    OpmLog::info(os.str());
                }
                if(write_system_){
                    this->writeMechSystem();
                }
            }
            this->calculateOutputQuantitiesMech();
            return result;
        }
        int numMechUnconverged() const { return num_mech_unconverged_; }

        //! The operators and the preconditioner are not saved: the first solve after a
        //! load assembles them and sets up this run's linear solver.
        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            namespace ckpt = GeomechCheckpoint;
            ckpt::fieldVectors(serializer, pressure_);
            ckpt::fieldVectors(serializer, mechPotentialForce_);
            ckpt::fieldVectors(serializer, mechPotentialForceMech_);
            ckpt::fieldVectors(serializer, mechPotentialPressForce_);
            ckpt::fieldVectors(serializer, mechPotentialPressForceFracture_);
            ckpt::fieldVectors(serializer, mechPotentialTempForce_);
            ckpt::fieldVectors(serializer, celldisplacement_);
            ckpt::fieldVectors(serializer, displacement_);
            ckpt::fieldVectors(serializer, flowDisplacement_);
            ckpt::fieldVectors(serializer, mechState_.stress);
            ckpt::fieldVectors(serializer, outputstress_);
            ckpt::fieldVectors(serializer, mechState_.strain);
            ckpt::fieldVectors(serializer, mechState_.displacement);
            serializer(fracHost_);
        }

        // used in FlowProblemGeoMech
        void init(bool /*restart*/){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech init\n");
            }
            size_t numDof = simulator_.model().numGridDof();
            pressure_.resize(numDof);
            mechPotentialForce_.resize(numDof);
            mechPotentialTempForce_.resize(numDof);
            mechPotentialPressForce_.resize(numDof);
            mechPotentialPressForceFracture_.resize(numDof);
            // hopefully temperature and pressure initilized
            if(!mechCtx_){
                mechMap_ = MechFlowMap::identity(numDof);
            }
            checkMechGridIsGeometric_();
#if HAVE_MPI
            Opm::Elasticity::validateMechanicalPartition(this->mechGrid());
#endif
            mechDisc_ = std::make_unique<Opm::Elasticity::VemDiscretisation<Grid>>(this->mechGrid());
            mechLinSolver_ = std::make_unique<Opm::Elasticity::MechLinearSolver<Grid>>(this->mechGrid());

            // Displacement, stress and strain live on the mechanics grid; the
            // loads above stay on the flow grid and are restricted to it.
            const size_t numMech = mechMap_.numMechCells();
            celldisplacement_.resize(numMech);
            std::fill(celldisplacement_.begin(),celldisplacement_.end(),0.0);
            mechState_.stress.resize(numMech);
            std::fill(mechState_.stress.begin(),mechState_.stress.end(),0.0);
            outputstress_.resize(numDof);
            std::fill(outputstress_.begin(),outputstress_.end(),0.0);
            mechState_.strain.resize(numMech);
            std::fill(mechState_.strain.begin(),mechState_.strain.end(),0.0);
            const auto& gv = this->mechGrid().leafGridView();
            displacement_.resize(gv.indexSet().size(3));
            if (!mechMap_.isIdentity()) {
                flowDisplacement_.resize(simulator_.gridView().size(3));
                flowDisplacement_ = 0.0;
            }
        };

        //! The Cartesian arrays supply mechanics cells outside the flow grid.
        void setMaterial(const std::vector<double>& ymodule,const std::vector<double>& pratio,
                         const std::vector<double>* ymoduleCartesian = nullptr,
                         const std::vector<double>* pratioCartesian = nullptr){
            if(mechMap_.isIdentity()){
                mechDisc_->setMaterial(ymodule,pratio);
                return;
            }
            std::vector<double> ymodule_mech, pratio_mech;
            mechMap_.restrict(ymodule, ymodule_mech);
            mechMap_.restrict(pratio, pratio_mech);
            if (ymoduleCartesian && pratioCartesian) {
                const auto [fromDeck, fromNeighbours] =
                    mechCtx_->fillMaterial(*ymoduleCartesian, *pratioCartesian,
                                           ymodule_mech, pratio_mech);
                if (fromDeck + fromNeighbours > 0) {
                    OpmLog::info("Mechanics grid: cells outside the flow grid take YMODULE "
                                 "and PRATIO from the deck (" + std::to_string(fromDeck)
                                 + ") or, where the deck has none, from their neighbours ("
                                 + std::to_string(fromNeighbours) + ")");
                }
            }
            std::size_t noStiffness = 0;
            for (const auto e : ymodule_mech) {
                noStiffness += (e <= 0.0);
            }
            if (noStiffness > 0) {
                OpmLog::warning("Mechanics grid: " + std::to_string(noStiffness)
                                + " cells have no YMODULE (outside the flow grid, no deck value)");
            }
            mechCtx_->communicate(ymodule_mech);
            mechCtx_->communicate(pratio_mech);
            ymoduleMech_ = ymodule_mech;
            pratioMech_ = pratio_mech;
            mechDisc_->setMaterial(ymodule_mech,pratio_mech);
        }
        /// Displacement of a flow grid vertex.
        const Dune::FieldVector<double,3>& displacement(size_t vertexIndex) const{
            return mechMap_.isIdentity() ? displacement_[vertexIndex]
                                         : flowDisplacement_[vertexIndex];
        }
        const double& mechPotentialForce(unsigned globalDofIdx) const
        {
            return mechPotentialForce_[globalDofIdx];
        }
        const double& pressure(unsigned globalDofIdx) const
        {
            return pressure_[globalDofIdx];
        }
        const double& mechPotentialTempForce(unsigned globalDofIdx) const
        {
            return mechPotentialTempForce_[globalDofIdx];
        }
        const double& mechPotentialPressForce(unsigned globalDofIdx) const
        {
            return mechPotentialPressForce_[globalDofIdx];
        }

        // Output forms of the potential forces, without the poroelastic factor
        // fac = (1-nu)/(1-2nu) that the stored (equation) quantities carry.
        double mechPotentialPressForceOutput(unsigned globalDofIdx) const
        {
            const double pratio = simulator_.problem().pRatio(globalDofIdx);
            const double fac = (1 - pratio) / (1 - 2 * pratio);
            return mechPotentialPressForce_[globalDofIdx] / fac;
        }

        double mechPotentialTempForceOutput(unsigned globalDofIdx) const
        {
            const double pratio = simulator_.problem().pRatio(globalDofIdx);
            const double fac = (1 - pratio) / (1 - 2 * pratio);
            return mechPotentialTempForce_[globalDofIdx] / fac;
        }
        // The bool is the opm-simulators output call's include_fracture flag; ignored.
        const Dune::FieldVector<double,3> disp(size_t globalIdx, bool = false) const{
            if (!this->hasMech(globalIdx)) {
                return Dune::FieldVector<double,3>(0.0);
            }
            return celldisplacement_[this->mechIdx(globalIdx)];
        }

        const SymTensor delstress(size_t globalIdx) const{
	  Dune::FieldVector<double,6> delStress = this->effstress(globalIdx);
	  double effPress = this->mechPotentialForce(globalIdx);
	  for(int i=0; i < 3; ++i){
	    delStress[i] += effPress;

	  }
	  return delStress;
        }

        const SymTensor linstress(size_t globalIdx) const{
	  if (!this->hasMech(globalIdx)) {
	    return SymTensor(0.0);
	  }
	  return mechState_.stress[this->mechIdx(globalIdx)];
        }

        const SymTensor outputstress(size_t globalIdx) const{
	        return outputstress_[globalIdx];
        }

        const SymTensor effstress(size_t globalIdx) const{
	  // make stress in with positive with compression
	         if (!this->hasMech(globalIdx)) {
	           return SymTensor(0.0);
	         }
	         return -1.0*mechState_.stress[this->mechIdx(globalIdx)];
        }

        const SymTensor strain(size_t globalIdx, bool = false) const{
            if (!this->hasMech(globalIdx)) {
                return SymTensor(0.0);
            }
            return mechState_.strain[this->mechIdx(globalIdx)];
        }

        const SymTensor stress(size_t globalIdx) const{
            Dune::FieldVector<double,6> effStress = this->effstress(globalIdx);
            effStress += simulator_.problem().initStress(globalIdx);
            double effPress = this->mechPotentialForce(globalIdx);
            for(int i=0; i < 3; ++i){
                effStress[i] += effPress;
            }
            return effStress;
         }

         const SymTensor fractureStress(size_t globalIdx) const{
            // need to know effect of absolute pressure
	   const auto& iq = simulator_.model().intensiveQuantities(globalIdx,0);
	   const auto& fs = iq.fluidState();
       const auto& press = Toolbox::value(fs.pressure(referencePhaseIdx(fs)));
            Dune::FieldVector<double,6> fracStress = this->stress(globalIdx);
            for(int i=0; i < 3; ++i){
                fracStress[i] -= press;
            }
            return fracStress;
         }

        // NB used in output should be eliminated

        double pressureDiff(unsigned dofIx) const{
            return mechPotentialForce_[dofIx];
        }

        // A flow vertex that is a mechanics vertex copies it; any other is
        // trilinear in the corners of its flow cell's mechanics cell.
        void buildFlowVertexSource_()
        {
            const auto& mgv = this->mechGrid().leafGridView();
            const auto key = [](const auto& p) {
                return std::array<long long,3>{std::llround(p[0]*1e6), std::llround(p[1]*1e6),
                                               std::llround(p[2]*1e6)};
            };
            std::map<std::array<long long,3>, int> mechVertexAt;
            for (const auto& v : Dune::vertices(mgv)) {
                mechVertexAt.emplace(key(v.geometry().center()), mgv.indexSet().index(v));
            }
            using MechElement = typename std::decay_t<decltype(mgv)>::template Codim<0>::Entity;
            std::vector<MechElement> mechCells;
            mechCells.reserve(mgv.size(0));
            for (const auto& e : Dune::elements(mgv)) {
                mechCells.push_back(e);
            }

            const auto& fgv = simulator_.gridView();
            flowVertexSource_.assign(fgv.size(3), FlowVertexSource{});
            int copied = 0, interpolated = 0;
            double misfit = 0.0;
            for (const auto& cell : Dune::elements(fgv)) {
                const int m = mechMap_.mechCell(fgv.indexSet().index(cell));
                for (const auto& v : Dune::subEntities(cell, Dune::Codim<3>{})) {
                    auto& src = flowVertexSource_[fgv.indexSet().index(v)];
                    if (src.count > 0) {
                        continue;
                    }
                    const auto x = v.geometry().center();
                    if (const auto it = mechVertexAt.find(key(x)); it != mechVertexAt.end()) {
                        src.vertex[0] = it->second;
                        src.weight[0] = 1.0;
                        src.count = 1;
                        ++copied;
                        continue;
                    }
                    if (m < 0) {
                        continue;
                    }
                    const auto& mcell = mechCells[m];
                    const auto u = localInCell_(mcell.geometry(), x);
                    for (int c = 0; c < 8; ++c) {
                        src.vertex[c] = mgv.indexSet().index(mcell.template subEntity<3>(c));
                        src.weight[c] = ((c & 1) ? u[0] : 1.0 - u[0])
                                      * ((c & 2) ? u[1] : 1.0 - u[1])
                                      * ((c & 4) ? u[2] : 1.0 - u[2]);
                    }
                    src.count = 8;
                    ++interpolated;
                    Dune::FieldVector<double,3> y(0.0);
                    for (int c = 0; c < 8; ++c) {
                        auto p = mcell.geometry().corner(c);
                        p *= src.weight[c];
                        y += p;
                    }
                    y -= x;
                    misfit = std::max(misfit, y.two_norm());
                }
            }
            std::ostringstream os;
            os << "Mechanics grid: flow vertices " << copied << " copied, " << interpolated
               << " interpolated (worst position misfit " << misfit << " m)";
            OpmLog::info(os.str());
        }

        // Newton on the trilinear map, capped: collapsed cells need not converge.
        template <class Geometry>
        static Dune::FieldVector<double,3>
        localInCell_(const Geometry& geo, const Dune::FieldVector<double,3>& x)
        {
            // Newton from the centre can stall in strongly distorted cells; the
            // corners are tried too and the best fit is kept.
            Dune::FieldVector<double,3> best(0.5);
            double bestMisfit = std::numeric_limits<double>::max();
            for (int start = -1; start < 8 && bestMisfit > 1e-8; ++start) {
                Dune::FieldVector<double,3> u(0.5);
                if (start >= 0) {
                    for (int a = 0; a < 3; ++a) {
                        u[a] = (start >> a) & 1;
                    }
                }
                for (int it = 0; it < 30; ++it) {
                    auto r = geo.global(u);
                    r -= x;
                    const auto jt = geo.jacobianTransposed(u);
                    Dune::FieldMatrix<double,3,3> j;
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            j[a][b] = jt[b][a];
                        }
                    }
                    Dune::FieldVector<double,3> du;
                    try {
                        j.solve(du, r);
                    } catch (const Dune::FMatrixError&) {
                        break;
                    }
                    u -= du;
                    if (du.two_norm() < 1e-10) {
                        break;
                    }
                }
                for (auto& c : u) {
                    c = std::clamp(c, 0.0, 1.0);
                }
                auto r = geo.global(u);
                r -= x;
                if (r.two_norm() < bestMisfit) {
                    bestMisfit = r.two_norm();
                    best = u;
                }
            }
            return best;
        }

        void prolongDisplacement_()
        {
            if (flowVertexSource_.empty()) {
                this->buildFlowVertexSource_();
            }
            flowDisplacement_.resize(flowVertexSource_.size());
            for (std::size_t i = 0; i < flowVertexSource_.size(); ++i) {
                const auto& src = flowVertexSource_[i];
                flowDisplacement_[i] = 0.0;
                for (int c = 0; c < src.count; ++c) {
                    auto d = displacement_[src.vertex[c]];
                    d *= src.weight[c];
                    flowDisplacement_[i] += d;
                }
            }
        }

        void makeDisplacement(const Opm::Elasticity::Vector& field) {
            // make displacement on all nodes used for output to vtk
            const auto& grid = this->mechGrid();
            const auto& gv = grid.leafGridView();
            int dim = 3;
            for (const auto& vertex : Dune::vertices(gv)){
                auto index = gv.indexSet().index(vertex);
                for(int k=0; k < dim; ++k){
                    displacement_[index][k] = field[index*dim+k];
                }
            }
            if (!mechMap_.isIdentity()) {
                this->prolongDisplacement_();
            }
            for (const auto& cell: elements(gv)){
                auto cellindex = gv.indexSet().index(cell);
                celldisplacement_[cellindex] = 0.0;
                const auto& vertices = Dune::subEntities(cell, Dune::Codim<Grid::dimension>{});
                for (const auto& vertex : vertices){
                    auto nodeidex = gv.indexSet().index(vertex);
                    for(int k=0; k < dim; ++k){
                        celldisplacement_[cellindex][k] += displacement_[nodeidex][k] ;
                    }
                }
                celldisplacement_[cellindex] /= vertices.size();
            }
        }

        std::string finalTimingSummary() const
        {
            std::ostringstream os;
            auto report_time = [&os](const char* label, const double value)
            {
                os << std::left << std::setw(28) << label
                   << std::right << std::setw(7) << std::fixed << std::setprecision(2) << value
                   << " s\n";
            };
            auto report_count = [&os](const char* label, const int value)
            {
                os << std::left << std::setw(28) << label
                   << std::right << std::setw(7) << value << "\n";
            };

            report_time("  Mechanical solve time:", total_mechanical_solve_time_seconds_);
            report_count("Overall Mechanical Solves:", num_mech_solves_);
            report_count("Overall Mech Lin Iters:", total_mech_linear_iterations_);
            report_count("Mech Solves Unconverged:", num_mech_unconverged_);

            if (simulator_.problem().hasFractures()) {
                const auto total_stats = fracHost_.fractureModelPtr() ? fracHost_.fractureModelPtr()->totalSolveStats() : FractureSolveStats{};
                report_time("  Fracture solve time:", total_stats.solve_time_seconds);
                report_count("Overall Fracture Solves:", total_stats.fractures_solved);
                report_count("Overall Fracture Nl Iters:", total_stats.nonlinear_iterations);
                report_count("Overall Fracture Lin Solves:", total_stats.linear_solves);
                report_count("Overall Fracture Lin Iters:", total_stats.linear_iterations);
            }

            return os.str();
        }

        const MechOptions& mechOptions(){
            if (!mechOptions_) {
                mechOptions_ = MechOptions::read(simulator_.problem().getFractureParam());
            }
            return *mechOptions_;
        }

    private:
        std::optional<double> assembledGravity_; // body force of the assembled matrix
        int num_mech_solves_{0};
        int num_mech_unconverged_{0};
        bool solved_in_step_{false};
        int total_mech_linear_iterations_{0};
        bool write_system_{false};
        std::optional<MechOptions> mechOptions_;
        double last_mechanical_solve_time_seconds_{0.0};
        double total_mechanical_solve_time_seconds_{0.0};
        Simulator& simulator_;

        Dune::BlockVector<Dune::FieldVector<double,1>> pressure_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForce_;
        // Only filled when the mechanics grid is a coarsening of the flow grid.
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForceMech_;
        MechFlowMap mechMap_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForceFracture_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialTempForce_;
        Dune::BlockVector<Dune::FieldVector<double,3> > celldisplacement_;
        std::vector<double> ymoduleMech_, pratioMech_;
        Dune::BlockVector<Dune::FieldVector<double,3> > displacement_;
        // Flow vertices, when the mechanics grid is coarser: each takes the
        // weighted displacement of up to eight mechanics vertices.
        struct FlowVertexSource
        {
            std::array<int,8> vertex{};
            std::array<double,8> weight{};
            int count{0};
        };
        std::vector<FlowVertexSource> flowVertexSource_;
        Dune::BlockVector<Dune::FieldVector<double,3> > flowDisplacement_;
        Dune::BlockVector<Dune::FieldVector<double,6> > outputstress_;// used in to avoid trouble with initialstress_
        Opm::Elasticity::MechState mechState_; // stress is the linear (mechanics-cell) stress
        Opm::Elasticity::Vector mechRhs_;
        std::unique_ptr<Opm::Elasticity::VemDiscretisation<Grid>> mechDisc_;
        std::unique_ptr<Opm::Elasticity::MechLinearSolver<Grid>> mechLinSolver_;
        const MechGridContext* mechCtx_{nullptr};
        FractureMechHost<TypeTag> fracHost_{simulator_};
    };
}



#endif
