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

// Backend-agnostic fracture <-> flow coupling apparatus, extracted verbatim
// from NonlinearSystemBlackOilReservoirGeoMech so that nonlinear systems
// built on other mechanics backends (e.g. TPSA) can drive the same fracture
// outer loop.  CRTP: the Derived nonlinear system declares this class a
// friend and provides simulator_, flowNonlinearIteration(timer, solver),
// onConnectionsUpdated(), wellFracturePicard(timer) and
// wellFractureJointSolve(timer).

#include <cmath>
#include <limits>
#include <sstream>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

namespace Opm
{
    template <class TypeTag, class Derived>
    class FractureCouplingOuterLoop
    {
        public:
        using EqVector = GetPropType<TypeTag, Properties::EqVector>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

        protected:
        Derived& derived()
        { return static_cast<Derived&>(*this); }

        const Derived& derived() const
        { return static_cast<const Derived&>(*this); }

        struct CouplingChangeMetrics
        {
            bool structure_changed{false};
            bool ctf_changed{false};
            bool perf_pressure_changed{false};
            double max_ctf_change_abs{0.0};
            double max_ctf_change_rel{0.0};
            double max_perf_change_abs{0.0};
            double max_perf_change_rel{0.0};
            double composite_norm{0.0};

            bool changed() const
            {
                return structure_changed || ctf_changed || perf_pressure_changed;
            }
        };

        enum class KeepIterating { No, CouplingChanged, FractureUnconverged, EmbeddedRestructured };

        // Every option fractureOuterBlock reads, read once per call.
        struct CouplingConfig
        {
            bool legacyChangeLogic{false};
            bool legacyParentSetup{true};
            bool requireConvergedFractureForWi{true};
            bool addConnections{false};
            bool valueOnlyWiUpdate{false};
            bool residualMode{false};
            bool rebindAfterGrowth{false};
            int rebindAfterGrowthRounds{2};
            double couplingTolerance{1e-3};
            bool hysteresis{false};
            double hysteresisOn{1.0};
            double hysteresisOff{0.3};
            double hysteresisEmergency{10.0};
            int hysteresisConfirm{1};
            int hysteresisCooldown{0};

            static CouplingConfig read(const PropertyTree& prm)
            {
                CouplingConfig c;
                c.legacyChangeLogic = prm.get<bool>("fractureparam.solver.legacy_coupling_change_logic", false);
                c.legacyParentSetup = prm.get<bool>("fractureparam.solver.legacy_parent_setup_iteration", true);
                c.requireConvergedFractureForWi =
                    prm.get<bool>("fractureparam.require_converged_fracture_for_wi_update", true);
                c.addConnections = prm.get<bool>("fractureparam.addconnections");
                c.valueOnlyWiUpdate = prm.get<bool>("fractureparam.solver.value_only_wi_update", false);
                c.residualMode = prm.get<std::string>("fractureparam.solver.coupling_convergence_mode",
                                                      "threshold") == "residual";
                c.rebindAfterGrowth = prm.get<bool>("fractureparam.solver.rebind_after_growth", false);
                c.rebindAfterGrowthRounds = prm.get<int>("fractureparam.solver.rebind_after_growth_rounds", 2);
                c.couplingTolerance = prm.get<double>("fractureparam.solver.coupling_tolerance", 1e-3);
                c.hysteresis = prm.get<bool>("fractureparam.topology_hysteresis.enable", false);
                c.hysteresisOn = prm.get<double>("fractureparam.topology_hysteresis.on_threshold", 1.0);
                c.hysteresisOff = prm.get<double>("fractureparam.topology_hysteresis.off_threshold", 0.3);
                c.hysteresisConfirm = prm.get<int>("fractureparam.topology_hysteresis.confirm_iterations", 1);
                c.hysteresisCooldown = prm.get<int>("fractureparam.topology_hysteresis.cooldown_iterations", 0);
                c.hysteresisEmergency =
                    prm.get<double>("fractureparam.topology_hysteresis.emergency_threshold", 10.0);
                return c;
            }
        };

        bool shouldUpdateConnections(const CouplingChangeMetrics& metrics, const CouplingConfig& cfg)
        {
            if (cfg.legacyChangeLogic || !cfg.hysteresis) {
                return true;
            }

            const double on_threshold = cfg.hysteresisOn;
            const double off_threshold = cfg.hysteresisOff;
            const int confirm_iterations = cfg.hysteresisConfirm;
            const int cooldown_iterations = cfg.hysteresisCooldown;
            const double emergency_threshold = cfg.hysteresisEmergency;

            if (topology_cooldown_counter_ > 0 && metrics.composite_norm < emergency_threshold && !metrics.structure_changed) {
                --topology_cooldown_counter_;
                return false;
            }

            const bool trigger_on = metrics.structure_changed || metrics.composite_norm >= on_threshold;
            const bool trigger_off = (!metrics.structure_changed) && metrics.composite_norm <= off_threshold;

            if (trigger_on) {
                ++topology_pending_counter_;
            }
            else if (trigger_off) {
                topology_pending_counter_ = 0;
            }

            if (metrics.composite_norm >= emergency_threshold || topology_pending_counter_ >= confirm_iterations) {
                topology_pending_counter_ = 0;
                topology_cooldown_counter_ = std::max(0, cooldown_iterations);
                return true;
            }

            return false;
        }

        mutable int topology_pending_counter_{0};
        mutable int topology_cooldown_counter_{0};
        // Mech-skip (opt-in performance lever): composite coupling-change norm from
        // the previous outer iteration's fracture solve. Mechanics depends on the
        // fracture only indirectly via the change in well flow, so when this norm is
        // small the mechanics solution is stable and the mech solve can be skipped.
        // Defaults large so nothing is skipped until a positive threshold is set.
        mutable double last_coupling_composite_norm_{std::numeric_limits<double>::max()};
        // Diagnostics: mechanics solves actually performed vs skipped in the outer loop.
        mutable int mech_solves_performed_{0};
        mutable int mech_solves_skipped_{0};
        // Mechanics solves performed within the current timestep; used by the
        // opt-in per-step cap (solver.max_mech_solves_per_step). Reset at the
        // first outer iteration of each step.
        mutable int mech_solves_this_step_{0};
        // Growth sub-iterations spent within the current timestep
        // (solver.max_growth_iterations is a per-timestep budget: contact
        // chatter must not re-arm the loop on every outer iteration).
        mutable int growth_rounds_this_step_{0};
        // The legacy setup iteration left an unconverged Newton update in the flow state.
        mutable bool flow_state_modified_by_setup_{false};
        // Coupling accepted; flow must re-converge before the step may end.
        mutable bool reconverge_flow_pending_{false};

        struct StorageCacheBackup
        {
            bool enabled{false};
            std::vector<std::vector<EqVector>> values;
            std::vector<std::vector<unsigned char>> valid;
        };

        StorageCacheBackup backupStorageCache(const unsigned num_time_indices = 2) const
        {
            StorageCacheBackup backup;
            const auto& flow_model = derived().simulator_.model();
            if (!flow_model.enableStorageCache()) {
                return backup;
            }

            backup.enabled = true;
            // Over every degree of freedom, not just the grid ones: restoring invalidates
            // the whole cache before putting the saved entries back, so anything left out
            // here is not preserved but destroyed.  That matters as soon as the fracture
            // itself has degrees of freedom in this system.
            const unsigned num_dof = flow_model.numTotalDof();
            backup.values.assign(num_time_indices, std::vector<EqVector>(num_dof, EqVector(0.0)));
            backup.valid.assign(num_time_indices, std::vector<unsigned char>(num_dof, 0));

            for (unsigned time_idx = 0; time_idx < num_time_indices; ++time_idx) {
                for (unsigned dof_idx = 0; dof_idx < num_dof; ++dof_idx) {
                    if (flow_model.storageCacheIsUpToDate(dof_idx, time_idx)) {
                        backup.values[time_idx][dof_idx] = flow_model.cachedStorage(dof_idx, time_idx);
                        backup.valid[time_idx][dof_idx] = 1;
                    }
                }
            }

            return backup;
        }

        void restoreStorageCache(const StorageCacheBackup& backup) const
        {
            if (!backup.enabled) {
                return;
            }

            const auto& flow_model = derived().simulator_.model();
            for (unsigned time_idx = 0; time_idx < backup.values.size(); ++time_idx) {
                flow_model.invalidateStorageCache(time_idx);
                for (unsigned dof_idx = 0; dof_idx < backup.valid[time_idx].size(); ++dof_idx) {
                    if (backup.valid[time_idx][dof_idx]) {
                        flow_model.updateCachedStorage(dof_idx, time_idx, backup.values[time_idx][dof_idx]);
                    }
                }
            }
        }

        SolutionVector backupSolution0() const
        {
            return derived().simulator_.model().solution(0);
        }

        void restoreSolution0(const SolutionVector& solution_backup)
        {
            auto& flow_model = derived().simulator_.model();
            flow_model.solution(0) = solution_backup;
            flow_model.invalidateAndUpdateIntensiveQuantities(0);
        }

        // Run one parent nonlinear iteration as a pure setup pass after the well
        // structure changed mid-step: the parent must see iteration 0 to redo its
        // once-per-step initialization, so the NewtonIterationContext is reset for
        // the call and restored afterwards.
        template <class NonlinearSolverType>
        SimulatorReportSingle runParentSetupIteration(
            const SimulatorTimerInterface& timer,
            NonlinearSolverType& nonlinear_solver)
        {
            const SetupIterationContextGuard guard{derived().simulator_.problem()};
            return derived().flowNonlinearIteration(timer, nonlinear_solver);
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle runParentFirstIterationPreservingState(
            const SimulatorTimerInterface& timer,
            NonlinearSolverType& nonlinear_solver)
        {
            OpmLog::debug("Running parent first iteration to set up the system and caches, while preserving state for geomechanics solve");
            // this is a hack to call all setup need for the system without changing states which matter for the calculations
            const auto storage_cache_backup = this->backupStorageCache();
            const auto solution_backup = this->backupSolution0();
            // The setup iteration also solves the wells, mutating the well/group
            // state.  solution(0) + storage cache alone do not cover that: for a
            // MultisegmentWell the segment state is left inconsistent with the
            // restored solution and the reported rate collapses to zero while the
            // perforations keep injecting.  Snapshot and restore the WGState too so
            // the setup iteration is genuinely state-neutral.
            auto wgstate_backup = derived().simulator_.problem().wellModel().snapshotWGState();
            // Restore on every exit path.  A throwing setup iteration (a failed
            // linear solve) would otherwise leave its half-applied iterate and an
            // invalidated storage cache behind as the step's state.
            auto restore_state = [&] {
                derived().simulator_.problem().wellModel().restoreWGState(std::move(wgstate_backup));
                this->restoreSolution0(solution_backup);
                this->restoreStorageCache(storage_cache_backup);
            };
            SimulatorReportSingle report;
            try {
                report = this->runParentSetupIteration(timer, nonlinear_solver);
            } catch (...) {
                restore_state();
                throw;
            }
            restore_state();
            OpmLog::debug("Finished parent first iteration");
            return report;
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle runParentFirstIterationLegacy(
            const SimulatorTimerInterface& timer,
            NonlinearSolverType& nonlinear_solver)
        {
            OpmLog::debug("Running parent first iteration in legacy mode (no state restore)");
            flow_state_modified_by_setup_ = true;
            auto report = this->runParentSetupIteration(timer, nonlinear_solver);
            OpmLog::debug("Finished parent first iteration");
            return report;
        }

        public:

        CouplingChangeMetrics evaluateFractureCouplingChange(
            const std::vector<std::vector<RuntimePerforation>>& allWellIndices) const
        {
            const auto& prm = derived().simulator_.problem().getFractureParam();
            const double ctf_threshold_abs = prm.get("solver.ctf_change_threshold", 1e-15);
            const double ctf_threshold_rel = prm.get("solver.ctf_change_threshold_rel", 0.0);
            const double ctf_rel_scale = prm.get("solver.ctf_relative_scale", 1.0);
            const double perf_threshold_abs = prm.get("solver.perf_pressure_change_threshold", 1e30);
            const double perf_threshold_rel = prm.get("solver.perf_pressure_change_threshold_rel", 0.0);
            const double perf_rel_scale = prm.get("solver.perf_pressure_relative_scale", 1.0);
            const int verbosity = prm.get("solver.verbosity", 1);

            CouplingChangeMetrics metrics;
            const auto allWellIndices_new = derived().simulator_.problem().getAllExtraWellIndices();
            if (allWellIndices.size() != allWellIndices_new.size()) {
                metrics.structure_changed = true;
            }

            const size_t num_wells = std::min(allWellIndices.size(), allWellIndices_new.size());
            for (size_t k = 0; k < num_wells; ++k) {
                const auto& wellIndices = allWellIndices[k];
                const auto& wellIndices_new = allWellIndices_new[k];
                if (wellIndices.size() != wellIndices_new.size()) {
                    metrics.structure_changed = true;
                }

                const size_t num_perfs = std::min(wellIndices.size(), wellIndices_new.size());
                for (size_t i = 0; i < num_perfs; ++i) {
                    if (wellIndices[i].cell != wellIndices_new[i].cell) {
                        metrics.structure_changed = true;
                    }

                    const double ctf_old = wellIndices[i].ctf;
                    const double ctf_new = wellIndices_new[i].ctf;
                    const double ctf_change_abs = std::abs(ctf_old - ctf_new);
                    const double ctf_scale = std::max({std::abs(ctf_old), std::abs(ctf_new), ctf_rel_scale});
                    const double ctf_change_rel = ctf_change_abs / ctf_scale;

                    metrics.max_ctf_change_abs = std::max(metrics.max_ctf_change_abs, ctf_change_abs);
                    metrics.max_ctf_change_rel = std::max(metrics.max_ctf_change_rel, ctf_change_rel);

                    const double perf_old = wellIndices[i].pressure;
                    const double perf_new = wellIndices_new[i].pressure;
                    const double perf_change_abs = std::abs(perf_old - perf_new);
                    const double perf_scale = std::max({std::abs(perf_old), std::abs(perf_new), perf_rel_scale});
                    const double perf_change_rel = perf_change_abs / perf_scale;

                    metrics.max_perf_change_abs = std::max(metrics.max_perf_change_abs, perf_change_abs);
                    metrics.max_perf_change_rel = std::max(metrics.max_perf_change_rel, perf_change_rel);
                }
            }

            if (ctf_threshold_rel > 0.0) {
                metrics.ctf_changed = metrics.max_ctf_change_abs > ctf_threshold_abs
                                   && metrics.max_ctf_change_rel > ctf_threshold_rel;
            }
            else {
                metrics.ctf_changed = metrics.max_ctf_change_abs > ctf_threshold_abs;
            }

            if (perf_threshold_rel > 0.0) {
                metrics.perf_pressure_changed = metrics.max_perf_change_abs > perf_threshold_abs
                                             && metrics.max_perf_change_rel > perf_threshold_rel;
            }
            else {
                metrics.perf_pressure_changed = metrics.max_perf_change_abs > perf_threshold_abs;
            }

            const double ctf_abs_norm = (ctf_threshold_abs > 0.0)
                ? metrics.max_ctf_change_abs / ctf_threshold_abs
                : 0.0;
            const double ctf_rel_norm = (ctf_threshold_rel > 0.0)
                ? metrics.max_ctf_change_rel / ctf_threshold_rel
                : 0.0;
            const double perf_abs_norm = (perf_threshold_abs > 0.0)
                ? metrics.max_perf_change_abs / perf_threshold_abs
                : 0.0;
            const double perf_rel_norm = (perf_threshold_rel > 0.0)
                ? metrics.max_perf_change_rel / perf_threshold_rel
                : 0.0;
            metrics.composite_norm = std::max({ctf_abs_norm,
                                               ctf_rel_norm,
                                               perf_abs_norm,
                                               perf_rel_norm,
                                               metrics.structure_changed ? 1.0 : 0.0});

            if (verbosity > 0) {
                std::stringstream os;
                os << "Fracture coupling change: structure_changed=" << (metrics.structure_changed ? "true" : "false")
                   << ", ctf_changed=" << (metrics.ctf_changed ? "true" : "false")
                   << ", perf_changed=" << (metrics.perf_pressure_changed ? "true" : "false")
                   << ", max_ctf_change_abs=" << metrics.max_ctf_change_abs
                   << ", max_ctf_change_rel=" << metrics.max_ctf_change_rel
                   << ", max_perf_change_abs=" << metrics.max_perf_change_abs
                   << ", max_perf_change_rel=" << metrics.max_perf_change_rel
                   << ", composite_norm=" << metrics.composite_norm
                   << ", abs_threshold=" << ctf_threshold_abs;
                if (ctf_threshold_rel > 0.0) {
                    os << ", rel_threshold=" << ctf_threshold_rel;
                }
                if (perf_threshold_abs < 1e29) {
                    os << ", perf_abs_threshold=" << perf_threshold_abs;
                }
                if (perf_threshold_rel > 0.0) {
                    os << ", perf_rel_threshold=" << perf_threshold_rel;
                }
                OpmLog::info(os.str());
            }

            return metrics;
        }

        bool fractureChanged(
            const std::vector<std::vector<RuntimePerforation>>& allWellIndices) const
        {
            return evaluateFractureCouplingChange(allWellIndices).changed();
        }

        protected:
        // Why the last fractureOuterBlock kept the outer loop active.
        KeepIterating lastKeepIterating_{KeepIterating::No};

        // The fracture part of one outer nonlinear iteration: solve the
        // fractures against the current mechanics state, couple the result to
        // the flow and keep the outer loop active until the coupling has
        // converged. The caller has already decided that fractures should be
        // solved this iteration (do_fracture && problem.hasFractures()).
        template <class NonlinearSolverType>
        void fractureOuterBlock(SimulatorReportSingle& report,
                                const SimulatorTimerInterface& timer,
                                NonlinearSolverType& nonlinear_solver)
        {
            auto& problem = derived().simulator_.problem();
            OpmLog::info("Solve Fractures:");
            flow_state_modified_by_setup_ = false;
            const auto wiBefore = problem.getAllExtraWellIndices();
            this->solveFracturesThisIteration(timer);
            // A FRAC deck with no WSEED has no fracture: nothing to couple.
            if (!problem.fractureHost().fractureModelActive()) {
                return;
            }
            const auto cfg = CouplingConfig::read(problem.getGeoMechParam());
            const bool fractureConvergedLocal =
                problem.fractureHost().fractureModel().lastSolveStats().converged;
            auto& comm = derived().simulator_.gridView().comm();
            const bool fractureConverged = comm.min(fractureConvergedLocal);
            // Computed on embedded runs too: the mech skip reads its composite norm.
            const auto metrics = this->globalCouplingMetrics(wiBefore, cfg);

            const KeepIterating why = problem.fractureFlowIsEmbedded()
                ? this->coupleEmbedded(cfg, fractureConverged)
                : this->coupleThroughWellIndices(metrics, cfg, fractureConverged,
                                                 fractureConvergedLocal, timer, nonlinear_solver);
            lastKeepIterating_ = why;
            if (why != KeepIterating::No) {
                report.converged = false;
            }
            comm.barrier();
        }

        private:
        void solveFracturesThisIteration(const SimulatorTimerInterface& timer)
        {
            auto& problem = derived().simulator_.problem();
            if (!problem.fractureFlowIsEmbedded()) {
                problem.fractureHost().solveFractures();
                return;
            }
            // Per-iteration embedded flags; the Picard and joint solves bind and set them.
            problem.clearEmbeddedRestructured();
            const PropertyTree& prm = problem.getGeoMechParam();
            // The fracture model is built by the first solveFractures().
            const bool pressureFromFlow = prm.get<bool>("fractureparam.solver.pressure_from_flow", false)
                && problem.fractureHost().fractureModelActive();
            if (!pressureFromFlow) {
                // embedded: plain solve
                problem.fractureHost().solveFractures();
            } else if (prm.get<bool>("fractureparam.solver.well_joint_solve", false)) {
                // embedded: joint Newton
                derived().wellFractureJointSolve(timer);
            } else {
                // embedded: well-fracture Picard
                derived().wellFracturePicard(timer);
            }
        }

        // Change of the upscaled well indices over this iteration's fracture
        // solve, reduced over all ranks; also feeds the mech skip.
        CouplingChangeMetrics globalCouplingMetrics(
            const std::vector<std::vector<RuntimePerforation>>& wiBefore,
            const CouplingConfig& cfg)
        {
            auto& comm = derived().simulator_.gridView().comm();
            const auto local = this->evaluateFractureCouplingChange(wiBefore);
            CouplingChangeMetrics metrics = local;
            metrics.structure_changed = comm.max(local.structure_changed);
            metrics.ctf_changed = comm.max(local.ctf_changed);
            metrics.perf_pressure_changed = comm.max(local.perf_pressure_changed);
            metrics.composite_norm = comm.max(local.composite_norm);
            if (cfg.residualMode) {
                metrics.max_ctf_change_rel = comm.max(local.max_ctf_change_rel);
                metrics.max_perf_change_rel = comm.max(local.max_perf_change_rel);
            }
            last_coupling_composite_norm_ = metrics.composite_norm;
            return metrics;
        }

        // The fracture flows through degrees of freedom of its own: re-describe
        // them to the reservoir. No schedule rebuild and no upscaled well index;
        // within the step the well's perforations of the fracture keep the
        // indices they got at the step boundary.
        KeepIterating coupleEmbedded(const CouplingConfig& cfg, const bool fractureConverged)
        {
            auto& problem = derived().simulator_.problem();
            auto& comm = derived().simulator_.gridView().comm();
            // A grown fracture needs a restructure, which a value-only refresh cannot express.
            problem.bindFractureAuxCells(
                /*allowTopologyChange=*/cfg.rebindAfterGrowth && fractureConverged,
                /*requireStableLayout=*/true);

            // The residual is the change in what the flow was fed (total conductance and
            // pore volume of the binding); contact chatter cancels in these aggregates.
            const double embedded_change = problem.embeddedCouplingChange();
            if (!problem.embeddedRestructuredMidStep() && problem.embeddedLayoutPending()
                && problem.takeEmbeddedPendingRound(cfg.rebindAfterGrowthRounds)) {
                // Confirm the grown layout by re-solving on it without growing.
                auto& fractureModel = problem.fractureHost().fractureModel();
                auto setSuppressed = [&](const bool on) {
                    for (auto& wellFractures : fractureModel.wellFractures()) {
                        for (auto& fracture : wellFractures) {
                            fracture.suppressPropagation(on);
                        }
                    }
                };
                setSuppressed(true);
                problem.fractureHost().solveFractures();
                setSuppressed(false);
                const bool confirmed = comm.min(fractureModel.lastSolveStats().converged);
                problem.bindFractureAuxCells(/*allowTopologyChange=*/confirmed,
                                             /*requireStableLayout=*/true);
            }

            if (problem.embeddedRestructuredMidStep()) {
                // The flow has not been solved on the grown fracture yet.
                OpmLog::info("Keeping outer loop active: fracture binding restructured "
                             "after growth (rebind_after_growth)");
                problem.holdEmbeddedOuterLoop();
                return KeepIterating::EmbeddedRestructured;
            }
            if (!fractureConverged) {
                OpmLog::info("Keeping outer loop active: fracture solve did not converge");
                return KeepIterating::FractureUnconverged;
            }
            if (embedded_change > cfg.couplingTolerance) {
                OpmLog::info("Keeping outer loop active: embedded coupling change "
                             + std::to_string(embedded_change) + " > "
                             + std::to_string(cfg.couplingTolerance));
                return KeepIterating::CouplingChanged;
            }
            return KeepIterating::No;
        }

        template <class NonlinearSolverType>
        KeepIterating coupleThroughWellIndices(const CouplingChangeMetrics& metrics,
                                               const CouplingConfig& cfg,
                                               const bool fractureConverged,
                                               const bool fractureConvergedLocal,
                                               const SimulatorTimerInterface& timer,
                                               NonlinearSolverType& nonlinear_solver)
        {
            const bool may_update_connections = cfg.addConnections
                && (cfg.legacyChangeLogic || !cfg.requireConvergedFractureForWi || fractureConverged);
            if (may_update_connections && this->shouldUpdateConnections(metrics, cfg)) {
                this->updateWellConnections(metrics, cfg, timer, nonlinear_solver);
            }
            else if (may_update_connections && !cfg.legacyChangeLogic) {
                std::stringstream hos;
                hos << "Deferring fracture-driven connection update due to topology hysteresis"
                    << ", composite_norm=" << metrics.composite_norm
                    << ", pending=" << topology_pending_counter_
                    << ", cooldown=" << topology_cooldown_counter_;
                OpmLog::info(hos.str());
            }
            else if (cfg.addConnections && cfg.requireConvergedFractureForWi && !fractureConvergedLocal) {
                OpmLog::info("Skipping fracture-driven connection update because fracture solve did not converge");
            }

            if (cfg.residualMode) {
                // Largest relative change of any fracture CTF or perforation pressure;
                // a structure change or an unconverged fracture solve always keeps the loop.
                const double outer_residual = std::max(metrics.max_ctf_change_rel,
                                                       metrics.max_perf_change_rel);
                const bool outer_converged = fractureConverged
                    && !metrics.structure_changed
                    && outer_residual <= cfg.couplingTolerance;
                std::stringstream ros;
                ros << "Outer coupling residual: " << outer_residual
                    << " (tol " << cfg.couplingTolerance << ")"
                    << ", structure_changed=" << (metrics.structure_changed ? "true" : "false")
                    << ", fracture_converged=" << (fractureConverged ? "true" : "false")
                    << ", outer_converged=" << (outer_converged ? "true" : "false");
                OpmLog::info(ros.str());
                if (outer_converged) {
                    return KeepIterating::No;
                }
                return fractureConverged ? KeepIterating::CouplingChanged
                                         : KeepIterating::FractureUnconverged;
            }

            const bool fracture_changed = metrics.changed();
            OpmLog::debug("Fracture changed: " + std::to_string(fracture_changed));
            // TODO check convergence properly
            const bool unconverged = !cfg.legacyChangeLogic && !fractureConverged;
            if (!unconverged && !fracture_changed) {
                return KeepIterating::No;
            }
            if (!fractureConverged) {
                OpmLog::info("Keeping outer nonlinear loop active because fracture solve did not converge");
            }
            return unconverged ? KeepIterating::FractureUnconverged : KeepIterating::CouplingChanged;
        }

        // Puts the new fracture well indices into the well model mid-step.
        template <class NonlinearSolverType>
        void updateWellConnections(const CouplingChangeMetrics& metrics,
                                   const CouplingConfig& cfg,
                                   const SimulatorTimerInterface& timer,
                                   NonlinearSolverType& nonlinear_solver)
        {
            auto& problem = derived().simulator_.problem();
            if (cfg.valueOnlyWiUpdate && !metrics.structure_changed) {
                // Same connections: overwrite the CTFs in place, no well-structure rebuild.
                OpmLog::info("Updating fracture CTFs in place (value_only_wi_update, no structure change)");
                problem.addConnectionsToWell();
            } else {
                OpmLog::debug("Add connections in iterations");
                problem.addConnectionsToSchedual();
                problem.wellModel().beginTimeStep(); // reinitialize well structure
                problem.addConnectionsToWell();
                problem.emptyFractureLogger();
                {
                    // prepareTimeStep logs through the group-state helper's deferred
                    // logger, which the caller must establish.
                    auto& group_state_helper = problem.wellModel().groupStateHelper();
                    auto logger_guard = group_state_helper.pushLogger();
                    problem.wellModel().prepareTimeStep(group_state_helper.deferredLogger());
                }
                // Required after a well-structure change: skipping it kills the
                // well/fracture on fine grids.
                [[maybe_unused]] auto tmp_report = cfg.legacyParentSetup
                    ? this->runParentFirstIterationLegacy(timer, nonlinear_solver)
                    : this->runParentFirstIterationPreservingState(timer, nonlinear_solver);
                OpmLog::debug("End connections in iterations");
            }
            derived().onConnectionsUpdated();
        }
    };
}
