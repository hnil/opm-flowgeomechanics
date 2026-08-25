#pragma once

// Backend-agnostic fracture <-> flow coupling apparatus, extracted verbatim
// from NonlinearSystemBlackOilReservoirGeoMech so that nonlinear systems
// built on other mechanics backends (e.g. TPSA) can drive the same fracture
// outer loop.  CRTP: the Derived nonlinear system declares this class a
// friend and provides flowNonlinearIteration(timer, solver) (one nonlinear
// iteration of its flow parent).

#include <cmath>
#include <iostream>
#include <limits>
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

        bool shouldUpdateConnections(const CouplingChangeMetrics& metrics,
                                     const bool legacy_coupling_change_logic,
                                     const PropertyTree& prm)
        {
            if (legacy_coupling_change_logic) {
                return true;
            }

            const bool hysteresis_enable = prm.get<bool>("fractureparam.topology_hysteresis.enable", false);
            if (!hysteresis_enable) {
                return true;
            }

            const double on_threshold = prm.get<double>("fractureparam.topology_hysteresis.on_threshold", 1.0);
            const double off_threshold = prm.get<double>("fractureparam.topology_hysteresis.off_threshold", 0.3);
            const int confirm_iterations = prm.get<int>("fractureparam.topology_hysteresis.confirm_iterations", 1);
            const int cooldown_iterations = prm.get<int>("fractureparam.topology_hysteresis.cooldown_iterations", 0);
            const double emergency_threshold = prm.get<double>("fractureparam.topology_hysteresis.emergency_threshold", 10.0);

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
        // structure changed mid-step.  The pre-rebase code called
        // Parent::nonlinearIteration(0, ...): the parent saw iteration 0 and
        // redid its once-per-step initialization, and the outer iteration count
        // was unaffected.  With the iteration state moved into the problem's
        // NewtonIterationContext, reproduce that by resetting the context for
        // the parent call and restoring it afterwards.
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
            std::cout << "Running parent first iteration to set up the system and caches, while preserving state for geomechanics solve" << std::endl;
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
            auto report = this->runParentSetupIteration(timer, nonlinear_solver);
            derived().simulator_.problem().wellModel().restoreWGState(std::move(wgstate_backup));
            this->restoreSolution0(solution_backup);
            this->restoreStorageCache(storage_cache_backup);
            std::cout << "Finished parent first iteration" << std::endl;
            return report;
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle runParentFirstIterationLegacy(
            const SimulatorTimerInterface& timer,
            NonlinearSolverType& nonlinear_solver)
        {
            std::cout << "Running parent first iteration in legacy mode (no state restore)" << std::endl;
            auto report = this->runParentSetupIteration(timer, nonlinear_solver);
            std::cout << "Finished parent first iteration" << std::endl;
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
        // The fracture part of one outer nonlinear iteration: solve the
        // fractures against the current mechanics state, evaluate the
        // coupling change, decide on and perform well-connection updates,
        // and keep the outer loop active until the coupling has converged.
        // The caller has already decided that fractures should be solved
        // this iteration (do_fracture && problem.hasFractures()).
        template <class NonlinearSolverType>
        void fractureOuterBlock(SimulatorReportSingle& report,
                                const SimulatorTimerInterface& timer,
                                NonlinearSolverType& nonlinear_solver)
        {
            const PropertyTree& prm = derived().simulator_.problem().getGeoMechParam();
            {
                std::stringstream os;
                os << "Solve Fractures:";
                OpmLog::info(os.str());
            }
            const bool legacy_coupling_change_logic =
                prm.get<bool>("fractureparam.solver.legacy_coupling_change_logic", false);
            const bool legacy_parent_setup_iteration =
                prm.get<bool>("fractureparam.solver.legacy_parent_setup_iteration", false);
            const auto allwellIndices = derived().simulator_.problem().getAllExtraWellIndices();
            derived().simulator_.problem().fractureHost().solveFractures();
            const bool fracture_converged = derived().simulator_.problem().fractureHost().fractureModel().lastSolveStats().converged;
            const bool require_converged_fracture_for_wi_update =
                prm.get<bool>("fractureparam.require_converged_fracture_for_wi_update", true);

            const auto coupling_metrics_local = this->evaluateFractureCouplingChange(allwellIndices);
            auto& comm = derived().simulator_.gridView().comm();
            CouplingChangeMetrics coupling_metrics = coupling_metrics_local;
            coupling_metrics.structure_changed = comm.max(coupling_metrics_local.structure_changed);
            coupling_metrics.ctf_changed = comm.max(coupling_metrics_local.ctf_changed);
            coupling_metrics.perf_pressure_changed = comm.max(coupling_metrics_local.perf_pressure_changed);
            coupling_metrics.composite_norm = comm.max(coupling_metrics_local.composite_norm);
            // Opt-in outer-residual convergence (coupling_convergence_mode=residual):
            // declare outer convergence from a dimensionless residual instead of the
            // threshold booleans. Needs the rel-change extrema globally consistent.
            const std::string coupling_convergence_mode = prm.get<std::string>(
                "fractureparam.solver.coupling_convergence_mode", "threshold");
            const bool residual_mode = (coupling_convergence_mode == "residual");
            if (residual_mode) {
                coupling_metrics.max_ctf_change_rel = comm.max(coupling_metrics_local.max_ctf_change_rel);
                coupling_metrics.max_perf_change_rel = comm.max(coupling_metrics_local.max_perf_change_rel);
            }
            // Remember this iteration's fracture->well coupling change so the next
            // iteration can decide whether the mechanics solve may be skipped.
            last_coupling_composite_norm_ = coupling_metrics.composite_norm;
            const bool fracture_converged_global = comm.min(fracture_converged);

            bool addconnections = prm.get<bool>("fractureparam.addconnections");
            const bool may_update_connections = addconnections && (legacy_coupling_change_logic
                || !require_converged_fracture_for_wi_update
                || fracture_converged_global);
            const bool do_update_connections = may_update_connections
                && this->shouldUpdateConnections(coupling_metrics, legacy_coupling_change_logic, prm);

            // Opt-in cheap update path: when only the CTF *values* changed (no new
            // or removed connections), write them in place on the existing well
            // perforations (addFracturePerforations overwrites well_index_fracture_)
            // and skip the schedule rebuild + well-model re-init + prepareTimeStep
            // + parent setup iteration. Structure changes always take the full path.
            const bool value_only_wi_update =
                prm.get<bool>("fractureparam.solver.value_only_wi_update", false);
            if (derived().simulator_.problem().fractureFlowIsEmbedded()) {
                // The fracture flows through degrees of freedom of its own, so a fracture
                // solve changed apertures and possibly opened cells: re-describe them to
                // the reservoir and refresh the well's perforations of them.  No schedule
                // rebuild and no upscaled well index -- re-adding those here is how the
                // wells would silently fall back to the representation being replaced.
                derived().simulator_.problem().bindFractureAuxCells(/*allowTopologyChange=*/false);
                // The well's perforations of the fracture keep the indices they were
                // given at the step boundary: within the step both the topology and the
                // perforation factors are held fixed -- the sequentially implicit lag --
                // and re-registering here would demand a well-structure rebuild the
                // step is not going to perform.

                // The coupling residual is the change in what the flow was just fed --
                // total conductance and pore volume of the binding -- not the well-index
                // change list, which is computed above but never applied here.  Contact
                // chatter moves individual cells' indices every iteration; the
                // aggregates cancel it the way the upscaled index summed over it.
                const double coupling_tolerance =
                    prm.get<double>("fractureparam.solver.coupling_tolerance", 1e-3);
                const double embedded_change =
                    derived().simulator_.problem().embeddedCouplingChange();
                if (!fracture_converged_global || embedded_change > coupling_tolerance) {
                    if (!fracture_converged_global) {
                        OpmLog::info("Keeping outer loop active: fracture solve did not converge");
                    } else {
                        OpmLog::info("Keeping outer loop active: embedded coupling change "
                                     + std::to_string(embedded_change) + " > "
                                     + std::to_string(coupling_tolerance));
                    }
                    report.converged = false;
                }
                comm.barrier();
                return;
            }
            else if(do_update_connections){
                if (value_only_wi_update && !coupling_metrics.structure_changed) {
                    OpmLog::info("Updating fracture CTFs in place (value_only_wi_update, no structure change)");
                    derived().simulator_.problem().addConnectionsToWell();
                } else {
                std::cout << "Add connections in iterations" << std::endl;
                derived().simulator_.problem().addConnectionsToSchedual();// add new connections in the schedual
                derived().simulator_.problem().wellModel().beginTimeStep(); // reinitialize well structure
                derived().simulator_.problem().addConnectionsToWell(); // set the new well indices
                derived().simulator_.problem().emptyFractureLogger();
                {
                    // prepareTimeStep logs through the group-state
                    // helper's deferred logger, which the caller must
                    // establish (same pattern as assemble()).
                    auto& group_state_helper = derived().simulator_.problem().wellModel().groupStateHelper();
                    auto logger_guard = group_state_helper.pushLogger();
                    derived().simulator_.problem().wellModel().prepareTimeStep(group_state_helper.deferredLogger());
                }

                // The parent setup iteration is REQUIRED when the well structure
                // actually changed (new fracture connections): skipping it kills
                // the well/fracture on fine grids (verified 2026-07-24 on
                // CASE_REFINE_nx_21..._THERMAL_MSW).  It is redundant only for
                // value-only CTF changes, which take the value_only path above.
                // Keep it; runParentFirstIterationPreservingState is now
                // state-neutral for the wells (WGState snapshot/restore).
                [[maybe_unused]] auto tmp_report =
                    legacy_parent_setup_iteration
                        ? this->runParentFirstIterationLegacy(timer, nonlinear_solver)
                        : this->runParentFirstIterationPreservingState(timer, nonlinear_solver);
                std::cout << "End connections in iterations" << std::endl;
                }
                derived().onConnectionsUpdated();
            }
            else if (may_update_connections && !legacy_coupling_change_logic) {
                std::stringstream hos;
                hos << "Deferring fracture-driven connection update due to topology hysteresis"
                    << ", composite_norm=" << coupling_metrics.composite_norm
                    << ", pending=" << topology_pending_counter_
                    << ", cooldown=" << topology_cooldown_counter_;
                OpmLog::info(hos.str());
            }
            else if (addconnections && require_converged_fracture_for_wi_update && !fracture_converged) {
                OpmLog::info("Skipping fracture-driven connection update because fracture solve did not converge");
            }
            if (residual_mode) {
                // Outer convergence from a physical residual: the largest relative
                // change of any fracture CTF or perforation pressure produced by
                // this iteration's fracture solve. Structure changes and a
                // non-converged inner solve always keep the loop active.
                const double coupling_tolerance =
                    prm.get<double>("fractureparam.solver.coupling_tolerance", 1e-3);
                const double outer_residual = std::max(coupling_metrics.max_ctf_change_rel,
                                                       coupling_metrics.max_perf_change_rel);
                const bool outer_converged = fracture_converged_global
                    && !coupling_metrics.structure_changed
                    && outer_residual <= coupling_tolerance;
                std::stringstream ros;
                ros << "Outer coupling residual: " << outer_residual
                    << " (tol " << coupling_tolerance << ")"
                    << ", structure_changed=" << (coupling_metrics.structure_changed ? "true" : "false")
                    << ", fracture_converged=" << (fracture_converged_global ? "true" : "false")
                    << ", outer_converged=" << (outer_converged ? "true" : "false");
                OpmLog::info(ros.str());
                if (!outer_converged) {
                    report.converged = false;
                }
            } else {
            const bool fracture_changed = coupling_metrics.changed();
            std::cout << "Fracture changed: " << fracture_changed << std::endl;
            // TODO check convergence properly
            if((!legacy_coupling_change_logic && !fracture_converged_global) || fracture_changed){
                if (!fracture_converged_global) {
                    OpmLog::info("Keeping outer nonlinear loop active because fracture solve did not converge");
                }
                report.converged = false;
             }
            }
            comm.barrier();
        }
    };
}
