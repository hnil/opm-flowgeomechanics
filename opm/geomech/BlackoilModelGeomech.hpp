#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <opm/simulators/flow/BlackoilModel.hpp>
namespace Opm
{
    template<class TypeTag>
    class BlackoilModelGeomech : public BlackoilModel<TypeTag>
    {
        public:
        using Parent = BlackoilModel<TypeTag>;
        using Simulator = typename Parent::Simulator;
        using Scalar = typename Parent::Scalar;
        using ModelParameters = typename Parent::ModelParameters;
        using EqVector = GetPropType<TypeTag, Properties::EqVector>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        //using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        //using ModelParameters = BlackoilModelParameters<Scalar>;
        BlackoilModelGeomech(Simulator& simulator,
                  const ModelParameters& param,
                  BlackoilWellModel<TypeTag>& well_model,
                  const bool terminal_output):
                    Parent(simulator, param, well_model, terminal_output)
                  {
                    
                  }

        private:
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

        struct StorageCacheBackup
        {
            bool enabled{false};
            std::vector<std::vector<EqVector>> values;
            std::vector<std::vector<unsigned char>> valid;
        };

        StorageCacheBackup backupStorageCache(const unsigned num_time_indices = 2) const
        {
            StorageCacheBackup backup;
            const auto& flow_model = this->simulator_.model();
            if (!flow_model.enableStorageCache()) {
                return backup;
            }

            backup.enabled = true;
            const unsigned num_dof = flow_model.numGridDof();
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

            const auto& flow_model = this->simulator_.model();
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
            return this->simulator_.model().solution(0);
        }

        void restoreSolution0(const SolutionVector& solution_backup)
        {
            auto& flow_model = this->simulator_.model();
            flow_model.solution(0) = solution_backup;
            flow_model.invalidateAndUpdateIntensiveQuantities(0);
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
            // keep state0 = state1 while settin up the system again
            auto& flow_model = this->simulator_.model();
            //flow_model.solution(0) = flow_model.solution(1); // well solve but will make excplict quantities change
            //flow_model.invalidateAndUpdateIntensiveQuantities(0);
            auto report = Parent::nonlinearIteration(0, timer, nonlinear_solver);
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
            auto report = Parent::nonlinearIteration(0, timer, nonlinear_solver);
            std::cout << "Finished parent first iteration" << std::endl;
            return report;
        }

        public:


        CouplingChangeMetrics evaluateFractureCouplingChange(
            const std::vector<std::vector<RuntimePerforation>>& allWellIndices) const
        {
            const auto& prm = this->simulator_.problem().getFractureParam();
            const double ctf_threshold_abs = prm.get("solver.ctf_change_threshold", 1e-15);
            const double ctf_threshold_rel = prm.get("solver.ctf_change_threshold_rel", 0.0);
            const double ctf_rel_scale = prm.get("solver.ctf_relative_scale", 1.0);
            const double perf_threshold_abs = prm.get("solver.perf_pressure_change_threshold", 1e30);
            const double perf_threshold_rel = prm.get("solver.perf_pressure_change_threshold_rel", 0.0);
            const double perf_rel_scale = prm.get("solver.perf_pressure_relative_scale", 1.0);
            const int verbosity = prm.get("solver.verbosity", 1);

            CouplingChangeMetrics metrics;
            const auto allWellIndices_new = this->simulator_.problem().getAllExtraWellIndices();
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



        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const int iteration,
                                                   const SimulatorTimerInterface& timer,
                                                  NonlinearSolverType& nonlinear_solver){
                SimulatorReportSingle report;
                const PropertyTree& prm = this->simulator_.problem().getGeomechParam();
                std::string method = prm.get<std::string>("solver.method");
                if (method == "PostSolve") {
                    report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
                } else if (method == "SeqMechFrac") {
                    report = this->nonlinearIterationSeqMechFrac(iteration, timer, nonlinear_solver);
                } else if (method == "SeqMech") {
                    report = this->nonlinearIterationSeqMech(iteration, timer, nonlinear_solver);
                } else if (method == "FullyImplicitMech") {
                    assert(false);
                } else {
                    assert(false);
                    std::cout << "Geomech nonlinearIterationNewton with mechanical solve:" << iteration;// << std::endl;
                    Parent::nonlinearIterationNewton(iteration, timer, nonlinear_solver);
                }
                return report;       
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationSeqMechFrac(const int iteration,
                                            const SimulatorTimerInterface& timer,
                                            NonlinearSolverType& nonlinear_solver){
            const PropertyTree& prm = this->simulator_.problem().getGeomechParam();
            if (iteration == 0) {
                topology_pending_counter_ = 0;
                topology_cooldown_counter_ = 0;
                // First outer iteration of a step always solves mechanics.
                last_coupling_composite_norm_ = std::numeric_limits<double>::max();
            }
            bool implicit_flow = prm.get<bool>("solver.implicit_flow");
            SimulatorReportSingle report;
            {
                std::stringstream os;
                os << "Geomech nonlinearIterationSeqMechFrac with mechanical and fracture solve: " << iteration << std::endl;
                //std::cout << "Nonlinear itration Flow Solve:" << std::endl;
                report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
                os << "Flow solve report converged: " << report.converged;// << std::endl;         
                OpmLog::info(os.str());
            }
            bool do_mech = true;
            bool do_fracture = true;
            if(implicit_flow){
                do_mech = report.converged;
                do_fracture = report.converged;
                int mech_max_it = prm.get<int>("solver.max_mech_it",10);
                if(iteration >= mech_max_it){
                    do_mech = false;
                }
                int frac_max_it = prm.get<int>("solver.max_frac_it",10);
                if(iteration >= frac_max_it){
                    do_fracture = false;
                }
            }
            
            //const PropertyTree& prm_frac = this->simulator_.problem().getFractureParam();
            

            // Mech-skip lever (opt-in, default off): skip the mechanics solve when
            // the previous outer iteration's fracture->well coupling change was below
            // threshold, i.e. the well flow (the only path by which the fracture
            // affects mechanics) is stable, so the mech solution would not move. The
            // fracture then reuses the previous stress. Self-correcting: if a skip was
            // wrong, the fracture changes the wells next iteration, the norm rises and
            // skipping stops. NOTE: this is a fracture-coupling proxy; a purely
            // thermal stress transient that does not move the wells is not detected,
            // so keep the threshold conservative on strongly thermal cases.
            const double mech_skip_threshold =
                prm.get<double>("solver.mech_skip_coupling_threshold", 0.0);
            const bool mech_skippable = (mech_skip_threshold > 0.0) && (iteration > 0)
                && (last_coupling_composite_norm_ < mech_skip_threshold);
            if(do_mech || do_fracture){
                if (mech_skippable) {
                    ++mech_solves_skipped_;
                    OpmLog::info("Skipping mechanics solve (prev coupling change "
                                 + std::to_string(last_coupling_composite_norm_) + " < "
                                 + std::to_string(mech_skip_threshold) + "); reusing stress. "
                                 + "mech solves performed/skipped: "
                                 + std::to_string(mech_solves_performed_) + "/"
                                 + std::to_string(mech_solves_skipped_));
                } else {
                    OpmLog::info("Solve Geomechanics:");
                    this->simulator_.problem().geomechModel().solveGeomechanics();
                    ++mech_solves_performed_;
                }
            }
            if(do_fracture && this->simulator_.problem().hasFractures()){
                //std::cout << "Solve Fractures:" << std::endl;
                std::stringstream os;
                //os << "Geomech nonlinearIteration with mechanical and fracture solve:" << iteration << std::endl;
                os << "Solve Fractures:";
                OpmLog::info(os.str());
                const bool legacy_coupling_change_logic =
                    prm.get<bool>("fractureparam.solver.legacy_coupling_change_logic", false);
                const bool legacy_parent_setup_iteration =
                    prm.get<bool>("fractureparam.solver.legacy_parent_setup_iteration", false);
                const auto allwellIndices = this->simulator_.problem().getAllExtraWellIndices();
                this->simulator_.problem().geomechModel().solveFractures();
                const bool fracture_converged = this->simulator_.problem().geomechModel().fractureModel().lastSolveStats().converged;
                const bool require_converged_fracture_for_wi_update =
                    prm.get<bool>("fractureparam.require_converged_fracture_for_wi_update", true);

                const auto coupling_metrics_local = this->evaluateFractureCouplingChange(allwellIndices);
                auto& comm = this->simulator_.gridView().comm();
                CouplingChangeMetrics coupling_metrics = coupling_metrics_local;
                coupling_metrics.structure_changed = comm.max(coupling_metrics_local.structure_changed);
                coupling_metrics.ctf_changed = comm.max(coupling_metrics_local.ctf_changed);
                coupling_metrics.perf_pressure_changed = comm.max(coupling_metrics_local.perf_pressure_changed);
                coupling_metrics.composite_norm = comm.max(coupling_metrics_local.composite_norm);
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

                if(do_update_connections){
                    std::cout << "Add connections in iterations" << std::endl;
                    this->simulator_.problem().addConnectionsToSchedual();// add new connections in the schedual
                    this->simulator_.problem().wellModel().beginTimeStep(); // reinitialize well structure
                    this->simulator_.problem().addConnectionsToWell(); // set the new well indices
                    this->simulator_.problem().emptyFractureLogger();
                    auto& local_deferredLogger = FractureModel::fractureLogger;
                    //this->simulator_.problem().wellModel().calculateExplicitQuantities(local_deferredLogger); //calcualte new explicite quantities TOFIX hould us old values
                    //this->simulator_.problem().wellModel().prepareTimeStep(local_deferredLogger);// hopefully set up all well realated stuff correctly

                    [[maybe_unused]] auto tmp_report =
                        legacy_parent_setup_iteration
                            ? this->runParentFirstIterationLegacy(timer, nonlinear_solver)
                            : this->runParentFirstIterationPreservingState(timer, nonlinear_solver);
                    std::cout << "End connections in iterations" << std::endl;
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
                const bool fracture_changed = coupling_metrics.changed();
                std::cout << "Fracture changed: " << fracture_changed << std::endl;
                // TODO check convergence properly
                if((!legacy_coupling_change_logic && !fracture_converged_global) || fracture_changed){
                    if (!fracture_converged_global) {
                        OpmLog::info("Keeping outer nonlinear loop active because fracture solve did not converge");
                    }
                    report.converged = false;
                 }
                // if(implicit_flow){
                //     // may do some extra things and + updating explicite quantities
                //     auto tmp_report = Parent::nonlinearIteration(0, timer, nonlinear_solver);
                // }
                comm.barrier();                
            }
            
            return report;
            
        }
        // template <class NonlinearSolverType>
        // SimulatorReportSingle nonlinearIterationFullyImplicit(const int iteration,
        //                                                       const SimulatorTimerInterface& timer,
        //                                                       NonlinearSolverType& nonlinear_solver)
        // {
        //     // -----------   Set up reports and timer   -----------
        //     SimulatorReportSingle report;
        //     Dune::Timer perfTimer;

        //     this->initialLinearization(
        //         report, iteration, nonlinear_solver.minIter(), nonlinear_solver.maxIter(), timer);

        //     this->geomecModel().setupAndUpdateGemechanics();

        //     // check convergence
        //     // -----------   If not converged, solve linear system and do Newton update  -----------
        //     if (!report.converged) {
        //         perfTimer.reset();
        //         perfTimer.start();
        //         report.total_newton_iterations = 1;

        //         // Compute the nonlinear update.
        //         unsigned nc = this->simulator_.model().numGridDof();
        //         BVector x(nc);

        //         // Solve the linear system.
        //         linear_solve_setup_time_ = 0.0;
        //         try {
        //             // Apply the Schur complement of the well model to
        //             // the reservoir linearized equations.
        //             // Note that linearize may throw for MSwells.
        //             this->wellModel().linearize(simulator().model().linearizer().jacobian(),
        //                                         simulator().model().linearizer().residual());
        //             this->geomecModel().setupAndUpdateGemechanics();

        //             // mechmatrix
        //             // reservoir matrix/operator

        //             // ---- Solve linear system ----
        //             solveJacobianSystem(x);

        //             report.linear_solve_setup_time += linear_solve_setup_time_;
        //             report.linear_solve_time += perfTimer.stop();
        //             report.total_linear_iterations += this->linearIterationsLastSolve();
        //         } catch (...) {
        //             report.linear_solve_setup_time += linear_solve_setup_time_;
        //             report.linear_solve_time += perfTimer.stop();
        //             report.total_linear_iterations += this->linearIterationsLastSolve();

        //             failureReport_ += report;
        //             throw; // re-throw up
        //         }

        //         perfTimer.reset();
        //         perfTimer.start();

        //         // handling well state update before oscillation treatment is a decision based
        //         // on observation to avoid some big performance degeneration under some circumstances.
        //         // there is no theorectical explanation which way is better for sure.
        //         this->wellModel().postSolve(x);

        //         if (param_.use_update_stabilization_) {
        //             // Stabilize the nonlinear update.
        //             bool isOscillate = false;
        //             bool isStagnate = false;
        //             nonlinear_solver.detectOscillations(
        //                 residual_norms_history_, residual_norms_history_.size() - 1, isOscillate, isStagnate);
        //             if (isOscillate) {
        //                 current_relaxation_ -= nonlinear_solver.relaxIncrement();
        //                 current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
        //                 if (terminalOutputEnabled()) {
        //                     std::string msg = "    Oscillating behavior detected: Relaxation set to "
        //                         + std::to_string(current_relaxation_);
        //                     OpmLog::info(msg);
        //                 }
        //             }
        //             nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
        //         }

        //         // ---- Newton update ----
        //         // Apply the update, with considering model-dependent limitations and
        //         // chopping of the update.
        //         this->updateSolution(x);

        //         report.update_time += perfTimer.stop();
        //     }

        //     return report;
        // }


        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationSeqMech(const int iteration,
                                             const SimulatorTimerInterface& timer,
                                            NonlinearSolverType& nonlinear_solver){
            const PropertyTree& prm = this->simulator_.problem().getGeomechParam();
            bool implicit_flow = prm.get<bool>("solver.implicit_flow");
            SimulatorReportSingle report;
            if(implicit_flow){
                assert(false);
            }else{
                report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
            }
            //const PropertyTree& prm_frac = this->simulator_.problem().getFractureParam();
            int mech_max_it = prm.get<int>("solver.max_mech_it");
            if(iteration < mech_max_it){
                //simulator_.problem().geomechModel().solveFracture();
                this->simulator_.problem().geomechModel().solveGeomechanics();
                std::cout << "Geomech nonlinearIteration with mechanical solve:";// << iteration << std::endl;
                // TODO check convergence properly
                report.converged = false;
            }
            return report;
            
        }
        // SimulatorReportSingle nonlinearIterationFullyImplicit(const int iteration,
        //                                      const SimulatorTimerInterface& timer,
        //                                     NonlinearSolverType& nonlinear_solver){
        //                                     }

        //void updateFailed(){
        //    Parent::updateFailed();
        //}

        private:
        

    };   
}
