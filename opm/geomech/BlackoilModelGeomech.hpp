#pragma once
#include <iostream>
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

        public:


                  bool
                  fractureChanged(const std::vector<std::vector<RuntimePerforation>>& allWellIndices)
                  {                    
                       const auto& prm = this->simulator_.problem().getFractureParam();
                       double ctf_threshold = prm.get("solver.ctf_change_threshold",1e-15);
                       int verbosity = prm.get("solver.verbosity",1);
                      const auto allWellIndices_new
                          = this->simulator_.problem().getAllExtraWellIndices();
                      bool structure_changed = false;    
                      double max_ctf_change = 0.0;
                      if ( !( allWellIndices.size() == allWellIndices_new.size()) ) {
                           structure_changed = true;
                          //return true;
                      } else {
                          //bool changed = false;
                          for (size_t k = 0; k < allWellIndices.size(); ++k) {
                              const auto& wellIndices = allWellIndices[k];
                              const auto& wellIndices_new = allWellIndices_new[k];
                              if (!(wellIndices.size() == wellIndices_new.size())) {
                                  structure_changed = true;
                                  return true;
                              } else {
                                  for (size_t i = 0; i < wellIndices.size(); ++i) {
                                      if (!(wellIndices[i].cell == wellIndices_new[i].cell)) {
                                          //changed = true;
                                          structure_changed = true;
                                        //    std::cout << "Fracture structure changed for well cell " << wellIndices[i].cell << std::endl;        
                                        //  return true;
                                      };
                                      double ctf_change = std::abs(wellIndices[i].ctf - wellIndices_new[i].ctf);
                                      max_ctf_change = std::max(max_ctf_change, ctf_change);
                                      if (max_ctf_change  > ctf_threshold) {
                                        if (verbosity > 1){
                                          std::stringstream os;
                                          os << "Fracture ctf changed for well  cell " << wellIndices[i].cell << std::endl;
                                          os << "ctf changed from " << wellIndices[i].ctf << " to " << wellIndices_new[i].ctf << std::endl;
                                          OpmLog::info(os.str());
                                        }
                                        //  changed = true;
                                        //  return true;
                                      }
                                  }
                              }
                          }
                        //   if (changed) {
                        //       assert(false);
                        //       return true;
                        //   }else{
                        //       return false;
                        //   }
                      }
                      if(structure_changed){
                        std::stringstream os;                
                        os << "Fracture structure changed, max ctf change: " << max_ctf_change << std::endl;
                        OpmLog::info(os.str());
                        return true;
                      }
                      std::stringstream os;                
                      os << "Fracture ctf changed, max ctf change: " << max_ctf_change << std::endl;
                      OpmLog::info(os.str());
                      if(max_ctf_change > ctf_threshold){
                        
                        return true;
                      }
                      return false;
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
            

            if(do_mech || do_fracture){
                std::stringstream os;
                os << "Solve Geomechanics:";// << std::endl;
                OpmLog::info(os.str());
                this->simulator_.problem().geomechModel().solveGeomechanics();
            }
            if(do_fracture && this->simulator_.problem().hasFractures()){
                //std::cout << "Solve Fractures:" << std::endl;
                std::stringstream os;
                //os << "Geomech nonlinearIteration with mechanical and fracture solve:" << iteration << std::endl;
                os << "Solve Fractures:";
                OpmLog::info(os.str());
                const auto allwellIndices = this->simulator_.problem().getAllExtraWellIndices();
                this->simulator_.problem().geomechModel().solveFractures();
                bool addconnections = prm.get<bool>("fractureparam.addconnections");
                if(addconnections){
                    std::cout << "Add connections in iterations" << std::endl;
                    this->simulator_.problem().addConnectionsToSchedual();// add new connections in the schedual
                    this->simulator_.problem().wellModel().beginTimeStep(); // reinitialize well structure
                    this->simulator_.problem().addConnectionsToWell(); // set the new well indices
                    this->simulator_.problem().emptyFractureLogger();
                    auto& local_deferredLogger = FractureModel::fractureLogger;
                    //this->simulator_.problem().wellModel().calculateExplicitQuantities(local_deferredLogger); //calcualte new explicite quantities TOFIX hould us old values
                    //this->simulator_.problem().wellModel().prepareTimeStep(local_deferredLogger);// hopefully set up all well realated stuff correctly

                    [[maybe_unused]] auto tmp_report =
                        this->runParentFirstIterationPreservingState(timer, nonlinear_solver);
                    std::cout << "End connections in iterations" << std::endl;
                }
                auto& comm = this->simulator_.gridView().comm();
                bool fracture_changed_local = this->fractureChanged(allwellIndices);
                bool fracture_changed = comm.max(fracture_changed_local);
                std::cout << "Fracture changed: " << fracture_changed << std::endl;
                // TODO check convergence properly
                if(fracture_changed){
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
