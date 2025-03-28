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
        //using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        //using ModelParameters = BlackoilModelParameters<Scalar>;
        BlackoilModelGeomech(Simulator& simulator,
                  const ModelParameters& param,
                  BlackoilWellModel<TypeTag>& well_model,
                  const bool terminal_output):
                    Parent(simulator, param, well_model, terminal_output)
                  {
                    
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
                    std::cout << "Geomech nonlinearIterationNewton with mechanical solve:" << iteration << std::endl;
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
            if(implicit_flow){
                assert(false);
            }else{
                report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
            }
            const PropertyTree& prm_frac = this->simulator_.problem().getFractureParam();
            int max_it = prm.get<int>("solver.max_mech_it");
            if(iteration < max_it){
                this->simulator_.problem().geomechModel().solveGeomechanics();
                if(this->simulator_.problem().hasFractures()){
                    this->simulator_.problem().geomechModel().solveFractures();
                }
                std::cout << "Geomech nonlinearIteration with mechanical and fracture solve:" << iteration << std::endl;
                bool addconnections = prm.get<bool>("fractureparam.addconnections");
                if(addconnections){
                    std::cout << "Add connections in iterations" << std::endl;
                    this->simulator_.problem().addConnectionsToSchedual();
                    this->simulator_.problem().wellModel().beginTimeStep();
                    this->simulator_.problem().addConnectionsToWell();
                }
                // TODO check convergence properly
                report.converged = false;
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
            bool implicit_flow = prm.get<bool>("method.implicit_flow");
            SimulatorReportSingle report;
            if(implicit_flow){
                assert(false);
            }else{
                report = Parent::nonlinearIteration(iteration, timer, nonlinear_solver);
            }
            const PropertyTree& prm_frac = this->simulator_.problem().getFractureParam();
            int max_it = prm.get<int>("method.max_mech_it");
            if(iteration < max_it){
                //simulator_.problem().geomechModel().solveFracture();
                this->simulator_.problem().geomechModel().solveGeomechanics();
                std::cout << "Geomech nonlinearIteration with mechanical solve:" << iteration << std::endl;
                // TODO check convergence properly
                report.converged = false;
            }
            return report;
            
        }
        // SimulatorReportSingle nonlinearIterationFullyImplicit(const int iteration,
        //                                      const SimulatorTimerInterface& timer,
        //                                     NonlinearSolverType& nonlinear_solver){
        //                                     }

        private:
        

    };   
}
