#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <opm/geomech/FractureCouplingOuterLoop.hpp>
#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/flow/NonlinearSystemBlackOilReservoir.hpp>
namespace Opm
{
    template<class TypeTag>
    class NonlinearSystemBlackOilReservoirGeoMech
        : public NonlinearSystemBlackOilReservoir<TypeTag>
        , public FractureCouplingOuterLoop<TypeTag, NonlinearSystemBlackOilReservoirGeoMech<TypeTag>>
    {
        public:
        using Parent = NonlinearSystemBlackOilReservoir<TypeTag>;
        using CouplingLoop = FractureCouplingOuterLoop<TypeTag, NonlinearSystemBlackOilReservoirGeoMech<TypeTag>>;
        friend CouplingLoop;
        using Simulator = typename Parent::Simulator;
        using Scalar = typename Parent::Scalar;
        using ModelParameters = typename Parent::ModelParameters;
        using EqVector = GetPropType<TypeTag, Properties::EqVector>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

        NonlinearSystemBlackOilReservoirGeoMech(Simulator& simulator,
                  const ModelParameters& param,
                  GetPropType<TypeTag, Properties::WellModel>& well_model,
                  const bool terminal_output):
                    Parent(simulator, param, well_model, terminal_output)
                  {

                  }

        // One nonlinear iteration of the flow parent, used by the coupling
        // loop for its setup passes.
        template <class NonlinearSolverType>
        SimulatorReportSingle flowNonlinearIteration(const SimulatorTimerInterface& timer,
                                                     NonlinearSolverType& nonlinear_solver)
        {
            return Parent::nonlinearIteration(timer, nonlinear_solver);
        }

        // Hook called by the coupling loop after a fracture-driven well
        // connection update; nothing to do for the end-of-step VEM solve.
        void onConnectionsUpdated()
        {
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver){
                SimulatorReportSingle report;
                const PropertyTree& prm = this->simulator_.problem().getGeoMechParam();
                std::string method = prm.get<std::string>("solver.method");
                if (method == "PostSolve") {
                    report = Parent::nonlinearIteration(timer, nonlinear_solver);
                } else if (method == "SeqMechFrac") {
                    report = this->nonlinearIterationSeqMechFrac(timer, nonlinear_solver);
                } else if (method == "SeqMech") {
                    report = this->nonlinearIterationSeqMech(timer, nonlinear_solver);
                } else if (method == "FullyImplicitMech") {
                    assert(false);
                } else {
                    assert(false);
                    Parent::nonlinearIterationNewton(timer, nonlinear_solver);
                }
                return report;
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationSeqMechFrac(const SimulatorTimerInterface& timer,
                                            NonlinearSolverType& nonlinear_solver){
            const PropertyTree& prm = this->simulator_.problem().getGeoMechParam();
            // The parent iteration advances the problem's iteration context, so
            // capture the current outer iteration number up front.
            const int iteration = this->simulator_.problem().iterationContext().iteration();
            if (iteration == 0) {
                this->topology_pending_counter_ = 0;
                this->topology_cooldown_counter_ = 0;
                // First outer iteration of a step always solves mechanics.
                this->last_coupling_composite_norm_ = std::numeric_limits<double>::max();
                this->mech_solves_this_step_ = 0;
            }
            bool implicit_flow = prm.get<bool>("solver.implicit_flow");
            SimulatorReportSingle report;
            {
                std::stringstream os;
                os << "Geomech nonlinearIterationSeqMechFrac with mechanical and fracture solve: " << iteration << std::endl;
                //std::cout << "Nonlinear itration Flow Solve:" << std::endl;
                report = Parent::nonlinearIteration(timer, nonlinear_solver);
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
            // Opt-in per-step cap (default 0 = unlimited = legacy behavior): once this
            // many mechanics solves have been performed within the current timestep,
            // further outer iterations reuse the stress and only re-solve the
            // fracture. The next timestep always starts with a fresh mech solve, so
            // the stress lags at most one step.
            const int max_mech_solves_per_step =
                prm.get<int>("solver.max_mech_solves_per_step", 0);
            const bool mech_cap_reached = (max_mech_solves_per_step > 0)
                && (this->mech_solves_this_step_ >= max_mech_solves_per_step);
            const bool mech_skippable = ((mech_skip_threshold > 0.0) && (iteration > 0)
                && (this->last_coupling_composite_norm_ < mech_skip_threshold))
                || mech_cap_reached;
            if(do_mech || do_fracture){
                if (mech_skippable) {
                    ++this->mech_solves_skipped_;
                    if (mech_cap_reached) {
                        OpmLog::info("Skipping mechanics solve (per-step cap "
                                     + std::to_string(max_mech_solves_per_step)
                                     + " reached); reusing stress. "
                                     + "mech solves performed/skipped: "
                                     + std::to_string(this->mech_solves_performed_) + "/"
                                     + std::to_string(this->mech_solves_skipped_));
                    } else {
                        OpmLog::info("Skipping mechanics solve (prev coupling change "
                                     + std::to_string(this->last_coupling_composite_norm_) + " < "
                                     + std::to_string(mech_skip_threshold) + "); reusing stress. "
                                     + "mech solves performed/skipped: "
                                     + std::to_string(this->mech_solves_performed_) + "/"
                                     + std::to_string(this->mech_solves_skipped_));
                    }
                } else {
                    OpmLog::info("Solve Geomechanics:");
                    this->simulator_.problem().geoMechModel().solveGeomechanics();
                    ++this->mech_solves_performed_;
                    ++this->mech_solves_this_step_;
                }
            }
            if(do_fracture && this->simulator_.problem().hasFractures()){
                this->fractureOuterBlock(report, timer, nonlinear_solver);
            }

            return report;

        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationSeqMech(const SimulatorTimerInterface& timer,
                                            NonlinearSolverType& nonlinear_solver){
            const PropertyTree& prm = this->simulator_.problem().getGeoMechParam();
            const int iteration = this->simulator_.problem().iterationContext().iteration();
            bool implicit_flow = prm.get<bool>("solver.implicit_flow");
            SimulatorReportSingle report;
            if(implicit_flow){
                assert(false);
            }else{
                report = Parent::nonlinearIteration(timer, nonlinear_solver);
            }
            //const PropertyTree& prm_frac = this->simulator_.problem().getFractureParam();
            int mech_max_it = prm.get<int>("solver.max_mech_it");
            if(iteration < mech_max_it){
                //simulator_.problem().geoMechModel().solveFracture();
                this->simulator_.problem().geoMechModel().solveGeomechanics();
                std::cout << "Geomech nonlinearIteration with mechanical solve:";// << iteration << std::endl;
                // TODO check convergence properly
                report.converged = false;
            }
            return report;

        }
    };
}
