#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
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
                this->growth_rounds_this_step_ = 0;
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
            // hasFractures() only says the deck has FRAC; a deck with no WSEED
            // has no fracture model, and then this whole block is a no-op.
            const bool have_fracture = this->simulator_.problem().hasFractures()
                && this->simulator_.problem().fractureHost().fractureModelActive();
            if(do_fracture && this->simulator_.problem().hasFractures()){
                this->fractureOuterBlock(report, timer, nonlinear_solver);
                if (have_fracture)
                    this->simulator_.problem().geoMechModel().fractureModel()
                        .writeIterationSnapshots(timer.currentStepNum(), iteration, "outer");

                // Opt-in growth sub-iterations (default 0 = off): during fracture
                // establishment the coupling legitimately changes every round, so
                // iterate flow+mech+fracture here, inside one outer Newton
                // iteration, instead of burning the flow Newton budget (whose
                // exhaustion chops dt, which cannot help: the required growth is
                // a ratio per solve, not per unit time).
                const int max_growth_it = prm.get<int>("solver.max_growth_iterations", 0);
                // Growth-driven switch (opt-in): keep iterating only while a
                // fracture is growing "faster than the timestep" — i.e. it hit
                // its per-solve growth cap (finite maxFlowTimeStep). When growth
                // is small the lagged (PostSolve-style) coupling of this round is
                // accepted, which is what PostSolve does on every step; the
                // sequential-implicit rounds are spent only where they matter.
                const bool growth_switch = prm.get<bool>("solver.growth_driven_switch", false);
                auto fracture_at_growth_cap = [&]() {
                    if (!have_fracture)
                        return false;
                    return this->simulator_.problem().geoMechModel().fractureModel()
                        .template maxFlowTimeStep<TypeTag, Simulator>(this->simulator_)
                        < std::numeric_limits<double>::max();
                };
                if (growth_switch && !report.converged && !fracture_at_growth_cap()) {
                    OpmLog::info("Fracture coupling change accepted lagged (growth_driven_switch: "
                                 "no fracture at its growth cap this round)");
                    report.converged = true;
                }
                // Opt-in propagation-aware acceptance (default off). The step is
                // otherwise accepted on the coupling change alone, which can settle
                // while the front is far past its propagation criterion - the state
                // behind the establishment pressure spike ends the step at
                // K1/K1c ~ 30. Refusing to converge while K1 exceeds the limit sends
                // the growth loop round again instead.
                const bool k1_accept = prm.get<bool>("solver.require_k1_below_k1c", false);
                const double k1_accept_factor = prm.get<double>("solver.k1_accept_factor", 1.0);
                // K1 > K1c means the front wants to advance, but the cells it would
                // advance into are retracted again by removeNewZeroWithCells when
                // they end up zero-width - the K1 (front-profile) and the local
                // force/width criteria do not have to agree.  Blocking acceptance
                // on K1 alone then grows and retracts forever, so only block while
                // the fracture is actually getting bigger; a blocked front is
                // accepted with a warning instead of cycling.
                double k1_prev_area = -1.0;
                auto k1_blocks_acceptance = [&]() {
                    if (!k1_accept || !have_fracture)
                        return false;
                    const auto& fmod = this->simulator_.problem().geoMechModel().fractureModel();
                    auto& comm = this->simulator_.vanguard().grid().comm();
                    const double worst = comm.max(fmod.maxK1Ratio());
                    if (worst <= k1_accept_factor) {
                        k1_prev_area = -1.0;
                        return false;
                    }
                    const double area = comm.sum(fmod.totalFractureArea());
                    const bool advancing = (k1_prev_area < 0.0) || (area > k1_prev_area * 1.001);
                    k1_prev_area = area;
                    if (!advancing) {
                        OpmLog::warning("Fracture propagation criterion still unmet (max K1/K1c = "
                                        + std::to_string(worst)
                                        + ") but the front stopped advancing at "
                                        + std::to_string(area) + " m2; accepting");
                        return false;
                    }
                    OpmLog::info("Fracture propagation criterion not met (max K1/K1c = "
                                 + std::to_string(worst) + " > "
                                 + std::to_string(k1_accept_factor) + "); not accepting yet");
                    return true;
                };
                if (k1_blocks_acceptance())
                    report.converged = false;
                while (max_growth_it > 0 && implicit_flow && !report.converged
                       && this->growth_rounds_this_step_ < max_growth_it) {
                    const int growth_round = ++this->growth_rounds_this_step_;
                    // Re-solve flow to convergence without advancing the outer
                    // Newton counter; after a connection update this needs
                    // several iterations, exactly like a fresh timestep.
                    SimulatorReportSingle flow_report;
                    bool flow_ok = false;
                    {
                        const SetupIterationContextGuard guard{this->simulator_.problem()};
                        const int max_flow_it =
                            prm.get<int>("solver.growth_flow_iterations", 12);
                        for (int fit = 0; fit < max_flow_it; ++fit) {
                            flow_report = Parent::nonlinearIteration(timer, nonlinear_solver);
                            report += flow_report;
                            if (flow_report.converged) {
                                flow_ok = true;
                                break;
                            }
                        }
                    }
                    if (!flow_ok) {
                        // genuine flow trouble: hand back to the outer Newton loop
                        report.converged = false;
                        OpmLog::info("Fracture growth sub-iteration "
                                     + std::to_string(growth_round)
                                     + ": flow did not reconverge, returning to outer loop");
                        break;
                    }
                    const bool growth_mech_cap_reached = (max_mech_solves_per_step > 0)
                        && (this->mech_solves_this_step_ >= max_mech_solves_per_step);
                    if (growth_mech_cap_reached) {
                        ++this->mech_solves_skipped_;
                    } else {
                        this->simulator_.problem().geoMechModel().solveGeomechanics();
                        ++this->mech_solves_performed_;
                        ++this->mech_solves_this_step_;
                    }
                    SimulatorReportSingle round_report = flow_report;
                    this->fractureOuterBlock(round_report, timer, nonlinear_solver);
                    if (have_fracture)
                        this->simulator_.problem().geoMechModel().fractureModel()
                            .writeIterationSnapshots(timer.currentStepNum(), 1000 + growth_round, "growth");
                    report.converged = round_report.converged;
                    if (report.converged && k1_blocks_acceptance())
                        report.converged = false;
                    if (growth_switch && !report.converged && !fracture_at_growth_cap()) {
                        OpmLog::info("Fracture growth stopped; accepting lagged coupling "
                                     "(growth_driven_switch)");
                        report.converged = true;
                    }
                    OpmLog::info("Fracture growth sub-iteration "
                                 + std::to_string(growth_round) + "/"
                                 + std::to_string(max_growth_it)
                                 + (round_report.converged ? ": coupling settled"
                                                           : ": still growing"));
                }
                // Per-step growth guard (opt-in): a fracture that grew more than
                // allowed within this step propagated on a pressure the flow had
                // not yet responded to. Raised here, inside the retry scope, so
                // the timestepper chops dt and re-solves from the checkpoint.
                if (report.converged && have_fracture) {
                    const std::string violation = this->simulator_.problem()
                        .geoMechModel().fractureModel().growthGuardViolation();
                    const int any_violation =
                        this->simulator_.vanguard().grid().comm().max(violation.empty() ? 0 : 1);
                    if (any_violation) {
                        OpmLog::warning("Fracture growth guard: "
                                        + (violation.empty() ? std::string("violation on another rank")
                                                             : violation));
                        OPM_THROW_NOLOG(NumericalProblem, "fracture growth guard: " + violation);
                    }
                }
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
                // Never implemented; assert(false) meant a release build walked
                // straight into undefined behaviour instead of saying so.
                OPM_THROW(std::runtime_error,
                          "solver.method=SeqMech with solver.implicit_flow=true is not "
                          "implemented - use SeqMechFrac, PostSolve, or set implicit_flow=false");
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
