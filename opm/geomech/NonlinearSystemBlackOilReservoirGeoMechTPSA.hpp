#pragma once

// Nonlinear system combining the TPSA fixed-stress/lagged flow-mechanics
// coupling (upstream NonlinearSystemBlackOilReservoirTPSA) with the fracture
// coupling outer loop shared with the VEM path (FractureCouplingOuterLoop).
//
// Composition: the fracture solves ONLY after the fixed-stress sequence has
// converged -- the fracture must see a mechanics state consistent with the
// flow state, and running it mid-sequence would both feed it a stale stress
// and invalidate the numLinearizations()==1 convergence signal.  When the
// fracture block changes well connections (or has not converged) it sets
// report.converged=false, the outer Newton loop calls nonlinearIteration
// again, and the fixed-stress sub-loop restarts cleanly (its history roll is
// once-per-timestep, see TpsaModel::prepareTPSA).

#include <cmath>
#include <iostream>
#include <limits>
#include <opm/geomech/FractureCouplingOuterLoop.hpp>
#include <opm/simulators/flow/NonlinearSystemBlackOilReservoirTPSA.hpp>

namespace Opm
{
    template<class TypeTag>
    class NonlinearSystemBlackOilReservoirGeoMechTPSA
        : public NonlinearSystemBlackOilReservoirTPSA<TypeTag>
        , public FractureCouplingOuterLoop<TypeTag, NonlinearSystemBlackOilReservoirGeoMechTPSA<TypeTag>>
    {
        public:
        using Parent = NonlinearSystemBlackOilReservoirTPSA<TypeTag>;
        using FlowParent = NonlinearSystemBlackOilReservoir<TypeTag>;
        using CouplingLoop = FractureCouplingOuterLoop<TypeTag, NonlinearSystemBlackOilReservoirGeoMechTPSA<TypeTag>>;
        friend CouplingLoop;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using ModelParameters = typename Parent::ModelParameters;

        NonlinearSystemBlackOilReservoirGeoMechTPSA(Simulator& simulator,
                  const ModelParameters& param,
                  GetPropType<TypeTag, Properties::WellModel>& well_model,
                  const bool terminal_output):
                    Parent(simulator, param, well_model, terminal_output)
                  {
                  }

        //! One nonlinear iteration of the plain FLOW parent (not the TPSA
        //! sequence), used by the coupling loop for its setup passes after
        //! a well-structure change.
        template <class NonlinearSolverType>
        SimulatorReportSingle flowNonlinearIteration(const SimulatorTimerInterface& timer,
                                                     NonlinearSolverType& nonlinear_solver)
        {
            return FlowParent::nonlinearIteration(timer, nonlinear_solver);
        }

        //! A fracture-driven well-structure change invalidates the current
        //! fixed-stress sequence; restart it defensively.
        void onConnectionsUpdated()
        {
            this->resetFixedStressIterations();
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
        {
            const PropertyTree& prm = this->simulator_.problem().getGeoMechParam();
            const std::string method = prm.get<std::string>("solver.method");
            if (method == "PostSolve") {
                // Mechanics inside the TPSA sequence; fractures from
                // FlowProblemGeoMechTPSA::endTimeStep.
                return Parent::nonlinearIteration(timer, nonlinear_solver);
            }
            if (method == "SeqMech") {
                OPM_THROW(std::runtime_error,
                          "solver.method SeqMech is meaningless with the TPSA "
                          "backend: mechanics is already solved inside the "
                          "nonlinear iteration.  Use PostSolve or SeqMechFrac.");
            }
            if (method != "SeqMechFrac") {
                OPM_THROW(std::runtime_error,
                          "Unknown solver.method '" + method +
                          "' for the TPSA geomech nonlinear system");
            }

            const int iteration = this->simulator_.problem().iterationContext().iteration();
            if (iteration == 0) {
                this->topology_pending_counter_ = 0;
                this->topology_cooldown_counter_ = 0;
                this->last_coupling_composite_norm_ = std::numeric_limits<double>::max();
            }

            // The whole fixed-stress (or lagged) flow<->mechanics machinery,
            // untouched.  Not converged => still inside the sequence.
            SimulatorReportSingle report = Parent::nonlinearIteration(timer, nonlinear_solver);
            if (!report.converged) {
                return report;
            }

            if (this->simulator_.problem().hasFractures()) {
                // The fracture samples the mechanics state through the
                // problem's stress(); refresh the cached reconstruction for
                // the converged state first.
                this->simulator_.problem().geoMechModel().updateStressCache();
                this->fractureOuterBlock(report, timer, nonlinear_solver);
            }

            return report;
        }
    };
}
