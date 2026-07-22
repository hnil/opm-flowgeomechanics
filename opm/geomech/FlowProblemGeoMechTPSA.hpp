#ifndef OPM_FLOW_PROBLEM_GEOMECH_TPSA_HPP
#define OPM_FLOW_PROBLEM_GEOMECH_TPSA_HPP

// Fracture-enabled flow problem with the TPSA mechanics backend: the shared
// FlowProblemMech layer (field properties, STRESSEQUIL initial stress,
// fracture configuration and connection plumbing) over FlowProblemTPSA, with
// a FractureMechHost driving the fracture machinery against the TPSA stress
// field.
//
// Total stress mirrors GeoMechModel::stress for the VEM backend:
//   stress(cell) = initstress_[cell]              (STRESSEQUIL, shared layer)
//                + TPSA effective-stress delta    (cached reconstruction)
//                + diag(mechPotentialForce(cell)) (biot*(p-p0) + tcoeff*(T-T0))

#include <opm/geomech/FlowProblemMech.hpp>
#include <opm/geomech/FractureMechHost.hpp>

#include <opm/simulators/flow/FlowProblemTPSA.hpp>

#include <string>
#include <vector>

namespace Opm{

    template<typename TypeTag>
    class FlowProblemGeoMechTPSA
        : public FlowProblemMech<TypeTag, FlowProblemTPSA<TypeTag>, FlowProblemGeoMechTPSA<TypeTag>>
    {
    public:
        using MechParent = FlowProblemMech<TypeTag, FlowProblemTPSA<TypeTag>, FlowProblemGeoMechTPSA<TypeTag>>;
        using Parent = FlowProblemTPSA<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using SymTensor = Dune::FieldVector<double,6>;

        explicit FlowProblemGeoMechTPSA(Simulator& simulator):
            MechParent(simulator),
            fracHost_(simulator)
        {
        }

        static void registerParameters(){
            MechParent::registerParameters();
        }

        //! The fracture-driving host used by the shared problem layer and
        //! the coupling loop.
        FractureMechHost<TypeTag>& fractureHost()
        { return fracHost_; }

        const FractureMechHost<TypeTag>& fractureHost() const
        { return fracHost_; }

        //! Total in-situ stress at a cell, as consumed by the fracture
        //! model.
        //!
        //! Unlike the VEM GeoMechModel::stress, no potential (Biot/thermal)
        //! diagonal is added here: the TPSA solid pressure is a total-pressure
        //! variable (p_s = lambda*div(u) - biot*dp - tcoeff*dT), so the
        //! reconstructed stress delta is already the TOTAL stress change.
        //! Adding the potential force again would double-count the coupling.
        SymTensor stress(size_t globalIdx) const
        {
            // outputstress carries the initial-stress offset once
            // applyInitialOutputStress has run (STRESSEQUIL path), so this
            // is initstress + TPSA total-stress delta.
            return this->geoMechModel().outputstress(globalIdx);
        }

        // ///
        // Backend hooks for the shared initial-stress handling
        // ///
        void prepareMechForInit() override
        {
        }

        void initializeStressFromMechSolve() override
        {
            OPM_THROW(std::runtime_error,
                      "TPSA mechanics requires STRESSEQUIL for now; gravity "
                      "equilibration of the initial stress is not yet "
                      "implemented for this backend");
        }

        void applyInitialOutputStress() override
        {
            // Seed the TPSA output stress with the initial stress so
            // outputstress()/BSTRSS report total stress; the TPSA solve
            // itself stays relative to the unstressed initial state.
            std::vector<SymTensor> offset(this->initstress_.size());
            for (std::size_t i = 0; i < offset.size(); ++i) {
                offset[i] = this->initstress_[i];
            }
            this->geoMechModel().setOutputStressOffset(std::move(offset));
        }

        void beginTimeStep() override
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start beginTimeStep (TPSA geomech)-------------------\n"
                << std::flush;
            }
            Parent::beginTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures()){
                  if(!(int(this->cstress_.size()) == int(this->gridView().size(0)))){
                        OPM_THROW(std::runtime_error,"CSTRESS not set but fractures exists");
                    }
                }
                fracHost_.beginTimeStep();
                if(this->hasFractures()){
                    this->wellModel().beginTimeStep();// just to be sure well conteiner is reinitialized
                    this->addConnectionsToWell(); // modify wells WI wiht fracture well
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("Begin time step geomech failed:",this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
        }

        void endTimeStep() override
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep (TPSA geomech)-------------------\n"
                << std::flush;
            }
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                // The TPSA mechanics is already converged inside the
                // nonlinear iteration (fixed-stress/lagged); refresh the
                // stress snapshot the fracture model samples and run the
                // post-step fracture solve.
                this->geoMechModel().updateStressCache();
                if(this->hasFractures()){
                    fracHost_.solveFractures();
                    if(fracHost_.fractureModelActive()){
                        this->addConnectionsToSchedual();
                        this->gridView().comm().barrier();
                        fracHost_.fractureModel().
                            assignGeoMechWellState(this->wellModel_.wellState());
                    }
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("End time step geomech failed: ", this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
            Parent::endTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures() ){
                    fracHost_.updateFilterCakePropertiesOnFractures();
                }
            }
            if(Parameters::Get<Parameters::EnableWriteAllSolutions>()){
                fracHost_.writeFractureSolution();
            }
        }

    private:
        FractureMechHost<TypeTag> fracHost_;
    };
}
#endif
