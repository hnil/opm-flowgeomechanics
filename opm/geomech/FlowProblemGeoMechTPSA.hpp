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

#include <opm/common/utility/SymmTensor.hpp>

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

        //! Total in-situ stress at a cell, as consumed by the fracture model.
        //! No Biot/thermal diagonal is added, unlike GeoMechModel::stress: the TPSA
        //! solid pressure (p_s = lambda*div(u) - biot*dp - tcoeff*dT) already makes
        //! the reconstructed delta the TOTAL stress change.
        SymTensor stress(size_t globalIdx) const
        {
            // outputstress carries the initial-stress offset once
            // applyInitialOutputStress has run (STRESSEQUIL path), so this
            // is initstress + TPSA total-stress delta.
            const auto t = this->geoMechModel().outputstress(globalIdx);
            SymTensor s;
            for (int k = 0; k < 6; ++k) {
                s[k] = t[static_cast<VoigtIndex>(k)];
            }
            return s;
        }

        // Embedded fracture flow (aux cells) is VEM-only; the coupling loop
        // still compiles against these.
        void registerAuxiliaryCellModules()
        {
            MechParent::registerAuxiliaryCellModules();
            if (this->hasFractures()
                && this->getFractureParam().template get<std::string>(
                       "solver.fracture_flow_mode", std::string {"wi_upscaling"}) == "embedded") {
                OPM_THROW(std::runtime_error,
                          "solver.fracture_flow_mode=embedded is not supported with TPSA mechanics");
            }
        }
        bool fractureFlowIsEmbedded() const { return false; }
        void clearEmbeddedRestructured() {}
        void bindFractureAuxCells(bool = true, bool = false) {}
        double embeddedCouplingChange() const { return 0.0; }
        bool embeddedRestructuredMidStep() const { return false; }
        bool embeddedLayoutPending() const { return false; }
        bool takeEmbeddedPendingRound(int) { return false; }
        void holdEmbeddedOuterLoop() {}

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
            std::vector<SymmTensor<double>> offset(this->initstress_.size());
            for (std::size_t i = 0; i < offset.size(); ++i) {
                for (int k = 0; k < 6; ++k) {
                    offset[i][static_cast<VoigtIndex>(k)] = this->initstress_[i][k];
                }
            }
            this->geoMechModel().setOutputStressOffset(std::move(offset));
        }

        void beginTimeStep() override
        {
            if (this->gridView().comm().rank() == 0){
                OpmLog::debug("----------------------Start beginTimeStep (TPSA geomech)-------------------");
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
                OpmLog::debug("----------------------Start endTimeStep (TPSA geomech)-------------------");
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
