/*
  Copyright 2020, NORCE AS

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
#include "config.h"
#include <exception>
#include <ebos/eclproblem.hh>
#include <ebos/eclnewtonmethod.hh>
#include <ebos/ebos.hh>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm{
    template<typename TypeTag>
    class MonitoringAuxModule : public BaseAuxiliaryModule<TypeTag>
    {      
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    public:
        void postSolve(GlobalEqVector& deltaX){
            std::cout << "Dummy PostSolve Aux" << std::endl;
        }
        void addNeighbors(std::vector<NeighborSet>& neighbors) const{};
        void applyInitial(){};
        unsigned numDofs() const{return 0;};
        void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual){
            std::cout << "Dummy Linearize Aux" << std::endl;
        };
    };

    
    template<typename TypeTag>
    class EclProblemFlow: public EclProblem<TypeTag>{
    public:
        using Parent = EclProblem<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper =  AdaptiveTimeSteppingEbos<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        EclProblemFlow(Simulator& simulator): EclProblem<TypeTag>(simulator){
        }
        void timeIntegration()
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start TimeIntegration-------------------\n"
                << std::flush;                                                     
            }
            Parent::timeIntegration();
        }
        void beginTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start beginTimeStep-------------------\n"
                << std::flush;                                                     
            }
            Parent::beginTimeStep();
        }
        void endTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;                                                     
            }
            Parent::endTimeStep();
        }
        void finishInit(){
            Parent::finishInit();
            this->model().addAuxiliaryModule(&monitorAux_);
            //this->model().addOutputModule(new VtkMechModule<TypeTag>(simulator));
        }
    private:
        using MonitorAuxType = MonitoringAuxModule<TypeTag>;
        MonitorAuxType monitorAux_;    
    
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };

    template<typename TypeTag>
    class BlackoilWellModelFvExtra: public BlackoilWellModel<TypeTag>{
        using Parent = BlackoilWellModel<TypeTag>;
    public:
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        BlackoilWellModelFvExtra(Simulator& ebosSimulator): Parent(ebosSimulator) {};
        void beginIteration(){
            Parent::beginIteration();
            std::cout << "EclWellModelFvExtra begin iteration" << std::endl;
        }
        void endIteration(){
            Parent::endIteration();
            std::cout << "EclWellModelFvExtra end iteration" << std::endl;
        }
    };
        
    template<typename TypeTag>
    class BlackOilModelFvLocal: public BlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>; 
    public:
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        BlackOilModelFvLocal(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }
        const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            const auto& primaryVars = this->solution(timeIdx);
            const auto& problem = this->simulator_.problem();
            const auto intquant = this->cachedIntensiveQuantities(globalIdx, timeIdx);
            if (!this->enableIntensiveQuantityCache_){
                OPM_THROW(std::logic_error, "Run without intentive quantites not enabled: Use --enable-intensive-quantity=true");
            }
            if(!intquant){
                OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
            }    
            return *intquant;    
        }
//         void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx){
//             std::cout << "----------------------Update quantities-------------------\n"
//                 << std::flush;
// //            Parent::invalidateAndUpdateIntensiveQuantities(timeIdx);
// //             Parent::invalidateAndUpdateIntensiveQuantitiesSimple(*this,solution,/*timeIdx*/0);
//             const auto& primaryVars = this->solution(timeIdx);
//             const auto& problem = this->simulator_.problem();
//             size_t numGridDof = primaryVars.size();
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif                
//             for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
//                 const auto& primaryVar = primaryVars[dofIdx];
//                 auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
//                 intquant.update(problem, primaryVar, dofIdx, timeIdx);
//             }
            
//             std::fill(this->intensiveQuantityCacheUpToDate_[timeIdx].begin(),
//                       this->intensiveQuantityCacheUpToDate_[timeIdx].end(),
//                       /*value=*/true);
            
//         }        
    };    
    
} // namespace Opm

// the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemMech {    
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemMech> {
    using type = EclProblemFlow<TypeTag>;
};

template<class TypeTag>
struct Model<TypeTag, TTag::EclFlowProblemMech> {
    using type = BlackOilModelFvLocal<TypeTag>;
};    

template<class TypeTag>
struct EclWellModel<TypeTag, TTag::EclFlowProblemMech> {
    using type = BlackoilWellModelFvExtra<TypeTag>;
};
    
// template<class TypeTag>
// struct NewtonMethod<TypeTag, TTag::EclFlowProblemMech> {
//     using type = EclNewtonMethodLinesearch<TypeTag>;
// };    
    
template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr int value = 1;
};

template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EclNewtonSumTolerance<TypeTag, TTag::EclFlowProblemMech> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-5;
};
    
// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EclFlowProblemMech> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class TypeTag>
struct EclNewtonRelaxedVolumeFraction<TypeTag, TTag::EclFlowProblemMech> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, TTag::EclFlowProblemMech> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10*getPropValue<TypeTag, Properties::NewtonTolerance>();
};

// template<class TypeTag>
// struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemMech> {
//     //using type = EclBlackOilIntensiveQuantities<TypeTag>;
//     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
//     //using type = BlackOilIntensiveQuantities<TypeTag>;
//     //using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
// };

template<class TypeTag>
struct Linearizer<TypeTag, TTag::EclFlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::EclFlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
    
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemMech> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemMech> { static constexpr bool value = false; };    

//static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();    
    
template<class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemMech> { using type = Opm::Simulator<TypeTag>; };    
    
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemMech> {
     static constexpr bool value = false;
};
}
}
int main(int argc, char** argv)
{
    
    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemMech;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
    //return Opm::start<TypeTag>(argc, argv);
}
