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
#if USE_TRACY
#define DETAILED_PROFILING 1
#endif

#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>
#include <opm/grid/polyhedralgrid.hh>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/PolyhedralGridVanguard.hpp>

#include <opm/simulators/aquifers/SupportsFaceTag.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>
#include <opm/simulators/flow/CollectDataOnIORank_impl.hpp>
#include <opm/simulators/flow/EclGenericWriter_impl.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>
#include <opm/simulators/flow/GenericThresholdPressure_impl.hpp>
#include <opm/simulators/flow/GenericTracerModel_impl.hpp>
#include <opm/simulators/flow/Transmissibility_impl.hpp>
#include <opm/simulators/utils/GridDataOutput_impl.hpp>

// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemMech {
    using InheritsFrom = std::tuple<FlowProblem,VtkGeoMech,FlowGeomechIstlSolverParams>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemMech> {
    using type = EclProblemGeoMech<TypeTag>;
};

// template<class TypeTag>
// struct Model<TypeTag, TTag::EclFlowProblemMech> {
//     using type = BlackOilModelFvLocal<TypeTag>;
// };


// template<class TypeTag>
// struct EclWellModel<TypeTag, TTag::EclFlowProblemMech> {
//     using type = BlackoilWellModelFvExtra<TypeTag>;
// };

// template<class TypeTag>
// struct NewtonMethod<TypeTag, TTag::EclFlowProblemMech> {
//     using type = EclNewtonMethodLinesearch<TypeTag>;
// };
template<class TypeTag>
struct EnableMech<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = true;
};


template<class TypeTag>
struct VtkWriteMoleFractions<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableVtkOutput<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableOpmRstFile<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr int value = 1;
};

template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, TTag::EclFlowProblemMech> {
    static constexpr bool value = false;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EclFlowProblemMech> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};

// template<class TypeTag>
// struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemMech> {
//     //using type = EclBlackOilIntensiveQuantities<TypeTag>;
//     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
//     //using type = BlackOilIntensiveQuantities<TypeTag>;
//     //using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
// };

// template<class TypeTag>
// struct Linearizer<TypeTag, TTag::EclFlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

// template<class TypeTag>
// struct LocalResidual<TypeTag, TTag::EclFlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemMech> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemMech> { static constexpr bool value = false; };

//static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();

template<class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemMech> { using type = Opm::Simulator<TypeTag>; };

// set grid to polygrid
    template<class TypeTag>
    struct Grid<TypeTag, TTag::EclFlowProblemMech> {
        using type = Dune::PolyhedralGrid<3, 3>;
    };
    template<class TypeTag>
    struct EquilGrid<TypeTag, TTag::EclFlowProblemMech> {
        //using type = Dune::CpGrid;
        using type = GetPropType<TypeTag, Properties::Grid>;
    };

    template<class TypeTag>
    struct Vanguard<TypeTag, TTag::EclFlowProblemMech> {
        using type = Opm::PolyhedralGridVanguard<TypeTag>;
    };

}

template<>
class SupportsFaceTag<Dune::PolyhedralGrid<3, 3>>
    : public std::bool_constant<true>
{};

}

int main(int argc, char** argv)
{

    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemMech;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
    //return Opm::start<TypeTag>(argc, argv);
}
