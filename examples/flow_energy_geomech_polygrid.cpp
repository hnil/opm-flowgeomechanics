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
//#if HAVE_TRACY
#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"
#define OPM_TIMEBLOCK(blockname) ZoneNamedN(blockname, #blockname, true);
//#define OPM_TIMEBLOCK_LOCAL(blockname) ZoneNamedN(blockname, #blockname, true);
//#endif

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
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>

#include <ebos/eclpolyhedralgridvanguard.hh>

#include <opm/grid/polyhedralgrid.hh>
#include <ebos/collecttoiorank_impl.hh>
#include <ebos/eclgenericproblem_impl.hh>
#include <ebos/eclgenericthresholdpressure_impl.hh>
#include <ebos/eclgenerictracermodel_impl.hh>
#include <ebos/ecltransmissibility_impl.hh>
#include <ebos/eclgenericwriter_impl.hh>
#include <ebos/equil/initstateequil_impl.hh>

// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemMech {
    using InheritsFrom = std::tuple<EclFlowProblem,VtkGeoMech,FlowGeomechIstlSolverParams>;
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

template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemMech> {
     static constexpr bool value = false;
};
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
        using type = Opm::EclPolyhedralGridVanguard<TypeTag>;
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
