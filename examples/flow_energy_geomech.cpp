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

#include <exception>
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>
#include "MechTypeTag.hpp"
#include <opm/geomech/BlackoilGeomechWellModel.hpp>
// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm
{
    namespace Properties
    {
        namespace TTag
        {
            struct EclFlowProblemMechTemp {
                using InheritsFrom = std::tuple<EclFlowProblemMech>;
            };
        }

        // template <class TypeTag>
        // struct WellModel<TypeTag, TTag::EclFlowProblemMechTemp> {
        //      using type = BlackoilGeomechWellModel<TypeTag>;
        // };
    }
}

int main(int argc, char** argv)
{

    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemMechTemp;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
    //return Opm::start<TypeTag>(argc, argv);
}
