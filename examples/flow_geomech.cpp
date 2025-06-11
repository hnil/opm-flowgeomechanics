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
#define DETAILED_PROFILING 0
#endif

#include <exception>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>

#include <opm/grid/polyhedralgrid.hh>
#ifdef HAVE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <opm/simulators/flow/AluGridVanguard.hpp>
#endif
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/PolyhedralGridVanguard.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers_impl.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>
//#include <ebos/eclpolyhedralgridvanguard.hh>
// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
#include "MechTypeTag.hpp"
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemMechNoTemp {
    using InheritsFrom = std::tuple<EclFlowProblemMech>;
};
}
        template <class TypeTag>
        struct EnableEnergy<TypeTag, TTag::EclFlowProblemMechNoTemp> {
            static constexpr bool value = false;
        };


}
}


int main(int argc, char** argv)
{
    // Opm::Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
    // Opm::Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
    // Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
    // Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
    // Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemMechNoTemp;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
    //return Opm::start<TypeTag>(argc, argv);
}
