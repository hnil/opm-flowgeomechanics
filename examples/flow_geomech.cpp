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
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/MechGridPartition.hpp>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <opm/grid/polyhedralgrid.hh>
#ifdef HAVE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <opm/simulators/flow/AluGridVanguard.hpp>
#endif
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/PolyhedralGridVanguard.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers_impl.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>
#include "MechTypeTag.hpp"
namespace Opm
{
namespace Properties
{
    namespace TTag
    {
        struct FlowProblemMechNoTemp
        {
            using InheritsFrom = std::tuple<FlowProblemMech>;
        };
    } // namespace TTag
    template <class TypeTag>
    struct EnableEnergy<TypeTag, TTag::FlowProblemMechNoTemp>
    {
        static constexpr bool value = false;
    };
    template<class TypeTag>
        struct EnergyModuleType<TypeTag, TTag::FlowProblemMechNoTemp>
        { static constexpr EnergyModules value = EnergyModules::NoTemperature; };

} // namespace Properties
} // namespace Opm


int
main(int argc, char** argv)
{
    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::FlowProblemMechNoTemp;
    // A mechanics grid of its own decides the flow partition too.
    Opm::MechPartition::installFlowPartition<TypeTag>(argc, argv);
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
