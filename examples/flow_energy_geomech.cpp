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
#endif

#include <exception>
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include "MechTypeTag.hpp"
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/MechGridPartition.hpp>
#include <opm/simulators/flow/flux/FluxActivation.hpp>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
namespace Opm
{
namespace Properties
{
    namespace TTag
    {
        struct FlowProblemMechTemp
        {
            using InheritsFrom = std::tuple<FlowProblemMech>;
        };
    } // namespace TTag

} // namespace Properties
} // namespace Opm

namespace {
// Main::runStatic, with the schedule's wells handed to the mechanics sector
// partition, which runs while the vanguard is built.
struct GeoMechMain : public Opm::Main
{
    using Opm::Main::Main;

    template <class TypeTag>
    int run()
    {
        int exitCode = EXIT_SUCCESS;
        if (this->initialize_<TypeTag>(exitCode)) {
            this->setupVanguard();
            if (const auto& schedule = Opm::FlowGenericVanguard::modelParams_.eclSchedule_) {
                Opm::MechPartition::SectorWells::capture(*schedule);
            }
            return Opm::flowMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
        }
        return exitCode;
    }
};
} // namespace

int
main(int argc, char** argv)
{

    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::FlowProblemMechTemp;
    // A mechanics grid of its own decides the flow partition too.
    Opm::MechPartition::installFlowPartition<TypeTag>(argc, argv);
    // The mechanics grid may extend past a USEFLUX sector.
    Opm::retainBeforeUseFlux({"YMODULE", "PRATIO"});
    auto mainObject = GeoMechMain(argc, argv);
    return mainObject.run<TypeTag>();
}
