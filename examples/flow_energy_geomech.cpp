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

namespace Opm
{
template <typename TypeTag>
class BlackoilGeomechWellModel : public BlackoilWellModel<TypeTag>
{
    using Parent = BlackoilWellModel<TypeTag>;
    using Simulator = typename Parent::Simulator;
public:
    BlackoilGeomechWellModel(Simulator& simulator):
    Parent(simulator)
    {

    };
    //using BlackoilWellModel::BlackoilWellModel;

    using NeighborSet = typename Parent::NeighborSet;
    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        if (!this->param_.matrix_add_well_contributions_) {
            return;
        }
        OPM_THROW(std::runtime_error, "Not implemented");
    };
    void createWellContainer(const int reportStepIdx)
    {
        Parent::createWellContainer(reportStepIdx);
        // only add effect of fracture after one report step
        // NB everything is not explicit and ministeps are not considered
        if (reportStepIdx > 0) {
            const auto& problem = this->simulator_.problem();
            if (problem.hasFractures()) {
                for (auto& wellPtr : this->well_container_) {
                    auto wellName = wellPtr->name();

                    const auto& geomechmodel = problem.geomechModel();
                    const auto& fracturemodel = geomechmodel.fractureModel();
                    auto wellIndices = fracturemodel.getExtraWellIndices(wellName);
                    wellPtr->addPerforations(wellIndices);
                }
            }
        }
    };
};
} // namespace Opm

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

        template <class TypeTag>
        struct WellModel<TypeTag, TTag::EclFlowProblemMechTemp> {
            using type = BlackoilGeomechWellModel<TypeTag>;
        };
    }
}

int main(int argc, char** argv)
{

    OPM_TIMEBLOCK(fullSimulation);
    Opm::Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
    Opm::Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemMechTemp;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
    //return Opm::start<TypeTag>(argc, argv);
}
