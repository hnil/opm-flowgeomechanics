/*
  Copyright 2025 Equinor ASA.

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

#pragma once
#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <sstream>
namespace Opm
{
template <typename TypeTag>
class BlackoilGeoMechWellModel : public BlackoilWellModel<TypeTag>
{
    using Parent = BlackoilWellModel<TypeTag>;
    using Simulator = typename Parent::Simulator;
public:
    BlackoilGeoMechWellModel(Simulator& simulator,
                             const NewtonIterationContext& iter_ctx):
    Parent(simulator, iter_ctx)
    {

    };
    //using BlackoilWellModel::BlackoilWellModel;

    using NeighborSet = typename Parent::NeighborSet;
    // Wells-in-matrix (--matrix-add-well-contributions): the parent adds the
    // clique of every schedule connection; auxiliary (fracture-DOF)
    // perforations are not schedule connections, so their clique is added here.
    void addNeighbors(std::vector<NeighborSet>& neighbors) const override
    {
        if (!this->param_.matrix_add_well_contributions_) {
            return;
        }
        Parent::addNeighbors(neighbors);
        for (const auto& [wname, auxPerfs] : this->auxiliary_perforations_) {
            if (auxPerfs.empty() || !this->schedule().hasWell(wname)) {
                continue;
            }
            auto cells = this->getCellsForConnections(this->schedule().back().wells(wname));
            for (const auto& aux : auxPerfs) {
                cells.push_back(aux.cell_index);
            }
            for (const int c : cells) {
                neighbors[c].insert(cells.begin(), cells.end());
            }
        }
    }
    // Fracture-created connections are added to the schedule at runtime and
    // may not (yet) exist on every rank that hosts part of the well; log and
    // continue instead of aborting the run.
    bool continueOnMissingWellConnections() const override
    {
        return true;
    }

    void createWellContainer(const int reportStepIdx)
    {
        Parent::createWellContainer(reportStepIdx);
        // only add effect of fracture after one report step
        // NB everything is not explicit and ministeps are not considered
        if (reportStepIdx > 0) {
            const auto& problem = this->simulator_.problem();
            const auto& fractureHost = problem.fractureHost();
            if (problem.hasFractures() && fractureHost.fractureModelActive()) {
                for (auto& wellPtr : this->well_container_) {
                    auto wellName = wellPtr->name();
                    const auto& fracturemodel = fractureHost.fractureModel();
                    auto wellIndices = fracturemodel.getExtraWellIndices(wellName);
                    wellPtr->addFracturePerforations(wellIndices);
                }
            }
        }
    };
};
} // namespace Opm
