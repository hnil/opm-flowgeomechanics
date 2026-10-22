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
#include <cstdlib>
#include <iomanip>
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

    // Opt-in (OPM_GEOMECH_MASS_TRACE=1): per well, the booked surface water
    // rate against the water actually sent to the reservoir residual, to
    // localise an unbooked source.
    void traceWellSources() const
    {
        static const bool on = (std::getenv("OPM_GEOMECH_MASS_TRACE") != nullptr);
        if (!on) {
            return;
        }
        using FluidSystem = typename Parent::FluidSystem;
        const unsigned wcomp = FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx);
        // well-state arrays are ordered by phase position, not component index
        // well-state arrays are indexed by active phase position
        const int wpos = this->phaseUsage().canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
        for (const auto& well : this->well_container_) {
            const auto& ws = this->wellState().well(well->indexOfWell());
            double booked = ws.surface_rates[wpos] * 86400.0;
            double perfsum = 0.0;
            const int np = static_cast<int>(ws.perf_data.phase_rates.size() / std::max<std::size_t>(1, ws.perf_data.cell_index.size()));
            for (std::size_t p = 0; p < ws.perf_data.cell_index.size(); ++p) {
                perfsum += ws.perf_data.phase_rates[p * np + wpos] * 86400.0;
            }
            double residual = 0.0;
            for (const auto& cr : well->connectionRates()) {
                residual += cr[wcomp].value() * 86400.0;
            }
            std::stringstream os;
            os << "WELLTRACE " << well->name() << " t=" << this->simulator_.time()
               << " booked_wat=" << std::setprecision(12) << booked
               << " perfsum_wat=" << perfsum << " residual_wat=" << residual
               << " nperf=" << ws.perf_data.cell_index.size();
            OpmLog::info(os.str());
        }
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
