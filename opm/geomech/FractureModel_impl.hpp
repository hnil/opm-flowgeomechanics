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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/PerforationData.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include <fmt/format.h>

namespace Opm {
    inline void
    FractureModel::updateFractureReservoirCells(const Dune::CpGrid& cpgrid)
    {
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cell_search_tree_, cpgrid, cell_seeds_);
            }
        }
    }

    template<class Grid>
    FractureModel::FractureModel(const Grid&              grid,
                                 const std::vector<Well>& wells,
                                 const PropertyTree&      param,
                                 const ConnCellResolver&  connCellResolver)
        : prm_(param)
    {
        GeometryHelper geomhelp(grid);
	// NB: need to be carefull in parallel
        vtkwritewells_= prm_.get<bool>("vtkwritewells", false);

        for (const auto& well : wells) {
            std::vector<Point3D> vertices;

            if (true) {
                auto conns = std::vector<FractureWell::Connection>{};

                // should be done property with welltraj
                for (const auto& connection : well.getConnections()) {
                    const int cell_index = connCellResolver
                        ? connCellResolver(well, connection)
                        : geomhelp.compressedIndex(connection.global_index());

                    if (cell_index < 0) {
                        // Not an interior cell on this rank (parallel), or
                        // unresolvable connection: skip.
                        continue;
                    }

                    vertices.push_back(geomhelp.centroid(cell_index));

                    auto& conn = conns.emplace_back();
                    conn.cell = cell_index;
                    conn.segment = connection.segment();
                    conn.perf_range = connection.perf_range();
                }

                if (conns.empty()) {
                    OpmLog::warning(fmt::format("No connections found for well {}", well.name()));
                    continue;
                }

                Point3D refpoint(vertices[0]);
                if(well.hasRefDepth()){
                    if(std::abs(well.getRefDepth()-refpoint[2])< 10){
                        refpoint[2] = well.getRefDepth()-10;//avoid zero
                    }else{
                        refpoint[2] = well.getRefDepth();
                    }
                }else{
                    refpoint[2] -= 10;
                }
                vertices.insert(vertices.begin(),refpoint);

                std::vector<Segment> segments;
                //NB unsinged so grid can not be to big
                for (unsigned i=0; i < unsigned(conns.size()); ++i) {
                    segments.push_back(Segment({i,i+1}));
                }

                // NB should add gri cells
                this->addWell(well.name(), conns, vertices, segments);
            }
        }

        external::buildBoundingBoxTree(cell_search_tree_, cell_seeds_, grid);
    }

    template<class wseed_collection>
    void
    FractureModel::updateActive(const wseed_collection& current_wseed)
    {
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.setActive(false);
            }
        }
        if (current_wseed().empty()) {
            return;
        }
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                auto wellInfo = fracture.wellInfo();
                if (current_wseed.has(wellInfo.name)) {
                    const auto& well_wseed = current_wseed(wellInfo.name);
                    for (const auto& seedcell : well_wseed.seedCells()) {
                        if (seedcell == wellInfo.global_index) {
                            fracture.setActive(true);
                            OpmLog::info(fmt::format("Activating fracture {}", fracture.name()));
                        }
                    }
                }
            }
        }
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::solve(const Simulator& simulator)
    {
        last_solve_stats_ = {};
        // Opt-in: outside the timestep retry scope (endTimeStep) a fracture solve
        // failure is fatal; "warn_keep_state" logs and continues on whatever state
        // the solve left. Serial-tested only.
        const bool warn_keep_state =
            prm_.get<std::string>("solver.solve_failure_policy", "throw") == "warn_keep_state";
        for (size_t i = 0; i < wells_.size(); ++i) {
            std::vector<Fracture>& fractures = well_fractures_[i];
            for (size_t j = 0; j < fractures.size(); ++j) {
                if (fractures[j].isActive()) {
                    OpmLog::debug(fmt::format("Solving fracture {}", fractures[j].name()));
                    if (warn_keep_state) {
                        try {
                            fractures[j].solve<TypeTag, Simulator>(cell_search_tree_, cell_seeds_, simulator);
                        } catch (const std::exception& e) {
                            OpmLog::warning("Fracture solve FAILED for " + fractures[j].name()
                                            + "; solve_failure_policy=warn_keep_state, continuing: "
                                            + e.what());
                        }
                    } else {
                        fractures[j].solve<TypeTag, Simulator>(cell_search_tree_, cell_seeds_, simulator);
                    }
                    last_solve_stats_ += fractures[j].lastSolveStats();
                    fractures[j].warnIfOpenButNotConducting();
                }
            }
        }
        total_solve_stats_ += last_solve_stats_;
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::updateFilterCakeProperties(const Simulator& simulator)
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            std::vector<Fracture>& fractures = well_fractures_[i];
            for (size_t j = 0; j < fractures.size(); ++j) {
                fractures[j].updateFilterCakePropertiesPost<TypeTag, Simulator>(simulator);
            }
        }
    }

    template <class TypeTag, class Simulator>
    double
    FractureModel::maxFlowTimeStep(const Simulator& simulator) const
    {
        double dt_min = std::numeric_limits<double>::max();
        for (size_t i = 0; i < wells_.size(); ++i) {
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for (size_t j = 0; j < fractures.size(); ++j) {
                if (fractures[j].isActive()) {
                    double dt_max = fractures[j].maxFlowTimeStep();
                    if (dt_max > 0) {
                        dt_min = std::min(dt_max, dt_min);
                    }
                }
            }
        }
        dt_min = simulator.vanguard().grid().comm().min(dt_min);
        return dt_min;
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::initReservoirProperties(const Simulator& simulator)
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]) {
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator, true);
            }
        }
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::updateReservoirAndWellProperties(const Simulator& simulator)
    {
        this->updateReservoirProperties<TypeTag, Simulator>(simulator);
        this->updateWellProperties<TypeTag, Simulator>(simulator);
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::resetFractures(const Simulator& simulator)
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]) {
                fracture.resetFracture();
            }
        }
        const auto& cpgrid = simulator.vanguard().grid();
        updateFractureReservoirCells(cpgrid);
        this->initReservoirProperties<TypeTag, Simulator>(simulator);
        this->updateReservoirProperties<TypeTag, Simulator>(simulator);
        this->updateWellProperties<TypeTag, Simulator>(simulator);
    }

    inline void
    FractureModel::writeIterationSnapshots(int step, int round, const std::string& tag) const
    {
        for (size_t i = 0; i < wells_.size(); ++i)
            for (const auto& fracture : well_fractures_[i])
                fracture.writeIterationSnapshot(step, round, tag);
    }

    inline bool
    FractureModel::anyLastSolveUnconverged() const
    {
        for (size_t i = 0; i < wells_.size(); ++i)
            for (const auto& fracture : well_fractures_[i])
                if (fracture.isActive() && !fracture.lastSolveStats().converged)
                    return true;
        return false;
    }

    inline std::string
    FractureModel::growthGuardViolation() const
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            for (const auto& fracture : well_fractures_[i]) {
                const std::string msg = fracture.growthGuardViolation();
                if (!msg.empty()) return msg;
            }
        }
        return {};
    }

    inline void
    FractureModel::moveForwardInTime(double dt_last)
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]) {
                fracture.moveForwardInTime(dt_last);
            }
        }
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::updateReservoirWellProperties(const Simulator& simulator)
    {
        for (auto& well : wells_) {
            well.updateReservoirProperties<TypeTag, Simulator>(simulator);
        }
    }

    inline bool
    FractureModel::addPertsToSchedule()
    {
        return prm_.get<bool>("addperfs_to_schedule");
    }

    inline Opm::PropertyTree&
    FractureModel::getParam()
    {
        return prm_;
    }

    inline const FractureSolveStats&
    FractureModel::lastSolveStats() const
    {
        return last_solve_stats_;
    }

    inline const FractureSolveStats&
    FractureModel::totalSolveStats() const
    {
        return total_solve_stats_;
    }

    template <class TypeTag, class Simulator>
    void
    FractureModel::updateReservoirProperties(const Simulator& simulator)
    {
        for (size_t i = 0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]) {
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator);
            }
        }
    }

  template<class WellState>
    int findPerf(const WellState& wellstate, int cell_index_frac){
        const auto& perf_data = wellstate.perf_data;
        auto it = std::find(perf_data.cell_index.begin(),
                            perf_data.cell_index.end(),
                            cell_index_frac);
                            // check if perforation exists
        if(it == perf_data.cell_index.end()) {
            OpmLog::warning(fmt::format("Could not find perforation for well {} in cell index {}",
                                        wellstate.name, cell_index_frac));
            return -1;
        }// skip if not found
                            
        int perf_index = it - perf_data.cell_index.begin();
        return perf_index;
    }

    template <class TypeTag, class Simulator>
    void FractureModel::updateWellProperties(const Simulator& simulator)
    {   
       int verbosity = prm_.get<int>("verbosity",0);
       using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        for (size_t i=0; i < wells_.size(); ++i) {
            // TO DO set wells to active even without fractures
            double injection_rate = 0.0;
            bool well_control_is_rate = true; // active well constraint (for control.type=="well")
            double well_bhp = -1.0;           // active BHP (bhp_well constraint target)
            std::vector<int> perf_cell_indices;
            double well_depth = 0.0;
            double total_wellindex = 0.0;
            int n_aux_perfs = 0;
            double aux_wellindex = 0.0;
            double wi_dz = 0.0;
            double wi_respress = 0.0;
            std::vector<double> perf_depths;
            {
              auto well_index = simulator.problem().wellModel().wellState().index(wells_[i].name());
                if(well_index.has_value()){
                    const auto& wellmodel = simulator.problem().wellModel();
                    const auto& wellstate = wellmodel.wellState().well(*well_index);
                    // get cell index from well perforation data
                    // check if well is open
                    if(wellstate.status == ::Opm::WellStatus::OPEN) {
                      // the ACTIVE constraint decides the fracture BC when
                      // control.type == "well": pressure-constrained -> perf_pressure,
                      // anything rate-like (RATE/RESV/GRUP/...) -> rate_well
                      bool is_thp = false;
                      if (wellstate.producer) {
                          const auto c = wellstate.production_cmode;
                          is_thp = (c == WellProducerCMode::THP);
                          well_control_is_rate = !(c == WellProducerCMode::BHP || is_thp);
                      } else {
                          const auto c = wellstate.injection_cmode;
                          is_thp = (c == WellInjectorCMode::THP);
                          well_control_is_rate = !(c == WellInjectorCMode::BHP || is_thp);
                      }
                      // THP is not supported by the automatic BC selection yet
                      // (needs the generalized well-constraint row with a VFP
                      // linearization) - fail loudly rather than mis-couple.
                      well_bhp = wellstate.bhp;
                      if (is_thp && prm_.get<std::string>("control.type", "") == "well")
                          OPM_THROW(std::runtime_error,
                                    "control.type='well': THP-constrained well '"
                                        + wells_[i].name()
                                        + "' is not supported - set an explicit fracture control");
                      // wellState() indexes all wells, localNonshutWells() only
                      // the non-shut ones - with shut wells in the deck the two
                      // indices diverge, so look up by name.
                      const auto& wells = wellmodel.localNonshutWells();
                      const auto well_it = std::find_if(wells.begin(), wells.end(),
                          [&](const auto& w) { return w->name() == wells_[i].name(); });
                      if (well_it != wells.end()) {
                      const auto& well = *well_it;

                      for (int perf_index=0; perf_index < wellstate.perf_data.cell_index.size(); ++perf_index) {    
                            const int cell_idx = wellstate.perf_data.cell_index[perf_index];
                            // Auxiliary perforations are the embedded fracture cells.
                            // The fracture's own solve already carries that same
                            // connection in perfinj_, so folding them in here counts
                            // the fracture twice -- and at fractureWI scale they swamp
                            // the matrix connections the reference index is meant to be.
                            if (perf_index < static_cast<int>(wellstate.perf_data.ecl_index.size())
                                && wellstate.perf_data.ecl_index[perf_index]
                                    == AUX_PERFORATION_ECL_INDEX) {
                                ++n_aux_perfs;
                                aux_wellindex += well->wellIndex()[perf_index];
                                continue;
                            }
                            perf_cell_indices.push_back(cell_idx);
                            const auto& intQuants = simulator.model()
                                .intensiveQuantities(cell_idx, /*timeIdx=*/0);
                            using Scalar = double;
                            auto obtain = [](const auto& value) { return getValue(value); };
                            const auto trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(intQuants, cell_idx, obtain);
                            // Matrix-only connection factor: upstream getTw() would add the
                            // fracture contribution. The WINJDAM cake multiplier is applied to
                            // Tw in the well model, not via wellTransMultiplier; without it the
                            // duplicated well row sees caked perfs as open.
                            double fc_mult = 1.0;
                            const auto& fc = well->filterCakeMultipliers();
                            if (!fc.empty() && static_cast<size_t>(perf_index) < fc.size()) {
                                const auto perf_ecl_index = well->perforationData()[perf_index].ecl_index;
                                const auto& connection = well->wellEcl().getConnections()[perf_ecl_index];
                                if (connection.filterCakeActive())
                                    fc_mult = fc[perf_index];
                            }
                            const auto effective_well_index_perf =
                                well->wellIndex()[perf_index] * trans_mult * fc_mult;

                                                 
                            double lambda = 0.0;
                            int numPhases = 3;
                            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                                            // assume sum should only be water;
                                            auto val = intQuants.mobility(phaseIdx);
                                            lambda += val.value();
                                    }
                            }   
                            const auto effective_well_index   =  effective_well_index_perf*lambda;
                            double density = intQuants.fluidState().density(FluidSystem::waterPhaseIdx).value();
                            total_wellindex += effective_well_index;
                            double dzwell = (well->perfDepth()[perf_index]- well->refDepth());  
                            wi_dz += effective_well_index*dzwell;
                            double pressure = intQuants.fluidState().pressure(FluidSystem::waterPhaseIdx).value();
                            wi_respress += effective_well_index*pressure;

                        }
                        injection_rate = wellstate.reservoir_rates[FluidSystem::waterPhaseIdx];
                        well_depth = well->refDepth();
                        perf_depths = well->perfDepth();
                      } // well found in localNonshutWells
                    }
                 }
            }
            assert(well_fractures_[i].size() < 2);// for now asser this if not better "well model is need"
            for (auto& fracture : well_fractures_[i]){
              WellInfo wellinfo = fracture.wellInfo();
              if(verbosity > 1){
                std::stringstream os;
                os << " Well " << wellinfo.name << " injection "
                        << injection_rate << " WI " << total_wellindex
                        << " (matrix perforations only; " << n_aux_perfs
                        << " auxiliary fracture perforations excluded, their CTF sum "
                        << aux_wellindex << ")" << std::endl;
                OpmLog::info(os.str());
              }  
              fracture.setTimeStep(simulator.timeStepSize());
              fracture.setWellControl(well_control_is_rate, well_bhp);
              // Consistency check: the fracture's duplicated well row should
              // track the reservoir well model. A large BHP gap means the row's
              // frozen WI/mobility no longer represents the well (the missing
              // filter-cake multiplier hid exactly this way). Warn only.
              {
                  const double warn_bar =
                      prm_.get<double>("solver.well_consistency_warn_bhp", 50.0);
                  const std::string ctrl = fracture.effectiveControlType();
                  if (warn_bar > 0.0 && well_bhp > 0.0
                      && fracture.hasPressureState()
                      && (ctrl == "rate_well" || ctrl == "bhp_well")) {
                      const double frac_bhp = fracture.injectionBhp();
                      if (std::isfinite(frac_bhp) && frac_bhp > 1.0e5
                          && std::abs(frac_bhp - well_bhp) > warn_bar * 1.0e5) {
                          OpmLog::warning(
                              "Fracture " + fracture.name() + ": well DOF "
                              + std::to_string(frac_bhp / 1.0e5)
                              + " bar deviates from the reservoir well model BHP "
                              + std::to_string(well_bhp / 1.0e5)
                              + " bar by more than "
                              + std::to_string(warn_bar)
                              + " bar - the duplicated well row (frozen WI/mobility) "
                                "is out of sync with the well");
                      }
                  }
              }
              fracture.setWellProps(injection_rate,  total_wellindex,  wi_dz,  wi_respress,  well_depth);
              fracture.setWellPerfCells(perf_cell_indices);
                // do update wells
                // set well properties

                // need to be double checked how to assosiate correct perforation/segment
                int perf_index_frac = wellinfo.perf;
                int cell_index_frac = wellinfo.well_cell;
                // see if well exist
                fracture.setActive(false);
                auto well_index = simulator.problem().wellModel().wellState().index(wellinfo.name);
                if(!well_index.has_value()){
                    wells_[i].setActive(false);
                    fracture.setActive(false);
                    continue;
                }
                const auto& wellstate = simulator.problem().wellModel().wellState().well(*well_index);
                // check if well is open
                if(wellstate.status != ::Opm::WellStatus::OPEN) {
                    wells_[i].setActive(false);
                    OpmLog::debug(fmt::format("Well {} is not open, skipping update.", wellinfo.name));
                    fracture.setActive(false);
                    continue; // skip if not open
                }
                // well is active
                wells_[i].setActive(true);
                // get well perforation
                const auto& perf_data = wellstate.perf_data;
                int  perf_index =  findPerf(wellstate,cell_index_frac);
                if(perf_index == -1) {
                    fracture.setActive(false);
                    wells_[i].setPerfActive(perf_index_frac,false);          
                    continue; // skip if not found
                }
                double perf_pressure = perf_data.pressure[perf_index];
                double perf_rate = perf_data.rates[perf_index];
                if(verbosity > 1){
                    std::stringstream os;
                    os << "Perf index perf " << perf_index << " fracture " << perf_index_frac << " pressure " << perf_pressure << " rate " << perf_rate << std::endl;
                    OpmLog::info(os.str());
                }
                double perf_depth = perf_depths[perf_index];
                fracture.setPerfProps(perf_pressure, perf_depth, perf_rate);
                wells_[i].setPerfPressure(perf_index_frac, perf_pressure);
                fracture.setActive(true);
                wells_[i].setPerfActive(perf_index_frac,true);
                //NB do we need some rates? need to be summed over "peforations of the fractures"
            }
        }
    }
    

}
