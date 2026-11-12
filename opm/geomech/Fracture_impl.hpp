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

#include <opm/geomech/DiscreteDisplacement.hpp>
#include "RegularTrimesh.hpp"
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/vem/vemutils.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh>
#include <fmt/format.h>
namespace Opm
{

template <class State>
struct BlackOilFluidStateConfig;

template <class ScalarT,
          class FluidSystemT,
          bool enableTemperature,
          bool enableEnergy,
          bool enableDissolution,
          bool enableVapwat,
          bool enableBrine,
          bool enableSaltPrecipitation,
          bool enableDissolutionInWater,
          bool enableSolvent,
          unsigned numStoragePhases>
struct BlackOilFluidStateConfig<BlackOilFluidState<ScalarT,
                                                   FluidSystemT,
                                                   enableTemperature,
                                                   enableEnergy,
                                                   enableDissolution,
                                                   enableVapwat,
                                                   enableBrine,
                                                   enableSaltPrecipitation,
                                                   enableDissolutionInWater,
                                                   enableSolvent,
                                                   numStoragePhases>>
{
    static constexpr bool EnableTemperature = enableTemperature;
    static constexpr bool EnableEnergy = enableEnergy;
    static constexpr bool EnableDissolution = enableDissolution;
    static constexpr bool EnableVapwat = enableVapwat;
    static constexpr bool EnableBrine = enableBrine;
    static constexpr bool EnableSaltPrecipitation = enableSaltPrecipitation;
    static constexpr bool EnableDissolutionInWater = enableDissolutionInWater;
    static constexpr bool EnableSolvent = enableSolvent;
};

template <class TargetState, class SourceState, class ValueFactory>
void copyBlackOilFluidState(TargetState& target, const SourceState& source, ValueFactory&& valueFactory)
{
    using TargetConfig = BlackOilFluidStateConfig<TargetState>;
    using SourceConfig = BlackOilFluidStateConfig<SourceState>;

    target.setPvtRegionIndex(source.pvtRegionIndex());

    for (unsigned phaseIdx = 0; phaseIdx < TargetState::numPhases; ++phaseIdx) {
        if (!TargetState::FluidSystem::phaseIsActive(phaseIdx))
            continue;

        target.setSaturation(phaseIdx, valueFactory(source.saturation(phaseIdx)));
        target.setPressure(phaseIdx, valueFactory(source.pressure(phaseIdx)));
        target.setDensity(phaseIdx, valueFactory(source.density(phaseIdx)));
        target.setInvB(phaseIdx, valueFactory(source.invB(phaseIdx)));

        if constexpr (TargetConfig::EnableEnergy && SourceConfig::EnableEnergy)
            target.setEnthalpy(phaseIdx, valueFactory(source.enthalpy(phaseIdx)));
    }

    if constexpr (TargetConfig::EnableTemperature || TargetConfig::EnableEnergy)
        target.setTemperature(valueFactory(source.temperature(/*phaseIdx=*/0)));

    if constexpr (TargetConfig::EnableDissolution) {
        target.setRs(valueFactory(source.Rs()));
        target.setRv(valueFactory(source.Rv()));
    }

    if constexpr (TargetConfig::EnableVapwat)
        target.setRvw(valueFactory(source.Rvw()));

    if constexpr (TargetConfig::EnableDissolutionInWater)
        target.setRsw(valueFactory(source.Rsw()));

    if constexpr (TargetConfig::EnableBrine)
        target.setSaltConcentration(valueFactory(source.saltConcentration()));

    if constexpr (TargetConfig::EnableSaltPrecipitation)
        target.setSaltSaturation(valueFactory(source.saltSaturation()));
}

// ----------------------------------------------------------------------------

template <typename Vec>
std::vector<double>
make_vector(const Vec& data, size_t sz = 0)
{
    const size_t N = (sz == 0) ? data.size() : sz;
    std::vector<double> res(N);
    for (size_t i = 0; i != N; ++i)
        res[i] = data[i];
    return res;
}

inline WellInfo&
Fracture::wellInfo()
{
    return wellinfo_;
}

inline void
Fracture::setActive(bool active)
{
    active_ = active;
}

inline bool
Fracture::isActive() const
{
    return active_;
}

inline const FractureSolveStats&
Fracture::lastSolveStats() const
{
    return last_solve_stats_;
}

inline double
Fracture::maxFlowTimeStep() const
{
    // onset dt-hold: keep the step at solver.onset_dt_limit while the seed is
    // still closed and the well pressure is rising, so the seed is evaluated
    // near its opening pressure instead of after the well has run away
    const double onset_dt = prm_.get<double>("solver.onset_dt_limit", 0.0);
    double cap = max_flow_time_step_;
    if (coupling_dt_cap_ > 0.0) {
        cap = (cap > 0.0) ? std::min(cap, coupling_dt_cap_) : coupling_dt_cap_;
    }
    if (stress_dt_cap_ > 0.0) {
        cap = (cap > 0.0) ? std::min(cap, stress_dt_cap_) : stress_dt_cap_;
    }
    if (onset_dt > 0.0 && onsetHoldActive()) {
        return (cap > 0.0) ? std::min(cap, onset_dt * 86400.0) : onset_dt * 86400.0;
    }
    return cap;
}

inline bool
Fracture::onsetHoldActive() const
{
    if (!active_) return false;
    // no checkpoint yet = first active step: hold conservatively, otherwise the
    // first step can run far past the seed's opening pressure (a 2-day step-1
    // overshoot makes establishment a round-off coin flip)
    if (area_step_start_ < 0.0) return true;
    const double seed_area = prm_.get<double>("solver.onset_seed_area", 2.0); // m2
    if (area_step_start_ > seed_area) return false; // established: no hold
    const double rise = prm_.get<double>("solver.onset_pressure_rise", 1.0e5); // Pa
    return perf_pressure_ - perf_pressure_step_start_ > rise
        || perf_pressure_step_start_ < 0.0; // first step: hold
}

inline std::string
Fracture::growthGuardViolation() const
{
    const double max_fac = prm_.get<double>("solver.max_area_growth_per_step", 0.0);
    if (max_fac <= 0.0 || !active_ || area_step_start_ < 0.0) return {};
    const double area_now = calculateFractureProperties().area;
    const double abs_floor = prm_.get<double>("solver.min_area_growth_guard", 200.0); // m2
    const double allowed = std::max(max_fac * area_step_start_, area_step_start_ + abs_floor);
    if (area_now > allowed) {
        return "Fracture " + name() + " grew " + std::to_string(area_step_start_) + " -> "
            + std::to_string(area_now) + " m2 in one step (limit "
            + std::to_string(allowed) + "); rejecting step";
    }
    return {};
}

inline int
Fracture::numWellEquations() const
{
    const std::string t = effectiveControlType();
    return (t == "rate_well" || t == "bhp_well") ? 1 : 0;
}

inline Fracture::DynamicMatrix&
Fracture::fractureMatrix() const
{
    if (fracture_matrix_ == nullptr)
        assembleFractureMatrix();
    return *fracture_matrix_;
}

template <class TypeTag, class Simulator>
void Fracture::initReservoirProperties(const Simulator& simulator)
{
    const auto& problem = simulator.problem();
    size_t ncf = reservoir_cells_.size();
    reservoir_perm_.resize(ncf);
    reservoir_cstress_.resize(ncf);
    reservoir_density_.resize(ncf);
    double dist = prm_.get<double>("reservoir.dist");
    reservoir_dist_.resize(ncf, dist);
    double numax = -1e99;
    double Emax = -1e99;
    assert(ncf > 0);
    for (size_t i = 0; i < ncf; ++i) {
        int cell = reservoir_cells_[i];
        if (cell < 0) {
            reservoir_perm_[i] = 0.0;
            reservoir_cstress_[i] = 1e20;
            reservoir_density_[i] = 1000.0;
            reservoir_dist_[i] = dist;
            continue;
        }
        auto normal = this->cell_normals_[i];
        {
            auto permmat = problem.intrinsicPermeability(cell);
            auto cstress = problem.cStress(cell);
            auto np = normal;
            permmat.mv(normal, np);
            double value = np.dot(normal);
            reservoir_perm_[i] = value;
            reservoir_cstress_[i] = cstress;
        }
        Emax = std::max(Emax, problem.yModule(cell));
        numax = std::max(numax, problem.pRatio(cell));
    }
    E_ = Emax;
    nu_ = numax;
    assert(E_ > 0);
    assert(nu_ > 0 && nu_ < 1);
}

template <class TypeTag, class Simulator>
void Fracture::updateReservoirProperties(const Simulator& simulator, bool init_constant_vals,bool /*update_filtercake*/)
    {
        // if `init_contant_vals` is true, the fields that should normally not
        // change, i.e.  E_, nu_ and reservoir_perm_ will be updated. This should
        // normally only be needed the first time.
        if (init_constant_vals) {
            initReservoirProperties<TypeTag, Simulator>(simulator);
        }

        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using ScalarFluidState = Opm::BlackOilFluidState<double, FluidSystem>;
        const auto& problem = simulator.problem();
        const auto& grid = simulator.vanguard().grid();
        GeometryHelper ghelper(grid);
        const auto value_of = [](const auto& value) { return Opm::scalarValue(value); };
        // NB reservoir dist not calculated
        size_t ncf = reservoir_cells_.size();
        assert(ncf > 0);
        reservoir_pressure_.resize(ncf);
        reservoir_stress_.resize(ncf);
        reservoir_mobility_.resize(ncf);
        reservoir_perm_.resize(ncf);
        reservoir_cstress_.resize(ncf);
        reservoir_cell_z_.resize(ncf, 0.0);
        filtercake_thikness_.resize(ncf, 0.0);
        // should be calculated
        bool calculate_dist = prm_.get<bool>("reservoir.calculate_dist");
        const double dist_factor = leakoffDistanceFactor();
        if(!calculate_dist){
          double dist = prm_.get<double>("reservoir.dist");
          reservoir_dist_.resize(ncf, dist);
        }

        std::vector<ScalarFluidState> fracture_scalar_fluid_states(ncf);
        std::vector<bool> fracture_scalar_fluid_state_valid(ncf, false);
        {
        // well cell properties
        const auto& intQuants = simulator.model().intensiveQuantities(wellinfo_.well_cell, /*timeIdx*/ 0);
        const auto& fs = intQuants.fluidState();
        density_perf_ = fs.density(FluidSystem::waterPhaseIdx).value();
        mobility_water_perf_ = intQuants.mobility(FluidSystem::waterPhaseIdx).value();
        }




        // optionally evaluate the reservoir stress at the fracture cell centroids by
        // trilinear interpolation of patch-recovered nodal stress, instead of the
        // piecewise-constant stress of the containing cell.  With large reservoir
        // cells the latter gives large jumps in the traction along the fracture.
        const bool interpolate_stress = prm_.get<bool>("reservoir.interpolate_stress", false);
        std::vector<Dune::FieldVector<double,6>> nodal_stress;
        std::vector<Dune::FieldVector<double,3>> frac_centroids;
        const auto stress_mode = centre_stress::parseMode(
            prm_.get<std::string>("reservoir.stress_interpolation", "none"));
        if (stress_mode != centre_stress::Mode::None) {
            if (interpolate_stress) {
                OPM_THROW(std::runtime_error,
                          "reservoir.stress_interpolation and reservoir.interpolate_stress exclude each other");
            }
            frac_centroids.resize(ncf);
            ElementMapper fracElemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
            for (const auto& element : elements(grid_->leafGridView()))
                frac_centroids[fracElemMapper.index(element)] = element.geometry().center();
        }
        const std::size_t stress_sides_logged = stress_sides_.size();
        const auto centre_of_cell = [&grid](int c) {
            return Dune::cpgrid::Entity<0>(*grid.currentData().back(), c, true).geometry().center();
        };
        const auto cell_sides = [&](int c) -> const centre_stress::CellSides& {
            const auto hit = stress_sides_.find(c);
            if (hit != stress_sides_.end()) {
                return hit->second;
            }
            if (!stress_faults_) {
                stress_faults_ = centre_stress::DeclaredFaults(simulator.vanguard().eclState().getFaults());
            }
            const auto& mapper = simulator.vanguard().cartesianIndexMapper();
            const auto logical = [&mapper](int cc) {
                std::array<int, 3> ijk;
                mapper.cartesianCoordinate(cc, ijk);
                return ijk;
            };
            const auto cartesian = [&mapper](int cc) { return std::size_t(mapper.cartesianIndex(cc)); };
            const auto flows = [&problem](int a, int b) {
                const auto trans = problem.eclTransmissibilities().findTransmissibility(a, b);
                return trans.has_value() && *trans > 0.0;
            };
            const auto has_mech = [&problem](int cc) {
                if constexpr (requires { problem.geoMechModel().hasMech(std::size_t(cc)); }) {
                    return problem.geoMechModel().hasMech(std::size_t(cc));
                } else {
                    return true;
                }
            };
            const auto gv = grid.leafGridView();
            const Dune::cpgrid::Entity<0> element(*grid.currentData().back(), c, true);
            return stress_sides_.emplace(c, centre_stress::classifySides(
                       gv, element, logical, cartesian, *stress_faults_, flows, has_mech,
                       prm_.get<bool>("reservoir.stress_interpolation_across_faults", false)))
                .first->second;
        };
        if (interpolate_stress) {
            OPM_TIMEBLOCK(interpolateReservoirStress);
            const auto& gv = grid.leafGridView();
            Dune::BlockVector<Dune::FieldVector<double,6>> cell_stress(gv.size(0));
            for (const auto& elem : elements(gv)) {
                const int cidx = gv.indexSet().index(elem);
                cell_stress[cidx] = problem.stress(cidx);
            }
            nodal_stress = vem::patchRecoveryNodal6(grid, cell_stress);
            frac_centroids.resize(ncf);
            ElementMapper fracElemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
            for (const auto& element : elements(grid_->leafGridView()))
                frac_centroids[fracElemMapper.index(element)] = element.geometry().center();
        }

        std::vector<int> unique_reservoir_cells_;
        for(auto & cells : all_reservoir_cells_){
            for(auto & cell : cells){
                unique_reservoir_cells_.push_back(cell);
            }
        }
        std::sort(unique_reservoir_cells_.begin(),unique_reservoir_cells_.end());
        auto last = std::unique (unique_reservoir_cells_.begin(), unique_reservoir_cells_.end());
        unique_reservoir_cells_.erase(last,unique_reservoir_cells_.end());


        for (size_t i = 0; i < unique_reservoir_cells_.size(); ++i) {
            int cell = unique_reservoir_cells_[i];
            if (!(cell < 0)) {
                const auto& intQuants = simulator.model().intensiveQuantities(cell, /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                double reservoir_mobility = 0.0;
                enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        // assume sum should only be water;
                        auto val = intQuants.mobility(phaseIdx);
                        reservoir_mobility += val.value();
                    }
                }
                map_reservoir_mobility_[cell] = reservoir_mobility;
                map_reservoir_density_[cell] = fs.density(FluidSystem::waterPhaseIdx).value();
                auto pval = fs.pressure(FluidSystem::waterPhaseIdx);
                map_reservoir_pressure_[cell] = pval.value();
                map_reservoir_cell_z_[cell] = problem.dofCenterDepth(cell);
            }
        }

        reservoir_water_mobility_.resize(ncf, 0.0);
        for (size_t i = 0; i < ncf; ++i) {
            int cell = reservoir_cells_[i];
            if (!(cell < 0)) {
                auto normal = this->cell_normals_[i];
                const auto& intQuants = simulator.model().intensiveQuantities(cell, /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                                copyBlackOilFluidState(fracture_scalar_fluid_states[i], fs, value_of);
                                fracture_scalar_fluid_state_valid[i] = true;
                                reservoir_pressure_[i] = value_of(fs.pressure(FluidSystem::waterPhaseIdx));
                enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
                reservoir_mobility_[i] = 0.0;
                                reservoir_density_[i] = value_of(fs.density(FluidSystem::waterPhaseIdx));
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                  if (FluidSystem::phaseIsActive(phaseIdx)) {
                    // assume sum should only be water;
                    auto val = intQuants.mobility(phaseIdx);
                    reservoir_mobility_[i] += val.value();
                  }
                }
                // Kept apart from the total: the flow's water equation upwinds this
                // one, so leakoff_mobility=upwind needs it when the reservoir feeds
                // the fracture.
                reservoir_water_mobility_[i] =
                    value_of(intQuants.mobility(FluidSystem::waterPhaseIdx));
                for (int dim = 0; dim < 3; ++dim) {
                  reservoir_stress_[i] = problem.stress(cell);
                }
                if (interpolate_stress) {
                    // trilinear interpolation of the recovered nodal stress at the
                    // fracture cell centroid, inside the containing reservoir cell
                    const auto& cData = grid.currentData();
                    const auto& relem = Dune::cpgrid::Entity<0>(*cData.back(), cell, true);
                    auto loc = relem.geometry().local(frac_centroids[i]);
                    bool loc_ok = true;
                    for (int d = 0; d < 3; ++d) {
                        // guard against failed Newton inversion on degenerate cells;
                        // centroids slightly outside the cell are clamped
                        if (!std::isfinite(loc[d]) || loc[d] < -0.5 || loc[d] > 1.5) {
                            loc_ok = false;
                            break;
                        }
                        loc[d] = std::min(1.0, std::max(0.0, loc[d]));
                    }
                    if (loc_ok && relem.geometry().corners() == 8) {
                        const auto& gv = grid.leafGridView();
                        Dune::FieldVector<double,6> interp(0.0);
                        for (int k = 0; k < 8; ++k) {
                            double w = 1.0;
                            for (int d = 0; d < 3; ++d)
                                w *= ((k >> d) & 1) ? loc[d] : 1.0 - loc[d];
                            interp.axpy(w, nodal_stress[gv.indexSet().subIndex(relem, k, 3)]);
                        }
                        reservoir_stress_[i] = interp;
                    } // else: keep the containing-cell value
                }
                if (stress_mode != centre_stress::Mode::None) {
                    const auto samples = centre_stress::samplesFromSides(
                        cell_sides(cell), centre_of_cell,
                        [&problem](int c) { return problem.stress(c); });
                    reservoir_stress_[i] = centre_stress::interpolate(
                        frac_centroids[i], normal, centre_of_cell(cell), problem.stress(cell), samples, stress_mode,
                        [&normal](const auto& s) { return ddm::tractionVoigt(s, normal); });
                }
                const auto& perm =  problem.intrinsicPermeability(cell);
                const auto& cstress =  problem.cStress(cell);
                auto pn = normal;
                perm.mv(normal, pn);
                double npn = pn.dot(normal);
                reservoir_perm_[i] = npn;
                reservoir_cstress_[i] = cstress;
                const auto& currentData = grid.currentData();
                const auto& elem = Dune::cpgrid::Entity<0>(*currentData.back(), cell, true);
                const auto& geom = elem.geometry();
                auto cell_center = geom.center();
                if(calculate_dist){
                  double dist = 0;
                  int num_corners = geom.corners();
                  for(int li=0; li < geom.corners();++li){
                    auto cdist = cell_center;
                    cdist -= geom.corner(li);
                    dist += std::abs(normal.dot(cdist));
                  }
                  dist /= num_corners; // = L/2 for a box cell: mean |n.(centre-corner)|
                  reservoir_dist_[i] = 2.0*dist*dist_factor;
                }
                // depth at which Flow defines the cell pressure
                reservoir_cell_z_[i] = problem.dofCenterDepth(cell);
            } else {
                // probably outside reservoir set all to zero
                double stressval = 300*300e5;// large value to keep it closed
                //NB should maybe be interpolated out from reservoir
                reservoir_stress_[i][0] = stressval;
                reservoir_stress_[i][1] = stressval;
                reservoir_stress_[i][2] = stressval; //???
                reservoir_pressure_[i] = injectionPressure();
                reservoir_mobility_[i] = 1000.0;
                reservoir_water_mobility_[i] = 1000.0;
                reservoir_density_[i] = 1000.0;
                reservoir_perm_[i] = 0.0;
                reservoir_cstress_[i]  = 1e99;
                reservoir_dist_[i] = 1e99;
                reservoir_cell_z_[i] = origo_[2];// should maybe be cell center of triangle
            }
            // assume reservoir distance is calculated
        }

        if (stress_mode != centre_stress::Mode::None && stress_sides_.size() != stress_sides_logged) {
            std::array<int, centre_stress::numSideKinds> counts {};
            for (const auto& entry : stress_sides_) {
                for (const auto& side : entry.second) {
                    ++counts[static_cast<std::size_t>(side.kind)];
                }
            }
            OpmLog::debug(fmt::format("Stress interpolation sides of {} reservoir cells: clean {}, boundary {}, "
                                      "non-conforming {}, logical offset {}, declared fault {}, no flow {}, no mechanics {}",
                                      stress_sides_.size(), counts[0], counts[1], counts[2], counts[3],
                                      counts[4], counts[5], counts[6]));
        }

        fracture_water_property_evaluator_ =
            [fluid_states = std::move(fracture_scalar_fluid_states),
             valid_states = std::move(fracture_scalar_fluid_state_valid),
             fallback_density = reservoir_density_](const size_t cell_idx, const double pressure)
                -> std::pair<CellFluidProperty, CellFluidProperty>
        {
            using PressureEval = Opm::DenseAd::Evaluation<double, 1>;
            using AdFluidState = Opm::BlackOilFluidState<PressureEval, FluidSystem>;
            const auto make_constant = [](const auto& value) {
                return PressureEval::createConstant(Opm::scalarValue(value));
            };

            const auto fallback = [&, cell_idx]() {
                const double density = cell_idx < fallback_density.size() ? fallback_density[cell_idx] : 1000.0;
                return std::make_pair(CellFluidProperty{density, 0.0}, CellFluidProperty{1.0, 0.0});
            };

            if (cell_idx >= fluid_states.size() || !valid_states[cell_idx])
                return fallback();

            AdFluidState ad_state;
            copyBlackOilFluidState(ad_state, fluid_states[cell_idx], make_constant);
            ad_state.setPressure(FluidSystem::waterPhaseIdx,
                                 PressureEval::createVariable(pressure, 0));

            const auto make_property = [](const auto& property_eval) {
                return CellFluidProperty{property_eval.value(), property_eval.derivative(0)};
            };

            const auto density = make_property(
                FluidSystem::density(ad_state, FluidSystem::waterPhaseIdx, ad_state.pvtRegionIndex()));
            const auto viscosity = make_property(
                FluidSystem::viscosity(ad_state, FluidSystem::waterPhaseIdx, ad_state.pvtRegionIndex()));
            return std::make_pair(density, viscosity);
        };

        fracture_dgh_.resize(ncf);
        for (auto element : Dune::elements(grid_->leafGridView())) {
            int i = grid_->leafGridView().indexSet().index(element);
            double z = element.geometry().center()[2];
            fracture_dgh_[i] = gravity_ * reservoir_density_[i] * z;
        }
        int reportStepIdx = simulator.episodeIndex();

        const auto& schedule =  simulator.vanguard().schedule();
        // get properties from well connections in case of filter cake
        if(schedule[reportStepIdx].wells.has(wellinfo_.name) == false){
            OpmLog::warning(fmt::format("Well {} not found in schedule step {}", wellinfo_.name, reportStepIdx));
            return;
        }
        const auto& well = schedule[reportStepIdx].wells(wellinfo_.name);
        const auto& connections = well.getConnections();
        total_WI_well_ = 0.0;
        for(const auto& conn : connections){
          if(conn.state() == Connection::State::OPEN){
              continue;
          }
          assert(conn.CF() >= 0);
          total_WI_well_ += conn.CF();
        }
        const auto* connection = wellinfo_.connection_id.find(connections);
        if (connection == nullptr) {
            OpmLog::warning(fmt::format("Well connection with global index {} in grid {} for well {} not found in schedule step {}",
                                        wellinfo_.connection_id.global_index, wellinfo_.connection_id.lgr_level,
                                        wellinfo_.name, reportStepIdx));
            return;
        }
        const auto& wellstates = simulator.problem().wellModel().wellState();
        const auto& well_index = wellstates.index(wellinfo_.name);
        if (!well_index.has_value()) {
            OpmLog::warning(fmt::format("Well {} not found in well state at step {}", wellinfo_.name, reportStepIdx));
            has_filtercake_ = false;// prevois state did not have this well
            return;
        }
        const auto& wellstate = wellstates[*well_index];

    }



template <class TypeTag, class Simulator>
void Fracture::updateFilterCakePropertiesPost(const Simulator& simulator)
    {

        int reportStepIdx = simulator.episodeIndex();
        const auto& schedule =  simulator.vanguard().schedule();
        // get properties from well connections in case of filter cake
        if(schedule[reportStepIdx].wells.has(wellinfo_.name) == false){
            OpmLog::warning(fmt::format("Well {} not found in schedule step {}", wellinfo_.name, reportStepIdx));
            return;
        }
        const auto& well = schedule[reportStepIdx].wells(wellinfo_.name);
        const auto& connections = well.getConnections();

        const auto* connection = wellinfo_.connection_id.find(connections);
        if (connection == nullptr) {
            OpmLog::warning(fmt::format("Well connection with global index {} in grid {} for well {} not found in schedule step {}",
                                        wellinfo_.connection_id.global_index, wellinfo_.connection_id.lgr_level,
                                        wellinfo_.name, reportStepIdx));
            return;
        }
        const auto& wellstates = simulator.problem().wellModel().wellState();
        const auto& well_index = wellstates.index(wellinfo_.name);
        if (!well_index.has_value()) {
            OpmLog::warning(fmt::format("Well {} not found in well state at step {}", wellinfo_.name, reportStepIdx));
            has_filtercake_ = false;// prevois state did not have this well
            return;
        }
        const auto& wellstate = wellstates[*well_index];

        double dt = simulator.timeStepSize();
        updateFilterCakeProps(*connection, wellstate, dt);

    }

// ----------------------------------------------------------------------------
template <class TypeTag, class Simulator>
void Fracture::solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
                     const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
                     const Simulator& simulator)
// ----------------------------------------------------------------------------
{
  const auto& cpgrid = simulator.vanguard().grid();
    if(!active_){
                last_solve_stats_ = {};
        OpmLog::debug(fmt::format("Fracture {} is not active, skipping solve.", this->name()));
        return;
    }
        Dune::Timer solve_timer;
        last_solve_stats_ = {};
        // control.type=="well": latch the BC for this whole solve (the well DOF
        // count must not change between rounds - the remap asserts otherwise),
        // and resize the persistent pressure vector across a control transition.
        if (prm_.get_child("control").get<std::string>("type") == "well") {
            // generalized well row: both constraints keep the well DOF, so a
            // control switch changes one row's content only - no state resize
            const std::string ctrl_now = well_control_is_rate_ ? "rate_well" : "bhp_well";
            if (ctrl_now != effective_control_latched_) {
                OpmLog::info(fmt::format("Fracture {}: well constraint {}{}", name(),
                                         effective_control_latched_.empty()
                                             ? std::string("initialized to ")
                                             : effective_control_latched_ + " -> ",
                                         ctrl_now));
                effective_control_latched_ = ctrl_now;
            }
        }
        last_solve_stats_.fractures_solved = 1;
    OPM_TIMEBLOCK(SolveFracture);
    {
        OpmLog::debug("Solve Fracture Pressure using Iterative Fracture with Trimesh Propagation");
        const auto opt = PropagationOptions::read(prm_);
        VolumePacing pacing = VolumePacing::start(*this, simulator.timeStepSize());
        if (!opt.remapSolution) {
            reinitializeTrimeshState_();
        }
        // the grid and filter cake at entry, to map onto each evolved grid
        const FilterCakeSnapshot cake0 = snapshotFilterCake_(opt);

        // Regrid to a candidate mesh, solve on it and score its front cells.
        auto score_function = [&](const RegularTrimesh& trimesh, const int level) -> std::vector<double> {
            auto maps = regridTo_<TypeTag, Simulator>(trimesh, level, opt,
                                                      cell_search_tree, entity_seeds, simulator);
            remapState_(maps, level, opt, cake0);
            const auto solve = solveOnCurrentGrid_(opt);
            pacing.afterSolve(*this);
            reportGridSolve_(solve, opt);
            return propagationScores_(maps.boundary, solve, opt);
        };

        const FractureProperties& fprop_before = calculateFractureProperties();
        auto stop_criterion= [&]() -> bool
        {
          // a joint round scores the current state once; the well's solve comes next
          if (joint_growth_round_ && joint_scores_ > 0) {
              return true;
          }
          const bool area_stop = expantionMax(fprop_before);
          if (pacing.on && pacing.budgetSpent) {
              if (pacing.verbosity > 0) OpmLog::debug("VolumePacing: stopping expansion (budget spent)");
              return true;
          }
          return area_stop;
        };

        const double threshold = 1.0;
        const std::vector<CellRef> fixed_cells = well_source_cellref_;
        const int target_cellcount = prm_.get<int>("solver.target_cellcount");
        // expand_to_criterion refines back only after a re-solve per expansion, which a joint
        // round does not do: coarsen the seed once and hold that level.
        int cellcount_threshold = prm_.get<int>("solver.cellcount_threshold");
        if (joint_growth_round_) {
            cellcount_threshold = joint_level_held_ ? std::numeric_limits<int>::max()
                                                    : prm_.get<int>("solver.well_joint_seed_cellcount", 100);
            joint_level_held_ = true;
        }
        const int max_iter_on_same_level = prm_.get<int>("solver.max_iter_on_same_level");
        const auto initial_mesh= *trimesh_;
        auto [mesh, cur_level] =
            expand_to_criterion(*trimesh_, score_function, threshold,
                                fixed_cells,
                                target_cellcount,
                                cellcount_threshold,
                                max_iter_on_same_level, stop_criterion);




        // this will set trimesh_,change dune grid  and solve
        if (joint_growth_round_) {
          // new cells have no width yet; the well's solve decides whether they open
          joint_grew_ = mesh.numCells() != initial_mesh.numCells();
          if (joint_grew_) {
            score_function(mesh, cur_level);
          }
        } else {
        bool any_removed = removeNewZeroWithCells(mesh,cur_level, initial_mesh);
        if(any_removed){
          mesh.removeSawtooths();
          score_function(mesh, cur_level);
        }
        }


            // make current level become the reference (finest) level
            // note that the well_source_cellref_ is already set from the last call to the score function
            for (auto& cell : well_source_cellref_)
                cell = RegularTrimesh::fine_to_coarse(cell, cur_level);

            // volume pacing: the drawn-down BC was for the expansion rounds only;
            // the coupling continues from the flow's perf pressure
            pacing.restorePerfPressure(*this);

            // Opt-in damped fixed point on the width exposed to the coupling
            // (WI, volume, K1 all derive from it): prevents the contact
            // open/close bang-bang from snapping the coupling between rounds.
            // Skipped when the mesh changed size this round.
            const double omega = prm_.get<double>("solver.coupling_state_relaxation", 1.0);
            if (omega < 1.0) {
                if (width_round_valid_
                    && width_round_prev_.size() == fracture_width_.size()) {
                    for (size_t i = 0; i < fracture_width_.size(); ++i) {
                        fracture_width_[i] = omega * fracture_width_[i]
                            + (1.0 - omega) * width_round_prev_[i];
                    }
                }
                width_round_prev_ = fracture_width_;
                width_round_valid_ = true;
            }

            // ----------------------------------------------------------------------------
        }
        bool first_solve = well_indices_.size() == 0;
        if(first_solve){
            well_indices_.resize(2);
        }
        auto well_indices_new = wellIndices_();
        auto well_indices_old = wellIndices(); // as used by the simulator at the previous step
        if(first_solve){
            well_indices_[0] = well_indices_new;
            well_indices_[1] = well_indices_new;
        }else{
            well_indices_[1] = well_indices_old;
            well_indices_[0] = well_indices_new;
        }
        // Opt-in vector acceleration of the CTF coupling fixed point; the mixed
        // vector is what wellIndices() then hands to the well model.
        if (prm_.get<std::string>("solver.wi_vector_acceleration", "none") != "none") {
            well_indices_accel_ = applyVectorCouplingUpdate(well_indices_new);
        }
        last_solve_stats_.solve_time_seconds = solve_timer.stop();
        // summary of solve
        summary_of_solve();


}

template <class TypeTag, class Simulator>
Fracture::RegridMaps
Fracture::regridTo_(const RegularTrimesh& trimesh, const int level, const PropagationOptions& opt,
                    const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
                    const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
                    const Simulator& simulator)
{
    const auto& cpgrid = simulator.vanguard().grid();
    *trimesh_ = trimesh;
    // well sources at this level; their cells are kept intact by createDuneGrid
    std::vector<CellRef> wsources = well_source_cellref_;
    for (auto& cell : wsources)
        cell = RegularTrimesh::fine_to_coarse(cell, level);

    auto [grid, fsmap, bmap] = trimesh_->createDuneGrid(opt.maxNumCoarsening, wsources,
                                                        opt.smoothBoundary, opt.numcellThreshold);
    RegridMaps maps;
    maps.before = grid_mesh_map_;
    grid_mesh_map_ = fsmap;
    grid_fixed_cells_ = wsources;
    composeFlowStateDonor(maps.before, fsmap, level);
    setFractureGrid(std::move(grid));

    std::map<CellRef, size_t> fsmap_inv;
    for (int i = 0; i != fsmap.size(); ++i)
        if (size(fsmap[i]) == 1) // a fine-scale cell
            fsmap_inv[fsmap[i].front()] = i;
    well_source_.clear();
    for (const auto& cell : wsources)
        well_source_.push_back(fsmap_inv[cell]);

    updateReservoirCells(cell_search_tree, cpgrid, entity_seeds);
    updateReservoirProperties<TypeTag, Simulator>(simulator, true, false);
    initPressureMatrix();
    rhs_pressure_.resize(0);
    coupling_matrix_ = nullptr;

    maps.after = std::move(fsmap);
    maps.boundary = std::move(bmap);
    return maps;
}

template <class Serializer>
void
Fracture::serializeState(Serializer& serializer)
{
    namespace ckpt = GeomechCheckpoint;
    if (!trimesh_) {
        OPM_THROW(std::runtime_error, "Fracture " + name() + ": only trimesh fractures can be saved.");
    }
    const bool loading = !serializer.isSerializing();

    serializer(*trimesh_);
    serializer(grid_fixed_cells_);
    std::vector<std::vector<CellRef>> saved_map = grid_mesh_map_;
    serializer(saved_map);
    if (loading) {
        rebuildGridForLoad(saved_map);
    }
    bool has_prev = trimesh_prev_ != nullptr;
    serializer(has_prev);
    if (has_prev) {
        if (loading) {
            trimesh_prev_ = std::make_unique<RegularTrimesh>(*trimesh_);
        }
        serializer(*trimesh_prev_);
    }
    serializer(grid_mesh_map_prev_);
    serializer(filtercake_thikness_prev_);
    serializer(well_source_cellref_prev_);
    serializer(well_source_prev_);

    serializer(active_);
    serializer(filtercake_thikness_);
    serializer(filtercake_perm_);
    serializer(filtercake_poro_);
    serializer(has_filtercake_);
    serializer(well_source_cellref_);
    serializer(well_source_);

    serializer(reservoir_cells_);
    serializer(all_reservoir_cells_);
    serializer(all_reservoir_areas_);
    std::size_t num_center_lists = all_reservoir_centers_.size();
    serializer(num_center_lists);
    all_reservoir_centers_.resize(num_center_lists);
    for (auto& centers : all_reservoir_centers_) {
        ckpt::fieldVectors(serializer, centers);
    }
    serializer(reservoir_perm_);
    serializer(reservoir_cstress_);
    serializer(reservoir_mobility_);
    serializer(reservoir_water_mobility_);
    serializer(reservoir_density_);
    serializer(reservoir_cell_z_);
    serializer(map_reservoir_mobility_);
    serializer(map_reservoir_density_);
    serializer(map_reservoir_cell_z_);
    serializer(map_reservoir_pressure_);
    serializer(reservoir_dist_);
    serializer(reservoir_pressure_);
    ckpt::fieldVectors(serializer, reservoir_stress_);

    ckpt::fieldVectors(serializer, fracture_width_);
    ckpt::fieldVectors(serializer, fracture_pressure_);
    ckpt::fieldVectors(serializer, fracture_width_prev_);
    ckpt::fieldVectors(serializer, fracture_pressure_prev_);
    ckpt::fieldVectors(serializer, width_round_prev_);
    serializer(width_round_valid_);
    serializer(area_step_start_);
    serializer(perf_pressure_step_start_);
    serializer(volume_step_start_);
    serializer(coupling_dt_cap_);
    serializer(traction_step_start_);
    serializer(stress_dt_cap_);
    serializer(well_control_is_rate_);
    serializer(well_target_bhp_);
    serializer(effective_control_latched_);
    serializer(coupling_quiet_steps_);

    serializer(htrans_);
    serializer(perfinj_);
    serializer(well_perf_cells_);
    serializer(fb_cell_residual_);
    serializer(current_dt_);
    serializer(external_pressure_);
    serializer(external_cell_mask_);
    serializer(flow_state_donor_);
    serializer(flow_state_remapped_);
    serializer(cell_last_toggle_iter_);
    serializer(current_solve_iter_);
    serializer(perf_pressure_);
    serializer(leakof_);
    serializer(leakoff_mobility_);

    serializer(E_);
    serializer(nu_);
    serializer(min_width_);
    serializer(fracture_dgh_);
    serializer(total_WI_well_);

    const auto perforations = [&serializer](std::vector<RuntimePerforation>& perfs) {
        std::size_t num = perfs.size();
        serializer(num);
        perfs.resize(num);
        for (auto& perf : perfs) {
            serializer(perf.cell);
            serializer(perf.ctf);
            serializer(perf.depth);
            serializer(perf.pressure);
            serializer(perf.segment);
            serializer(perf.perf_range);
            serializer(perf.fracture_pressure_offset);
        }
    };
    std::size_t num_index_lists = well_indices_.size();
    serializer(num_index_lists);
    well_indices_.resize(num_index_lists);
    for (auto& perfs : well_indices_) {
        perforations(perfs);
    }
    perforations(well_indices_accel_);
    serializer(max_flow_time_step_);
    serializer(well_rate_);
    serializer(total_wellindex_);
    serializer(wi_dz_);
    serializer(wi_respress_);
    serializer(well_ref_depth_);
    serializer(perf_ref_depth_);
    serializer(well_perf_rate_);

    const auto mix = [&serializer](CouplingMixState& state) {
        serializer(state.initialized);
        serializer(state.has_previous_residual);
        serializer(state.previous_target);
        serializer(state.previous_residual);
        serializer(state.omega);
    };
    mix(perf_pressure_mix_);
    mix(well_rate_mix_);
    mix(total_wellindex_mix_);
    serializer(wi_vec_mix_.initialized);
    serializer(wi_vec_mix_.has_prev_residual);
    serializer(wi_vec_mix_.omega);
    serializer(wi_vec_mix_.prev_mixed);
    serializer(wi_vec_mix_.prev_target);
    serializer(wi_vec_mix_.prev_residual);

    serializer(density_);
    serializer(density_perf_);
    serializer(mobility_water_perf_);
    serializer(closed_cells_);
    serializer(cell_flip_counts_);
    serializer(last_solve_stats_.converged);
}

} // namespace Opm
