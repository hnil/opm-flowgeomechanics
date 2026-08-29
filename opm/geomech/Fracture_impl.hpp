#pragma once

#include "RegularTrimesh.hpp"
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/vem/vemutils.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> 
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

inline double
compute_target_expansion(const double K1_target,
                         const double aperture,
                         const double E, // young
                         const double nu) // poisson
{
    const double mu = E / (2 * (1 + nu)); // shear modulus
    const double fac = mu * std::sqrt(M_PI) / (2 * (1 - nu) * 1.834);
    return pow(fac * aperture / K1_target, 2);
}


// ----------------------------------------------------------------------------


  
 // template<class CoordType> inline void ensure_convexity(const CoordType& ipoint,
//                                                        std::vector<CoordType>& pts) {
//   while(true) {
//     bool modified = false;
//     for (int ix = 0; ix != pts.size(); ++ix) {
//         const int ixp = (ix + 1) % pts.size();
//         const int ixm = (ix + pts.size() - 1) % pts.size();

//         CoordType& p = pts[ix];
//         const CoordType& pp = pts[ixp];
//         const CoordType& pm = pts[ixm];

//         const CoordType v1 = pp - p;
//         const CoordType v2 = pm - p;
//         const CoordType vi = p - ipoint;

//         const CoordType cross1 { v1[1] * v2[2] - v1[2] * v2[1],
//                                  v1[2] * v2[0] - v1[0] * v2[2],
//                                  v1[0] * v2[1] - v1[1] * v2[0] };

//         const CoordType cross2 { v1[1] * vi[2] - v1[2] * vi[1],
//                                  v1[2] * vi[0] - v1[0] * vi[2],
//                                  v1[0] * vi[1] - v1[1] * vi[0] };

//         const double sprod = cross1[0] * cross2[0] +
//                              cross1[1] * cross2[1] +
//                              cross1[2] * cross2[2];

//         if (sprod < 0) {
//           // concave corner
//           p = 0.5 * (pp + pm);
//           modified = true;
//         }
//     }
//     if (!modified)
//       break;
//   }

// };


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
        // NB burde truleg interpolere
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


        //double numax = 0, Emax = 0;
        //auto& enitity_seeds = problem.elementEntitySeeds();//used for radom axes better way?
        // get properties from rservoir


        // optionally evaluate the reservoir stress at the fracture cell centroids by
        // trilinear interpolation of patch-recovered nodal stress, instead of the
        // piecewise-constant stress of the containing cell.  With large reservoir
        // cells the latter gives large jumps in the traction along the fracture.
        // ("NB burde truleg interpolere" -- now optional via reservoir.interpolate_stress)
        const bool interpolate_stress = prm_.get<bool>("reservoir.interpolate_stress", false);
        std::vector<Dune::FieldVector<double,6>> nodal_stress;
        std::vector<Dune::FieldVector<double,3>> frac_centroids;
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

        std::vector<int> unique_reservoir_cells_;//NB move to member?
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
                const auto& currentData = grid.currentData();
                const auto& elem = Dune::cpgrid::Entity<0>(*currentData.back(), cell, true);
                //const auto& geom = elem.geometry();
                const auto& cell_center = elem.geometry().center();
                map_reservoir_mobility_[cell] = reservoir_mobility;
                map_reservoir_density_[cell] = fs.density(FluidSystem::waterPhaseIdx).value();
                auto pval = fs.pressure(FluidSystem::waterPhaseIdx);
                map_reservoir_pressure_[cell] = pval.value();
                map_reservoir_cell_z_[cell] = cell_center[2];
            }
        }

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
                const auto& perm =  problem.intrinsicPermeability(cell);
                const auto& cstress =  problem.cStress(cell);
                auto pn = normal;
                perm.mv(normal, pn);
                double npn = pn.dot(normal); 
                reservoir_perm_[i] = npn;
                reservoir_cstress_[i] = cstress;
                const auto& currentData = grid.currentData();
                const auto& elem = Dune::cpgrid::Entity<0>(*currentData.back(), cell, true);
                //auto cell_center = ghelper.centroid(i);
                const auto& geom = elem.geometry();
                auto cell_center = geom.center();
                if(calculate_dist){
                  //assert(false);                  
                  double dist = 0;
                  int num_corners = geom.corners();
                  for(int li=0; li < geom.corners();++li){
                    auto cdist = cell_center;
                    cdist -= geom.corner(li);
                    dist += std::abs(normal.dot(cdist));
                  }
                  dist /= num_corners;
                  reservoir_dist_[i] = dist/2.0;
                }
                reservoir_cell_z_[i] = cell_center[2];
                   
                
            } else {
                // probably outside reservoir set all to zero
                double stressval = 300*300e5;// large value to keep it closed
                //NB should maybe be interpolated out from reservoir
                reservoir_stress_[i][0] = stressval;
                reservoir_stress_[i][1] = stressval;
                reservoir_stress_[i][2] = stressval; //???
                //
                reservoir_pressure_[i] = injectionPressure();
                reservoir_mobility_[i] = 1000.0;
                reservoir_density_[i] = 1000.0;
                reservoir_perm_[i] = 0.0;
                reservoir_cstress_[i]  = 1e99;
                reservoir_dist_[i] = 1e99;
                reservoir_cell_z_[i] = origo_[2];// should maybe be cell center of triangle
            }
            // assume reservoir distance is calculated
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
           // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
            // Dune::mcmgElementLayout>::index(element);
            int i = grid_->leafGridView().indexSet().index(element);
            double z = element.geometry().center()[2];
            fracture_dgh_[i] = gravity_ * reservoir_density_[i] * z;
        }
        int reportStepIdx = simulator.episodeIndex();
        
        const auto& schedule =  simulator.vanguard().schedule();
        //const Opm::Well& well = wells[wellinfo_.name];
        //const std::vector<Opm::Well>& wells = problem.wellModel().getLocalWells(reportStepIdx);
        // get properties from well connections in case of filter cake
        if(schedule[reportStepIdx].wells.has(wellinfo_.name) == false){
            std::cerr << "Warning: Well " << wellinfo_.name << " not found in schedule step " << reportStepIdx << std::endl;
            return;
        }
        const auto& well = schedule[reportStepIdx].wells(wellinfo_.name);
        //const auto& well = problem.wellModel().getWellEcl(wellinfo_.name);
        const auto& connections = well.getConnections();
        total_WI_well_ = 0.0;
        for(const auto& conn : connections){
          if(conn.state() == Connection::State::OPEN){
              continue;
          }
          assert(conn.CF() >= 0);
          total_WI_well_ += conn.CF();
        }
        //const auto& connection = connections[wellinfo_.perf];// probably wrong
        if(connections.hasGlobalIndex(wellinfo_.global_index) == false){
            std::cerr << "Warning: Well connection with global index " << wellinfo_.global_index << " not found in schedule step " << reportStepIdx << std::endl;
            return;
        }
        const auto& wellstates = simulator.problem().wellModel().wellState();
        const auto& well_index = wellstates.index(wellinfo_.name);
        if (!well_index.has_value()) {
            std::cerr << "Warning: Well " << wellinfo_.name << " not found in well state at step " << reportStepIdx << std::endl;
            has_filtercake_ = false;// prevois state did not have this well
            return;
        }
        const auto& wellstate = wellstates[*well_index];

        // if(update_filtercake){
        //   double dt = simulator.timeStepSize();
        //   updateFilterCakeProps(connections, wellstate, dt);
        // } 
        // else {
        //   // remap
        //   int ncf = reservoir_cells_.size();
        //   filtercake_thikness_.resize(ncf, 0.0);
        // }
    }



template <class TypeTag, class Simulator>
void Fracture::updateFilterCakePropertiesPost(const Simulator& simulator)//, bool init_constant_vals,bool update_filtercake)
    {
        // if `init_contant_vals` is true, the fields that should normally not
        // change, i.e.  E_, nu_ and reservoir_perm_ will be updated. This should
        // normally only be needed the first time.

        int reportStepIdx = simulator.episodeIndex();
        const auto& schedule =  simulator.vanguard().schedule();
        //const Opm::Well& well = wells[wellinfo_.name];
        //const std::vector<Opm::Well>& wells = problem.wellModel().getLocalWells(reportStepIdx);
        // get properties from well connections in case of filter cake
        if(schedule[reportStepIdx].wells.has(wellinfo_.name) == false){
            std::cerr << "Warning: Well " << wellinfo_.name << " not found in schedule step " << reportStepIdx << std::endl;
            return;
        }
        const auto& well = schedule[reportStepIdx].wells(wellinfo_.name);
        //const auto& well = problem.wellModel().getWellEcl(wellinfo_.name);
        const auto& connections = well.getConnections();
        
        //const auto& connection = connections[wellinfo_.perf];// probably wrong
        if(connections.hasGlobalIndex(wellinfo_.global_index) == false){
            std::cerr << "Warning: Well connection with global index " << wellinfo_.global_index << " not found in schedule step " << reportStepIdx << std::endl;
            return;
        }
        const auto& wellstates = simulator.problem().wellModel().wellState();
        const auto& well_index = wellstates.index(wellinfo_.name);
        if (!well_index.has_value()) {
            std::cerr << "Warning: Well " << wellinfo_.name << " not found in well state at step " << reportStepIdx << std::endl;
            has_filtercake_ = false;// prevois state did not have this well
            return;
        }
        const auto& wellstate = wellstates[*well_index];

        double dt = simulator.timeStepSize();
        updateFilterCakeProps(connections, wellstate, dt);

    }   

// ----------------------------------------------------------------------------
template <class TypeTag, class Simulator>
void Fracture::solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
                     //const Dune::CpGrid& cpgrid,
                     const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
                     const Simulator& simulator)
// ----------------------------------------------------------------------------
{
  const auto& cpgrid = simulator.vanguard().grid();
    if(!active_){
                last_solve_stats_ = {};
        std::cout << "Fracture " << this->name() << " is not active, skipping solve." << std::endl;
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
                std::cout << "Fracture " << name() << ": well constraint "
                          << (effective_control_latched_.empty()
                                  ? std::string("initialized to ")
                                  : effective_control_latched_ + " -> ")
                          << ctrl_now << std::endl;
                effective_control_latched_ = ctrl_now;
            }
        }
        last_solve_stats_.fractures_solved = 1;
    OPM_TIMEBLOCK(SolveFracture);
    //std::cout << "Solve Fracture Pressure" << std::endl;
    std::string method = prm_.get<std::string>("solver.method");
    if (method == "nothing") {
    } else if (method == "simple") {
        this->solveFractureWidth();
        this->solvePressure();
    } else if (method == "only_pressure") {
        this->solvePressure();
    } else if (method == "only_width") {
        this->solveFractureWidth();
    } else if (method == "iterative") {
        int max_it = prm_.get<int>("max_iter");
        int it = 0;
        bool changed = true;
        while (changed && (it < max_it)) {
            initFractureStates(); // ensure initial fracture_width and fracture_pressure
                                  // set to something reasonable
            auto fracture_width = fracture_width_;
            auto fracture_pressure = fracture_pressure_;
            this->solveFractureWidth();
            // grow fracture
            this->solvePressure();
            it += 1;
            double tol = prm_.get<double>("solver.max_change");
            double max_change = 0;
            for (int i = 0; fracture_width_.size(); ++i) {
                double diff_width = fracture_width_[i] - fracture_width[i];
                double diff_press = fracture_pressure_[i] - fracture_pressure[i];
                max_change = std::max(max_change, diff_width / 1e-2);
                max_change = std::max(max_change, diff_press / 1e5);
            }
            changed = (max_change < tol);
        }

        // ----------------------------------------------------------------------------
    } else if (method == "if") {
        // ----------------------------------------------------------------------------
        // iterate full nonlinear system until convergence
        std::cout << "Solve Fracture Pressure using Iterative Fracture" << std::endl;
        double min_width = prm_.get<double>("solver.min_width");
        for (auto& width : fracture_width_) {
            width[0] = std::max(width[0], min_width); // Ensure not completely closed
        }
        // start by assuming pressure equal to confining stress (will also set
        // fracture_pressure_ to its correct size
        normalFractureTraction(fracture_pressure_);
        if (numWellEquations() > 0) // @@ it is implicitly assumed for now that
            // there is just one well equation.  We initializze
            // it with an existing value.
            fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];

        const double tol = prm_.get<int>("solver.tolerance"); // 1e-5; // @@
        const int max_iter = prm_.get<int>("solver.max_iter");
        int iter = 0;

        // solve flow-mechanical system
        const int nlin_verbosity = prm_.get<double>("solver.verbosity");
        while (!fullSystemIteration(tol,iter) && iter++ < max_iter) {
            if (nlin_verbosity > 1) {
                std::cout << "Iteration: " << iter << std::endl;
            }
        };
        if(iter >= max_iter){
            last_solve_stats_.converged = false;
            bool failure_on_convergence = prm_.get<bool>("solver.failure_on_nonconvergence",false);
            if(failure_on_convergence){
                throw std::runtime_error("Fracture solver did not converge within " + std::to_string(max_iter) + " iterations.");
            }else{
                std::cout << "Warning: Fracture solver did not converge within " << max_iter << " iterations." << std::endl;
            }
        }

        // @@ debug
        const std::vector<double> K1_not_nan = Fracture::stressIntensityK1();
        std::vector<double> K1;
        for (size_t i = 0; i != K1_not_nan.size(); ++i)
            if (!std::isnan(K1_not_nan[i]))
                K1.push_back(K1_not_nan[i]);

        std::cout << "K1: ";
        std::cout << *std::min_element(K1.begin(), K1.end()) << ", " << *std::max_element(K1.begin(), K1.end())
                  << std::endl;
        std::cout << "Pressure: ";
        std::cout << *std::min_element(fracture_pressure_.begin(), fracture_pressure_.end()) << ", "
                  << *std::max_element(fracture_pressure_.begin(), fracture_pressure_.end()) << std::endl;
        std::cout << "Normal traction: ";
        Dune::BlockVector<Dune::FieldVector<double, 1>> krull(fracture_width_);
        normalFractureTraction(krull, false);
        std::cout << *std::min_element(krull.begin(), krull.end()) << ", "
                  << *std::max_element(krull.begin(), krull.end()) << std::endl;
        std::cout << "Aperture: ";
        std::cout << *std::min_element(fracture_width_.begin(), fracture_width_.end()) << ", "
                  << *std::max_element(fracture_width_.begin(), fracture_width_.end()) << std::endl;

        // ----------------------------------------------------------------------------
    } else if (method == "if_propagate_trimesh") {
        // ----------------------------------------------------------------------------
        std::cout << "Solve Fracture Pressure using Iterative Fracture with Trimesh Propagation" << std::endl;
        // Volume-paced propagation (opt-in): see VOLUME_PACED_PROPAGATION.md.
        // V0 = fracture volume at the step checkpoint; V_avail = fluid the well
        // can deliver into new fracture volume during this step (eq. 1).
        const bool vol_pacing = prm_.get<bool>("solver.volume_paced_propagation", false);
        const double vol_fac = prm_.get<double>("solver.volume_pacing_factor", 1.0);
        const bool vol_pressure = prm_.get<bool>("solver.volume_pacing_pressure", false);
        const int vol_verb = prm_.get<int>("solver.volume_pacing_verbosity", 0);
        const double vol_dt = simulator.timeStepSize();
        const double vol_V0 = vol_pacing ? std::max(volume_step_start_, 0.0) : 0.0;
        const double vol_perf_p0 = perf_pressure_;
        bool vol_budget_spent = false;
        bool reinitialize= !(prm_.get<bool>("solver.remap_solution"));
        if(reinitialize){
            std::cout << "Reinitializing fracture states before IF with trimesh propagation" << std::endl;
            //if(true){
                fracture_width_ = 1e-3; // Ensure not completely closed
            //fracture_pressure_ = perf_pressure_;
            //}
            // start by assuming pressure equal to confining stress (will also set
            // fracture_pressure_ to its correct size
            auto traction = fracture_pressure_;
            normalFractureTraction(traction);
            for(size_t i=0; i < traction.size();++i){
                fracture_pressure_[i] = std::max(traction[i], fracture_pressure_[i]);
            }
            if (numWellEquations() > 0){
                fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];
            }
        }
        // It is implicitly assumed for now that there is just one well equation.
        // We initialize with an existing value. @@
        

        // save original grid and filtercake, to allow us to map it onto evolved grids
        const auto filtercake_thickness_0 = filtercake_thikness_; // copy
        const auto grid_mesh_map_0 = grid_mesh_map_; 
        double total_filtercake_0 = filterCakeVolume();
        // local function taking a trimesh, updates the Fracture object with it and
        // runs a simulation.  Its return value should be a vector of doubles:
        
        auto score_function = [&](const RegularTrimesh& trimesh, const int level) -> std::vector<double> {
            const int max_iter = prm_.get<int>("solver.max_iter");
            const double tol = prm_.get<double>("solver.tolerance");//,1e-8);
            const int nlin_verbosity = prm_.get<double>("solver.verbosity");
            *trimesh_ = trimesh;
            // save well sources before grid change
            std::vector<CellRef> wsources = well_source_cellref_; 
            for (auto& cell : wsources) 
                cell = RegularTrimesh::fine_to_coarse(cell, level);
            
            // setup fracture with new grid
            const int MAX_NUM_COARSENING = prm_.get<int>("solver.max_num_coarsening"); // should be enough for all practical purposes
            const int numcell_threshold = prm_.get<int>("solver.numcell_threshold");
            bool smooth_boundary = prm_.get<bool>("solver.smooth_boundary",false);
            auto [grid, fsmap, bmap] =
            trimesh_->createDuneGrid(MAX_NUM_COARSENING, wsources,/* smoothed triangels */ smooth_boundary, numcell_threshold); // well cells kept intact!
            auto org_map = grid_mesh_map_; 
            grid_mesh_map_ = fsmap;
            setFractureGrid(std::move(grid)); // true -> coarsen interior
            // generate the inverse map of fsmap_ (needed below)
            std::map<CellRef, size_t> fsmap_inv;
            for (int i = 0; i != fsmap.size(); ++i)
                if (size(fsmap[i]) == 1) // a fine-scale cell
                    fsmap_inv[fsmap[i].front()] = i;

            // update indices for well sources to the correct cells in the new grid
            well_source_.clear();
            for (const auto& cell : wsources) 
                well_source_.push_back(fsmap_inv[cell]);

            // Update the rest of the fracture object to adapt to grid change
            updateReservoirCells(cell_search_tree, cpgrid, entity_seeds);
            updateReservoirProperties<TypeTag, Simulator>(simulator, true, false);
            initPressureMatrix();
            //initFractureWidth();
            //initFracturePressureFromReservoir();
            rhs_pressure_.resize(0);
            coupling_matrix_ = nullptr;



            // solve flow-mechanical system
            bool remap_solution= prm_.get<bool>("solver.remap_solution");
            if (remap_solution) {

                for (size_t i = 0; i < fracture_pressure_.size(); ++i) {
                    assert(std::abs(fracture_pressure_[i][0]) < 1e10);
                }
                for (size_t i = 0; i < fracture_width_.size(); ++i) {
                    assert(std::abs(fracture_width_[i][0]) < 0.6);
                }
                auto old_fracture_width_ = fracture_width_;
                redistribute_values(fracture_width_, org_map, fsmap, level, /*point_wise*/true);
                double well_value = 0.0;
                if (numWellEquations() > 0){
                    well_value = fracture_pressure_[fracture_pressure_.size() - 1];
                    fracture_pressure_.resize(fracture_pressure_.size() - 1);
                }
                redistribute_values(fracture_pressure_, org_map, fsmap, level, /*point_wise*/true);
                if (numWellEquations() > 0){
                    assert(fracture_pressure_.size() == fracture_width_.size());
                    fracture_pressure_.resize(fracture_pressure_.size() + 1);
                    fracture_pressure_[fracture_pressure_.size() - 1] = well_value;
                }
                for (size_t i = 0; i < fracture_pressure_.size(); ++i) {
                    assert(std::abs(fracture_pressure_[i][0]) < 1e10);
                }
                for (size_t i = 0; i < fracture_width_.size(); ++i) {
                    assert(std::abs(fracture_width_[i][0]) < 0.6);
                }
                // solve flow-mechanical system
                // cell entries only: reservoir_pressure_ has no well-equation entry
                const size_t num_press_cells = fracture_pressure_.size() - numWellEquations();
                for (size_t i = 0; i < num_press_cells; ++i) {
                    if(fracture_pressure_[i] == 0.0){
                      // not initialized values
                      fracture_pressure_[i] = reservoir_pressure_[i];
                    }
                }
                if (numWellEquations() > 0 && fracture_pressure_[num_press_cells] == 0.0) {
                    fracture_pressure_[num_press_cells] = fracture_pressure_[0];
                }
            } else {
                initFractureWidth();
                initFracturePressureFromReservoir();
            }

            // total amount of filtercake be
        

            filtercake_thikness_ = redistribute_values(filtercake_thickness_0,
                                                       grid_mesh_map_0,
                                                       fsmap, level, /*point_wise*/ true);
            double total_filtercake = filterCakeVolume();
            if(std::abs(total_filtercake - total_filtercake_0) > 1e-6 && (nlin_verbosity > 1)){
                std::cout << "Total filtercake thickness before remapping: " << total_filtercake_0 << std::endl;
                std::cout << "Total filtercake thickness after remapping: " << total_filtercake << std::endl;
                //assert(false);
            }
            int iter = 0;
            while (!fullSystemIteration(tol,iter) && iter++ < max_iter) {
            };
            // Contact-set stability vs the previous coupling round (evaluated on
            // the first solve of this round, i.e. before any expansion this
            // round; a grid-size mismatch means the mesh moved: treat as unstable).
            if (vol_pacing) {
                // eq. (1)-(2): V(k)-V0 vs (q_perf - Q_leak(k))*dt*f_vol
                const FractureProperties fp = calculateFractureProperties();
                // budget source: this perforation's current CTF share (default) or
                // the whole well (solver.volume_pacing_source=well) — the share is
                // circular for an establishing fracture (small fracture -> small
                // share -> halted), the whole well is what an open fracture takes
                const std::string vol_src = prm_.get<std::string>("solver.volume_pacing_source", "perf");
                const double q_perf = (vol_src == "well") ? std::abs(well_rate_) : well_perf_rate_;
                const double v_avail = std::max(0.0, (q_perf - fp.flux) * vol_dt) * vol_fac;
                const double dV = fp.volume - vol_V0;
                vol_budget_spent = dV > v_avail;
                if (vol_verb > 0) {
                    std::cout << "VolumePacing: V=" << fp.volume << " V0=" << vol_V0
                              << " dV=" << dV << " q_perf=" << q_perf << " Q_leak=" << fp.flux
                              << " V_avail=" << v_avail << (vol_budget_spent ? "  -> BUDGET SPENT" : "")
                              << std::endl;
                }
                if (vol_pressure && dV > v_avail && fp.volume > 0.0) {
                    // eq. (3)-(4): draw the BC pressure down by the over-draw over the
                    // secant storage compliance C_eff = V / (p_perf - p_res)
                    const int pc = perfinj_.empty() ? 0 : std::get<0>(perfinj_[0]);
                    const double p_res = (pc < static_cast<int>(reservoir_pressure_.size()))
                        ? reservoir_pressure_[pc] : 0.0;
                    const double dp0 = std::max(vol_perf_p0 - p_res, 1.0e5);
                    const double c_eff = fp.volume / dp0;
                    const double p_new = std::max(vol_perf_p0 - (dV - v_avail) / c_eff, p_res);
                    if (vol_verb > 0) {
                        std::cout << "VolumePacing: perf pressure " << perf_pressure_ / 1e5
                                  << " -> " << p_new / 1e5 << " bar (C_eff=" << c_eff << ")" << std::endl;
                    }
                    perf_pressure_ = p_new;
                }
            }
            // remove closed cells which was not open before ??
            if(iter >= max_iter){
                last_solve_stats_.converged = false;
                bool failure_on_convergence = prm_.get<bool>("solver.failure_on_nonconvergence",false);
                if(failure_on_convergence){
                    throw std::runtime_error("Fracture solver did not converge within " + std::to_string(max_iter) + " iterations.");
                }else{
                    std::cout << "Warning: Fracture solver did not converge within " << max_iter << " iterations." << std::endl;
                }
            }
            if(nlin_verbosity > 0){
                int num_closed_cells = std::accumulate(closed_cells_.begin(), closed_cells_.end(), 0);
                std::cout << "Nonlinear iterations needed with fixed expansion: " << iter << std::endl;
                std::cout << " fracture cells: " << fracture_width_.size();
                std::cout << " closed cells: " << num_closed_cells;
                std::cout << std::endl;
             }

            // compute K1 stress intensity
            const std::vector<double> K1_not_nan = Fracture::stressIntensityK1();
            const std::vector<CellRef> boundary_cells = trimesh_->boundaryCells();
            std::vector<double> result(boundary_cells.size());
            for (size_t i = 0; i != result.size(); ++i){
                const double KImax = reservoir_cstress_[bmap[boundary_cells[i]]];
                const double KI = K1_not_nan[bmap[boundary_cells[i]]];
                result[i] = KI/KImax;
            }
            //
            const double force_limit = prm_.get<double>("solver.force_limit", 0.0);
            for (size_t i = 0; i != result.size(); ++i){
                int bind = bmap[boundary_cells[i]];
                if(fractureForce(bind) > force_limit ){
                    result[i] = -1.0; // do not expand here
                }
            }

            // Conservative propagation (opt-in, default off): do not let the
            // fracture grow off a non-converged inner solve. The K1 scores come
            // from the current width/pressure state; if that solve did not
            // converge they are unreliable, and opening/propagating on them can
            // irreversibly alter the solution (an opened cell changes the global
            // mechanics). Suppress expansion this round and let growth happen on a
            // later, converged step instead of relying on an expensive roll-back.
            if (prm_.get<bool>("solver.conservative_propagation", false)
                && iter >= max_iter) {
                std::fill(result.begin(), result.end(), -1.0);
            }


            // Per-cell propagation veto (opt-in): refuse to grow from a front cell
            // that is closed or whose FB complementarity has not converged, instead
            // of vetoing ALL growth when the solve hit its iteration cap
            // (conservative_propagation) - a global veto starves the fracture while
            // the well keeps pressurising. Only front cells are scored at all, so
            // this leaves every trustworthy open front free to propagate.
            if (prm_.get<bool>("solver.propagate_converged_cells_only", false)) {
                const double fbtol = prm_.get<double>("solver.propagate_fb_tol", 1e-3);
                int vetoed = 0;
                for (size_t i = 0; i != result.size(); ++i) {
                    const int bind = bmap[boundary_cells[i]];
                    if (bind < 0) continue;
                    const bool closed = bind < static_cast<int>(closed_cells_.size())
                                        && closed_cells_[bind];
                    const bool unconverged = bind < static_cast<int>(fb_cell_residual_.size())
                                             && fb_cell_residual_[bind] > fbtol;
                    if (closed || unconverged) {
                        result[i] = -1.0;
                        ++vetoed;
                    }
                }
                if (vetoed > 0 && prm_.get<int>("verbosity", 0) > 1)
                    std::cout << "Propagation veto: " << vetoed << "/" << result.size()
                              << " front cells closed or unconverged" << std::endl;
            }

            // Fixed-topology mode: freeze the mesh at the seed (Reveal comparison).
            if (prm_.get<bool>("solver.disable_propagation", false)) {
                std::fill(result.begin(), result.end(), -1.0);
            }

            // Only propagate where the front is OPEN to the fresh rock (opt-in):
            // a closed front cell carries no mechanical opening load, so growth
            // from it would ride on the width-floor conductivity, not on K1.
            // Combine with solver.conservative_propagation so new cells are only
            // added off a converged open/closed set.
            if (prm_.get<bool>("solver.propagate_open_front_only", false)) {
                for (size_t i = 0; i != result.size(); ++i) {
                    const int bind = bmap[boundary_cells[i]];
                    if (bind >= 0 && bind < static_cast<int>(closed_cells_.size())
                        && closed_cells_[bind])
                        result[i] = -1.0;
                }
            }

            return result;
        };

        const FractureProperties& fprop_before = calculateFractureProperties();
        auto stop_criterion= [&]() -> bool
        { 
          const bool area_stop = expantionMax(fprop_before);
          if (vol_pacing && vol_budget_spent) {
              if (vol_verb > 0) std::cout << "VolumePacing: stopping expansion (budget spent)" << std::endl;
              return true;
          }
          return area_stop;
        };

        //const double K1max = prm_.get<double>("KMax");
        const double threshold = 1.0;
        const std::vector<CellRef> fixed_cells = well_source_cellref_;
        const int target_cellcount = prm_.get<int>("solver.target_cellcount"); 
        const int cellcount_threshold = prm_.get<int>("solver.cellcount_threshold");
        const int max_iter_on_same_level = prm_.get<int>("solver.max_iter_on_same_level");
        const auto initial_mesh= *trimesh_;
        auto [mesh, cur_level] =
            expand_to_criterion(*trimesh_, score_function, threshold,
                                fixed_cells,
                                target_cellcount,
                                cellcount_threshold,
                                max_iter_on_same_level, stop_criterion);

                

        
        bool any_removed = removeNewZeroWithCells(mesh,cur_level, initial_mesh);
        // this will set trimesh_,change dune grid  and solve
        if(any_removed){
          mesh.removeSawtooths();
          score_function(mesh, cur_level);
        }


            // make current level become the reference (finest) level
            // note that the well_source_cellref_ is already set from the last call to the score function
            for (auto& cell : well_source_cellref_)
                cell = RegularTrimesh::fine_to_coarse(cell, cur_level);

            // volume pacing: the drawn-down BC was for the expansion rounds only;
            // the coupling continues from the flow's perf pressure
            if (vol_pacing && vol_pressure) perf_pressure_ = vol_perf_p0;

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
        else if (method == "if_propagate")
        {
            // ----------------------------------------------------------------------------
            // iterate full nonlinear system until convergence, and expand fracture if necessary
            if (false) {
                fracture_width_ = 1e-2; // Ensure not completely closed
                fracture_pressure_ = 0.0;
            }

            // start by assuming pressure equal to confining stress (will also set
            // fracture_pressure_ to its correct size
            normalFractureTraction(fracture_pressure_);
            if (numWellEquations() > 0) { // @@ it is implicitly assumed for now that
                // there is just one well equation.  We initializze
                // it with an existing value.
                fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];
            }
            const int max_iter = prm_.get<int>("solver.max_iter");
            const double tol = prm_.get<double>("solver.tolerance"); //,1e-8);


            const double efac = prm_.get<double>("solver.efac"); // 2; // @@ heuristic
            const double rfac = prm_.get<double>("solver.rfac"); // 2; // @@ heuristic
            double K1max = prm_.get<double>("KMax"); // @@ for testing.  Should be added as a proper data member
            const std::vector<size_t> boundary_cells = grid_stretcher_->boundaryCellIndices();
            const size_t N = boundary_cells.size(); // number of boundary nodes and boundary cells

            std::vector<double> total_bnode_disp(N, 0), bnode_disp(N, 0), cell_disp(N, 0);
            // const std::vector<GridStretcher::CoordType>
            //   bnode_normals_orig = grid_stretcher_->bnodenormals();
            int max_expand_iter = prm_.get<int>("solver.max_expand_iter");
            std::vector<GridStretcher::CoordType> displacements(N, {0, 0, 0});
            int count = 0; // @@
            while (true && (count < max_expand_iter)) {

                std::cout << "Iteration: " << ++count << std::endl;
                // solve flow-mechanical system
                int iter = 0;
                // open file "width" for appending data
                while (!fullSystemIteration(tol,iter) && iter++ < max_iter) {
                    // @@@@
                    // auto fw = make_vector(fracture_width_);
                    // auto fp = make_vector(fracture_pressure_, fracture_pressure_.size()-1);
                    // grid_stretcher_->dumpToVTK("stretchedgrid", { fw, fp });

                    // if (iter > 20) {
                    //   std::ofstream width_debug("width", std::ios::app);
                    //   std::ofstream pressure_debug("pressure", std::ios::app);

                    //   std::copy(fw.begin(), fw.end(), std::ostream_iterator<double>(width_debug, " "));
                    //   std::copy(fp.begin(), fp.end(), std::ostream_iterator<double>(pressure_debug, " "));

                    //   int krull=0;
                    //}
                };
                if(iter >= max_iter){
                    last_solve_stats_.converged = false;
                    bool failure_on_convergence = prm_.get<bool>("solver.failure_on_nonconvergence",false);
                    if(failure_on_convergence){
                        throw std::runtime_error("Fracture solver did not converge within " + std::to_string(max_iter) + " iterations.");
                    }else{
                        std::cout << "Warning: Fracture solver did not converge within " << max_iter << " iterations." << std::endl;
                    }
                }
                std::cout << "Iterations needed: " << iter << std::endl;


                // identify where max stress intensity is exceeded and propagation is needed
                const auto dist = grid_stretcher_->centroidEdgeDist();
                fill(bnode_disp.begin(), bnode_disp.end(), 0.0);

                const std::vector<double> K1_not_nan = Fracture::stressIntensityK1();
                std::vector<double> K1;
                bool should_fracture = false;
                for (size_t i = 0; i != K1_not_nan.size(); ++i) {
                    if (!std::isnan(K1_not_nan[i])) {
                        K1max = reservoir_cstress_[i];
                        if (K1_not_nan[i] > K1max) {
                            should_fracture = true;
                        }
                    }
                }
                if (!should_fracture) {
                    break;
                }

                // loop over cells, determine how much they should be expanded or contracted
                const double maxgrow = rfac * grid_stretcher_->maxBoxLength();
                for (size_t i = 0; i != N; ++i) {
                    K1max = reservoir_cstress_[boundary_cells[i]];
                    cell_disp[i] = efac
                        * (compute_target_expansion(K1max, fracture_width_[boundary_cells[i]], E_, nu_) - dist[i]);
                    cell_disp[i] = std::max(std::min(cell_disp[i], maxgrow), -maxgrow);
                }

                bnode_disp = grid_stretcher_->computeBoundaryNodeDisplacements(cell_disp); //@@
                // bnode_disp =
                //   grid_stretcher_->computeBoundaryNodeDisplacements(cell_disp, bnode_normals_orig);
                for (size_t i = 0; i != N; ++i)
                    bnode_disp[i] = std::max(std::min(bnode_disp[i], maxgrow), -maxgrow);

                // // ensure no boundary node moved inwards further than starting point
                // for (size_t i = 0; i != N; ++i) {
                //   bnode_disp[i] = std::max(bnode_disp[i], -total_bnode_disp[i]);
                //   total_bnode_disp[i] += bnode_disp[i];
                // } // @@@ no longer works after rebalanceBoundary introduced()

                // ensure convexity
                grid_stretcher_->adjustToConvex(bnode_disp, total_bnode_disp, grid_stretcher_->bnodenormals());
                // bnode_normals_orig);

                for (size_t i = 0; i != N; ++i)
                    displacements[i] = grid_stretcher_->bnodenormals()[i] * bnode_disp[i];
                // displacements[i] = bnode_normals_orig[i] * bnode_disp[i];

                grid_stretcher_->applyBoundaryNodeDisplacements(displacements);
                grid_stretcher_->rebalanceBoundary();

                // debug stuff
                // grid_stretcher_->dumpToVTK("stretchedgrid");

                // grid has changed its geometry, so we have to recompute discretizations
                updateCellNormals();
                updateReservoirCells(cell_search_tree, cpgrid, entity_seeds);
                updateReservoirProperties<TypeTag, Simulator>(simulator, true, false);
                initPressureMatrix();
                fracture_matrix_ = nullptr;

                //const auto& pts = grid_stretcher_->nodecoords();
                //const auto bix = grid_stretcher_->boundaryNodeIndices();
                // for(auto i : bix)
                //   os << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
                // os.close();
            }
            if (count >= max_expand_iter) {
                last_solve_stats_.converged = false;
                std::cout << "Fracture expansion did not converge within the maximum number of iterations" << std::endl;
            }
        }
        else
        {
            OPM_THROW(std::runtime_error, "Unknowns solution method");
        }
        bool first_solve = well_indices_.size() == 0;
        if(first_solve){
            well_indices_.resize(2);
        }
        auto well_indices_new = wellIndices_();//calculate wellindices
        auto well_indices_old = wellIndices();// use the function used from simulator at prevois step
        if(first_solve){
            well_indices_[0] = well_indices_new;
            well_indices_[1] = well_indices_new;
        }else{
            //well_indices_[1] = well_indices_[0];
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

} // namespace Opm
