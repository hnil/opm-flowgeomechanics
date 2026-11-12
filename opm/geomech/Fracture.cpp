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

#include "config.h"

#include <set>

#include <iostream> // std::cout
#include <sstream>
#include <string> // std::string

#include <opm/geomech/Fracture.hpp>

#include <opm/geomech/Math.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/wells/ConnFracStatistics.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/common/fmatrixev.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/istl/io.hh> // needed for printSparseMatrix??
#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/FracturePressureAssemblerAD.hpp>

#include <opm/grid/polyhedralgrid.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <limits>
#include <numeric>
#include <cassert>



namespace
{


}; // end anonymous namespace

namespace Opm
{
void
Fracture::init(const std::string& well,
               const int perf,
               const int well_cell,
               int global_index,
               const int segment,
               const std::optional<std::pair<double, double>>& perf_range,
               const Point3D& origo,
               const Point3D& normal,
               const Opm::PropertyTree& prm,
               FractureConnectionId connection_id)
{
    OPM_TIMEFUNCTION();
    max_flow_time_step_ = std::numeric_limits<double>::max();
    prm_ = prm;
    if(prm_.get<bool>("config.gravity_off", true)){
        gravity_ = 0.0;
    }

    min_width_ = prm_.get<double>("config.min_width", 1e-3);
    wellinfo_ = WellInfo({well, perf, well_cell, global_index,
        connection_id,
        segment, perf_range});

    origo_ = origo;
    axis_[2] = normal;
    axis_[0] = Point3D({std::copysign(normal[2], normal[0]),
                        std::copysign(normal[2], normal[1]),
                        -std::copysign(std::abs(normal[0]) + std::abs(normal[1]), normal[2])});
    axis_[1] = crossProduct(axis_[2], axis_[0]);
    double init_scale = prm_.get<double>("config.axis_scale");
    for (int i = 0; i < 3; ++i) {
        axis_[i] /= axis_[i].two_norm();
        naxis_[i] = axis_[i];
        axis_[i] *= init_scale;
    }

    const std::string method = prm_.get<std::string>("solver.method");
    if (method != "if_propagate_trimesh") {
        OPM_THROW(std::runtime_error,
                  "Unknown fracture solver.method '" + method + "': use if_propagate_trimesh, "
                  "with solver.disable_propagation=true for a fixed mesh");
    }
    {
        double trires = prm_.get<double>("config.trires");
        const double edgelen = init_scale / trires;
        const double radius = init_scale;
        const double fac = std::sqrt(3) / 2;
        const std::array<double, 3> ax1 {axis_[0][0], axis_[0][1], axis_[0][2]};
        const std::array<double, 3> ax2 {0.5 * ax1[0] + fac * axis_[1][0],
                                         0.5 * ax1[1] + fac * axis_[1][1],
                                         0.5 * ax1[2] + fac * axis_[1][2]};

        OpmLog::info(fmt::format("Creating trimesh with radius: {}, edgelen: {}, ax1: {},{},{}, ax2: {},{},{}",
                                 radius, edgelen, ax1[0], ax1[1], ax1[2], ax2[0], ax2[1], ax2[2]));

        trimesh_ = std::make_unique<RegularTrimesh>(radius,
                                                    std::array {origo_[0], origo_[1], origo_[2]},
                                                    ax1,
                                                    ax2,
                                                    std::array {edgelen, edgelen});

        trimesh_->removeSawtooths();

        // identify well cells (since this is not done in setFractureGrid when providing
        // a user-defined grid)
        well_source_cellref_ = RegularTrimesh::inner_ring_cells();

        for (const auto& cell : well_source_cellref_)
            well_source_.push_back(trimesh_->linearCellIndex(cell));

        auto [grid, fsmap, bmap] = trimesh_->createDuneGrid(1, well_source_cellref_);
        grid_mesh_map_ = fsmap;
        grid_fixed_cells_ = well_source_cellref_;

        setFractureGrid(std::move(grid)); // create the physical grid from trimesh
    }


    // NB: The Fracture object is still not fully initialized, since there is
    // not yet any connection with a surrounding reservoir. The following
    // needs to be done before the fracture can be used:
    // 1) Call `updateReservoirCells` to establish the mapping between frature grid
    //    cells and the cells in the reservoir grid (sets `reservoir_map_.cells`)
    // 2) Then, call `updateReservoirProperties` to import relevant reservoir properties
    //    to the fracture (`reservoir_XXX_` vectors, as well as `E_` and `nu_`)
}

void
Fracture::resetWriters()
{
    // nead to be reseat if grid is changed ??
    vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(),
                                                                       Dune::VTK::nonconforming);

    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename") + this->name();
    std::string multiFileName = "";
    if (!vtkmultiwriter_) {
        vtkmultiwriter_ = std::make_unique<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>>(
            /*async*/ false, grid_->leafGridView(), outputdir, simName, multiFileName);
    } else {
        vtkmultiwriter_->gridViewChanged(grid_->leafGridView()); // need to be called if grid is changed
    }
}

/**
 * @brief Removes cells from the grid if they are out side reservoir.
 */
void
Fracture::removeCells()
{
    // copy all to presistent container
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout()); // used id sets interally
    Dune::PersistentContainer<Grid, double> fracture_width(*grid_, 0);
    fracture_width.resize();
    for (const auto& elem : elements(grid_->leafGridView())) {
        size_t eIdx = mapper.index(elem);
        fracture_width[elem] = fracture_width_[eIdx];
    }

    for (const auto& elem : elements(grid_->leafGridView())) {
        size_t eIdx = mapper.index(elem);
        if (reservoir_map_.cells[eIdx] < 0) {
            grid_->removeElement(elem);
        }
    }

    grid_->grow();
    grid_->postGrow();
    // resize array

    fracture_width_.resize(numFractureCells());
    // copy back from presistent contatiner
    for (const auto& elem : elements(grid_->leafGridView())) {
        size_t eIdx = mapper.index(elem);
        fracture_width_[eIdx] = fracture_width[elem];
    }
    this->resetWriters();
}

void
Fracture::updateFilterCakeProps(const Opm::Connection& connection,
                                const Opm::SingleWellState<double,Fracture::IndexTraits>& wellstate,double dt)
{
    OPM_TIMEFUNCTION();
    // potentially remap filtercake thikness if grid has changed
    assert(filtercake_thikness_.size() == reservoir_map_.cells.size());
    const auto& perfdata = wellstate.perf_data;
    std::map<int, double> WI_fluxes;
    int water_index = 0;
    int np = 3;

    double qw_sum = 0.0;
    double qwpos_sum = 0.0;
    // should we do this based on fracture cross flow?
    for (size_t i = 0; i < perfdata.cell_index.size(); ++i) {
        int cell_index = perfdata.cell_index[i];
        // similar as in WellFilterCake.cpp:148
        const auto& connection_rates = perfdata.phase_rates;
        const double water_rate = std::max(0.0, connection_rates[i * np + water_index]);
        qwpos_sum += std::max(0.0, water_rate);
        qw_sum += water_rate;
        WI_fluxes[cell_index] = water_rate;
    }
    double xfact = 1.0;
    if (qwpos_sum > 1.0e-12) {
        xfact = qw_sum / qwpos_sum;
    }else{
      if(qwpos_sum == qw_sum){
        xfact = 1.0;
      }else{
        xfact = 0.0;
      }
    }
    //NB need to be checked
    bool scale_filtrate = prm_.get<bool>("config.scale_filtrate", false);
    double scale_factor = 1.0;
    if(scale_filtrate){
        // overkill but simplest
        auto frac_prop = this->calculateFractureProperties();
        double frac_flux = frac_prop.flux; // m3/s
        double sim_frac_flux = wellstate.frac_rate;
        //NB maybe we should only use posive fluxes (cross flow problems)
        if(frac_flux > 1.0e-12){
            scale_factor = std::max(0.0,sim_frac_flux) / frac_flux;
        }
    }
    double filtrate_conn = wellstate.filtrate_conc*xfact;
    has_filtercake_ = connection.filterCakeActive();
    if (has_filtercake_) {
        auto& filter_cake = connection.getFilterCake();
        const double poro = filter_cake.poro;
        const double perm = filter_cake.perm;
        filtercake_poro_ = poro;
        filtercake_perm_ = perm;

        std::map<int, double> reservoir_areas;
        std::map<int, double> reservoir_flux;
        if ((leakof_.size()
             > 0)) { // should be first step with seed fracture is clean could have used rates from solve
            std::vector<double> leakof_rate = leakOfRate();
            assert(leakof_.size() == fracture_pressure_.size()-numWellEquations());
            ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
            for (auto& element : Dune::elements(grid_->leafGridView())) {
                int eIdx = mapper.index(element);
                auto geom = element.geometry();
                double area = geom.volume(); // is the area of this face
                int res_cell = reservoir_map_.cells[eIdx];
                reservoir_areas[res_cell] += area;
                double flux = leakof_rate[eIdx]*area; // m3/s
                reservoir_flux[res_cell] += flux;
                if (flux < 0) {
                    if(prm_.get<bool>("verbose",false)){
                        OpmLog::debug(fmt::format("Negative flux {} for element index {} with reservoir cell {}",
                                                  flux, eIdx, res_cell));
                    }
                    flux = 0.0;
                }
                double water_rate = flux;
                double water_rate_scaled = water_rate*scale_factor;// account for difference fracture and flow solution
                double filtrate_volume = water_rate_scaled*filtrate_conn* dt; // m3
                double dh = filtrate_volume / (area * (1 - filtercake_poro_));
                assert(dh >= 0);
                filtercake_thikness_[eIdx] += dh;
            }
        }

        /// for setting or controling the fluxes between multiphase and single phase
        ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
        for (const auto& [res_cell, flux] : reservoir_flux) {
            if (reservoir_areas.find(res_cell) == reservoir_areas.end()) {
                OpmLog::warning(fmt::format("Reservoir area not found for element index {}", res_cell));
                continue;
            }
            double frac_flux = reservoir_flux[res_cell];
            if (res_cell != wellinfo_.well_cell) {
                if (std::abs(flux - WI_fluxes[res_cell]) > 0) {
                    OpmLog::debug(fmt::format("Fracture flux differs from flow flux {}\nFlux: frac {} vs res {}",
                                              res_cell, frac_flux, WI_fluxes[res_cell]));
                }
            } else {
                OpmLog::debug(fmt::format("Total WI flux {} for cell {} matches fracture flux {}",
                                          WI_fluxes[res_cell], res_cell, frac_flux));
            }
        }
    }
}



Dune::BlockVector<Dune::FieldVector<double, 3>>
Fracture::all_slips() const
{
    Dune::BlockVector<Dune::FieldVector<double, 3>> slips(grid_->leafGridView().size(0));
    slips = 0;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& elem : Dune::elements(grid_->leafGridView())) {
        size_t eIdx = mapper.index(elem);
        // only normal slip for now
        slips[eIdx][0] = fracture_width_[eIdx];
    }
    return slips;
}

Dune::FieldVector<double, 3>
Fracture::disp(Dune::FieldVector<double, 3> obs) const
{
    auto slips = all_slips();
    Dune::FieldVector<double, 3> disp = ddm::disp(obs, slips, *grid_, E_, nu_);
    return disp;
}

Dune::FieldVector<double, 6>
Fracture::strain(Dune::FieldVector<double, 3> obs) const
{
    // for now use full slip in interface even we only calculate normal slip
    Dune::BlockVector<Dune::FieldVector<double, 3>> slips = all_slips();
    Dune::FieldVector<double, 6> strain = ddm::strain(obs, slips, *grid_, E_, nu_);
    return strain;
}

Dune::FieldVector<double, 6>
Fracture::stress(Dune::FieldVector<double, 3> obs) const
{
    // for now use full slip in interface even we only calculate normal slip
    Dune::BlockVector<Dune::FieldVector<double, 3>> slips = all_slips();
    Dune::FieldVector<double, 6> strain = ddm::strain(obs, slips, *grid_, E_, nu_);
    Dune::FieldVector<double, 6> stress = ddm::strainToStress(E_, nu_, strain);
    return stress;
}

std::string
Fracture::name() const
{
    std::string name
        = "Fracture_on_" + wellinfo_.name + "_perf_" + std::to_string(wellinfo_.perf) + "_nr";
    return name;
}

void
Fracture::updateCellNormals()
{
    cell_normals_.resize(numFractureCells());
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : elements(grid_->leafGridView()))
        cell_normals_[elemMapper.index(element)] = ddm::normalOfElement(element);
}

void
Fracture::setFractureGrid(std::unique_ptr<Fracture::Grid> gptr)
{
    assert(gptr != nullptr);
    grid_ = std::move(gptr);
    // compute the cell normals
    updateCellNormals();

    // set sparsity of pressure matrix, but do not compute its entries

    // Since the grid has been created/updated/changed, any previous mapping to
    // reservoir cells has been invalidated, and the fracture matrix (for
    // mechanics) is obsolete.
    reservoir_map_.cells.clear();
    fracture_matrix_ = nullptr;

    this->resetWriters();
}

std::vector<double>
Fracture::stressIntensityK1() const
{
    size_t nc = numFractureCells();
    std::vector<double> stressIntensityK1(nc, std::nan("0"));
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());

    // Opt-in displacement-correlation K1 extraction: instead of evaluating the
    // asymptotic opening w(s) = C*sqrt(s) from the single boundary cell (legacy
    // "aperture_cell", mesh sensitive), fit C by least squares over the cells in
    // a corridor of physical length L behind the front edge and convert with the
    // same calibrated constant (fractureK1(1,C) == C * mu*sqrt(pi)/(2(1-nu)*1.834)).
    // With one point the fit reduces exactly to the legacy formula.
    const std::string k1_mode = prm_.get<std::string>("solver.k1_extraction", "aperture_cell");
    const bool correlate = (k1_mode == "displacement_correlation");
    if (k1_mode != "aperture_cell" && !correlate) {
        OPM_THROW(std::runtime_error, "Unknown solver.k1_extraction: " + k1_mode);
    }
    const double corr_length = prm_.get<double>("solver.k1_correlation_length", 0.0);
    const double corr_cells = prm_.get<double>("solver.k1_correlation_cells", 3.0);

    // cell adjacency, centers and areas (only needed for the correlation mode)
    std::vector<std::vector<int>> nbs;
    std::vector<Dune::FieldVector<double, 3>> centers;
    std::vector<double> areas;
    if (correlate) {
        nbs.resize(nc);
        centers.resize(nc);
        areas.resize(nc);
        for (const auto& elem : elements(grid_->leafGridView())) {
            const int ci = mapper.index(elem);
            centers[ci] = elem.geometry().center();
            areas[ci] = elem.geometry().volume();
            for (auto& is : Dune::intersections(grid_->leafGridView(), elem)) {
                if (is.neighbor()) {
                    nbs[ci].push_back(mapper.index(is.outside()));
                }
            }
        }
    }

    for (const auto& elem : elements(grid_->leafGridView())) {
        for (auto& is : Dune::intersections(grid_->leafGridView(), elem)) {
            if (is.boundary()) {
                int nIdx = mapper.index(elem);
                auto isCenter = is.geometry().center();
                auto elCenter = elem.geometry().center();
                auto vecC = isCenter - elCenter;
                auto distC = vecC.two_norm();
                if (!correlate) {
                    double K1 = ddm::fractureK1(distC, fracture_width_[nIdx], E_, nu_);
                    stressIntensityK1[nIdx] = K1;
                    continue;
                }

                // inward normal of the front edge and corridor dimensions
                auto inward = elCenter - isCenter;
                inward /= inward.two_norm();
                const double hloc = std::sqrt(areas[nIdx]);
                const double L = (corr_length > 0.0) ? corr_length : corr_cells * hloc;
                const double half_width = 0.9 * hloc;

                // BFS over neighbors collecting cells inside the corridor
                std::vector<int> stack {nIdx};
                std::vector<char> seen(nc, 0);
                seen[nIdx] = 1;
                double sum_ws = 0.0; // sum w_i*sqrt(s_i)
                double sum_s = 0.0; // sum s_i
                int npts = 0;
                while (!stack.empty()) {
                    const int c = stack.back();
                    stack.pop_back();
                    auto d = centers[c] - isCenter;
                    const double s = d.dot(inward);
                    auto lat = d;
                    lat.axpy(-s, inward);
                    if (s > 1e-3 * hloc && s <= L && lat.two_norm() <= half_width) {
                        sum_ws += fracture_width_[c][0] * std::sqrt(s);
                        sum_s += s;
                        ++npts;
                    }
                    // expand while still within reach of the corridor
                    if (s <= L) {
                        for (const int nb : nbs[c]) {
                            if (!seen[nb]) {
                                seen[nb] = 1;
                                stack.push_back(nb);
                            }
                        }
                    }
                }

                if (npts >= 2 && sum_s > 0.0) {
                    const double C = sum_ws / sum_s; // LSQ fit of w = C*sqrt(s)
                    stressIntensityK1[nIdx] = ddm::fractureK1(1.0, C, E_, nu_);
                } else {
                    // not enough points (tiny fracture): legacy single-cell value
                    stressIntensityK1[nIdx] = ddm::fractureK1(distC, fracture_width_[nIdx], E_, nu_);
                }
            }
        }
    }
    return stressIntensityK1;
}

void
Fracture::summary_of_solve(){
  if(numWellEquations()>0){
    int nc = numFractureCells();
    assert( int(fracture_pressure_.size()) == nc + 1);
    double well_pressure = fracture_pressure_[nc];
    const int cell = std::get<0>(perfinj_[0]); // TODO: verify this is the correct index
    const double pres = reservoir_.pressure[cell];
    const double dp = well_pressure - pres;
    FractureProperties fprop = calculateFractureProperties();
    // fracture rate
    OpmLog::debug(fmt::format("Reservoir pressure at perf cell {} is {} dp {} frac flux {}\n"
                              "Fracture properties:\n"
                              " Height: {} Length: {} Area: {}\n"
                              " Volume: {} Avg width: {}\n"
                              " WI: {} Inj pressure: {} Inj bhp: {} Inj rate: {}",
                              cell, pres, dp, fprop.flux,
                              fprop.height, fprop.length, fprop.area,
                              fprop.volume, fprop.avg_width,
                              fprop.WI, fprop.inj_pressure, fprop.inj_bhp, fprop.inj_wellrate));
  }
}


void
Fracture::write(int reportStep) const
{
    OPM_TIMEFUNCTION();
    std::vector<double> K1; // need to be in scope until written
    auto tmp(fracture_pressure_);
    tmp.resize(numFractureCells());
    if (reservoir_map_.cells.size() > 0) {
        vtkwriter_->addCellData(reservoir_map_.cells, "ReservoirCell");
    }
    if (reservoir_.perm.size() > 0) {
        vtkwriter_->addCellData(reservoir_.perm, "ReservoirPerm");
    }
    if (reservoir_.cstress.size() > 0) {
        vtkwriter_->addCellData(reservoir_.cstress, "ReservoirCStress");
    }
    if (reservoir_map_.dist.size() > 0) {
        vtkwriter_->addCellData(reservoir_map_.dist, "ReservoirDist");
    }
    if (fracture_pressure_.size() > 0) {
        // fracture_pressure_ may carry an extra well-pressure entry, so pass a
        // copy truncated to the cell count
        vtkwriter_->addCellData(tmp, "FracturePressure");
    }
    if (reservoir_.pressure.size() > 0) {
        vtkwriter_->addCellData(reservoir_.pressure, "ReservoirPressure");
    }
    if (reservoir_.mobility.size() > 0) {
        vtkwriter_->addCellData(reservoir_.mobility, "ReservoirMobility");
    }
    if (fracture_width_.size() > 0) {
        vtkwriter_->addCellData(fracture_width_, "FractureWidth");
        K1 = this->stressIntensityK1();
        vtkwriter_->addCellData(K1, "stressIntensityK1");
    }
    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename");
    std::string filename = outputdir + "/Static_" + simName + this->name();
    if (reportStep > 0) {
        filename = filename + "_step_" + std::to_string(reportStep);
    }
    vtkwriter_->write(filename.c_str());
};

double Fracture::reservoirTraction(int i) const{
   return ddm::tractionVoigt(reservoir_.stress[i], cell_normals_[i]);
}
double Fracture::fractureForce(int i) const{
  return reservoirTraction(i) - fracture_pressure_[i];
}

double Fracture::leakofDp(int i) const{
   return fracture_pressure_[i][0] - reservoir_.pressure[i]
                - (fracture_dgh_[i] - reservoir_.cellZ[i] * gravity_ * reservoir_.density[i]);
}
void
Fracture::writeIterationSnapshot(int step, int round, const std::string& tag) const
{
    if (!prm_.get<bool>("solver.write_iteration_snapshots", false) || !active_)
        return;
    const std::string dir = prm_.get<std::string>("solver.iteration_snapshot_dir", ".");
    const std::string base = dir + "/" + name();
    const bool per_cell = prm_.get<bool>("solver.iteration_snapshot_cells", true);
    static std::set<std::string> headers_written;
    const std::vector<double> K1 = stressIntensityK1();
    const std::vector<double> leak = leakOfRate();
    std::map<int, double> wimap;
    for (const auto& w : wellIndices()) wimap.insert_or_assign(w.cell, w.ctf);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    double area = 0, open_area = 0, vol = 0, qleak = 0, wisum = 0, psum = 0, wmax = 0;
    int nclosed = 0;
    std::ofstream cf;
    if (per_cell) {
        const bool fresh = !headers_written.count(base + "_iter.csv");
        cf.open(base + "_iter.csv", std::ios::app);
        if (fresh) {
            cf << "step,round,tag,cell,x,y,z,area,p_frac,p_res,traction,width,closed,K1,leakoff_rate,wi\n";
            headers_written.insert(base + "_iter.csv");
        }
    }
    for (const auto& e : Dune::elements(grid_->leafGridView())) {
        const int i = mapper.index(e);
        const auto g = e.geometry();
        const double a = g.volume();
        const auto c = g.center();
        const double w = fracture_width_[i][0];
        const double pf = fracture_pressure_[i][0];
        const int closed = (i < static_cast<int>(closed_cells_.size())) ? closed_cells_[i] : 0;
        area += a; if (w > 0) { open_area += a; vol += a * w; }
        qleak += (i < static_cast<int>(leak.size()) ? leak[i] : 0.0) * a;
        wisum += wimap.count(reservoir_map_.cells[i]) ? wimap[reservoir_map_.cells[i]] : 0.0;
        psum += pf; wmax = std::max(wmax, w); nclosed += closed;
        if (per_cell)
            cf << step << "," << round << "," << tag << "," << i << "," << c[0] << "," << c[1] << "," << c[2]
               << "," << a << "," << pf << "," << reservoir_.pressure[i] << "," << reservoirTraction(i)
               << "," << w << "," << closed << "," << (i < static_cast<int>(K1.size()) ? K1[i] : 0.0)
               << "," << (i < static_cast<int>(leak.size()) ? leak[i] : 0.0)
               << "," << (wimap.count(reservoir_map_.cells[i]) ? wimap[reservoir_map_.cells[i]] : 0.0) << "\n";
    }
    const std::string sname = base + "_iter_summary.csv";
    const bool fresh = !headers_written.count(sname);
    std::ofstream sf(sname, std::ios::app);
    if (fresh) {
        sf << "step,round,tag,perf_pressure,well_rate,perf_rate,ncells,area,open_area,volume,Q_leak,WI_sum,p_frac_mean,width_max,n_closed\n";
        headers_written.insert(sname);
    }
    const size_t n = numFractureCells();
    sf << step << "," << round << "," << tag << "," << perf_pressure_ << "," << well_rate_ << "," << well_perf_rate_
       << "," << n << "," << area << "," << open_area << "," << vol << "," << qleak << "," << wisum
       << "," << (n ? psum / n : 0.0) << "," << wmax << "," << nclosed << "\n";
}

void
Fracture::writemulti(double time) const
{
    OPM_TIMEFUNCTION();
    // NB it time == 0.0 we should
    bool first = (time == 0.0);
    std::stringstream message;
    message << "Writing fracture data to VTK files at time: " << time << "grid_size"
            << numFractureCells() << std::endl;
    OpmLog::info(message.str());
    //  need to have copies in case of async outout (and interface to functions)

    std::vector<double> K1 = this->stressIntensityK1();
    std::vector<double> fracture_pressure(numFractureCells(), 0.0);
    std::vector<double> leakof_dp(numFractureCells(), 0.0);
    std::vector<double> filtercake_thikness = filtercake_thikness_;
    std::vector<double> reservoir_pressure = reservoir_.pressure;
    std::vector<double> reservoir_dist = reservoir_map_.dist;
    std::vector<double> reservoir_perm = reservoir_.perm;
    std::vector<double> reservoir_cstress = reservoir_.cstress;
    std::vector<double> reservoir_mobility = reservoir_.mobility;
    std::vector<double> reservoir_traction(reservoir_.stress.size(), 0);
    std::vector<double> fracture_force(reservoir_.stress.size(), 0);
    std::vector<double> reservoir_cells(reservoir_map_.cells.size(), 0.0);
    std::vector<double> fracture_width(fracture_width_.size(), 0.0);
    std::vector<double> fracture_head(fracture_width_.size(), 0.0);
    std::vector<double> flow_width(fracture_width_.size(), 0.0);
    std::vector<double> rhs_width(rhs_width_.size(), 0.0);
    std::vector<double> well_index(numFractureCells(), 0.0);
    std::map<int, double> wellIndMap;
    std::vector<double> leakofrate = leakOfRate();
    {
        // make map to do it easy
        for (const auto& wind : this->wellIndices()) {
            wellIndMap.insert_or_assign(wind.cell, wind.ctf);
        }
    }
    // loop for need things only with fracture

        for (size_t i = 0; i < rhs_width_.size(); ++i) {
            rhs_width[i] = rhs_width_[i][0];
        }

        for (size_t i = 0; i < fracture_width_.size(); ++i) {
            fracture_width[i] = fracture_width_[i][0];
            flow_width[i] = std::max(fracture_width_[i][0], min_width_);
            fracture_pressure[i] = fracture_pressure_[i][0];
            reservoir_cells[i] = reservoir_map_.cells[i]; // only converts to double
            reservoir_traction[i] = reservoirTraction(i);
            fracture_force[i] = fractureForce(i);
            leakof_dp[i] = leakofDp(i);
            // make relative to origo
            fracture_head[i] = (fracture_pressure[i] - fracture_dgh_[i])
                - (injectionPressure() - (perf_ref_depth_ * gravity_ * density_perf_));
            assert(filtercake_thikness[i] >= 0.0);
        }
    if (!first) {
        for (size_t i = 0; i < numFractureCells(); ++i) {
                // maybe slow
                well_index[i] = wellIndMap[reservoir_map_.cells[i]];
        }

    }
    vtkmultiwriter_->beginWrite(time);


    if (reservoir_cells.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_cells, "ReservoirCell");
    }
    if (reservoir_perm.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_perm, "ReservoirPerm");
    }
    if (reservoir_dist.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_dist, "ReservoirDist");
    }

    if (reservoir_cstress.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_cstress, "ReservoirCStres");
    }
    if (fracture_pressure.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_pressure, "FracturePressure");
    }
    if (fracture_head.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_head, "FractureHead");
    }
    if (reservoir_pressure.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_pressure, "ReservoirPressure");
    }
    if (leakof_dp.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(leakof_dp, "LeakofDP");
    }
    if (filtercake_thikness.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(filtercake_thikness, "FilterCakeThickness");
    }
    if (reservoir_traction.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_traction, "ReservoirTraction");
    }
    if (fracture_force.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_force, "FractureForce");
    }

    if (reservoir_mobility.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_mobility, "ReservoirMobility");
    }
    if (rhs_width.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(rhs_width, "ForceOnFracture");
    }
    if (leakofrate.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(leakofrate, "LeakOfRate");
    }

    if (well_index.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(well_index, "WellIndex");
    }
    if (fracture_width.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_width, "FractureWidth");
        vtkmultiwriter_->attachScalarElementData(flow_width, "FlowWidth");
        vtkmultiwriter_->attachScalarElementData(K1, "stressIntensityK1");
    }
    vtkmultiwriter_->endWrite(false);
};

void
Fracture::updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                               const Dune::CpGrid& grid,
                               const std::vector<Fracture::EntitySeed>& entity_seeds)
{
    OPM_TIMEFUNCTION();
    reservoir_map_.cells.resize(numFractureCells());
    using GridView = typename Grid::LeafGridView;
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    int tri_divide = 0;
    int tri_outside = 0;
    reservoir_map_.allAreas.resize(numFractureCells());
    reservoir_map_.allCells.resize(numFractureCells());
    reservoir_map_.allCenters.resize(numFractureCells());
    bool full_intersections = prm_.get<bool>("solver.divide_wellidx", false);
    for (auto& element : elements(grid_->leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        auto geom = element.geometry();
        using Vec3d = external::cvf::Vec3d;
        auto vertex = geom.center();
        Vec3d point(vertex[0], vertex[1], vertex[2]);
        reservoir_map_.cells[elemIdx]  = external::cellOfPoint(cellSearchTree, grid, entity_seeds, point);
        // calculate area
        std::vector<std::array<double,3>> tri_corners;
        {
            for (size_t i = 0; i < geom.corners(); ++i) {
                auto corner = geom.corner(i);
                std::array<double,3> tmp({corner[0], corner[1], corner[2]});
                tri_corners.push_back(tmp);
            }
        }

        if(full_intersections){
            std::vector<int> rcells;
            std::vector<double> areas;
            std::vector<Dune::FieldVector<double,3>> centers;
            std::vector<std::array<double,3>> hex_corners;
            external::cellsOfTri(rcells, areas, centers, cellSearchTree, grid, entity_seeds, tri_corners);
            reservoir_map_.allAreas[elemIdx] = areas;
            reservoir_map_.allCells[elemIdx] = rcells;
            reservoir_map_.allCenters[elemIdx] = centers;
        }
    }
    auto it = std::find(reservoir_map_.cells.begin(), reservoir_map_.cells.end(), -1);
    auto extended_fractures = prm_.get<bool>("extended_fractures");
    if (!(it == reservoir_map_.cells.end()) && !(extended_fractures)) {
        OpmLog::warning("Remove fracture outside of model");
        // remove fracture outside of model
        this->removeCells();
        this->updateReservoirCells(cellSearchTree, grid, entity_seeds);
    }else{
      if(!(it == reservoir_map_.cells.end())){
          OpmLog::warning("Fracture outside of model");
      }
    }
}

void
Fracture::updateReservoirProperties()
{
    // updater for standalone test
    double perm = prm_.get<double>("reservoir.perm");
    double dist = prm_.get<double>("reservoir.dist");
    double cstress = prm_.get<double>("KMax");
    double mobility = prm_.get<double>("reservoir.mobility");
    size_t nc = numFractureCells();
    reservoir_.perm.resize(nc, perm);
    reservoir_map_.dist.resize(nc, dist);
    reservoir_.mobility.resize(nc, 1000);
    reservoir_.pressure.resize(nc, 100.0e5);
    reservoir_.stress.resize(nc);
    reservoir_.cstress.resize(nc, cstress);
    fracture_water_property_evaluator_ = {};

    for (size_t i = 0; i != nc; ++i)
        reservoir_.stress[i] = Dune::FieldVector<double, 6> {0, 0, 0, 0, 0, 0};

    nu_ = 0.25;
    E_ = 1e9;
    this->initFractureWidth();
}


// g * rho_face * (z_i - z_j) for each entry of htrans_; independent of the depth datum
std::vector<double>
Fracture::faceGravityHeads() const
{
    const auto& gv = grid_->leafGridView();
    std::vector<double> z(gv.size(0), 0.0);
    for (const auto& element : Dune::elements(gv)) {
        z[gv.indexSet().index(element)] = element.geometry().center()[2];
    }
    std::vector<double> dgh(htrans_.size(), 0.0);
    for (size_t k = 0; k < htrans_.size(); ++k) {
        const size_t i = std::get<0>(htrans_[k]);
        const size_t j = std::get<1>(htrans_[k]);
        const double rho = 0.5 * (reservoir_.density[i] + reservoir_.density[j]);
        dgh[k] = gravity_ * rho * (z[i] - z[j]);
    }
    return dgh;
}

void
Fracture::addSource()
{
    if (rhs_pressure_.size() == 0) {
        size_t nc = numFractureCells();
        rhs_pressure_.resize(nc + numWellEquations());
    }
    assert(rhs_pressure_.size() == numFractureCells() + numWellEquations());
    rhs_pressure_ = 0;
    // add contributions for fracture fracture gravity contributions
    int nc = grid_->leafGridView().size(0);
    fracture_dgh_.resize(nc, 0.0);
    for (auto element : Dune::elements(grid_->leafGridView())) {
        int i = grid_->leafGridView().indexSet().index(element);
        double z = element.geometry().center()[2];
        fracture_dgh_[i] = gravity_ * reservoir_.density[i] * z;
    }
    // could be put into the assemble loop
    const std::vector<double> face_dgh = faceGravityHeads();
    for (size_t k = 0; k < htrans_.size(); ++k) {
        const auto& matel = htrans_[k];
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double t1 = std::get<2>(matel);
        double t2 = std::get<3>(matel);
        const double h1 = std::max(fracture_width_[i][0],min_width_);
        const double h2 = std::max(fracture_width_[j][0],min_width_);
        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);

        const double mobility = 0.5 * (reservoir_.mobility[i] + reservoir_.mobility[j]);//NB shoul chane to mobility_water_perf_;
        value = 1 / value;
        value *= mobility;
        rhs_pressure_[i] += value * face_dgh[k];
        rhs_pressure_[j] -= value * face_dgh[k];
    }

    for (size_t i = 0; i < reservoir_.pressure.size(); ++i) {
        rhs_pressure_[i] += leakof_[i] * reservoir_.pressure[i];// gravity i handeld in next loop
    }
    // gravity contribution from fracture to reservoir
    for (auto element : Dune::elements(grid_->leafGridView())) {
        int i = grid_->leafGridView().indexSet().index(element);
        double z = element.geometry().center()[2];
        rhs_pressure_[i] += leakof_[i] * (z - reservoir_.cellZ[i]) * gravity_ * reservoir_.density[i];
    }
    // gravity contributions between fracture cells

    auto control = prm_.get_child("control");
    std::string control_type = effectiveControlType();
    if (control_type == "rate") {
        double scale = well_source_.size();
        double rate = control.get<double>("rate");
        for (auto cell : well_source_) {
            rhs_pressure_[cell] += rate / scale;
        }
    } else if (control_type == "pressure") {
        double pressure = control.get<double>("pressure");
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            double dh_perf = origo_[2] * gravity_ * reservoir_.density[cell];
            double dh_cell = fracture_dgh_[cell];
            ;
            double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (pressure - dh);
        }
    } else if (control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            double dh_perf = perf_ref_depth_ * gravity_ * reservoir_.density[cell];
            double dh_cell = fracture_dgh_[cell];
            ;
            double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (perf_pressure_ - dh);
        }
    } else if (control_type == "rate_well") {
        // Single-DOF well row: a static single-phase stand-in for the well model
        // (matrix WI/p_res frozen from the last flow solve, unit mobility).
        // Limitations: one fracture per well, lambda=1, no THP.
        assert(numWellEquations() == 1);
        const double well_rate = well_rate_;
        const double WI = total_wellindex_;
        const double density = reservoir_.density[std::get<0>(perfinj_[0])];
        const double lambda = 1.0;
        double wi_dzfac = wi_respress_;
        wi_dzfac -= wi_dz_*density*gravity_; // sum wi dz
        wi_dzfac += (WI* (perf_ref_depth_- well_ref_depth_))*gravity_*density; // convert from bhp ref depth to perf ref depth
        rhs_pressure_[rhs_pressure_.size() - 1] = well_rate + lambda * wi_dzfac;
    } else if (control_type == "bhp_well") {
        // constraint row: scale*(p_w) = scale*(target at the perf-ref datum);
        // scale MUST match bhpWellRowScale() in the assembler
        assert(numWellEquations() == 1);
        const double scale = std::max(total_wellindex_, 1e-14);
        const double density = reservoir_.density[std::get<0>(perfinj_[0])];
        const double p_target = well_target_bhp_
            + density * gravity_ * (perf_ref_depth_ - well_ref_depth_);
        rhs_pressure_[rhs_pressure_.size() - 1] = scale * p_target;
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
    if (control_type == "rate_well" || control_type == "bhp_well") {
        // the well DOF is at the perf-ref datum; fed cells see its hydrostatic column
        const auto input = makePressureAssemblyInput();
        double well_dh = 0.0;
        for (const auto& [cell, wi] : perfinj_) {
            const double v = wi * wellConnectionMobilityValue(input, cell);
            const double dh = fracture_dgh_[cell] - gravity_ * reservoir_.density[cell] * perf_ref_depth_;
            rhs_pressure_[cell] += v * dh;
            well_dh += v * dh;
        }
        if (control_type == "rate_well") {
            rhs_pressure_[rhs_pressure_.size() - 1] -= well_dh;
        }
    }
}

double
Fracture::injectionBhp() const
{

  double dz = perf_ref_depth_ - well_ref_depth_;;
  double density = density_perf_;
  return injectionPressure() - dz*gravity_*density;
}

double
Fracture::injectionPressure() const
{
    std::string control_type = effectiveControlType();
    if (control_type == "rate") {
        double bhp = 0.0;
        double scale = well_source_.size();
        // could have corrected for WI
        for (auto cell : well_source_) {
            bhp += fracture_pressure_[cell] / scale;
        }
        return bhp;
    } else if (control_type == "pressure") {
        double bhp = prm_.get<double>("control.pressure");
        return bhp;
    } else if (control_type == "perf_pressure") {
        double bhp = perf_pressure_;
        return bhp;
    } else if (control_type == "rate_well" || control_type == "bhp_well") {
        // the solved well DOF
        return fracture_pressure_[fracture_pressure_.size() - 1][0];
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
    return 0.0;
}

std::vector<double>
Fracture::leakOfRate() const
{
    if (leakof_.size() == 0) {
        std::vector<double> rate;
        return rate;
    }
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    std::vector<double> leakofrate(numFractureCells(), 0);
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        double dh_res = reservoir_.cellZ[eIdx] * gravity_ * reservoir_.density[eIdx];
        double dh_frac = fracture_dgh_[eIdx];
        double dp = ((fracture_pressure_[eIdx]-dh_frac) - (reservoir_.pressure[eIdx]-dh_res));
        double q = leakof_[eIdx] * dp;
        leakofrate[eIdx] = q / geom.volume();
    }
    return leakofrate;
}

std::vector<RuntimePerforation>
Fracture::wellIndices() const{
   std::vector<RuntimePerforation> wellindices;
   const std::string wi_accel = prm_.get<std::string>("solver.wi_vector_acceleration", "none");
   if (wi_accel != "none" && !well_indices_accel_.empty()) {
         wellindices = well_indices_accel_;
   } else if(well_indices_.size() == 2){
         wellindices = wellIndicesAvrg(well_indices_);
   }
   return wellindices;
}

double Fracture::applyCouplingUpdate(double current,
                                     double target,
                                     double damping_factor,
                                     const std::string& mode,
                                     const std::string& channel,
                                     CouplingMixState& state)
{
    auto clamp = [](const double value, const double low, const double high)
    {
        return std::max(low, std::min(value, high));
    };

    const double relax_min = prm_.get<double>("solver.coupling_relax_min", 0.05);
    const double relax_max = prm_.get<double>("solver.coupling_relax_max", 1.25);
    const double omega_init_default = 1.0 / (1.0 + std::max(0.0, damping_factor));
    const double omega_init = prm_.get<double>("solver.coupling_relax_init", omega_init_default);
    const double acceptance_factor = prm_.get<double>("solver.coupling_acceptance_factor", 1.0);
    const bool fallback_legacy = prm_.get<bool>("solver.coupling_fallback_legacy", true);
    const double anderson_alpha_min = prm_.get<double>("solver.coupling_anderson_alpha_min", -1.0);
    const double anderson_alpha_max = prm_.get<double>("solver.coupling_anderson_alpha_max", 2.0);
    const double eps = 1e-14;

    const auto legacyUpdate = [&](const double cur, const double tgt)
    {
        return (cur * damping_factor + tgt) / (1.0 + damping_factor);
    };

    if (!state.initialized) {
        state.initialized = true;
        state.previous_target = target;
        state.previous_residual = 0.0;
        state.has_previous_residual = false;
        state.omega = clamp(omega_init, relax_min, relax_max);
        return target;
    }

    if (mode == "legacy" || mode == "current") {
        const double mixed = legacyUpdate(current, target);
        state.previous_target = target;
        state.previous_residual = target - mixed;
        state.has_previous_residual = true;
        return mixed;
    }

    const double residual = target - current;
    const double residual_abs = std::abs(residual);

    if (mode == "anderson") {
        bool accepted = false;
        if (state.has_previous_residual) {
            const double denominator = residual - state.previous_residual;
            if (std::abs(denominator) > eps) {
                const double alpha = clamp(-residual / denominator, anderson_alpha_min, anderson_alpha_max);
                const double candidate = alpha * state.previous_target + (1.0 - alpha) * target;
                const double candidate_residual_abs = std::abs(target - candidate);
                if (candidate_residual_abs <= acceptance_factor * residual_abs) {
                    current = candidate;
                    accepted = true;
                }
            }
        }
        if (!accepted) {
            const double aitken_like = clamp(state.omega, relax_min, relax_max);
            const double candidate = current + aitken_like * residual;
            const double candidate_residual_abs = std::abs(target - candidate);
            if (candidate_residual_abs <= acceptance_factor * residual_abs || !fallback_legacy) {
                current = candidate;
            }
            else {
                current = legacyUpdate(current, target);
            }
        }
    }
    else {
        double omega = clamp(state.omega, relax_min, relax_max);
        if (state.has_previous_residual) {
            const double denominator = residual - state.previous_residual;
            if (std::abs(denominator) > eps) {
                omega = clamp(-omega * state.previous_residual / denominator, relax_min, relax_max);
            }
        }
        const double candidate = current + omega * residual;
        const double candidate_residual_abs = std::abs(target - candidate);
        if (candidate_residual_abs <= acceptance_factor * residual_abs || !fallback_legacy) {
            current = candidate;
            state.omega = omega;
        }
        else {
            current = legacyUpdate(current, target);
        }
    }

    if (prm_.get<int>("solver.verbosity", 0) > 2) {
        std::stringstream os;
        os << "Coupling update " << channel
           << ": mode=" << mode
           << ", current=" << current
           << ", target=" << target
           << ", residual=" << (target - current)
           << ", omega=" << state.omega;
        OpmLog::info(os.str());
    }

    state.previous_target = target;
    state.previous_residual = target - current;
    state.has_previous_residual = true;
    return current;
}

std::vector<RuntimePerforation>
Fracture::applyVectorCouplingUpdate(const std::vector<RuntimePerforation>& target)
{
    // Vector Aitken / depth-1 Anderson acceleration of the per-cell CTF fixed
    // point g(x) (opt-in via solver.wi_vector_acceleration). The legacy path
    // (wellIndicesAvrg) is the damped Picard x_k = (g_k + d*x_{k-1})/(1+d);
    // here the relaxation/mixing coefficient is computed from the residual
    // history of the whole CTF vector instead of being fixed. Safeguards: CTFs
    // are kept non-negative, and an accelerated step is only accepted when it
    // does not increase the residual (coupling_acceptance_factor), otherwise
    // the legacy damped update is used.
    const std::string mode = prm_.get<std::string>("solver.wi_vector_acceleration", "none");
    const double relax_min = prm_.get<double>("solver.coupling_relax_min", 0.05);
    const double relax_max = prm_.get<double>("solver.coupling_relax_max", 1.25);
    const double damping_factor = prm_.get<double>("solver.damping_factor_wi", 2.0);
    const double alpha_min = prm_.get<double>("solver.coupling_anderson_alpha_min", -1.0);
    const double alpha_max = prm_.get<double>("solver.coupling_anderson_alpha_max", 2.0);
    const double acceptance = prm_.get<double>("solver.coupling_acceptance_factor", 1.0);
    const bool fallback_legacy = prm_.get<bool>("solver.coupling_fallback_legacy", true);
    const int verbosity = prm_.get<int>("solver.verbosity", 0);

    auto clamp = [](const double v, const double lo, const double hi)
    { return std::max(lo, std::min(v, hi)); };

    std::map<int, double> tgt;
    for (const auto& p : target) {
        tgt[p.cell] = p.ctf;
    }

    auto& st = wi_vec_mix_;
    std::vector<RuntimePerforation> out = target;

    if (!st.initialized) {
        // first solve: take the target directly (matches the scalar channels)
        st.initialized = true;
        st.omega = clamp(1.0 / (1.0 + std::max(0.0, damping_factor)), relax_min, relax_max);
        st.prev_mixed = tgt;
        st.prev_target = tgt;
        st.prev_residual.clear();
        st.has_prev_residual = false;
        return out;
    }

    // residual r_k = g_k - x_{k-1} on cells with history; brand-new cells (the
    // fracture grew) have no history and take the target directly below.
    std::map<int, double> resid;
    double r_norm2 = 0.0, dr_norm2 = 0.0, r_dot_dr = 0.0, rprev_dot_dr = 0.0;
    for (const auto& [cell, g] : tgt) {
        const auto it = st.prev_mixed.find(cell);
        if (it == st.prev_mixed.end()) {
            continue;
        }
        const double r = g - it->second;
        resid[cell] = r;
        r_norm2 += r * r;
        if (st.has_prev_residual) {
            const auto rp = st.prev_residual.find(cell);
            if (rp != st.prev_residual.end()) {
                const double dr = r - rp->second;
                dr_norm2 += dr * dr;
                r_dot_dr += r * dr;
                rprev_dot_dr += rp->second * dr;
            }
        }
    }

    std::map<int, double> mixed;
    bool accelerated = false;
    std::string detail;
    if (st.has_prev_residual && dr_norm2 > 1e-28) {
        if (mode == "anderson") {
            const double alpha = clamp(-r_dot_dr / dr_norm2, alpha_min, alpha_max);
            for (const auto& [cell, g] : tgt) {
                const auto gp = st.prev_target.find(cell);
                const double gprev = (gp != st.prev_target.end()) ? gp->second : g;
                mixed[cell] = std::max(0.0, alpha * gprev + (1.0 - alpha) * g);
            }
            accelerated = true;
            detail = "alpha=" + std::to_string(alpha);
        } else { // aitken
            const double omega = clamp(-st.omega * rprev_dot_dr / dr_norm2, relax_min, relax_max);
            for (const auto& [cell, g] : tgt) {
                const auto xp = st.prev_mixed.find(cell);
                mixed[cell] = (xp != st.prev_mixed.end())
                    ? std::max(0.0, xp->second + omega * resid[cell])
                    : g;
            }
            st.omega = omega;
            accelerated = true;
            detail = "omega=" + std::to_string(omega);
        }
        if (acceptance > 0.0 && fallback_legacy) {
            double cand2 = 0.0;
            for (const auto& [cell, r] : resid) {
                const double c = tgt[cell] - mixed[cell];
                cand2 += c * c;
            }
            if (cand2 > acceptance * acceptance * r_norm2) {
                accelerated = false;
                detail += " rejected";
            }
        }
    }
    if (!accelerated) {
        // legacy damped Picard (same weights as wellIndicesAvrg)
        for (const auto& [cell, g] : tgt) {
            const auto xp = st.prev_mixed.find(cell);
            const double xprev = (xp != st.prev_mixed.end()) ? xp->second : g;
            mixed[cell] = (g + damping_factor * xprev) / (1.0 + damping_factor);
        }
    }

    for (auto& p : out) {
        const auto mx = mixed.find(p.cell);
        if (mx != mixed.end()) {
            p.ctf = mx->second;
        }
    }

    if (verbosity > 1) {
        std::stringstream os;
        os << "CTF vector coupling update: mode=" << mode
           << ", cells=" << tgt.size()
           << ", |r|=" << std::sqrt(r_norm2)
           << ", accelerated=" << (accelerated ? "true" : "false")
           << (detail.empty() ? "" : (", " + detail));
        OpmLog::info(os.str());
    }

    st.prev_residual = std::move(resid);
    st.has_prev_residual = true;
    st.prev_target = std::move(tgt);
    st.prev_mixed = std::move(mixed);
    return out;
}

void Fracture::setPerfProps(double perfpressure, double depth, double perfrate){
    const double damping_factor_perf = prm_.get<double>("solver.damping_factor_perf",2.0);
    const std::string coupling_mode = prm_.get<std::string>("solver.coupling_update_mode", "legacy");
    // When true, apply damping even on the first call (pre-0e6bb48 behaviour).
    const bool legacy_first_damp = prm_.get<bool>("solver.coupling_legacy_first_damp", false);
    if(well_indices_.size() >0){
        perf_pressure_ = perfpressure;
        perf_pressure_mix_.initialized = false;
        perf_pressure_mix_.has_previous_residual = false;
    }else{
        if (legacy_first_damp && !perf_pressure_mix_.initialized) {
            // Emulate the pre-0e6bb48 behaviour: always damp, even on the first call.
            perf_pressure_ = (perf_pressure_ * damping_factor_perf + perfpressure)
                             / (1.0 + damping_factor_perf);
            perf_pressure_mix_.initialized = true;
        } else {
            perf_pressure_ = applyCouplingUpdate(perf_pressure_,
                                                 perfpressure,
                                                 damping_factor_perf,
                                                 coupling_mode,
                                                 "perf_pressure",
                                                 perf_pressure_mix_);
        }
    }
    well_perf_rate_ = perfrate;
    perf_ref_depth_ = depth;
}

std::string Fracture::effectiveControlType() const
{
    // control.type == "well": the fracture BC is the well's currently ACTIVE
    // constraint - rate-like (RATE/RESV/GRUP) -> rate_well, BHP -> bhp_well
    // (THP throws in the sampler). Makes one config safe across well types and
    // mid-run control switches.
    const std::string t = prm_.get_child("control").get<std::string>("type");
    if (t == "well") {
        // generalized well row: the well DOF is kept for BOTH constraints, only
        // the well row's content differs (rate balance vs p_w = BHP target)
        if (!effective_control_latched_.empty())
            return effective_control_latched_;
        return well_control_is_rate_ ? "rate_well" : "bhp_well";
    }
    return t;
}

void Fracture::setWellProps(double wellrate,
                            double totalwi,
                            double wi_dz,
                            double wi_respress,
                            double ref_depth){
    const bool enable_wi_update = prm_.get<bool>("solver.enable_wi_coupling_update", false);
    const std::string coupling_mode = prm_.get<std::string>("solver.coupling_update_mode", "legacy");
    const double damping_factor_wi = prm_.get<double>("solver.damping_factor_wi",2.0);
    if (enable_wi_update) {
        well_rate_ = applyCouplingUpdate(well_rate_,
                                         wellrate,
                                         damping_factor_wi,
                                         coupling_mode,
                                         "well_rate",
                                         well_rate_mix_);
        total_wellindex_ = applyCouplingUpdate(total_wellindex_,
                                               totalwi,
                                               damping_factor_wi,
                                               coupling_mode,
                                               "total_wellindex",
                                               total_wellindex_mix_);
    }
    else {
        well_rate_ = wellrate;
        total_wellindex_ = totalwi;
        well_rate_mix_.initialized = false;
        total_wellindex_mix_.initialized = false;
        well_rate_mix_.has_previous_residual = false;
        total_wellindex_mix_.has_previous_residual = false;
    }
        wi_dz_ = wi_dz;
        wi_respress_ = wi_respress;
        well_ref_depth_ = ref_depth;
}

std::vector<RuntimePerforation>
Fracture::wellIndicesAvrg(const std::vector<std::vector<RuntimePerforation>>& well_indices) const{
    // average int time between to well indices.
    OPM_TIMEFUNCTION();
  // find union of cell indices
  assert(well_indices.size() == 2);
  std::set<int> cells;
  for(const auto& winds : well_indices){
    for (const auto& wind : winds){
      cells.insert(wind.cell);
    }
  }
  // map
  int count = 0;
  std::map<int,int> cell_pind;
  for(const auto& cell : cells) {
    cell_pind[cell] = count;
    count += 1;
  }

  std::vector<RuntimePerforation> wellindices(cells.size());
  std::vector<double> weight_ctf(well_indices.size(),0.0);
  double damping_factor_wi = prm_.get<double>("solver.damping_factor_wi",2.0);
  weight_ctf[0] = 1.0;
  weight_ctf[1] = damping_factor_wi;
  double total_weight = 0.0;
  for(const auto& winds : well_indices){
    for (const auto& wind : winds){
      int ind = cell_pind[wind.cell];
      wellindices[ind] = wind;
      wellindices[ind].ctf = 0.0;
    }
  }
  // the fracture pressure offset must match the newest fracture state
  for (const auto& wind : well_indices[0]) {
    wellindices[cell_pind[wind.cell]].fracture_pressure_offset = wind.fracture_pressure_offset;
  }
  for(size_t tind=0; tind < well_indices.size(); ++tind){
    for (const auto& wind : well_indices[tind]){
      int ind = cell_pind[wind.cell];
      wellindices[ind].ctf += weight_ctf[tind]*wind.ctf;
    }
    total_weight += weight_ctf[tind];
  }

  for(size_t i=0; i < wellindices.size(); ++i){
        wellindices[i].ctf /= total_weight;
  }
  return wellindices;
}

Fracture::WiOptions
Fracture::WiOptions::read(const PropertyTree& prm)
{
    WiOptions opt;
    const std::string wi_upscaling = prm.get<std::string>("solver.wi_upscaling", "legacy");
    if (wi_upscaling != "legacy" && wi_upscaling != "conductivity") {
        OPM_THROW(std::runtime_error, "Unknown solver.wi_upscaling: " + wi_upscaling);
    }
    // opt-in: Flow drives the fracture CTF with the fracture's own pressure
    const bool fracture_pressure = prm.get<bool>("solver.wi_fracture_pressure", false);
    if (fracture_pressure && wi_upscaling == "conductivity") {
        OPM_THROW(std::runtime_error,
                  "solver.wi_fracture_pressure and solver.wi_upscaling=conductivity exclude each other");
    }
    opt.mode = fracture_pressure ? WiMode::FracturePressure
        : (wi_upscaling == "conductivity" ? WiMode::Conductivity : WiMode::Legacy);
    opt.divideWellIdx = prm.get<bool>("solver.divide_wellidx", false);
    opt.fluxNorm = prm.get<bool>("solver.wi_flux_normalization", true);
    // Must stay well below near-closure leak-off differences (~1e2-1e3 Pa), or alpha is underestimated.
    opt.dpFloor = prm.get<double>("solver.wi_pressure_floor", 1.0e2);
    opt.alphaMax = prm.get<double>("solver.wi_normalization_max", 2.0);
    opt.signGate = prm.get<bool>("solver.wi_sign_gate", false);
    return opt;
}

std::vector<Fracture::CellExchange>
Fracture::gatherCellExchange_(const WiOptions& opt) const
{
    const bool only_to_one = !opt.divideWellIdx;
    std::vector<int> res_cells = reservoir_map_.cells;
    std::vector<double> tot_areas(reservoir_map_.allCells.size(), 0.0);
    if (!only_to_one) {
        for (size_t ind = 0; ind < reservoir_map_.allCells.size(); ++ind) {
            double total_area = 0.0;
            const auto& cells = reservoir_map_.allCells[ind];
            for (size_t lind = 0; lind < cells.size(); ++lind) {
                res_cells.push_back(cells[lind]);
                total_area += reservoir_map_.allAreas[ind][lind];
            }
            tot_areas[ind] = total_area;
        }
    }
    std::sort(res_cells.begin(), res_cells.end());
    res_cells.erase(std::unique(res_cells.begin(), res_cells.end()), res_cells.end());

    std::vector<CellExchange> ex(res_cells.size());
    for (size_t i = 0; i < res_cells.size(); ++i) {
        ex[i].cell = res_cells[i];
    }
    const bool fracture_pressure_drive = (opt.mode == WiMode::FracturePressure);
    const std::vector<double> leakofrate = this->leakOfRate();
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        const double q = leakofrate[eIdx] * geom.volume();

        std::vector<int> local_res_cells;
        if (only_to_one) {
            local_res_cells.push_back(reservoir_map_.cells[eIdx]);
        } else {
            local_res_cells = reservoir_map_.allCells[eIdx];
        }
        for (int lind = 0; lind < int(local_res_cells.size()); ++lind) {
            const int res_cell = local_res_cells[lind];
            double area_frac = 1.0;
            if (!only_to_one) {
                area_frac = reservoir_map_.allAreas[eIdx][lind] / tot_areas[eIdx];
            }
            const auto it = std::lower_bound(res_cells.begin(), res_cells.end(), res_cell);
            assert(it != res_cells.end() && *it == res_cell);
            auto& c = ex[it - res_cells.begin()];
            const double loc_area = area_frac * element.geometry().volume();
            c.area += loc_area;
            c.q += area_frac * q;
            // leakof_ = mob*perm*area/dist >= 0, so this conductivity is never negative
            if (eIdx < static_cast<int>(leakof_.size()) && reservoir_.mobility[eIdx] > 0.0) {
                c.trans += area_frac * leakof_[eIdx] / reservoir_.mobility[eIdx];
            }
            if (reservoir_map_.cells[eIdx] == res_cell) {
                c.p = reservoir_.pressure[eIdx];
                c.mob = reservoir_.mobility[eIdx];
                c.dens = reservoir_.density[eIdx];
                c.z = reservoir_.cellZ[eIdx];
            } else {
                // no fracture cell has its centroid in this reservoir cell
                const auto search = off_centroid_cells_.find(res_cell);
                assert(search != off_centroid_cells_.end());
                const ReservoirCellValues v = (search != off_centroid_cells_.end())
                    ? search->second : ReservoirCellValues{};
                c.p = v.pressure;
                c.mob = v.mobility;
                c.dens = v.density;
                c.z = v.cellZ;
            }
            if (fracture_pressure_drive) {
                // fracture pressure moved hydrostatically to the connection depth
                const double p_at_conn = fracture_pressure_[eIdx][0]
                    + gravity_ * reservoir_.density[eIdx] * (c.z - geom.center()[2]);
                const double w_leak = (eIdx < static_cast<int>(leakof_.size()))
                    ? area_frac * leakof_[eIdx] : 0.0;
                c.pfracLeak += w_leak * p_at_conn;
                c.leak += w_leak;
                c.pfracArea += loc_area * p_at_conn;
            }
        }
    }
    return ex;
}

RuntimePerforation
Fracture::perforation_(const CellExchange& c, const double inj_press) const
{
    RuntimePerforation perf;
    perf.cell = c.cell;
    perf.ctf = 0.0;
    // the CTF is derived for a connection pressure at the reservoir cell depth
    perf.depth = c.z;
    perf.segment = this->wellinfo_.segment;
    perf.perf_range = this->wellinfo_.perf_range;
    perf.pressure = inj_press;
    return perf;
}

// CTF back-calculated from the leak-off rate at the current pressure difference.
std::vector<RuntimePerforation>
Fracture::wiLegacy_(const std::vector<CellExchange>& cells, const double inj_press) const
{
    std::vector<RuntimePerforation> perfs;
    perfs.reserve(cells.size());
    for (const auto& c : cells) {
        auto& perf = perfs.emplace_back(perforation_(c, inj_press));
        const double dh_res = c.z * gravity_ * c.dens;
        const double dh_perf = gravity_ * c.dens * origo_[2];
        double WI = 0.0;
        double ctf = 0.0;
        if (c.mob <= 0.0) {
            OpmLog::debug(fmt::format("Zero mobility for perf cell: {}", c.cell));
        } else {
            WI = c.q / ((inj_press - dh_perf) - (c.p - dh_res));
            ctf = WI / c.mob;
        }
        if (WI < 0.0) {
            OpmLog::debug(fmt::format("Negative WI: {} for cell: {}", WI, c.cell));
            ctf = 0.0;
        }
        perf.ctf = ctf;
    }
    return perfs;
}

// Aggregated exchange conductivity, scaled by one global factor alpha that
// preserves the total exchange rate at the current operating point.
std::vector<RuntimePerforation>
Fracture::wiConductivity_(const std::vector<CellExchange>& cells,
                          const double inj_press,
                          const WiOptions& opt) const
{
    std::vector<RuntimePerforation> perfs;
    perfs.reserve(cells.size());
    double sum_q = 0.0;
    double sum_ctf_dp = 0.0;
    int legacy_negative_count = 0;
    for (const auto& c : cells) {
        auto& perf = perfs.emplace_back(perforation_(c, inj_press));
        const double dh_res = c.z * gravity_ * c.dens;
        const double dh_perf = gravity_ * c.dens * origo_[2];
        double ctf = c.trans;
        const double dp_well = (inj_press - dh_perf) - (c.p - dh_res);
        const bool not_fed = (c.q < 0.0 || (dp_well <= 0.0 && c.q > 0.0));
        if (not_fed) {
            ++legacy_negative_count; // legacy would have zeroed this perf
        }
        // A cell the fracture does not feed must get no CTF: the flow applies
        // it against its own well-cell dp and would drive the cell backwards.
        if (opt.signGate && not_fed) {
            ctf = 0.0;
        }
        sum_q += c.q;
        sum_ctf_dp += ctf * std::max(c.mob, 0.0) * std::max(dp_well, opt.dpFloor);
        perf.ctf = ctf;
    }
    // alpha stays 1 unless both sums are meaningfully positive
    double alpha = 1.0;
    if (opt.fluxNorm && sum_q > 0.0 && sum_ctf_dp > 0.0) {
        alpha = std::min(sum_q / sum_ctf_dp, opt.alphaMax);
    }
    if (prm_.get<int>("solver.verbosity", 0) > 0) {
        std::stringstream os;
        os << "WI upscaling (conductivity): alpha=" << alpha
           << ", sum_q=" << sum_q
           << ", legacy_negative_perfs=" << legacy_negative_count
           << "/" << cells.size();
        OpmLog::info(os.str());
    }
    if (alpha != 1.0) {
        for (auto& perf : perfs) {
            perf.ctf *= alpha;
        }
    }
    return perfs;
}

// Leak-off weighted conductivity; Flow drives it with the fracture pressure
// passed as an offset to the anchor pressure.
std::vector<RuntimePerforation>
Fracture::wiFracturePressure_(const std::vector<CellExchange>& cells, const double inj_press) const
{
    // Flow adds the offset to BHP, or to the seed segment pressure for multi-segment wells
    const double anchor_pressure = this->wellinfo_.perf_range.has_value()
        ? inj_press : injectionBhp();
    std::vector<RuntimePerforation> perfs;
    perfs.reserve(cells.size());
    for (const auto& c : cells) {
        auto& perf = perfs.emplace_back(perforation_(c, inj_press));
        perf.ctf = c.trans;
        const double p_frac = (c.leak > 0.0) ? c.pfracLeak / c.leak
            : (c.area > 0.0 ? c.pfracArea / c.area : inj_press);
        perf.fracture_pressure_offset = p_frac - anchor_pressure;
    }
    return perfs;
}

std::vector<RuntimePerforation>
Fracture::wellIndices_() const
{
    OPM_TIMEFUNCTION();
    if (leakof_.size() == 0) {
        // pressure not assembled yet
        return {};
    }
    const auto opt = WiOptions::read(prm_);
    auto cells = gatherCellExchange_(opt);
    const auto outside = std::remove_if(cells.begin(), cells.end(),
                                        [](const auto& c) { return c.cell < 0; });
    if (outside != cells.end()) {
        OpmLog::warning("Cells of fracture outside removed from perforation index");
        cells.erase(outside, cells.end());
    }
    const double inj_press = injectionPressure();
    switch (opt.mode) {
    case WiMode::Conductivity:
        return wiConductivity_(cells, inj_press, opt);
    case WiMode::FracturePressure:
        return wiFracturePressure_(cells, inj_press);
    case WiMode::Legacy:
        break;
    }
    return wiLegacy_(cells, inj_press);
}

template <typename Scalar>
void
Fracture::assignGeoMechWellState(PerfData<Scalar>& perfData) const
{
    auto perfPos = std::find(
    perfData.cell_index.begin(), perfData.cell_index.end(), this->wellInfo().well_cell);
    const auto perfIx = std::distance(perfData.cell_index.begin(), perfPos);
    ConnFracStatistics<Scalar>& stats = perfData.connFracStatistics[perfIx];

    using Quantity = typename ConnFracStatistics<Scalar>::Quantity;

    constexpr auto pressIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Pressure);
    constexpr auto rateIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::FlowRate);
    constexpr auto widthIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Width);

    const auto nCells = this->reservoir_map_.cells.size();

    stats.reset();

    for (auto cellIx = 0 * nCells; cellIx < nCells; ++cellIx) {
        auto samplePoint = typename ConnFracStatistics<Scalar>::SamplePoint {};

        samplePoint[pressIx] = this->fracture_pressure_[cellIx][0];
        double dh_res = reservoir_.cellZ[cellIx] * gravity_ * reservoir_.density[cellIx];
        double dh_frac = fracture_dgh_[cellIx];
        samplePoint[rateIx] = this->leakof_[cellIx]
          * ((this->fracture_pressure_[cellIx][0]-dh_frac)- (this->reservoir_.pressure[cellIx]-dh_res));

        samplePoint[widthIx] = this->fracture_width_[cellIx][0];

        stats.addSamplePoint(samplePoint);
    }
    //
    ConnFractureData<Scalar>& frac_prop = perfData.fracture_data;
    FractureProperties fprop = calculateFractureProperties();
    // set properties of fracture to structure used in output
    frac_prop.height[perfIx] = fprop.height;
    frac_prop.length[perfIx] = fprop.length;
    frac_prop.area[perfIx] = fprop.area;
    frac_prop.flux[perfIx] = fprop.flux;
    frac_prop.WI[perfIx] = fprop.WI;
    frac_prop.volume[perfIx] = fprop.volume;
    frac_prop.filter_volume[perfIx] = fprop.filter_volume;
    frac_prop.avg_width[perfIx] = fprop.avg_width;
    frac_prop.avg_filter_width[perfIx] = fprop.avg_filter_width;
    frac_prop.inj_pressure[perfIx] = fprop.inj_pressure;
    frac_prop.inj_bhp[perfIx] = fprop.inj_bhp;
    frac_prop.inj_wellrate[perfIx] = fprop.inj_wellrate;
}

bool
Fracture::expantionMax(const FractureProperties& fprop)
{
    FractureProperties fprop_current = calculateFractureProperties();
    const double area_change_fac = prm_.get<double>("solver.area_change_fac");
    const double dt_lim = prm_.get<double>("solver.dt_limit")*86400.0;
    if (fprop_current.area > area_change_fac * fprop.area) {
        OpmLog::debug("Fracture area doubled, stopping expansion");
        max_flow_time_step_  = dt_lim;
        return true;
    } else {
      max_flow_time_step_  = std::numeric_limits<double>::max();// no limit
        return false;
    }
    return false;
}

Fracture::PropagationOptions
Fracture::PropagationOptions::read(const PropertyTree& prm)
{
    PropagationOptions o;
    o.maxIter = prm.get<int>("solver.max_iter");
    o.tol = prm.get<double>("solver.tolerance");
    o.nlinVerbosity = prm.get<double>("solver.verbosity");
    o.maxNumCoarsening = prm.get<int>("solver.max_num_coarsening");
    o.numcellThreshold = prm.get<int>("solver.numcell_threshold");
    o.smoothBoundary = prm.get<bool>("solver.smooth_boundary", false);
    o.remapSolution = prm.get<bool>("solver.remap_solution");
    const std::string cake_remap = prm.get<std::string>("solver.filtercake_remap", "legacy");
    if (cake_remap != "legacy" && cake_remap != "conservative") {
        OPM_THROW(std::runtime_error, "Unknown solver.filtercake_remap: " + cake_remap);
    }
    o.conservativeCake = (cake_remap == "conservative");
    o.forceLimit = prm.get<double>("solver.force_limit", 0.0);
    o.conservative = prm.get<bool>("solver.conservative_propagation", false);
    o.convergedCellsOnly = prm.get<bool>("solver.propagate_converged_cells_only", false);
    o.fbTol = prm.get<double>("solver.propagate_fb_tol", 1e-3);
    o.stableIters = prm.get<int>("solver.propagate_stable_iters", 2);
    o.verbosity = prm.get<int>("verbosity", 0);
    o.disable = prm.get<bool>("solver.disable_propagation", false);
    o.openFrontOnly = prm.get<bool>("solver.propagate_open_front_only", false);
    return o;
}

Fracture::VolumePacing
Fracture::VolumePacing::start(const Fracture& f, const double dt)
{
    VolumePacing v;
    v.on = f.prm_.get<bool>("solver.volume_paced_propagation", false);
    v.fac = f.prm_.get<double>("solver.volume_pacing_factor", 1.0);
    v.pressure = f.prm_.get<bool>("solver.volume_pacing_pressure", false);
    v.verbosity = f.prm_.get<int>("solver.volume_pacing_verbosity", 0);
    v.dt = dt;
    v.v0 = v.on ? std::max(f.volume_step_start_, 0.0) : 0.0;
    v.perfP0 = f.perf_pressure_;
    return v;
}

void
Fracture::VolumePacing::afterSolve(Fracture& f)
{
    if (!on) {
        return;
    }
    // eq. (1)-(2): V(k)-V0 vs (q_perf - Q_leak(k))*dt*f_vol
    const FractureProperties fp = f.calculateFractureProperties();
    // The perforation's CTF share starves an establishing fracture; "well" uses the whole well rate.
    wholeWell = f.prm_.get<std::string>("solver.volume_pacing_source", "perf") == "well";
    const double q_perf = wholeWell ? std::abs(f.well_rate_) : f.well_perf_rate_;
    const double v_avail = std::max(0.0, (q_perf - fp.flux) * dt) * fac;
    const double dV = fp.volume - v0;
    budgetSpent = dV > v_avail;
    if (verbosity > 0) {
        OpmLog::debug(fmt::format("VolumePacing: V={} V0={} dV={} q_perf={} Q_leak={} V_avail={}{}",
                                  fp.volume, v0, dV, q_perf, fp.flux, v_avail,
                                  budgetSpent ? "  -> BUDGET SPENT" : ""));
    }
    if (pressure && dV > v_avail && fp.volume > 0.0) {
        // eq. (3)-(4): draw the BC pressure down by the over-draw over the
        // secant storage compliance C_eff = V / (p_perf - p_res)
        const int pc = f.perfinj_.empty() ? 0 : std::get<0>(f.perfinj_[0]);
        const double p_res = (pc < static_cast<int>(f.reservoir_.pressure.size()))
            ? f.reservoir_.pressure[pc] : 0.0;
        const double dp0 = std::max(perfP0 - p_res, 1.0e5);
        const double c_eff = fp.volume / dp0;
        const double p_new = std::max(perfP0 - (dV - v_avail) / c_eff, p_res);
        if (verbosity > 0) {
            OpmLog::debug(fmt::format("VolumePacing: perf pressure {} -> {} bar (C_eff={})",
                                      f.perf_pressure_ / 1e5, p_new / 1e5, c_eff));
        }
        f.perf_pressure_ = p_new;
    }
}

void
Fracture::VolumePacing::restorePerfPressure(Fracture& f) const
{
    if (on && pressure) {
        f.perf_pressure_ = perfP0;
    }
}

void
Fracture::reinitializeTrimeshState_()
{
    OpmLog::debug("Reinitializing fracture states before IF with trimesh propagation");
    fracture_width_ = 1e-3; // Ensure not completely closed
    // start from the confining stress; this also sizes fracture_pressure_
    auto traction = fracture_pressure_;
    normalFractureTraction(traction);
    for (size_t i = 0; i < traction.size(); ++i) {
        fracture_pressure_[i] = std::max(traction[i], fracture_pressure_[i]);
    }
    if (numWellEquations() > 0) {
        fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];
    }
}

Fracture::FilterCakeSnapshot
Fracture::snapshotFilterCake_(const PropagationOptions& opt) const
{
    FilterCakeSnapshot s;
    s.thickness = filtercake_thikness_;
    s.gridMeshMap = grid_mesh_map_;
    s.total = filterCakeVolume();
    if (opt.conservativeCake) {
        s.volumes = filterCakeVolumesByCell();
    }
    return s;
}

void
Fracture::remapState_(const RegridMaps& maps, const int level, const PropagationOptions& opt,
                      const FilterCakeSnapshot& cake0)
{
    if (opt.remapSolution) {
        for (size_t i = 0; i < fracture_pressure_.size(); ++i) {
            assert(std::abs(fracture_pressure_[i][0]) < 1e10);
        }
        for (size_t i = 0; i < fracture_width_.size(); ++i) {
            assert(std::abs(fracture_width_[i][0]) < 0.6);
        }
        const auto old_size = fracture_width_.size();
        redistribute_values(fracture_width_, maps.before, maps.after, level, /*point_wise*/ true);
        double well_value = 0.0;
        if (numWellEquations() > 0) {
            well_value = fracture_pressure_[fracture_pressure_.size() - 1];
            fracture_pressure_.resize(fracture_pressure_.size() - 1);
        }
        if (external_cell_mask_.size() == old_size) {
            Dune::BlockVector<Dune::FieldVector<double, 1>> m(external_cell_mask_.size());
            for (size_t i = 0; i < m.size(); ++i) {
                m[i] = external_cell_mask_[i] ? 1.0 : 0.0;
            }
            redistribute_values(m, maps.before, maps.after, level, /*point_wise*/ true);
            external_cell_mask_.assign(m.size(), 0);
            for (size_t i = 0; i < m.size(); ++i) {
                external_cell_mask_[i] = (m[i][0] > 0.5) ? 1 : 0;
            }
        } else {
            external_cell_mask_.clear();
        }
        redistribute_values(fracture_pressure_, maps.before, maps.after, level, /*point_wise*/ true);
        if (numWellEquations() > 0) {
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
        // cell entries only: reservoir_.pressure has no well-equation entry
        const size_t num_press_cells = fracture_pressure_.size() - numWellEquations();
        for (size_t i = 0; i < num_press_cells; ++i) {
            if (fracture_pressure_[i] == 0.0) {
                fracture_pressure_[i] = hydrostaticReservoirPressure(i);
            }
        }
        if (numWellEquations() > 0 && fracture_pressure_[num_press_cells] == 0.0) {
            fracture_pressure_[num_press_cells] = fracture_pressure_[0];
        }
    } else {
        initFractureWidth();
        initFracturePressureFromReservoir();
    }

    filtercake_thikness_ = opt.conservativeCake
        ? remapFilterCakeConservative(cake0.volumes, cake0.gridMeshMap, level)
        : redistribute_values(cake0.thickness, cake0.gridMeshMap, maps.after, level,
                              /*point_wise*/ true);
    const double total_filtercake = filterCakeVolume();
    if (std::abs(total_filtercake - cake0.total) > 1e-6 && (opt.nlinVerbosity > 1)) {
        OpmLog::debug(fmt::format("Total filtercake thickness before remapping: {}\n"
                                  "Total filtercake thickness after remapping: {}",
                                  cake0.total, total_filtercake));
    }
}

Fracture::GridSolve
Fracture::solveOnCurrentGrid_(const PropagationOptions& opt)
{
    int iter = 0;
    ++joint_scores_;
    // a grow-only round scores the current state without solving
    bool converged = joint_growth_round_;
    if (!joint_growth_round_) {
        inner_ = InnerSolveOptions::read(prm_);
        coupledSolver_.setOptions(inner_.linear);
        while (!(converged = fullSystemIteration(opt.tol, iter)) && iter++ < opt.maxIter) {
        }
    }
    return {converged, iter};
}

void
Fracture::reportGridSolve_(const GridSolve& solve, const PropagationOptions& opt)
{
    if (!solve.converged) {
        last_solve_stats_.converged = false;
        if (prm_.get<bool>("solver.failure_on_nonconvergence", false)) {
            throw std::runtime_error("Fracture solver did not converge within "
                                     + std::to_string(opt.maxIter) + " iterations.");
        } else {
            OpmLog::warning(fmt::format("Fracture solver did not converge within {} iterations.",
                                        opt.maxIter));
        }
    }
    if (opt.nlinVerbosity > 0) {
        const int num_closed_cells = std::accumulate(closed_cells_.begin(), closed_cells_.end(), 0);
        OpmLog::debug(fmt::format("Nonlinear iterations needed with fixed expansion: {}\n"
                                  " fracture cells: {} closed cells: {}",
                                  solve.iterations, fracture_width_.size(), num_closed_cells));
    }
}

std::vector<double>
Fracture::propagationScores_(std::map<CellRef, int>& boundary, const GridSolve& solve,
                             const PropagationOptions& opt) const
{
    const std::vector<double> K1 = stressIntensityK1();
    const std::vector<CellRef> boundary_cells = trimesh_->boundaryCells();
    // operator[] as before: a front cell missing from the map scores grid cell 0
    std::vector<int> front(boundary_cells.size());
    for (size_t i = 0; i != front.size(); ++i) {
        front[i] = boundary[boundary_cells[i]];
    }
    std::vector<double> score(front.size());
    for (size_t i = 0; i != score.size(); ++i) {
        score[i] = K1[front[i]] / reservoir_.cstress[front[i]];
    }
    vetoForceLimit_(score, front, opt);
    // K1 from an unconverged solve is unreliable, and growth on it cannot be undone.
    if (opt.conservative && !solve.converged) {
        std::fill(score.begin(), score.end(), -1.0);
    }
    if (opt.convergedCellsOnly) {
        vetoUnconvergedCells_(score, front, opt);
    }
    if (suppress_propagation_ || opt.disable) {
        std::fill(score.begin(), score.end(), -1.0);
    }
    if (opt.openFrontOnly) {
        vetoClosedFront_(score, front);
    }
    return score;
}

void
Fracture::vetoForceLimit_(std::vector<double>& score, const std::vector<int>& front,
                          const PropagationOptions& opt) const
{
    for (size_t i = 0; i != score.size(); ++i) {
        if (fractureForce(front[i]) > opt.forceLimit) {
            score[i] = -1.0;
        }
    }
}

// A closed or still-toggling front cell may not grow; unlike conservative
// propagation this leaves the trustworthy part of the front free.
void
Fracture::vetoUnconvergedCells_(std::vector<double>& score, const std::vector<int>& front,
                                const PropagationOptions& opt) const
{
    int vetoed = 0;
    for (size_t i = 0; i != score.size(); ++i) {
        const int bind = front[i];
        if (bind < 0)
            continue;
        const bool closed = bind < static_cast<int>(closed_cells_.size()) && closed_cells_[bind];
        const bool unconverged
            = (bind < static_cast<int>(fb_cell_residual_.size()) && fb_cell_residual_[bind] > opt.fbTol)
            || (bind < static_cast<int>(cell_last_toggle_iter_.size())
                && current_solve_iter_ - cell_last_toggle_iter_[bind] < opt.stableIters);
        if (closed || unconverged) {
            score[i] = -1.0;
            ++vetoed;
        }
    }
    if (vetoed > 0 && opt.verbosity > 1)
        OpmLog::debug(fmt::format("Propagation veto: {}/{} front cells closed or unconverged",
                                  vetoed, score.size()));
}

// A closed front cell carries no opening load; growth from it would ride on the width floor.
void
Fracture::vetoClosedFront_(std::vector<double>& score, const std::vector<int>& front) const
{
    for (size_t i = 0; i != score.size(); ++i) {
        const int bind = front[i];
        if (bind >= 0 && bind < static_cast<int>(closed_cells_.size()) && closed_cells_[bind])
            score[i] = -1.0;
    }
}

void
Fracture::initFractureStates()
{
    this->initFractureWidth();
    this->initFracturePressureFromReservoir();
}

void
Fracture::initFractureWidth()
{
    size_t nc = numFractureCells();
    fracture_width_.resize(nc);
    fracture_width_ = prm_.get<double>("config.initial_fracture_width");
}
void
Fracture::initFracturePressureFromReservoir()
{
    size_t nc = reservoir_map_.cells.size();
    fracture_pressure_.resize(nc + numWellEquations());
    fracture_pressure_ = 0;
    for (size_t i = 0; i < nc; ++i) {
        fracture_pressure_[i] = hydrostaticReservoirPressure(i);
    }
}

double
Fracture::hydrostaticReservoirPressure(size_t i) const
{
    // reservoir pressure moved hydrostatically from the reservoir cell depth to the fracture cell depth
    if (i >= fracture_dgh_.size()) {
        return reservoir_.pressure[i];
    }
    return reservoir_.pressure[i] + fracture_dgh_[i]
        - gravity_ * reservoir_.density[i] * reservoir_.cellZ[i];
}




bool
Fracture::ensureLeakoffCurrent()
{
    const std::size_t nc = numFractureCells();
    if (leakof_.size() == nc) {
        return true;
    }

    const bool inputsCurrent = (reservoir_.mobility.size() == nc)
        && (reservoir_.perm.size() == nc)
        && (reservoir_map_.dist.size() == nc)
        && (!has_filtercake_ || (filtercake_thikness_.size() == nc));

    if (!inputsCurrent) {
        return false;
    }

    updateLeakoff();

    return leakof_.size() == nc;
}

// Presets: legacy = 2 faces, d = L/4 (8 lambda k A/L, frozen default);
// onesided = 1 face, d = L/4 (4 lambda k A/L); reveal = 1 face, d = 0.632 L
// (1.58 lambda k A/L, fitted to the dz20 reference — no physical argument).
// solver.leakoff_sides / solver.leakoff_distance_factor override the preset.
int
Fracture::leakoffSides() const
{
    const std::string model = prm_.get<std::string>("solver.leakoff_model", "legacy");
    int sides = 2;
    if (model == "onesided" || model == "reveal") {
        sides = 1;
    } else if (model != "legacy") {
        OPM_THROW(std::runtime_error, "Unknown solver.leakoff_model: " + model);
    }
    return prm_.get<int>("solver.leakoff_sides", sides);
}

double
Fracture::leakoffDistanceFactor() const
{
    const std::string model = prm_.get<std::string>("solver.leakoff_model", "legacy");
    const double factor = (model == "reveal") ? 0.632 : 0.25;
    return prm_.get<double>("solver.leakoff_distance_factor", factor);
}

// Which mobility multiplies k*A/d on the leak-off connection. Legacy sums every
// phase mobility of the reservoir cell; upwind mirrors the flow's water equation --
// the fracture cell when it is feeding the formation, the reservoir cell's water
// mobility when the flow reverses.
double
Fracture::leakoffMobilityFor(const int eIdx, const bool upwind) const
{
    const double total = reservoir_.mobility[eIdx];
    if (!upwind) {
        return total;
    }

    const double dpot = leakoffPotentialDifference(eIdx);
    if (std::isnan(dpot)) {
        return total;
    }

    if (dpot > 0.0) {
        // Out of the fracture: the fracture cell is water filled by construction, so
        // its relative permeability is one and the mobility is the reciprocal
        // viscosity. Without the fluid system's water properties there is nothing
        // better than the legacy value.
        if (!fracture_water_property_evaluator_) {
            return total;
        }
        const auto props = fracture_water_property_evaluator_(eIdx, fracture_pressure_[eIdx][0]);
        const double mu = props.second.value;
        return (mu > 0.0) ? 1.0 / mu : total;
    }

    return (static_cast<size_t>(eIdx) < reservoir_.waterMobility.size())
        ? reservoir_.waterMobility[eIdx] : total;
}

// Fracture minus reservoir potential across the leak-off connection; NaN without a state.
double
Fracture::leakoffPotentialDifference(const int eIdx) const
{
    const bool have_state = (static_cast<size_t>(eIdx) < fracture_pressure_.size())
        && (static_cast<size_t>(eIdx) < reservoir_.pressure.size());
    if (!have_state) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double dh_frac = (static_cast<size_t>(eIdx) < fracture_dgh_.size())
        ? fracture_dgh_[eIdx] : 0.0;
    const double dh_res = reservoir_.cellZ[eIdx] * gravity_ * reservoir_.density[eIdx];
    return (fracture_pressure_[eIdx][0] - dh_frac) - (reservoir_.pressure[eIdx] - dh_res);
}

// Outflow: the water-filled fracture cell's mobility; inflow: all reservoir phases.
double
Fracture::flowLeakoffMobility(const int eIdx) const
{
    const double dpot = leakoffPotentialDifference(eIdx);
    return (dpot > 0.0) ? leakoffMobilityFor(eIdx, /*upwind=*/true) : reservoir_.mobility[eIdx];
}

void
Fracture::updateLeakoff()
{

    // find boundary nodes.
    constexpr int dim = 2;
    constexpr Dune::Codim<2> codimVertex;
    constexpr Dune::Codim<1> codimFace;
    constexpr Dune::Codim<0> codimCell;
    auto gv = grid_->leafGridView();
    auto& indexSet = gv.indexSet();
    std::vector<bool> isBoundaryNode(indexSet.size(codimVertex), false);
    for (const auto& element : elements(gv)) {
        for (const auto& is : intersections(gv, element)) {
            if (is.boundary()) {
               auto lf = is.indexInInside();
               auto face = element.subEntity<codimFace>(lf);
               int nv = face.subEntities(codimVertex);
               assert(nv == 2);
               const auto refElement = Dune::ReferenceElements<double,dim>::general(element.type());
               for (int v = 0; v < nv; v++) {
                 auto&& vertex = element.template subEntity<dim>(refElement.subEntity(lf, dim-1, v, dim));
                 auto vid = indexSet.index(vertex);
                 isBoundaryNode[vid] = true;
               }
            }
        }
    }
    // find all cells with  aboundary node
    std::vector<bool> elementHasBoundaryNode(gv.size(codimCell), false);
    for (const auto& element : elements(gv)) {
        auto eidx = indexSet.index(element);
        if(true){
        int nv = element.subEntities(codimVertex);
        for (int i = 0; i < nv; i++) {
          auto vertex = element.subEntity<codimVertex>(i);
            auto vid = indexSet.index(vertex);
            if (isBoundaryNode[vid]) {
                elementHasBoundaryNode[eidx] = true;
                break; // no need to check more vertices
            }
        }
        }
    }
    bool no_leakof_outercells = prm_.get<bool>("solver.no_leakof_outercells",false);
    const size_t nc = numFractureCells();
    leakof_.resize(nc, 0.0);
    leakoff_mobility_.resize(nc, 0.0);
    // "reservoir_total" (legacy): the reservoir cell's summed phase mobilities.
    // "upwind": what the flow's water equation actually applies across the same
    // connection, so that the embedded representation and the fracture's own solve
    // form the same conductance instead of differing by a mobility ratio.
    const std::string mob_mode =
        prm_.get<std::string>("solver.leakoff_mobility", "reservoir_total");
    const bool upwind_mobility = (mob_mode == "upwind");
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    // per leaking face: reservoir path in series with its share of the cake
    // (filtercake_thikness_ is the total over all faces); faces act in parallel
    const double sides = leakoffSides();
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        double area = geom.volume();
        double res_mob = leakoffMobilityFor(eIdx, upwind_mobility);
        leakoff_mobility_[eIdx] = res_mob;
        leakof_[eIdx] = res_mob * reservoir_.perm[eIdx] * area / reservoir_map_.dist[eIdx];
        double invtrans = 1 / leakof_[eIdx];
        if (has_filtercake_) {
            assert(filtercake_thikness_[eIdx] >= 0.0);
            if (filtercake_thikness_[eIdx] > 0.0) {
                double fitercaketrans = res_mob * filtercake_perm_ * area / (filtercake_thikness_[eIdx]/sides);
                invtrans += 1 / fitercaketrans;
            }
        }
        leakof_[eIdx] = sides / invtrans;
        if(elementHasBoundaryNode[eIdx] && no_leakof_outercells){
          leakof_[eIdx] = 0.0;
        }
    }
}
std::vector<Fracture::Htrans>
Fracture::currentHalfTrans() const
{
    std::vector<Htrans> htrans;

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        for (const auto& is : Dune::intersections(grid_->leafGridView(), element)) {
            if (is.boundary()) {
                continue;
            }

            const int nIdx = mapper.index(is.outside());
            if (eIdx >= nIdx) {
                continue; // each interior face once
            }

            const auto eCenter = geom.center();
            const auto nCenter = is.outside().geometry().center();
            const auto isCenter = is.geometry().center();
            const auto d_inside = eCenter - isCenter;
            const auto d_outside = nCenter - isCenter;
            const double area = is.geometry().volume();

            htrans.emplace_back(nIdx, eIdx,
                                area / d_inside.two_norm(),
                                area / d_outside.two_norm());
        }
    }

    return htrans;
}

std::vector<double>
Fracture::cellAreas() const
{
    std::vector<double> areas(numFractureCells(), 0.0);

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid_->leafGridView())) {
        areas[mapper.index(element)] = element.geometry().volume();
    }

    return areas;
}

std::vector<double>
Fracture::cellDepths() const
{
    std::vector<double> depths(numFractureCells(), 0.0);

    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid_->leafGridView())) {
        depths[mapper.index(element)] = element.geometry().center()[2];
    }

    return depths;
}

double
Fracture::filterCakeVolume() const
{
    double total_filtercake_0 = 0.0;
    for (const auto& elment : Dune::elements(grid_->leafGridView())) {
        int i = grid_->leafGridView().indexSet().index(elment);
        auto geom = elment.geometry();
        double area = geom.volume();
        total_filtercake_0 += filtercake_thikness_[i] * area;
    }
    return total_filtercake_0;
}

RegularTrimesh::CellVolumes
Fracture::filterCakeVolumesByCell() const
{
    RegularTrimesh::CellVolumes cake;
    const size_t nc = numFractureCells();
    cake.volume.resize(nc, 0.0);
    cake.group.assign(nc, -1);
    cake.centre.resize(nc);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid_->leafGridView())) {
        const auto eIdx = mapper.index(element);
        const auto geom = element.geometry();
        const auto c = geom.center();
        cake.centre[eIdx] = {c[0], c[1], c[2]};
        if (eIdx < filtercake_thikness_.size()) {
            cake.volume[eIdx] = filtercake_thikness_[eIdx] * geom.volume();
        }
        if (eIdx < reservoir_map_.cells.size()) {
            cake.group[eIdx] = reservoir_map_.cells[eIdx];
        }
    }
    return cake;
}

std::vector<double>
Fracture::remapFilterCakeConservative(const RegularTrimesh::CellVolumes& source,
                                      const std::vector<std::vector<CellRef>>& source_map,
                                      const int level) const
{
    const size_t nc = numFractureCells();
    std::vector<double> area(nc, 0.0);
    std::vector<int> group(nc, -1);
    std::vector<Coord3D> centre(nc);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid_->leafGridView())) {
        const auto eIdx = mapper.index(element);
        const auto geom = element.geometry();
        const auto c = geom.center();
        area[eIdx] = geom.volume();
        group[eIdx] = reservoir_map_.cells[eIdx];
        centre[eIdx] = {c[0], c[1], c[2]};
    }
    return RegularTrimesh::conservativeRemap(source, source_map, grid_mesh_map_, level, area, group, centre);
}

void Fracture::resetFracture()
{
    if (!trimesh_prev_) {
        // No checkpoint available: moveForwardInTime was never called before this
        // step, which happens when the fracture model was first initialized during
        // a Newton iteration that subsequently failed. Nothing to restore.
        return;
    }
    // Restore trimesh from checkpoint.  Use the saved well-source CellRefs so
    // that the Dune grid and grid_mesh_map_ are consistent with the checkpoint.
    trimesh_ = std::make_unique<Opm::RegularTrimesh>(*trimesh_prev_);
    const int MAX_NUM_COARSENING = prm_.get<int>("solver.max_num_coarsening");
    const int numcell_threshold  = prm_.get<int>("solver.numcell_threshold");
    bool smooth_boundary = prm_.get<bool>("solver.smooth_boundary", false);
    std::vector<CellRef> wsources = well_source_cellref_prev_;
    auto [grid, fsmap, bmap] =
        trimesh_->createDuneGrid(MAX_NUM_COARSENING, wsources, smooth_boundary, numcell_threshold);
    grid_mesh_map_ = fsmap;
    grid_fixed_cells_ = wsources;
    // Restore grid and cell normals (setFractureGrid clears reservoir_map_.cells;
    // they will be repopulated by FractureModel::resetFractures via
    // updateFractureReservoirCells + initReservoirProperties).
    setFractureGrid(std::move(grid));
    // Restore well source indices from checkpoint
    well_source_cellref_ = well_source_cellref_prev_;
    well_source_         = well_source_prev_;
    // Restore state variables from checkpoint
    fracture_width_      = fracture_width_prev_;
    fracture_pressure_   = fracture_pressure_prev_;
    filtercake_thikness_ = filtercake_thikness_prev_;
    // Invalidate derived/computed state
    coupling_matrix_ = nullptr;
    fracture_matrix_ = nullptr;
    rhs_pressure_.resize(0);
    // drop CTF-acceleration history; it belongs to the rolled-back solves
    wi_vec_mix_.has_prev_residual = false;
    wi_vec_mix_.prev_residual.clear();
    width_round_valid_ = false; // relaxation history belongs to the rolled-back solves too
    // Note: reservoir_map_.cells, reservoir_map_.allCells, reservoir_properties_, etc.
    // are recalculated by FractureModel::resetFractures after this call returns.
}
void Fracture::rebuildGridForLoad(const std::vector<std::vector<CellRef>>& saved_map)
{
    const int max_coarsening = prm_.get<int>("solver.max_num_coarsening");
    const int numcell_threshold = prm_.get<int>("solver.numcell_threshold");
    const bool smooth_boundary = prm_.get<bool>("solver.smooth_boundary", false);
    auto built = trimesh_->createDuneGrid(max_coarsening, grid_fixed_cells_, smooth_boundary, numcell_threshold);
    if (std::get<1>(built) != saved_map) {
        // a fracture that has not grown still has its seed grid, see init()
        built = trimesh_->createDuneGrid(1, grid_fixed_cells_);
    }
    if (std::get<1>(built) != saved_map) {
        OPM_THROW(std::runtime_error,
                  "Fracture " + name() + ": the saved mesh gives another grid with this run's "
                  "solver.max_num_coarsening, solver.numcell_threshold or solver.smooth_boundary; "
                  "these cannot change after a load.");
    }
    grid_mesh_map_ = saved_map;
    setFractureGrid(std::move(std::get<0>(built)));
}

void Fracture::moveForwardInTime(double dt_last)
{
    // CTF-vector acceleration: the secant information is only meaningful within
    // one timestep's outer iterations — the target jumps physically across the
    // step boundary. Keep prev_mixed (damping continuity) but drop the history.
    wi_vec_mix_.has_prev_residual = false;
    wi_vec_mix_.prev_residual.clear();
    width_round_valid_ = false; // state relaxation history is per-timestep too
    // checkpoint values for the per-step growth guard / onset dt-hold
    if (active_) {
        const auto fp = calculateFractureProperties();
        // Coupling-rate dt controller (opt-in via solver.coupling_dt_pressure_target,
        // bar per step): suggest the next flow step from how fast the coupling state
        // moved in the completed step - small steps arise exactly during transients
        // (establishment, cake re-pressurisation) and release when the state is quiet.
        const double dp_target = prm_.get<double>("solver.coupling_dt_pressure_target", 0.0) * 1e5;
        if (dp_target > 0.0 && dt_last > 0.0 && perf_pressure_step_start_ >= 0.0) {
            const double relA_target = prm_.get<double>("solver.coupling_dt_area_target", 0.2);
            const double dp = std::abs(perf_pressure_ - perf_pressure_step_start_);
            const double relA = (area_step_start_ > 0.0)
                ? std::abs(fp.area - area_step_start_) / area_step_start_ : 0.0;
            const bool violated = (dp > dp_target)
                || (relA_target > 0.0 && relA > relA_target);
            // Floor at the proven establishment step: per-step coupling changes do
            // not vanish with dt during the bang-bang transient, so an unbounded
            // multiplicative cap would ratchet dt to useless values.
            const double dt_min = prm_.get<double>("solver.coupling_dt_min", 0.1) * 86400.0;
            const int quiet_needed = prm_.get<int>("solver.coupling_dt_quiet_steps", 5);
            if (violated) {
                // shrink immediately, proportionally to the excursion
                double fac = 1.0;
                if (dp > dp_target) fac = std::min(fac, dp_target / dp);
                if (relA_target > 0.0 && relA > relA_target)
                    fac = std::min(fac, relA_target / relA);
                coupling_quiet_steps_ = 0;
                coupling_dt_cap_ = std::max(dt_last * std::max(fac, 0.5), dt_min);
            } else if (++coupling_quiet_steps_ >= quiet_needed) {
                // transients settle in episodes: release only after a run of
                // quiet steps, and gently
                coupling_dt_cap_ = std::max(dt_last * 1.5, dt_min);
            } else {
                coupling_dt_cap_ = std::max(dt_last, dt_min); // hold
            }
        }
        // Stress-rate dt controller (opt-in via solver.stress_dt_target, bar per
        // step): pace dt so the mean fracture traction moves at most the target
        // per step. The thermal transient is smooth, so a plain proportional
        // controller suffices (no quiet-step episodes needed).
        const double ds_target = prm_.get<double>("solver.stress_dt_target", 0.0) * 1e5;
        const size_t ncf = numFractureCells();
        if (ds_target > 0.0 && dt_last > 0.0 && reservoir_.stress.size() >= ncf
            && cell_normals_.size() >= ncf) {
            double tsum = 0.0;
            for (size_t i = 0; i < ncf; ++i)
                tsum += normalFractureTraction(i);
            const double tmean = (ncf > 0) ? tsum / ncf : 0.0;
            if (std::isfinite(traction_step_start_)) {
                const double ds = std::abs(tmean - traction_step_start_);
                const double dt_min = prm_.get<double>("solver.stress_dt_min", 0.1) * 86400.0;
                const double fac = std::clamp(ds_target / std::max(ds, 1.0), 0.5, 2.0);
                stress_dt_cap_ = std::max(dt_last * fac, dt_min);
            }
            traction_step_start_ = tmean;
        }
        area_step_start_ = fp.area;
        volume_step_start_ = fp.volume;
    } else {
        area_step_start_ = -1.0;
        volume_step_start_ = -1.0;
    }
    perf_pressure_step_start_ = perf_pressure_;
    // copy current state to "previous" state
    trimesh_prev_ = trimesh_ ? std::make_unique<Opm::RegularTrimesh>(*trimesh_) : nullptr;
    fracture_width_prev_ = fracture_width_;
    fracture_pressure_prev_ = fracture_pressure_;
    filtercake_thikness_prev_ = filtercake_thikness_; // properties of filter cake, if any
    // save well source indices (needed to rebuild the dune grid on reset)
    well_source_cellref_prev_ = well_source_cellref_;
    well_source_prev_ = well_source_;
    // reservoir cell mappings are NOT stored here; they are recomputed on reset
    // via FractureModel::resetFractures → updateFractureReservoirCells + initReservoirProperties.
}



void
Fracture::composeFlowStateDonor(const std::vector<std::vector<CellRef>>& map1,
                                const std::vector<std::vector<CellRef>>& map2,
                                const int level)
{
    // dominant overlap wins: the flow state is a set of primary variables, not averaged
    const auto g2g = RegularTrimesh::createGridToGridMap(map1, map2, level);
    std::vector<int> donor(map2.size(), -1);
    std::vector<double> best(map2.size(), 0.0);
    for (const auto& [i, j, w] : g2g) {
        if (w > best[j]) {
            best[j] = w;
            if (!flow_state_remapped_) {
                donor[j] = static_cast<int>(i);
            } else {
                donor[j] = (i < flow_state_donor_.size()) ? flow_state_donor_[i] : -1;
            }
        }
    }
    flow_state_donor_ = std::move(donor);
    flow_state_remapped_ = true;
}

void
Fracture::redistribute_values(Dune::BlockVector<Dune::FieldVector<double, 1>>& values,
                              const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level,
                              bool point_wise)
{
    std::vector<double> tmp_values(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        tmp_values[i] = values[i][0]; // assuming values is a BlockVector with one component
    }
    tmp_values = redistribute_values(tmp_values, map1, map2, level, point_wise);
    values.resize(tmp_values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        values[i][0] = tmp_values[i]; // assuming values is a BlockVector with one component
    }
}


std::vector<double>
Fracture::redistribute_values(const std::vector<double>& values,
                              const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level,
                              bool point_wise)
// ----------------------------------------------------------------------------
{
    const auto g2gmap = RegularTrimesh::createGridToGridMap(map1, map2, level);

    std::vector<double> redistributed_values(map2.size(), 0.0);
    std::vector<double> weight(map2.size(), 0.0);

    for (const auto& e : g2gmap) {
        assert(std::get<0>(e) < values.size());
        assert(std::get<1>(e) < redistributed_values.size());
        redistributed_values[std::get<1>(e)] += values[std::get<0>(e)] * std::get<2>(e);
        weight[std::get<1>(e)] += std::get<2>(e);
    }
    if (point_wise) {
        for (size_t i = 0; i < redistributed_values.size(); ++i) {
            assert(weight[i] >= 0.0);
            if (weight[i] > 0.0) {
                redistributed_values[i] /= weight[i];
            } else {
                redistributed_values[i] = 0.0; // or some other default value
            }
        }
    }

    return redistributed_values;
}

void
updateMaxMinVector(std::array<double, 2>& dvalvec, double val)
{
    if (val < 0) {
        dvalvec[0] = std::min(val, dvalvec[0]);
    } else {
        dvalvec[1] = std::max(val, dvalvec[1]);
    }
}


FractureProperties
Fracture::calculateFractureProperties() const
{
    std::array<double, 2> dhvec({0.0, 0.0});
    std::array<double, 2> dwvec({0.0, 0.0});
    double total_flux(0);
    double area(0.0);
    double volume(0.0);
    double filter_volume(0.0);
    std::vector<double> leak_of_rate = leakOfRate();
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
      //NB we do not filter on with
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        area += geom.volume();
        volume += geom.volume() * (fracture_width_[eIdx]);
        filter_volume += geom.volume() * filtercake_thikness_[eIdx];
        auto dist = geom.center() - origo_;
        double dh = dist.dot(naxis_[0]);
        double dw = dist.dot(naxis_[1]);
        updateMaxMinVector(dhvec, dh);
        updateMaxMinVector(dwvec, dw);
        if(leak_of_rate.size() != numFractureCells()){
            continue;
        }
        total_flux += leak_of_rate[eIdx] * geom.volume();
    }
    double height = (dhvec[1] - dhvec[0])/2.0;
    double length = (dwvec[1] - dwvec[0])/2.0;
    auto WIs = wellIndices();
    double WI=0.0;
    for(auto wi : WIs){
      WI += wi.ctf;
    }
    double avgh = volume/area;
    double avgfilter_h = filter_volume/area;
    double inj_pressure = injectionPressure();
    double inj_bhp = injectionBhp();
    double surface_factor = 1.0; //NB to do convert to surface volume
    double inj_wellrate = well_rate_*surface_factor;//only copy of wellrate from wells
    double solid_filter_volume = filter_volume*(1 - filtercake_poro_);// return the solid volume
    FractureProperties fracprop(height,
                                length,
                                total_flux,
                                area,
                                WI,
                                volume,
                                solid_filter_volume,
                                avgh,
                                avgfilter_h,
                                inj_pressure,
                                inj_bhp,
                                inj_wellrate
                            );
    return fracprop;
}

void
Fracture::initPressureMatrix()
{
    // Flow from wells to fracture cells
    const double fWI_cfg = prm_.get<double>("fractureWI");
    // fractureWI_mode: "constant" (config fractureWI, legacy) or "estimate"
    // (radial flow in the fracture plane towards the wellbore, cubic-law
    // permeability: WI = 2*pi*(w^3/12)/ln(re/rw), same estimate as the embedded
    // aux-cell path; w = max(current width, flow floor)).
    const std::string fwi_mode = prm_.get<std::string>("solver.fractureWI_mode", "constant");
    const double fwi_rw = prm_.get<double>("solver.fractureWI_rw", 0.1);
    perfinj_.clear();
    htrans_.clear();
    ElementMapper wsmapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    std::vector<double> cell_area(numFractureCells(), 0.0);
    for (const auto& e : Dune::elements(grid_->leafGridView()))
        cell_area[wsmapper.index(e)] = e.geometry().volume();
    auto fwi_for = [&](int cell) {
        if (fwi_mode != "estimate") return fWI_cfg;
        const double w = std::max(cell < static_cast<int>(fracture_width_.size()) ? fracture_width_[cell][0] : 0.0,
                                  min_width_);
        const double re = std::sqrt(std::max(cell_area[cell], 1e-12) / M_PI);
        const double lnTerm = std::log(std::max(re / fwi_rw, 1.1));
        return 2.0 * M_PI * (w * w * w / 12.0) / lnTerm;
    };
    // Well sources: the seed ring (legacy) and, with solver.well_source_all_perfs,
    // every fracture cell that lies in a reservoir cell the well perforates — so a
    // fracture that grows across the perforated interval is fed along it (as in a
    // rate-fed fracture model) instead of only at the seed.
    std::set<int> sources(well_source_.begin(), well_source_.end());
    // Embedded flow makes every source cell a well perforation, so feeding the
    // whole fracture is costly for the same answer; default there is the seed
    // ring, solver.embedded_well_source=config follows the well_source_* settings.
    const bool embedded_ring =
        prm_.get<std::string>("solver.fracture_flow_mode", "wi_upscaling") == "embedded"
        && prm_.get<std::string>("solver.embedded_well_source", "ring") == "ring";
    {
        const auto ews = prm_.get<std::string>("solver.embedded_well_source", "ring");
        if (ews != "ring" && ews != "config") {
            OPM_THROW(std::runtime_error, "Unknown solver.embedded_well_source: " + ews);
        }
    }
    if (!embedded_ring && prm_.get<bool>("solver.well_source_all_perfs", false) && !well_perf_cells_.empty()) {
        // The wellbore cuts the fracture along a line, so only cells within about a
        // wellbore radius of it are fed directly; on a coarse grid "shares a reservoir
        // cell with a perforation" is the whole fracture and short-circuits the cubic
        // law. solver.well_source_radius > 0 applies the limit (vertical wells: the
        // vertical projected into the fracture plane through origo_).
        const double wsr = prm_.get<double>("solver.well_source_radius", 0.0);
        Point3D wdir({0.0, 0.0, 1.0});
        if (wsr > 0.0) {
            const double dn = wdir * naxis_[2];
            for (int d = 0; d < 3; ++d) wdir[d] -= dn * naxis_[2][d];
            const double wn = wdir.two_norm();
            if (wn > 1e-12) wdir /= wn;
        }
        std::set<int> perfcells(well_perf_cells_.begin(), well_perf_cells_.end());
        for (const auto& e : Dune::elements(grid_->leafGridView())) {
            const int i = wsmapper.index(e);
            if (i >= static_cast<int>(reservoir_map_.cells.size())) continue;
            if (!perfcells.count(reservoir_map_.cells[i])) continue;
            if (wsr > 0.0) {
                const auto c = e.geometry().center();
                Point3D r({c[0] - origo_[0], c[1] - origo_[1], c[2] - origo_[2]});
                const double along = r * wdir;
                for (int d = 0; d < 3; ++d) r[d] -= along * wdir[d];
                if (r.two_norm() > wsr) continue; // too far from the wellbore line
            }
            sources.insert(i);
        }
    }
    for (int cell : sources) {
        perfinj_.push_back({cell, fwi_for(cell)});
    }
    // flow between fracture cells
    const size_t nc = numFractureCells() + numWellEquations();
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        // iterator over all intersections
        for (auto& is : Dune::intersections(grid_->leafGridView(), element)) {

            if (!is.boundary()) {
                auto eCenter = geom.center();
                int nIdx = mapper.index(is.outside());
                if (eIdx < nIdx) {
                    // calculate distance between the midpoints
                    auto nCenter = is.outside().geometry().center();
                    auto isCenter = is.geometry().center();
                    auto d_inside = eCenter - isCenter;
                    auto d_outside = nCenter - isCenter;
                    // probably should use projected distance
                    auto igeom = is.geometry();
                    double area = igeom.volume();
                    double h1 = area / d_inside.two_norm();
                    double h2 = area / d_outside.two_norm();

                    {
                        Htrans matel(nIdx, eIdx, h1, h2);
                        htrans_.push_back(matel);
                    }
                }
            }
        }
    }
    // not need if build mode is implicit
    //  build matrix
    // overflow budget must hold the well row (one entry per well-source cell)
    const double overflow = 0.4 + (nc > 0 ? 2.0 * perfinj_.size() / static_cast<double>(nc) : 0.0);
    pressure_matrix_ = std::make_unique<Matrix>(nc, nc, 4, overflow, Matrix::implicit);
    auto& matrix = *pressure_matrix_;
    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double zero_entry = 0.0;
        matrix.entry(i, j) = zero_entry;
        matrix.entry(j, i) = zero_entry;
        matrix.entry(j, j) = zero_entry;
        matrix.entry(i, i) = zero_entry;
    }
    // add in case a cell have only no internal boundaries
    for(int i=0; i < nc; ++i){
        matrix.entry(i, i) = 0.0;
    }

    if (numWellEquations() > 0) {
        assert(numWellEquations() == 1); // assumes a single well equation
        // add elements for last row and column
        const int weqix = nc - 1; // index of well equation
        for (const auto& pi : perfinj_) {
            const int i = std::get<0>(pi);
            matrix.entry(i, weqix) = 0.0;
            matrix.entry(weqix, i) = 0.0;
            matrix.entry(i, i) = 0.0;
        }
        matrix.entry(weqix, weqix) = 0.0;
    }
    matrix.compress();
}
void
Fracture::assemblePressure()
{
    // Delegate to the standalone implementation shared with unit tests and the
    // AD assembler path, eliminating the duplicate inline code that lived here.
    // updateLeakoff() is called by the caller before this function.
    auto result = assemblePressureOriginal(makePressureAssemblyInput());
    if (!pressure_matrix_) {
        // First call: take ownership of the newly created matrix.
        pressure_matrix_ = std::move(result);
    } else {
        // Subsequent calls: copy values so references to the matrix stay valid.
        *pressure_matrix_ = *result;
    }
}

bool
Fracture::removeNewZeroWithCells(RegularTrimesh& mesh,
                                 int cur_level,
                                 const RegularTrimesh& initial_mesh) const
{
    const double force_limit = prm_.get<double>("solver.force_limit", 0.0);
    // Opt-in (default off): also retract *established* cells (already active in the
    // initial mesh) that ended up zero-width or compressed, not just newly-grown
    // ones. Physically: a region that opened then, after later flow/mech solves,
    // carries net compression should not remain part of the fracture — currently
    // such a cell is only *closed* (contact, width pinned 0), a "half-open" state
    // that is hard for the nonlinear solve to resolve. Retracting it removes that
    // topology instead. WARNING: enabling this lets the fracture shrink, so it can
    // oscillate against growth across outer iterations; keep it gated/conservative.
    const bool retract_established =
        prm_.get<bool>("solver.retract_established_compressed", false);
    bool any_removed = false;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& cell : elements(grid_->leafGridView())) {
        const auto index = elemMapper.index(cell);
        if ((fracture_width_[index][0] == 0) || (fractureForce(index) > force_limit)) { // maybe to strong.
            auto cellrefs = grid_mesh_map_[index];
            for (const auto& cellref : cellrefs) {
                if (cur_level == 0) {
                    if (retract_established || !initial_mesh.isActive(cellref)) {
                        mesh.setInactive(cellref);
                        any_removed = true;
                    }
                } else {
                    // go through all fine cells of the initial mesh to see if they any is active
                    assert(cur_level > 0);
                    const auto& fine_cells = RegularTrimesh::coarse_to_fine(cellref, cur_level);
                    bool any_active = false;
                    bool all_active = true;
                    for (const auto& fine_cell : fine_cells) {
                        if (initial_mesh.isActive(fine_cell)) {
                            any_active = true;
                        }else{
                            all_active = false;
                        }
                    }
                    if (retract_established || !all_active) {
                        mesh.setInactive(cellref);
                        any_removed = true;
                    }
                }
            }
        }
    }
    return any_removed;
}



void
Fracture::warnIfOpenButNotConducting() const
{
    // An open, conductive fracture cannot carry zero rate; when it reports one the
    // well connection has been lost (e.g. legacy WI upscaling zeroing a connection
    // whose sign flips). Opt out with warn_open_not_conducting=0.
    if (!prm_.get<bool>("solver.warn_open_not_conducting", true) || !active_)
        return;
    if (fracture_width_.size() == 0)
        return;

    const auto fp = calculateFractureProperties();
    // "Open" means mechanically open, not merely meshed: cells carrying real
    // aperture above the flow floor.
    double open_area = 0.0;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& e : Dune::elements(grid_->leafGridView())) {
        const int i = mapper.index(e);
        if (i < static_cast<int>(fracture_width_.size())
            && fracture_width_[i][0] > min_width_)
            open_area += e.geometry().volume();
    }
    const double min_open = prm_.get<double>("solver.warn_open_area", 1.0); // m2
    const double q_tol = prm_.get<double>("solver.warn_flux_tol", 1e-12);   // m3/s
    if (open_area > min_open && std::abs(fp.flux) < q_tol) {
        OpmLog::warning("Fracture " + name() + ": " + std::to_string(open_area)
                        + " m2 is mechanically open but the fracture carries no flow ("
                        + std::to_string(fp.flux)
                        + " m3/s) - the well connection looks lost, not the fracture");
    }
}

double
Fracture::normalFractureTraction(size_t eIdx) const
{
    return ddm::tractionVoigt(reservoir_.stress[eIdx], cell_normals_[eIdx]);
}

void
Fracture::normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction,
                                 bool resize) const
{
    OPM_TIMEFUNCTION();
    const size_t nc = numFractureCells();
    if (resize)
        traction.resize(nc + numWellEquations());

    for (size_t eIdx = 0; eIdx < nc; ++eIdx)
        traction[eIdx] = normalFractureTraction(eIdx);
}
void
Fracture::assembleFractureMatrix() const
{
    OPM_TIMEFUNCTION();
    size_t nc = numFractureCells();
    if (!fracture_matrix_) {
        fracture_matrix_ = std::make_unique<DynamicMatrix>();
    }
    fracture_matrix_->resize(nc, nc);
    *fracture_matrix_ = 0.0;
    ddm::assembleMatrix_fast(*fracture_matrix_, E_, nu_, *grid_);
}


template void Fracture::assignGeoMechWellState(PerfData<float>& perfData) const;
template void Fracture::assignGeoMechWellState(PerfData<double>& perfData) const;

} // namespace Opm
