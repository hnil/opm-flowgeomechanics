#include "config.h"

#include <iostream> // std::cout
#include <sstream>
#include <string> // std::string

#include <opm/geomech/Fracture.hpp>

#include "opm/geomech/GridStretcher.hpp"
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

#include <opm/grid/polyhedralgrid.hh>

#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

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



namespace
{

/**
 * @brief Computes the target expansion for a fracture.
 *
 * This function calculates the target expansion based on the given parameters:
 * the target stress intensity factor (K1_target), the fracture aperture,
 * Young's modulus (E), and Poisson's ratio (nu).
 *
 * @param K1_target The target stress intensity factor.
 * @param aperture The fracture aperture.
 * @param E Young's modulus.
 * @param nu Poisson's ratio.
 * @return The computed target expansion.
 */
// double compute_target_expansion(const double K1_target,
//                                  const double aperture,
//                                  const double E, // young
//                                  const double nu) // poisson
// {
//   const double mu = E / (2 * (1+nu)); // shear modulus
//   const double fac = mu * std::sqrt(M_PI) /
//                      (2 * (1-nu) * 1.834);
//   return pow(fac * aperture / K1_target, 2);
// };

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
               const Opm::PropertyTree& prm)
{
    OPM_TIMEFUNCTION();
    max_flow_time_step_ = std::numeric_limits<double>::max();
    prm_ = prm;
    if(prm_.get<bool>("config.gravity_off", true)){
        gravity_ = 0.0;
    }

    min_width_ = prm_.get<double>("config.min_width", 1e-3);
    wellinfo_ = WellInfo({well, perf, well_cell, global_index, segment, perf_range});

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
    std::cout << "axis" << axis_[0][0] << "," << axis_[0][1] << "," << axis_[0][2] << " " << axis_[1][0]
              << "," << axis_[1][1] << "," << axis_[1][2] << " " << axis_[2][0] << "," << axis_[2][1]
              << "," << axis_[2][2] << std::endl;
    layers_ = 0;
    nlinear_ = 0;

    std::string method = prm_.get<std::string>("solver.method");

    if (method == "if_propagate_trimesh") {
        // const int trimeshlayers = 4;
        // const double init_scale = prm_.get<double>("config.axis_scale");
        double trires = prm_.get<double>("config.trires");
        const double edgelen = init_scale / trires; // 1;
        const double radius = init_scale;
        const double fac = std::sqrt(3) / 2;
        const std::array<double, 3> ax1 {axis_[0][0], axis_[0][1], axis_[0][2]};
        const std::array<double, 3> ax2 {0.5 * ax1[0] + fac * axis_[1][0],
                                         0.5 * ax1[1] + fac * axis_[1][1],
                                         0.5 * ax1[2] + fac * axis_[1][2]};

        std::cout << "Creating trimesh with radius: " << radius << ", edgelen: " << edgelen
                  << ", ax1: " << ax1[0] << "," << ax1[1] << "," << ax1[2] << ", ax2: " << ax2[0] << ","
                  << ax2[1] << "," << ax2[2] << std::endl;

        trimesh_ = std::make_unique<RegularTrimesh>(radius, // trimeshlayers,
                                                    std::array {origo_[0], origo_[1], origo_[2]},
                                                    ax1,
                                                    ax2,
                                                    std::array {edgelen, edgelen});

        trimesh_->removeSawtooths();

        // identify well cells (since this is not done in setFractureGrid when providing
        // a user-defined grid)
        // std::vector<CellRef> wellcells { {0, 0, 0}, {0, -1, 1}, {0, -1, 0}, {-1, -1, 1},
        //                                  {-1, 0, 0}, {-1, 0, 1} };
        well_source_cellref_ = RegularTrimesh::inner_ring_cells();

        for (const auto& cell : well_source_cellref_)
            well_source_.push_back(trimesh_->linearCellIndex(cell));

        auto [grid, fsmap, bmap] = trimesh_->createDuneGrid(1, well_source_cellref_);
        grid_mesh_map_ = fsmap;

        setFractureGrid(std::move(grid)); // create the physical grid from trimesh

    } else {
        setFractureGrid();
    }

    //setPerfPressure(0.0); // This can be changed by subsequently calling this
                          // function when the fracture is connected to the
                          // reservoir

    // NB: The Fracture object is still not fully initialized, since there is
    // not yet any connection with a surrounding reservoir. The following
    // needs to be done before the fracture can be used:
    // 1) Call `updateReservoirCells` to establish the mapping between frature grid
    //    cells and the cells in the reservoir grid (sets `reservoir_cells_`)
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
        // vtkmultiwriter_->gridChanged();// need to be called if grid is changed
        vtkmultiwriter_->gridViewChanged(grid_->leafGridView()); // need to be called if grid is changed
    }
}

void
Fracture::setupPressureSolver()
{
    OPM_TIMEFUNCTION();
    Opm::FlowLinearSolverParameters p;
    p.linsolver_ = prm_.get<std::string>("pressuresolver");
    prmpressure_ = Opm::setupPropertyTree(p, true, true);
    {
        std::size_t pressureIndex = 0; // Dummy
        const std::function<Vector()> weightsCalculator; // Dummy
        auto pressure_operator = std::make_unique<PressureOperatorType>(*pressure_matrix_);
        // using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
        pressure_operator_ = std::move(pressure_operator);
        auto psolver = std::make_unique<FlexibleSolverType>(
            *pressure_operator_, prmpressure_, weightsCalculator, pressureIndex);

        pressure_solver_ = std::move(psolver);
    }
}
/**
 * @brief Removes cells from the grid if they are out side reservoir.
 *
 * This function performs the following steps:
 * 1. Copies the current fracture width data to a persistent container.
 * 2. Iterates over all elements in the grid's leaf view and removes elements where the reservoir cell
 * index is negative.
 * 3. Grows the grid and performs post-growth operations.
 * 4. Resizes the fracture width array to match the new grid size.
 * 5. Copies the fracture width data back from the persistent container to the resized array.
 * 6. Resets the writers associated with the Fracture object.
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

    // const auto& indexSet = foamGridLeafView.indexSet();// for indices
    // const auto& indexSet = grid.localIdSet();// presitent numbering
    for (const auto& elem : elements(grid_->leafGridView())) {
        size_t eIdx = mapper.index(elem);
        if (reservoir_cells_[eIdx] < 0) {
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
Fracture::updateFilterCakeProps(const Opm::WellConnections& connections,
                                const Opm::SingleWellState<double,Fracture::IndexTraits>& wellstate,double dt)
{
    OPM_TIMEFUNCTION();
    // potentially remap filtercake thikness if grid has changed
    assert(filtercake_thikness_.size() == reservoir_cells_.size());
    const auto& perfdata = wellstate.perf_data;
    std::map<int, double> WI_fluxes;
    int water_index = 0;
    int np = 3;
    
    double qw_sum = 0.0;
    double qwpos_sum = 0.0;
    // should we do this based on fracture cross flow?        
    for (int i = 0; i < perfdata.cell_index.size(); ++i) {
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

    const auto& connection = connections.getFromGlobalIndex(wellinfo_.global_index); // probably wrong
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
                int res_cell = reservoir_cells_[eIdx];
                reservoir_areas[res_cell] += area;
                double flux = leakof_rate[eIdx]*area; // m3/s
                reservoir_flux[res_cell] += flux;
                if (flux < 0) {
                    if(prm_.get<bool>("verbose",false)){
                        std::cout << "Negative flux " << flux << " for element index " << eIdx
                                  << " with reservoir cell " << res_cell << std::endl;
                    }
                    flux = 0.0;
                }
                double water_rate = flux;
                double filtrate_volume = water_rate*filtrate_conn* dt; // m3
                double dh = filtrate_volume / (area * (1 - filtercake_poro_));
                assert(dh >= 0);
                filtercake_thikness_[eIdx] += dh;
            }
        }

        /// for setting or controling the fluxes between multiphase and single phase
        ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
        for (const auto& [res_cell, flux] : reservoir_flux) {
            // for (auto& element : Dune::elements(grid_->leafGridView())) {
            // int eIdx = mapper.index(element);
            // int res_cell = reservoir_cells_[eIdx];
            if (reservoir_areas.find(res_cell) == reservoir_areas.end()) {
                std::cout << "Reservoir area not found for element index " << res_cell << std::endl;
                continue;
            }
            double frac_flux = reservoir_flux[res_cell];
            if (res_cell != wellinfo_.well_cell) {
                if (std::abs(flux - WI_fluxes[res_cell]) > 0) {
                    std::cout << "Fracture flux differs from flow flux " << res_cell << std::endl;
                    std::cout << "Flux: frac " << frac_flux << " vs res " << WI_fluxes[res_cell]
                              << std::endl;
                }
            } else {
                std::cout << "Total WI flux " << WI_fluxes[res_cell] << " for cell " << res_cell
                          << " matches fracture flux " << frac_flux
                          << std::endl; //" for element index " << eIdx << std::endl;
            }
            // double res_area = reservoir_areas[res_area];
            //  To do use multiphase flux
            // double flux = reservoir_flux[res_area];
            // double dh =  flux/(res_area*(1-filtercake_poro_));
            // filtercake_thikness_[eIdx] += dh;
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
        = "Fracure_on_" + wellinfo_.name + "_perf_" + std::to_string(wellinfo_.perf) + "_nr";
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
    if (gptr == nullptr) {
        // create a new grid
        initFracture();
        int num_exp = prm_.get<int>("config.num_exp");
        if (num_exp > 0) {
            grow(num_exp, 0);
        }
        nlinear_ = layers_;
        int num_lin = prm_.get<int>("config.num_lin");
        if (num_lin > 0) {
            grow(num_lin, 1);
        }
        grid_->grow();
        grid_->postGrow();
    } else {
        // use grid provided by user
        grid_ = std::move(gptr);
        // @@ NB: we have not initialized the out_indices_ vector, which would be
        // necessary if we plan to grow the user-provided grid.
        // well_source_ is also not initialized.
    }
    // compute the cell normals
    updateCellNormals();

    // set object to allow stretching of the grid // @@ do not use with trimesh
    if (trimesh_ == nullptr)
        grid_stretcher_ = std::unique_ptr<GridStretcher>(new GridStretcher(*grid_));

    // set sparsity of pressure matrix, but do not compute its entries
    //initPressureMatrix(); need to be done after reservoir cells are known

    // Since the grid has been created/updated/changed, any previous mapping to
    // reservoir cells has been invalidated, and the fracture matrix (for
    // mechanics) is obsolete.
    reservoir_cells_.clear();
    fracture_matrix_ = nullptr;

    this->resetWriters();
}

void
Fracture::initFracture()
{
    Dune::GridFactory<Grid> factory; // Dune::FoamGrid<2,3>> factory;
    size_t N = 6;
    double radius = 1;
    std::vector<unsigned int> inner_indices;
    std::vector<Dune::FieldVector<double, 3>> vertices;
    {
        // Dune::FieldVector<double, 3> vertex({0, 0, 0});
        Point3D vertex = surfaceMap(0.0, 0.0);
        vertices.push_back(vertex);
    }

    for (size_t i = 0; i < N; ++i) {

        inner_indices.push_back(i + 1);
        double theta = (i * 2 * M_PI) / N;
        double x = radius * cos(theta);
        double y = radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        vertices.push_back(vertex);
        // assume the first 6 cells has source
        well_source_.push_back(i);
    }
    for (size_t i = 0; i < vertices.size(); i++) {
        factory.insertVertex(vertices[i]);
    }
    std::vector<std::vector<unsigned int>> cornerIDs;
    for (size_t i = 0; i < N; ++i) {
        unsigned int next = inner_indices[(i + 1) % N];
        std::vector<unsigned int> cornerID({unsigned(0), unsigned(inner_indices[i]), next});
        cornerIDs.push_back(cornerID);
    }
    for (size_t i = 0; i < N; ++i) {
        factory.insertElement(Dune::GeometryTypes::simplex(2), cornerIDs[i]); // std::move(mapping));
    }
    out_indices_ = inner_indices;
    grid_ = factory.createGrid();
    // grid_ = factory.createGrid();
}

std::vector<double>
Fracture::stressIntensityK1() const
{
    size_t nc = numFractureCells();
    std::vector<double> stressIntensityK1(nc, std::nan("0"));
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (const auto& elem : elements(grid_->leafGridView())) {
        // bool isboundary = false;
        for (auto& is : Dune::intersections(grid_->leafGridView(), elem)) {
            if (is.boundary()) {
                int nIdx = mapper.index(elem);
                auto isCenter = is.geometry().center();
                auto elCenter = elem.geometry().center();
                auto vecC = isCenter - elCenter;
                auto distC = vecC.two_norm();
                double K1 = ddm::fractureK1(distC, fracture_width_[nIdx], E_, nu_);
                stressIntensityK1[nIdx] = K1;
            }
        }
    }
    return stressIntensityK1;
}

void
Fracture::summary_of_solve(){
  if(numWellEquations()>0){
    int nc = numFractureCells();
    assert( fracture_pressure_.size() == nc + 1);
    double well_pressure = fracture_pressure_[nc];
    std::cout << "Perf pressure " << well_pressure << std::endl;
    const int cell = std::get<0>(perfinj_[0]); // @@ will this be the correct index?
    const double pres = reservoir_pressure_[cell];
    const double dp = well_pressure - pres;
    FractureProperties fprop = calculateFractureProperties();
    // set properties of fracture to structure used in output
    // frac_prop.height[perfIx] = fprop.height;
    // frac_prop.length[perfIx] = fprop.length;
    // frac_prop.area[perfIx] = fprop.area;
    // frac_prop.flux[perfIx] = fprop.flux;
    // frac_prop.WI[perfIx] = fprop.WI;
    // frac_prop.volume[perfIx] = fprop.volume;
    // frac_prop.filter_volume[perfIx] = fprop.filter_volume;
    // frac_prop.avg_width[perfIx] = fprop.avg_width;
    // frac_prop.avg_filter_width[perfIx] = fprop.avg_filter_width;
    // frac_prop.inj_pressure[perfIx] = fprop.inj_pressure;
    // frac_prop.inj_bhp[perfIx] = fprop.inj_bhp;
    // frac_prop.inj_wellrate[perfIx] = fprop.inj_wellrate;
    // fracture rate
    std::cout << "Reservoir pressure at perf cell " << cell << " is " << pres << " dp " <<   dp << " frac flux " << fprop.flux << std::endl;
    std::cout << "Fracture properties: " << std::endl;
    std::cout << " Height: " << fprop.height << " Length: " << fprop.length << " Area: " << fprop.area << std::endl;
    std::cout << " Volume: " << fprop.volume << " Avg width: " << fprop.avg_width << std::endl;
    std::cout << " WI: " << fprop.WI << " Inj pressure: " << fprop.inj_pressure << " Inj bhp: " << fprop.inj_bhp << " Inj rate: " << fprop.inj_wellrate << std::endl;       
  }
}

  
void
Fracture::write(int reportStep) const
{
    OPM_TIMEFUNCTION();
    std::vector<double> K1; // need to be in scope until written
    auto tmp(fracture_pressure_); //@@
    tmp.resize(numFractureCells());
    if (reservoir_cells_.size() > 0) {
        vtkwriter_->addCellData(reservoir_cells_, "ReservoirCell");
    }
    if (reservoir_perm_.size() > 0) {
        vtkwriter_->addCellData(reservoir_perm_, "ReservoirPerm");
    }
    if (reservoir_cstress_.size() > 0) {
        vtkwriter_->addCellData(reservoir_cstress_, "ReservoirCStress");
    }
    if (reservoir_dist_.size() > 0) {
        vtkwriter_->addCellData(reservoir_dist_, "ReservoirDist");
    }
    if (fracture_pressure_.size() > 0) {
        // since fracture_pressure_ may have an additional entry for well
        // pressure, we need to ensure we pass the correct size (@@ ideally this
        // should be fixed by moving the well pressure value into perf_pressure_
        // once it has been computed, and shrink fracture_pressure_ accordingly.

        // for (size_t i = 0; i != numFractureCells(); ++i) tmp[i] = fracture_pressure_[i];
        vtkwriter_->addCellData(tmp, "FracturePressure");
        // vtkwriter_->addCellData(fracture_pressure_, "FracturePressure");
    }
    if (reservoir_pressure_.size() > 0) {
        vtkwriter_->addCellData(reservoir_pressure_, "ReservoirPressure");
    }
    if (reservoir_mobility_.size() > 0) {
        vtkwriter_->addCellData(reservoir_mobility_, "ReservoirMobility");
    }
    if (fracture_width_.size() > 0) {
        vtkwriter_->addCellData(fracture_width_, "FractureWidth");
        K1 = this->stressIntensityK1();
        vtkwriter_->addCellData(K1, "stressIntensityK1");
    }
    // std::stringstream ss(this->name());
    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename");
    std::string filename = outputdir + "/Static_" + simName + this->name();
    if (reportStep > 0) {
        filename = filename + "_step_" + std::to_string(reportStep);
    }
    vtkwriter_->write(filename.c_str());
};

double Fracture::reservoirTraction(int i) const{  
   return ddm::tractionSymTensor(reservoir_stress_[i], cell_normals_[i]);
}
double Fracture::fractureForce(int i) const{
  return reservoirTraction(i) - fracture_pressure_[i];
}

double Fracture::leakofDp(int i) const{
   return fracture_pressure_[i][0] - reservoir_pressure_[i]
                - (fracture_dgh_[i] - reservoir_cell_z_[i] * gravity_ * reservoir_density_[i]);
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
    // vtkmultiwriter_->gridChanged();// need to be called if grid is changed
    //  need to have copies in case of async outout (and interface to functions)
    
    std::vector<double> K1 = this->stressIntensityK1();
    std::vector<double> fracture_pressure(numFractureCells(), 0.0);
    std::vector<double> leakof_dp(numFractureCells(), 0.0);
    std::vector<double> filtercake_thikness = filtercake_thikness_; //.size(),0.0);
    std::vector<double> reservoir_pressure = reservoir_pressure_; //.size(),0.0);
    std::vector<double> reservoir_dist = reservoir_dist_; //.size(),0.0);
    std::vector<double> reservoir_perm = reservoir_perm_; //.size(),0.0);
    std::vector<double> reservoir_cstress = reservoir_cstress_; //.size(),0.0);
    std::vector<double> reservoir_mobility = reservoir_mobility_; //.size(),0.0);
    std::vector<double> reservoir_traction(reservoir_stress_.size(), 0);
    std::vector<double> fracture_force(reservoir_stress_.size(), 0);
    std::vector<double> reservoir_cells(reservoir_cells_.size(), 0.0);
    std::vector<double> fracture_width(fracture_width_.size(), 0.0);
    std::vector<double> fracture_head(fracture_width_.size(), 0.0);
    std::vector<double> flow_width(fracture_width_.size(), 0.0);
    std::vector<double> rhs_width(rhs_width_.size(), 0.0);
    std::vector<double> well_index(numFractureCells(), 0.0);
    std::map<int, double> wellIndMap;
    std::vector<double> leakofrate = leakOfRate();
    // assert(leakofrate.size() == fracture_pressure_.size());
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
            reservoir_cells[i] = reservoir_cells_[i]; // only converts to double
            reservoir_traction[i] = reservoirTraction(i);//ddm::tractionSymTensor(reservoir_stress_[i], cell_normals_[i]);
            fracture_force[i] = fractureForce(i);//reservoir_traction[i] - fracture_pressure[i];
            leakof_dp[i] = leakofDp(i);//fracture_pressure[i] - reservoir_pressure[i]
                //- (fracture_dgh_[i] - reservoir_cell_z_[i] * gravity_ * reservoir_density_[i]);
            // make relative to origo
            fracture_head[i] = (fracture_pressure[i] - fracture_dgh_[i])
                - (injectionPressure() - (perf_ref_depth_ * gravity_ * density_perf_));
            // make it relative to origo
            assert(filtercake_thikness[i] >= 0.0);
        }
    if (!first) {    
        for (size_t i = 0; i < numFractureCells(); ++i) {
                // maybe slow
                well_index[i] = wellIndMap[reservoir_cells_[i]];
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
Fracture::grow(int layers, int method)
{
    while (layers_ < layers) {
        std::vector<unsigned int> inner_indices = out_indices_;
        out_indices_.resize(0);
        if (method == 0) {
            this->insertExp(inner_indices);
        } else {
            this->insertLinear(inner_indices);
            nlinear_ += 1;
        }
        ++layers_;
        inner_indices = out_indices_;
    }
}
void
Fracture::insertLinear(const std::vector<unsigned int>& inner_indices)
{
    double radius = 1.0;
    size_t N = inner_indices.size();
    for (size_t i = 0; i < N; ++i) {
        double theta = (i * 2 * M_PI) / (N);
        theta += (layers_ - nlinear_) * 0.5 * (2 * M_PI / N);
        double out_radius = radius + layers_ + 1;
        double x = out_radius * cos(theta);
        double y = out_radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        int new_ind = grid_->insertVertex(vertex);
        out_indices_.push_back(new_ind);
    }
    for (size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i % N], out_indices_[(i) % (N)], inner_indices[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }

        {
            std::vector<unsigned int> cornerID(
                {inner_indices[(i + 1) % N], out_indices_[(i) % (N)], out_indices_[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

void
Fracture::insertExp(const std::vector<unsigned int>& inner_indices)
{
    size_t N = inner_indices.size();
    double radius = 1.0;

    for (size_t i = 0; i < N * 2; ++i) {
        double theta = (i * 2 * M_PI) / (N * 2);
        double out_radius = radius + layers_ + 1;
        double x = out_radius * cos(theta);
        double y = out_radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        int new_ind = grid_->insertVertex(vertex);
        out_indices_.push_back(new_ind);
    }


    for (size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID({inner_indices[i],
                                                out_indices_[(2 * i) % (N * 2)],
                                                out_indices_[(2 * i + 1) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i], out_indices_[(2 * i + 1) % (N * 2)], inner_indices[(i + 1) % N]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
        {
            std::vector<unsigned int> cornerID({inner_indices[(i + 1) % N],
                                                out_indices_[(2 * i + 1) % (N * 2)],
                                                out_indices_[(2 * i + 2) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

Dune::FieldVector<double, 3>
Fracture::surfaceMap(double x, double y)
{
    Point3D vec(0.0);
    vec += x * axis_[0];
    vec += y * axis_[1];
    // for (int i = 0; i < 2; ++i) {
    // }
    vec += origo_;
    return vec;
}

void
Fracture::updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                               const Dune::CpGrid& grid,
                               const std::vector<Fracture::EntitySeed>& entity_seeds)
{
    OPM_TIMEFUNCTION();
    reservoir_cells_.resize(numFractureCells());
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    int tri_divide = 0;
    int tri_outside = 0;
    all_reservoir_areas_.resize(numFractureCells());
    all_reservoir_cells_.resize(numFractureCells());
    all_reservoir_centers_.resize(numFractureCells());
    bool full_intersections = prm_.get<bool>("solver.full_intersections", true);
    for (auto& element : elements(grid_->leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        auto geom = element.geometry();
        //external::cvf::BoundingBox bb;
        using Vec3d = external::cvf::Vec3d;
        auto vertex = geom.center();
        Vec3d point(vertex[0], vertex[1], vertex[2]);
        //bb.add(point);
        reservoir_cells_[elemIdx]  = external::cellOfPoint(cellSearchTree, grid, entity_seeds, point);
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
            all_reservoir_areas_[elemIdx] = areas;
            all_reservoir_cells_[elemIdx] = rcells;
            all_reservoir_centers_[elemIdx] = centers;
        }
    }
    auto it = std::find(reservoir_cells_.begin(), reservoir_cells_.end(), -1);
    auto extended_fractures = prm_.get<bool>("extended_fractures");
    if (!(it == reservoir_cells_.end()) && !(extended_fractures)) {
        std::cout << "Remove fracture outside of model" << std::endl;
        // remove fracture outside of model
        this->removeCells();
        this->updateReservoirCells(cellSearchTree, grid, entity_seeds);
    }else{
      if(!(it == reservoir_cells_.end())){
          std::cout << "Warning:: Fracture outside of model" << std::endl;
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
    // assert(reservoir_cells_.size() == nc);
    reservoir_perm_.resize(nc, perm);
    reservoir_dist_.resize(nc, dist);
    reservoir_mobility_.resize(nc, 1000);
    reservoir_pressure_.resize(nc, 100.0e5);
    reservoir_stress_.resize(nc);
    reservoir_cstress_.resize(nc, cstress);

    for (size_t i = 0; i != nc; ++i)
        reservoir_stress_[i] = Dune::FieldVector<double, 6> {0, 0, 0, 0, 0, 0};

    nu_ = 0.25;
    E_ = 1e9;
    this->initFractureWidth();
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
        // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
        // Dune::mcmgElementLayout>::index(element);
        int i = grid_->leafGridView().indexSet().index(element);
        double z = element.geometry().center()[2];
        fracture_dgh_[i] = gravity_ * reservoir_density_[i] * z;
    }
    // could be put into the assemble loop
    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double t1 = std::get<2>(matel);
        double t2 = std::get<3>(matel);
        //double h1 = fracture_width_[i] + min_width_;
        //double h2 = fracture_width_[j] + min_width_;
        const double h1 = std::max(fracture_width_[i][0],min_width_);
        const double h2 = std::max(fracture_width_[j][0],min_width_);
        //const double h1 = min_width_;
        //const double h2 = min_width_;
        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);

        const double mobility = 0.5 * (reservoir_mobility_[i] + reservoir_mobility_[j]);//NB shoul chane to mobility_water_perf_;
        value = 1 / value;
        value *= mobility;
        double dh = (fracture_dgh_[i] - fracture_dgh_[j])*(-1.0);//NB -1.0??
        rhs_pressure_[i] -= value * dh;
        rhs_pressure_[j] += value * dh;
    }

    for (size_t i = 0; i < reservoir_pressure_.size(); ++i) {
        rhs_pressure_[i] += leakof_[i] * reservoir_pressure_[i];// gravity i handeld in next loop
    }
    // gravity contribution from fracture to reservoir
    for (auto element : Dune::elements(grid_->leafGridView())) {
        // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
        // Dune::mcmgElementLayout>::index(element);
        int i = grid_->leafGridView().indexSet().index(element);
        double z = element.geometry().center()[2];
        // double gravity = 9.81;
        rhs_pressure_[i] += leakof_[i] * (z - reservoir_cell_z_[i]) * gravity_ * reservoir_density_[i];
    }
    // gravity contributions between fracture cells

    auto control = prm_.get_child("control");
    std::string control_type = control.get<std::string>("type");
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
            double dh_perf = origo_[2] * gravity_ * reservoir_density_[cell];
            double dh_cell = fracture_dgh_[cell];
            ;
            double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (pressure - dh);
        }
    } else if (control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            //double dh_perf = origo_[2] * gravity_ * reservoir_density_[cell];
            double dh_perf = perf_ref_depth_ * gravity_ * reservoir_density_[cell];
            double dh_cell = fracture_dgh_[cell];
            ;
            double dh = dh_perf - dh_cell;
            rhs_pressure_[cell] += value * (perf_pressure_ - dh);
        }
    } else if (control_type == "rate_well") {
      // this tries to make control equation for single phase standard wells with one phase and static reservoir
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        const double well_rate = well_rate_;//control.get<double>("rate") / 24 / 60 / 60; // convert to m3/sec
        const double WI = total_wellindex_; //eclude lambdacontrol.get<double>("WI");
        const int cell = std::get<0>(perfinj_[0]); // @@ will this be the correct index?
        const double pres = reservoir_pressure_[cell];
        const double density = reservoir_density_[cell];
        //const double lambda = reservoir_mobility_[cell]; // @@ only correct if mobility is constant!
        const double lambda = 1.0;
        if (false) {
            rhs_pressure_[rhs_pressure_.size() - 1] = well_rate + WI * lambda * pres; // well source term
        } else {
            double wi_dzfac = wi_respress_; 
            wi_dzfac -= wi_dz_*density*gravity_; // sum wi dz
            //wi_dzfac += wi_respress_; // sum wi res_press
            wi_dzfac += (WI* (perf_ref_depth_- well_ref_depth_))*gravity_*density; // counvert to perf_ref_depth from bhp ref_depth for primary
                                            // variable reservoir_pressure_[nc] i.e. perf pressure
            rhs_pressure_[rhs_pressure_.size() - 1] = well_rate + lambda * wi_dzfac;
        }
        if (false){
            for (const auto& perfinj : perfinj_) {
                int cell = std::get<0>(perfinj);
                double value = std::get<1>(perfinj)*mobility_water_perf_;
                assert(origo_[2] == perf_ref_depth_);
                double dh_perf = perf_ref_depth_ * gravity_ * reservoir_density_[cell];
                double dh_cell = fracture_dgh_[cell];
                ;
                double dh = dh_perf - dh_cell;
                rhs_pressure_[cell] += value * (-dh);
            }
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
}

double
Fracture::injectionBhp() const
{

  //if(!(perfinj_.size()>0)){
     // no solves done
  //   return 0.0;
 // }
  double dz = perf_ref_depth_ - well_ref_depth_;;
  //double density = reservoir_density_[std::get<0>(perfinj_[0])];// hack may have small inconsistency
  double density = density_perf_;
  return injectionPressure() - dz*gravity_*density;
}

double
Fracture::injectionPressure() const
{
    std::string control_type = prm_.get<std::string>("control.type");
    if (control_type == "rate") {
        double bhp = 0.0;
        double scale = well_source_.size();
        // could have corrected for WI
        for (auto cell : well_source_) {
            bhp += fracture_pressure_[cell] / scale;
        }
        return bhp;
        // } else if (control_type == "bhp") {
        //     double bhp = prm_.get<double>("control.bhp");
    } else if (control_type == "pressure") {
        double bhp = prm_.get<double>("control.pressure");
        return bhp;
    } else if (control_type == "perf_pressure") {
        double bhp = perf_pressure_;
        return bhp;
    } else if (control_type == "rate_well") {
        // @@ We should use perf_pressure_ here too, but ensure it is updated
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
        double dh_res = reservoir_cell_z_[eIdx] * gravity_ * reservoir_density_[eIdx];
        double dh_frac = fracture_dgh_[eIdx];
        double dp = ((fracture_pressure_[eIdx]-dh_frac) - (reservoir_pressure_[eIdx]-dh_res));
        double q = leakof_[eIdx] * dp;
        leakofrate[eIdx] = q / geom.volume();
    }
    return leakofrate;
}

std::vector<RuntimePerforation>
Fracture::wellIndices() const{
   std::vector<RuntimePerforation> wellindices;
   if(well_indices_.size() == 2){
         wellindices = wellIndicesAvrg(well_indices_);
   }
   for(int i=0; i < wellindices.size(); ++i){
     wellindices[i].ref_ctf = wellindices[i].ctf;
     wellindices[i].ref_pressure = wellindices[i].pressure;
   }
   return wellindices;
}

void Fracture::setPerfProps(double perfpressure, double depth, double perfrate){
    double damping_factor_perf = prm_.get<double>("solver.damping_factor_perf",2.0);
    if(well_indices_.size() >0){
        perf_pressure_ = perfpressure;
    }else{
        // dampted updating for stability of fracture growth
        perf_pressure_ = (perf_pressure_*damping_factor_perf+ perfpressure)/(1.0+damping_factor_perf);//(perfpressure-perf_pressure_)*damping_factor;
    }
    well_perf_rate_ = perfrate;
    perf_ref_depth_ = depth;
}

void Fracture::setWellProps(double wellrate,
                            double totalwi,
                            double wi_dz,
                            double wi_respress,
                            double ref_depth){
    //if(well_indices_.size() >0){
        well_rate_ = wellrate;
        total_wellindex_ = totalwi;
        wi_dz_ = wi_dz;
        wi_respress_ = wi_respress;
        well_ref_depth_ = ref_depth;
    // }else{
    //     // dampted updating for stability of fracture growth
    //     well_rate_ = (well_rate_ + damping_factor_wi*wellrate)/(1.0+damping_factor_wi);//(wellrate-well_rate_)*damping_factor;
    //     total_wellindex_ = (total_wellindex_ + damping_factor_wi*totalwi)/(1.0+damping_factor_wi);//(totalwi-total_wellindex_)*damping_factor;
    // }
}

std::vector<RuntimePerforation>
Fracture::wellIndicesAvrg(const std::vector<std::vector<RuntimePerforation>>& well_indices) const{
    // average int time between to well indices.
    OPM_TIMEFUNCTION();
  // find union of cell indices
  assert(well_indices.size() == 2);// || well_indices_.size() == 0);
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
  //std::vector<RuntimePerforation> wellindices(cells.size(),0);
  int tind=0;
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
      wellindices[ind].ref_ctf = 0.0;
      //wellindices[ind].pressure = 0.0;
    }
  }
  for(int tind=0; tind < well_indices.size(); ++tind){ 
    for (const auto& wind : well_indices[tind]){
      int ind = cell_pind[wind.cell];
      wellindices[ind].ctf += weight_ctf[tind]*wind.ctf;
      //wellindices[ind].pressure *= damping_factor*wind.pressure;
    }
    total_weight += weight_ctf[tind];
  }

  for(int i=0; i < wellindices.size(); ++i){
        wellindices[i].ctf /= total_weight;
        //wellindices[i].pressure /= total_weight;
  }
  return wellindices;
}

std::vector<RuntimePerforation>
Fracture::wellIndices_() const
{
    OPM_TIMEFUNCTION();
    // find unique reservoir cells
    if (leakof_.size() == 0) {
        // if pressure is not assembled return empty
        return {};
    }
    bool only_to_one = !(prm_.get<bool>("solver.divide_wellidx",false));
    std::vector<int> res_cells = reservoir_cells_;
    std::vector<double> tot_areas(all_reservoir_cells_.size(), 0.0);
    if(only_to_one){
        res_cells = reservoir_cells_;
    }else{
        for(size_t ind=0; ind < all_reservoir_cells_.size(); ++ind){
            double total_area = 0.0;
            auto cells = all_reservoir_cells_[ind];
            for(size_t lind=0; lind < cells.size(); ++lind){
                res_cells.push_back(cells[lind]);
                total_area += all_reservoir_areas_[ind][lind];
            }
            // this area should be the same as areal of cell if it is inside full model
            tot_areas[ind] = total_area;
        }
    }

    
    std::sort(res_cells.begin(), res_cells.end());
    auto last = std::unique(res_cells.begin(), res_cells.end());
    res_cells.erase(last, res_cells.end());
    std::vector<double> q_cells(res_cells.size(), 0.0);
    std::vector<double> p_cells(res_cells.size(), 0.0);
    std::vector<double> mob_cells(res_cells.size(), 0.0);
    std::vector<double> dens_cells(res_cells.size(), 0.0);
    std::vector<double> z_cells(res_cells.size(), 0.0);
    std::vector<double> traction(res_cells.size(), 0.0);
    std::vector<double> area(res_cells.size(), 0.0);
    std::vector<double> leakofrate = this->leakOfRate();
    double q_prev = 0;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        double q = leakofrate[eIdx]*geom.volume(); // just to avoid warning
        
        std::vector<int> local_res_cells;
        if(only_to_one){
            local_res_cells.push_back(reservoir_cells_[eIdx]);
        }else{
            local_res_cells = all_reservoir_cells_[eIdx];
        }    
        // just to search
        for(int lind=0; lind < local_res_cells.size(); ++lind){
            int res_cell = local_res_cells[lind];
            double area_frac = 1.0;
            if(!only_to_one){
                area_frac = all_reservoir_areas_[eIdx][lind]/tot_areas[eIdx];
            }
            auto it = std::find(res_cells.begin(), res_cells.end(), res_cell);
            int ind_wellIdx = it - res_cells.begin();
            assert(!(it == res_cells.end()));
            if (q_prev * q < 0) {
                OPM_THROW(std::runtime_error, "Cross flow in fracture ??");
            }
            double loc_area = area_frac*element.geometry().volume();
            area[ind_wellIdx] += loc_area;
            q_cells[ind_wellIdx] += area_frac*q;
            if(reservoir_cells_[eIdx] == res_cell){
                traction[ind_wellIdx] += ddm::tractionSymTensor(reservoir_stress_[eIdx], cell_normals_[eIdx])*loc_area; //normalFractureTraction(elIx);
                // assume all of this is the same so no need for interpolation
                p_cells[ind_wellIdx] = reservoir_pressure_[eIdx]; // is set multiple times
                mob_cells[ind_wellIdx] = reservoir_mobility_[eIdx]; // is set multiple times
                dens_cells[ind_wellIdx] = reservoir_density_[eIdx]; // is set multiple times
                z_cells[ind_wellIdx] = reservoir_cell_z_[eIdx]; // is set multiple times
            }else{
                // logic could be changed after introduction of map
              auto getFromMap = [](const std::map<int,double> map, int key){
                if (auto search = map.find(key); search != map.end()){
                  return search->second;
                }else{
                  assert(false);
                  return 0.0;
                }
              };

              p_cells[ind_wellIdx] = getFromMap(map_reservoir_pressure_,res_cell); // is set multiple times
              mob_cells[ind_wellIdx] = getFromMap(map_reservoir_mobility_,res_cell); // is set multiple time
              dens_cells[ind_wellIdx] = getFromMap(map_reservoir_density_,res_cell); // is set multiple time
              z_cells[ind_wellIdx] = getFromMap(map_reservoir_cell_z_,res_cell); // 
                // need to collect this values if no cell has centroid nearest to this
                // reservoir cell
            }
        }
    }
    std::vector<RuntimePerforation> wellIndices(res_cells.size());
    double inj_press = injectionPressure();
    bool cells_outside = false;
    for (size_t i = 0; i < res_cells.size(); ++i) {
        auto& perf = wellIndices[i];
        perf.cell = res_cells[i];
        if(res_cells[i]<0){
            // NB will be removed make dummy values
            cells_outside = true;
            perf.depth = this->origo_[2];
            perf.ctf = 0;
            perf.ref_ctf = 0;
            perf.pressure = inj_press;
            perf.ref_pressure = inj_press-10e5;// dummy 
            perf.segment = this->wellinfo_.segment;
            perf.perf_range = this->wellinfo_.perf_range;
            continue;
        }
        double dh_res = z_cells[i] * gravity_ * dens_cells[i];
        double perf_density = dens_cells[i];
        double dh_perf = gravity_ * perf_density * origo_[2];
        double WI = 0.0;
        if(mob_cells[i] <= 0.0){
            std::cout << "Warning zero mobility for perf cell: " << res_cells[i] << std::endl;
            WI = 0.0;
        }else{
            WI = q_cells[i] / ((inj_press - dh_perf) - (p_cells[i] - dh_res));//NB d_perf def
        }
        if (WI < 0.0) {
          // keep but could maybe be removed
            std::cout << "Negative WI: " << WI << " for cell: " << res_cells[i] << std::endl;
            WI = 0.0;
        }
        perf.ctf = WI/mob_cells[i]; // convert to ctf
        {
                // NBsould probably be removed
            perf.depth = this->origo_[2];
            perf.segment = this->wellinfo_.segment;
            perf.perf_range = this->wellinfo_.perf_range;
            perf.pressure = inj_press;
            double tmp_press =  traction[i]/area[i];
            //if(tmp_press < inj_press){
            //    perf.ref_pressure = tmp_press;        
            //    perf.ref_ctf = 0.0; // should be the closed value
            //}else{
            //perf.ref_pressure = tmp_press;//inj_press- 10e5; // dummy value
            //perf.ref_ctf = 0.0;//perf.ctf; 
            perf.ref_pressure = inj_press- 10e5; // dummy value
            perf.ref_ctf = perf.ctf; 
        }
    }
    // remove connections to outside of reservoir
    if(cells_outside){
        std::cout << "Cells of fracture outside removed from perforation index" << std::endl;
        wellIndices.erase(
                          std::remove_if(wellIndices.begin(),
                                         wellIndices.end(),
                          [](const auto& perf){ return perf.cell < 0;}),
                          wellIndices.end());
    }

    return wellIndices;
}

template <typename Scalar>
void
Fracture::assignGeomechWellState(PerfData<Scalar>& perfData) const
{
    auto perfPos = std::find(
    perfData.cell_index.begin(), perfData.cell_index.end(), this->wellInfo().well_cell);
    const auto perfIx = std::distance(perfData.cell_index.begin(), perfPos);
    ConnFracStatistics<Scalar>& stats = perfData.connFracStatistics[perfIx];

    using Quantity = typename ConnFracStatistics<Scalar>::Quantity;

    constexpr auto pressIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Pressure);
    constexpr auto rateIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::FlowRate);
    constexpr auto widthIx = static_cast<std::underlying_type_t<Quantity>>(Quantity::Width);

    const auto nCells = this->reservoir_cells_.size();

    stats.reset();

    for (auto cellIx = 0 * nCells; cellIx < nCells; ++cellIx) {
        auto samplePoint = typename ConnFracStatistics<Scalar>::SamplePoint {};

        samplePoint[pressIx] = this->fracture_pressure_[cellIx][0];
        double dh_res = reservoir_cell_z_[cellIx] * gravity_ * reservoir_density_[cellIx];
        double dh_frac = fracture_dgh_[cellIx];
        samplePoint[rateIx] = this->leakof_[cellIx]
          * ((this->fracture_pressure_[cellIx][0]-dh_frac)- (this->reservoir_pressure_[cellIx]-dh_res));

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
        std::cout << "Fracture area doubled, stopping expansion" << std::endl;
        max_flow_time_step_  = dt_lim;
        return true;
    } else {
      max_flow_time_step_  = std::numeric_limits<double>::max();// no limit
        return false;
    }
    return false;
}

void
Fracture::writePressureSystem() const
{
    if (prm_.get<bool>("write_pressure_system")) {
        Dune::storeMatrixMarket(*pressure_matrix_, "pressure_matrix");
        Dune::storeMatrixMarket(rhs_pressure_, "pressure_rhs");
    }
}

void
Fracture::writeFractureSystem() const
{
    if (prm_.get<bool>("write_fracture_system")) {
        // Dune::storeMatrixMarket(*fracture_matrix_, "fracture_matrix");
        Dune::storeMatrixMarket(rhs_width_, "rhs_width");
    }
}

void
Fracture::solvePressure()
{
    OPM_TIMEFUNCTION();
    size_t nc = numFractureCells();
    assert(numWellEquations() == 0); // @@ not implemented/tested for "rate_well" control
    fracture_pressure_.resize(nc);
    fracture_pressure_ = 1e5;
    assert(pressure_matrix_); // should always be constructed at this pointn
    // if(!pressure_matrix_){
    //     this->initPressureMatrix();
    // }
    this->assemblePressure();
    this->addSource(); // probably include reservoir pressure
    this->writePressureSystem();
    try {
        if (!pressure_solver_) {
            this->setupPressureSolver();
        }
        Dune::InverseOperatorResult r;
        fracture_pressure_.resize(rhs_pressure_.size());
        fracture_pressure_ = 0;
        pressure_solver_->apply(fracture_pressure_, rhs_pressure_, r);
    } catch (Dune::ISTLError& e) {
        std::cerr << "exception thrown " << e << std::endl;
    }
}

// based on the (given) values of fracture_pressure__, compute rhs_width_ and
// fracture_width_
void
Fracture::solveFractureWidth()
{
    fractureMatrix().solve(fracture_width_, rhs_width_);
    double max_width = prm_.get<double>("solver.max_width");
    double min_width = prm_.get<double>("solver.min_width");
    for (int i = 0; i < fracture_width_.size(); ++i) {
        assert(std::isfinite(fracture_width_[i]));
        if (fracture_width_[i] > max_width) {
            std::cout << "Limit Fracture width" << std::endl;
            fracture_width_[i] = max_width;
        }
        if (fracture_width_[i] < min_width) {
            std::cout << "Remove small Fracture width" << std::endl;
            fracture_width_[i] = min_width;
        }
        assert(std::isfinite(fracture_width_[i]));
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
    // ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    // for (auto& element : elements(grid_->leafGridView())) {
    //     const auto elemIdx = elemMapper.index(element);
    //     auto geom = element.geometry();
    //     auto vec_origo = geom.center() - origo_;
    //     double dist_origo = vec_origo.two_norm();
    //     fracture_width_[elemIdx] *= dist_origo;
    // }
}
void
Fracture::initFracturePressureFromReservoir()
{
    size_t nc = reservoir_cells_.size();
    fracture_pressure_.resize(nc + numWellEquations());
    fracture_pressure_ = 0;
    for (size_t i = 0; i < nc; ++i) {
        fracture_pressure_[i] = reservoir_pressure_[i];
    }
}




void
Fracture::updateLeakoff()
{

    // find boundary nodes.
    //Dune::Codim<> {}
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
              //auto geom = is.geometry();
               auto lf = is.indexInInside();
               auto face = element.subEntity<codimFace>(lf);
               //for (const auto& vertex : Dune::subEntities(face, codimVertex)){
               //}
               int nv = face.subEntities(codimVertex);
               assert(nv == 2);
               const auto refElement = Dune::ReferenceElements<double,dim>::general(element.type());
               for (int v = 0; v < nv; v++) {
                 auto&& vertex = element.template subEntity<dim>(refElement.subEntity(lf, dim-1, v, dim));
                 //auto vertex = face.subEntity(codimVertex,v);
                   //boundaryVerts.insert(indexSet.index(vertex));
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
        // for (const auto& is : intersections(gv, element)) {
        //     if (is.boundary()) {
        //       elementHasBoundaryNode[eidx] = true;
        //     }
        // }
        if(true){
        int nv = element.subEntities(codimVertex);
        //for (const auto& f : Dune::subEntities(element, codimVertex))
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
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto geom = element.geometry();
        double area = geom.volume();
        double res_mob = reservoir_mobility_[eIdx];
        leakof_[eIdx] = res_mob * reservoir_perm_[eIdx] * area / reservoir_dist_[eIdx];
        if (has_filtercake_) {
            double invtrans = 1 / leakof_[eIdx];
            assert(filtercake_thikness_[eIdx] >= 0.0);
            if (filtercake_thikness_[eIdx] > 0.0) {
                double fitercaketrans = res_mob * filtercake_perm_ * area / filtercake_thikness_[eIdx];
                invtrans += 1 / fitercaketrans;
                // assert(filtercake_perm_ > 0.0);
                // assert(filtercake_thikness_[eIdx] > 0.0);
            }
            leakof_[eIdx] = 1 / invtrans;
            if(elementHasBoundaryNode[eIdx] && no_leakof_outercells){
              leakof_[eIdx] = 0.0;
            }
        }
    }
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
        volume += geom.volume() * (fracture_width_[eIdx]);// + min_width_);
        filter_volume += geom.volume() * filtercake_thikness_[eIdx];//
        auto dist = geom.center() - origo_;
        double dh = dist.dot(naxis_[0]);
        double dw = dist.dot(naxis_[1]);
        updateMaxMinVector(dhvec, dh);
        updateMaxMinVector(dwvec, dw);
        if(leak_of_rate.size() != numFractureCells()){
            //std::cout << "Error in leak_of_rate size" << std::endl;
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
    //double WI = std::accumulate(WIs.begin(), WIs.end(), 0)
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
    // size_t num_columns = 0;
    //  index of the neighbour    //
    // Flow from wells to fracture cells
    double fWI = prm_.get<double>("fractureWI");
    perfinj_.clear();
    htrans_.clear();
    for (int cell : well_source_) {
        perfinj_.push_back({cell, fWI});
    }
    // flow between fracture cells
    const size_t nc = numFractureCells() + numWellEquations();
    // leakof_.resize(nc,0.0);
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
    // std::sort(transes.begin(),transes.end(),sortcsr);
    //  build matrix
    // if (pressure_matrix_.bu == 0){
    // size_t nc = numFractureCells();
    pressure_matrix_ = std::make_unique<Matrix>(nc, nc, 4, 0.4, Matrix::implicit);
    auto& matrix = *pressure_matrix_;
    // matrix.setBuildMode(Matrix::implicit);
    //  map from dof=3*nodes at a cell (ca 3*3*3) to cell
    // matrix.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
    // size_t nc = numFractureCells();
    // matrix.setSize(nc, nc);
    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double zero_entry = 0.0; // 1e-11;
        matrix.entry(i, j) = zero_entry; // 0;
        matrix.entry(j, i) = zero_entry; // 0;
        matrix.entry(j, j) = zero_entry; // 0;
        matrix.entry(i, i) = zero_entry; // 0;
    }
    // add in case a cell have only no internal boundaries
    for(int i=0; i < nc; ++i){
        matrix.entry(i, i) = 0.0; 
    }

    if (numWellEquations() > 0) {
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        // add elements for last row and column
        const int weqix = nc - 1; // index of well equation
        //for (int i : well_source_) {
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
    updateLeakoff();

    auto& matrix = *pressure_matrix_;
    matrix = 0.0;
    // double mobility=1e4; //1e4; // @@ 1.0
    //  get head in all fracture cells
    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double t1 = std::get<2>(matel);
        double t2 = std::get<3>(matel);
        //double h1 = fracture_width_[i] + min_width_;
        //double h2 = fracture_width_[j] + min_width_;
        const double h1 = std::max(fracture_width_[i][0],min_width_);
        const double h2 = std::max(fracture_width_[j][0],min_width_);
        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);

        const double mobility = 0.5 * (reservoir_mobility_[i] + reservoir_mobility_[j]);
        value = 1 / value;
        value *= mobility;
        // matrix.entry(i, j) -= value;
        // matrix.entry(j, i) -= value;
        // //
        // matrix.entry(i, i) += value;
        // matrix.entry(j, j) += value;

        matrix[i][j] -= value;
        matrix[j][i] -= value;
        //
        matrix[i][i] += value;
        matrix[j][j] += value;
        // gravity contributions
    }
    auto control = prm_.get_child("control");
    std::string control_type = control.get<std::string>("type");
    for (size_t i = 0; i < leakof_.size(); ++i) {
        // matrix.entry(i, i) += leakof_[i];
        matrix[i][i] += leakof_[i];
    }
    if (control_type == "rate") {
        // no extra tings in matrix
    } else if (control_type == "pressure" || control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            matrix[cell][cell] += value;
        }
    } else if (control_type == "rate_well") {
        assert(numWellEquations() == 1); // @@ for now, we assume there is just one well equation
        const int nc = numFractureCells() + numWellEquations();
        const int cell = std::get<0>(perfinj_[0]); // @@ will this be the correct index? 
        const double lambda = 1.0;//reservoir_mobility_[cell]; // @@ If not constant, this might be wrong
        // control.get<double>("WI")
        const double WI_lambda = total_wellindex_ * lambda; // we could ahve used total mobility here
        matrix[nc - 1][nc - 1] = WI_lambda;
        // NB: well_source_[i] is assumed to be the same as get<0>(perfinj_[i])
        for (const auto& pi : perfinj_) {
            const int i = std::get<0>(pi);
            const double value = std::get<1>(pi) * mobility_water_perf_;// only flow in fracture not into reservoir
            matrix[nc - 1][i] = -value; // well equation
            matrix[nc - 1][nc - 1] += value; // well equation
            matrix[i][nc - 1] = -value;
            matrix[i][i] += value;
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknown control of injection into Fracture");
    }
}

bool
Fracture::removeNewZeroWithCells(RegularTrimesh& mesh,
                                 int cur_level,
                                 const RegularTrimesh& initial_mesh) const
{
    const double force_limit = prm_.get<double>("solver.force_limit", 0.0);
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
                    if (!initial_mesh.isActive(cellref)) {
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
                    //if (!any_active) {
                    if (!all_active) {
                        mesh.setInactive(cellref);
                        any_removed = true;
                    }
                }
            }
        }
    }
    return any_removed;
}



double
Fracture::normalFractureTraction(size_t eIdx) const
{
    return ddm::tractionSymTensor(reservoir_stress_[eIdx], cell_normals_[eIdx]);
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
Fracture::updateFractureRHS()
{
    assert(numWellEquations() == 0); // @@ not implemented/tested for rate-controlled systems
    rhs_width_ = fracture_pressure_;
    for (size_t i = 0; i < rhs_width_.size(); ++i) {
        std::cout << i << " " << rhs_width_[i] << std::endl;
        rhs_width_[i] = rhs_width_[i] - normalFractureTraction(i);
        if (rhs_width_[i] < 0.0)
            rhs_width_[i] = 0.0; // @@ not entirely accurate, but will avoid
        // unphysical negative normal displacements
    }
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

void
Fracture::printPressureMatrix() const // debug purposes
{
    Dune::printSparseMatrix(std::cout, *pressure_matrix_, "matname", "linameo");
}

void
Fracture::printMechMatrix() const // debug purposes
{
    Dune::printmatrix(std::cout, fractureMatrix(), "matname", "linameo");
}

// bool sortcsr(const tuple<size_t,size_t,double>& a,
//              const tuple<size_t,size_t,double>& b){
//     if(get<0>(a) < get<0>(b)){
//         return true;
//     }else{
//         if(get<1>(a) < get<1>(b)){
//             return true;
//         }else{
//             return false;
//         }
//     }
// }

// void buildRowWise(const std::tuple<size_t,size_t,double>& transes,Matrix& matrix){
//     matrix.setBuildMode(Matrx::row_wise);
//     matrix.setSize(1,nc,nc);
//     size_t row_index = 0;
//     auto trans_it = tranes.begin();
//     for (auto row = matrix.createbegin(),
//              end = matrix.createend(); row != end; ++row) {
//         while(row.index() == get<0>(*trans_it)){
//             row.insert(get<1>(*trans_it))
//                 ++trans_it;
//         }else{
//             ++row;
//         }
//     }
// }


// template void
// Fracture::updateReservoirCells<Dune::CpGrid>(const external::cvf::ref<external::cvf::BoundingBoxTree>&
// cellSearchTree,
//                                              const Dune::CpGrid& grid3D);
// template void
// Fracture::updateReservoirCells<Dune::PolyhedralGrid<3,3,double>>(const
// external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
//                                              const Dune::PolyhedralGrid<3,3,double>& grid3D);

template void Fracture::assignGeomechWellState(PerfData<float>& perfData) const;
template void Fracture::assignGeomechWellState(PerfData<double>& perfData) const;

} // namespace Opm
