#pragma once

#include "RegularTrimesh.hpp"
#include <opm/common/TimingMacros.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> 
namespace Opm
{

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

template <class TypeTag, class Simulator>
void Fracture::updateReservoirProperties(const Simulator& simulator, bool init_constant_vals,bool update_filtercake)
    {
        // if `init_contant_vals` is true, the fields that should normally not
        // change, i.e.  E_, nu_ and reservoir_perm_ will be updated. This should
        // normally only be needed the first time.
        if (init_constant_vals) {
            initReservoirProperties<TypeTag, Simulator>(simulator);
        }
      
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto& problem = simulator.problem();
        const auto& grid = simulator.vanguard().grid();
        GeometryHelper ghelper(grid);
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
        
        double numax = 0, Emax = 0;
        //auto& enitity_seeds = problem.elementEntitySeeds();//used for radom axes better way?
        // get properties from rservoir
        for (size_t i = 0; i < ncf; ++i) {
            int cell = reservoir_cells_[i];
            if (!(cell < 0)) {
                auto normal = this->cell_normals_[i];
                const auto& intQuants = simulator.model().intensiveQuantities(cell, /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                {
                  auto val = fs.pressure(FluidSystem::waterPhaseIdx);
                  reservoir_pressure_[i] = val.value();
                }
                enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
                reservoir_mobility_[i] = 0.0;
                reservoir_density_[i] = intQuants.fluidState().density(FluidSystem::waterPhaseIdx).value();
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
                const auto& perm =  problem.intrinsicPermeability(cell);
                const auto& cstress =  problem.cStress(cell);
                auto pn = normal;
                perm.mv(normal, pn);
                double npn = pn.dot(normal); 
                reservoir_perm_[i] = npn;
                reservoir_cstress_[i] = cstress;
                if(calculate_dist){
                  //assert(false);
                  const auto& currentData = grid.currentData();
                  const auto& elem = Dune::cpgrid::Entity<0>(*currentData.back(), cell, true);
                  //auto cell_center = ghelper.centroid(i);
                  const auto& geom = elem.geometry();
                  auto cell_center = geom.center();
                  reservoir_cell_z_[i] = cell_center[2];
                  double dist = 0;
                  int num_corners = geom.corners();
                  for(int li=0; li < geom.corners();++li){
                    auto cdist = cell_center;
                    cdist -= geom.corner(li);
                    dist += std::abs(normal.dot(cdist));
                  }
                  dist /= num_corners;
                  reservoir_dist_[i] = dist;
                }
                
            } else {
                // probably outside reservoir set all to zero
                double stressval = 300e5;
                reservoir_stress_[i][0] = stressval;
                reservoir_stress_[i][1] = stressval;
                reservoir_stress_[i][2] = stressval; //???
                reservoir_pressure_[i] = injectionPressure();
                reservoir_mobility_[i] = 0.0;
                reservoir_density_[i] = 1000.0;
                reservoir_perm_[i] = 0.0;
            }
            // assume reservoir distance is calculated
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

        if(update_filtercake){
            updateFilterCakeProps(connections, wellstate);
        } 
        // else {
        //   // remap
        //   int ncf = reservoir_cells_.size();
        //   filtercake_thikness_.resize(ncf, 0.0);
        // }
    }

// ----------------------------------------------------------------------------
template <class TypeTag, class Simulator>
void Fracture::solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
                     const Simulator& simulator)
// ----------------------------------------------------------------------------
{
    if(!active_){
        std::cout << "Fracture " << this->name() << " is not active, skipping solve." << std::endl;
        return;
    }
    OPM_TIMEBLOCK(SolveFracture);
    std::cout << "Solve Fracture Pressure" << std::endl;
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
        while (!fullSystemIteration(tol) && iter++ < max_iter) {
            if (nlin_verbosity > 1) {
                std::cout << "Iteration: " << iter << std::endl;
            }
        };

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

        if(true){
            fracture_width_ = 1e-3; // Ensure not completely closed
            fracture_pressure_ = perf_pressure_;
        }
        // start by assuming pressure equal to confining stress (will also set
        // fracture_pressure_ to its correct size
        normalFractureTraction(fracture_pressure_);

        // It is implicitly assumed for now that there is just one well equation.
        // We initialize with an existing value. @@
        if (numWellEquations() > 0)
            fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];

        // save original grid and filtercake, to allow us to map it onto evolved grids
        const auto filtercake_thickness_0 = filtercake_thikness_; // copy
        const auto grid_mesh_map_0 = grid_mesh_map_; 
        
        // local function taking a trimesh, updates the Fracture object with it and
        // runs a simulation.  Its return value should be a vector of doubles:

        auto score_function = [&](const RegularTrimesh& trimesh, const int level) -> std::vector<double> {
            const int max_iter = prm_.get<int>("solver.max_iter");
            const double tol = prm_.get<double>("solver.tolerance");//,1e-8);
            *trimesh_ = trimesh;
            // save well sources before grid change
            std::vector<CellRef> wsources = well_source_cellref_; 
            for (auto& cell : wsources) 
                cell = RegularTrimesh::fine_to_coarse(cell, level);
            
            // setup fracture with new grid
            const int MAX_NUM_COARSENING = prm_.get<int>("solver.max_num_coarsening"); // should be enough for all practical purposes
            const int numcell_threshold = prm_.get<int>("solver.numcell_threshold");
            auto [grid, fsmap, bmap] =
            trimesh_->createDuneGrid(MAX_NUM_COARSENING, wsources); // well cells kept intact!
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
            updateReservoirCells(cell_search_tree);
            updateReservoirProperties<TypeTag, Simulator>(simulator, true, false);
            initPressureMatrix();
            //initFractureWidth();
            //initFracturePressureFromReservoir();
            rhs_pressure_.resize(0);
            coupling_matrix_ = nullptr;


            bool point_wise = true;
            // solve flow-mechanical system
            bool remap_solution=false;
            if (remap_solution) {
                for (size_t i = 0; i < fracture_pressure_.size(); ++i) {
                    assert(std::abs(fracture_pressure_[i][0]) < 1e10);
                    assert(std::abs(fracture_width_[i][0]) < 0.6);
                }
                auto old_fracture_width_ = fracture_width_;
                redistribute_values(fracture_width_, grid_mesh_map_, fsmap, level, point_wise);

                redistribute_values(fracture_pressure_, grid_mesh_map_, fsmap, level, point_wise);
                // solve flow-mechanical system
                for (size_t i = 0; i < fracture_pressure_.size(); ++i) {
                    assert(std::abs(fracture_pressure_[i][0]) < 1e10);
                    assert(std::abs(fracture_width_[i][0]) < 0.6);
                }
            } else {
                initFractureWidth();
                initFracturePressureFromReservoir();
            }

            // filtercake is explicite
            filtercake_thikness_ = redistribute_values(filtercake_thickness_0,
                                                       grid_mesh_map_0,
                                                       fsmap, level, point_wise);
            
            int iter = 0;
            while (!fullSystemIteration(tol) && iter++ < max_iter) {
            };
            std::cout << "Iterations needed: " << iter << std::endl;

            // compute K1 stress intensity
            const std::vector<double> K1_not_nan = Fracture::stressIntensityK1();
            const std::vector<CellRef> boundary_cells = trimesh_->boundaryCells();
            std::vector<double> result(boundary_cells.size());
            for (size_t i = 0; i != result.size(); ++i){
                const double KImax = reservoir_cstress_[bmap[boundary_cells[i]]];
                const double KI = K1_not_nan[bmap[boundary_cells[i]]];
                result[i] = KI/KImax;
            }
            return result;
        };

        //const double K1max = prm_.get<double>("KMax");
        const double threshold = 1.0;
        const std::vector<CellRef> fixed_cells = well_source_cellref_;
        const int target_cellcount = prm_.get<int>("solver.target_cellcount"); 
        const int cellcount_threshold = prm_.get<int>("solver.cellcount_threshold");
        const auto [mesh, cur_level] =
            expand_to_criterion(*trimesh_, score_function, threshold,
                                fixed_cells, target_cellcount, cellcount_threshold);
        // make current level become the reference (finest) level
        // note that the well_source_cellref_ is already set from the last call to the score function
        for (auto& cell : well_source_cellref_) 
            cell = RegularTrimesh::fine_to_coarse(cell, cur_level);
        
        // ----------------------------------------------------------------------------
    } else if (method == "if_propagate") {
        // ----------------------------------------------------------------------------
        // iterate full nonlinear system until convergence, and expand fracture if necessary
        if(false){
            fracture_width_ = 1e-2; // Ensure not completely closed
            fracture_pressure_ = 0.0;
        }

        // start by assuming pressure equal to confining stress (will also set
        // fracture_pressure_ to its correct size
        normalFractureTraction(fracture_pressure_);
        if (numWellEquations() > 0){ // @@ it is implicitly assumed for now that
            // there is just one well equation.  We initializze
            // it with an existing value.
            fracture_pressure_[fracture_pressure_.size() - 1] = fracture_pressure_[0];
        }
        const int max_iter = prm_.get<int>("solver.max_iter");
        const double tol = prm_.get<double>("solver.tolerance");//,1e-8);


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
            while (!fullSystemIteration(tol) && iter++ < max_iter) {
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
            std::cout << "Iterations needed: " << iter << std::endl;


            // identify where max stress intensity is exceeded and propagation is needed
            const auto dist = grid_stretcher_->centroidEdgeDist();
            fill(bnode_disp.begin(), bnode_disp.end(), 0.0);

            const std::vector<double> K1_not_nan = Fracture::stressIntensityK1();
            std::vector<double> K1;
            bool should_fracture = false;
            for (size_t i = 0; i != K1_not_nan.size(); ++i){
                if (!std::isnan(K1_not_nan[i])){
                    K1max = reservoir_cstress_[i];
                    if(K1_not_nan[i] > K1max){
                        should_fracture = true;
                    }
                }
            }
            if(!should_fracture){
                break;
            }

            // loop over cells, determine how much they should be expanded or contracted
            const double maxgrow = rfac * grid_stretcher_->maxBoxLength();
            for (size_t i = 0; i != N; ++i) {
                K1max = reservoir_cstress_[boundary_cells[i]];
                cell_disp[i]
                    = efac * (compute_target_expansion(K1max, fracture_width_[boundary_cells[i]], E_, nu_) - dist[i]);
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
            updateReservoirCells(cell_search_tree);
            updateReservoirProperties<TypeTag, Simulator>(simulator, true, false);
            initPressureMatrix();
            fracture_matrix_ = nullptr;

            const auto& pts = grid_stretcher_->nodecoords();
            const auto bix = grid_stretcher_->boundaryNodeIndices();
            // for(auto i : bix)
            //   os << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
            // os.close();
        }
        if(count >= max_expand_iter){
           std::cout << "Fracture expansion did not converge within the maximum number of iterations" << std::endl;
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknowns solution method");
    }
}

} // namespace Opm
