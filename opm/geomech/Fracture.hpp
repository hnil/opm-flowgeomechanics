/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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

#ifndef OPM_FRACTURE_HH
#define OPM_FRACTURE_HH



#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include "GeometryHelpers.hpp"
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>

// for system solve
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <fstream>

#include <dune/common/indices.hh>
#include <dune/common/timer.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>
// for linear solve
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <dune/istl/matrixmarket.hh>

#include <opm/geomech/GridStretcher.hpp>
#include <opm/geomech/RegularTrimesh.hpp>
#include <opm/geomech/FracturePressureAssemblerAD.hpp>




#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>


namespace Dune{
   // Support FieldTraits for reference-based vectors
  template <typename... Args>
  struct FieldTraits< MultiTypeBlockVector<Args&...> >
  {
    using field_type = typename MultiTypeBlockVector<Args&...>::field_type;
    using real_type  = typename FieldTraits<field_type>::real_type;
  };
  template<typename FirstRow, typename... Args>
  struct FieldTraits< MultiTypeBlockMatrix<FirstRow&, Args&...> >
  {
    using field_type = typename MultiTypeBlockMatrix<FirstRow&, Args&...>::field_type;
    using real_type  = typename FieldTraits<field_type>::real_type;
  };
  
  /**
     @brief A reference-only MultiTypeBlockMatrix type

     This variant stores only references to rows and does not own data.
   */
  template<typename FirstRow, typename... Args>
  using MultiTypeBlockMatrixRef = MultiTypeBlockMatrix<FirstRow&, Args&...>;

  /**
     @brief Create a reference-only MultiTypeBlockMatrix
   */
  template<typename FirstRow, typename... Args>
  auto makeMultiTypeBlockMatrixRef(FirstRow& firstRow, Args&... args)
  {
    return MultiTypeBlockMatrixRef<FirstRow, Args...>(firstRow, args...);
  }

    /**
     @brief A MultiTypeBlockVector that stores only references to elements

     This is a non-owning view of vector elements.
   */
  template <typename... Args>
  using MultiTypeBlockVectorRef = MultiTypeBlockVector<Args&...>;

  /**
     @brief Create a reference-only view for MultiTypeBlockVector elements
   */
  template <typename... Args>
  auto makeMultiTypeBlockVectorRef(Args&... args)
  {
    return MultiTypeBlockVectorRef<Args...>(args...);
  }
}

namespace Opm {
    template <typename Scalar>
    class ConnFracStatistics;
}

namespace Opm::Properties {
    template <class TypeTag, class MyType>
    struct FluidSystem;

    template <class TypeTag, class MyType>
    struct NumPhases;
} // namespace Opm::Properties

namespace Opm
{
struct FractureProperties
{
  double height;
  double length;
  double flux;
  double area;
  double WI;
  double volume;
  double filter_volume;
  double avg_width;
  double avg_filter_width;
  double inj_pressure;
  double inj_bhp;
  double inj_wellrate;
  FractureProperties(double height_, double length_, double flux_, double area_,
                     double WI_, double volume_, double filter_volume_,
                     double avg_width_, double avg_filter_width_,
                     double inj_pressure, double inj_bhp, double inj_wellrate):
    height(height_),
    length(length_),
    flux(flux_),
    area(area_),
    WI(WI_),
    volume(volume_),
    filter_volume(filter_volume_),
    avg_width(avg_width_),
    avg_filter_width(avg_filter_width_),
    inj_pressure(inj_pressure),
    inj_bhp(inj_bhp),
    inj_wellrate(inj_wellrate)
  {
  };
};

  
struct WellInfo
{
    std::string name;
    int perf;
    int well_cell;
    int global_index;
    int segment{};
    std::optional<std::pair<double, double>> perf_range{};
};

  struct FractureSolveStats
  {
    int fractures_solved = 0;
    int nonlinear_iterations = 0;
    int linear_solves = 0;
    int linear_iterations = 0;
    double solve_time_seconds = 0.0;
    bool converged = true;
    // Contact (open/close) chatter: number of cells that flipped open<->closed
    // state, summed over nonlinear iterations. High values indicate the contact
    // set is oscillating and driving up the iteration count (Phase 2 signal).
    int closed_cell_toggles = 0;
    // Robustness signals (Phase 1): coupled linear solves that did not converge,
    // and how many of those were recovered by the fallback ladder.
    int linear_solve_failures = 0;
    int ladder_rescues = 0;

    FractureSolveStats& operator+=(const FractureSolveStats& other)
    {
      fractures_solved += other.fractures_solved;
      nonlinear_iterations += other.nonlinear_iterations;
      linear_solves += other.linear_solves;
      linear_iterations += other.linear_iterations;
      solve_time_seconds += other.solve_time_seconds;
      converged = converged && other.converged;
      closed_cell_toggles += other.closed_cell_toggles;
      linear_solve_failures += other.linear_solve_failures;
      ladder_rescues += other.ladder_rescues;
      return *this;
    }
  };

struct RuntimePerforation;

/// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
class Fracture
{
public:
    using IndexTraits = Opm::BlackOilDefaultFluidSystemIndices;
    using Grid = Dune::FoamGrid<2, 3>;
    using Point3D = Dune::FieldVector<double, 3>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using DynamicMatrix = Dune::DynamicMatrix<double>;
    using EntitySeed = Dune::CpGrid::Codim<0>::Entity::EntitySeed;

    void init(const std::string& well,
              const int perf,
              const int well_cell,
              const int global_index,
              const int segment,
              const std::optional<std::pair<double, double>>& perf_range,
              const Point3D& origo,
              const Point3D& normal,
              const PropertyTree& prm);

    void grow(int layers, int method);
    std::string name() const;
    void write(int reportStep = -1) const;
    void writemulti(double time) const;
    void updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                            const Dune::CpGrid& grid,
                            const std::vector<EntitySeed>& entity_seeds);
    template <class TypeTag, class Simulator>
    void updateFilterCakePropertiesPost(const Simulator& simulator);
    // solver related
    void updateReservoirProperties();
    void removeCells();

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator, bool init_constant_vals=false,bool update_filtercake=true);
    
    void updateFilterCakeProps(const Opm::WellConnections& connections,
                               const Opm::SingleWellState<double,IndexTraits>& wellstate,double dt);    
    void initFracturePressureFromReservoir();
    void initFractureStates();
    void initFractureWidth();
    void solveFractureWidth();
    void solvePressure();

    template <class TypeTag, class Simulator>
    void solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
               //const Dune::CpGrid& grid,
               const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
               const Simulator& simulator);
    FractureProperties calculateFractureProperties() const;
    void printPressureMatrix() const; // debug purposes
    void printMechMatrix() const; // debug purposes
    void writeFractureSystem()  const;
    void writePressureSystem()  const;
    void setFractureGrid(std::unique_ptr<Fracture::Grid> gptr = nullptr); // a hack to allow use of another grid
    std::vector<RuntimePerforation> wellIndices() const;
    std::vector<RuntimePerforation> wellIndicesAvrg(const std::vector<std::vector<RuntimePerforation>>& well_indices) const;
    WellInfo& wellInfo();
    const WellInfo& wellInfo() const { return wellinfo_; }
    std::vector<double> leakOfRate() const;

    // ---- Access for the embedded (auxiliary-cell) fracture flow representation ----
    //
    // The upscaled well-index formulation consumes these quantities internally; the
    // embedded one hands them to the reservoir discretization instead, which needs to
    // read them from outside.

    //! Number of cells in the fracture grid.
    std::size_t numCells() const { return numFractureCells(); }

    //! Reservoir cell each fracture cell leaks into, by fracture cell index.
    const std::vector<int>& reservoirCells() const { return reservoir_cells_; }

    //! Leak-off coefficient per fracture cell.  Includes the reservoir mobility.
    const std::vector<double>& leakOf() const { return leakof_; }

    //! Reservoir mobility used to form leakOf(), so that it can be divided back out.
    const std::vector<double>& reservoirMobility() const { return reservoir_mobility_; }

    //! Well-to-fracture connections: (fracture cell, well index).  The well injects
    //! into these cells directly, which in the embedded representation is a perforation
    //! of the fracture's own degrees of freedom rather than of the reservoir.
    const std::vector<std::tuple<int, double>>& wellPerforations() const { return perfinj_; }

    //! Half transmissibilities (i, j, t_i, t_j) between neighbouring fracture cells.
    const std::vector<Htrans>& halfTrans() const { return htrans_; }

    //! Aperture per fracture cell.
    const Dune::BlockVector<Dune::FieldVector<double, 1>>& fractureWidth() const
    { return fracture_width_; }

    //! Area of each fracture cell.
    std::vector<double> cellAreas() const;

    //! Depth of each fracture cell's centre, positive downwards.
    std::vector<double> cellDepths() const;
    double injectionPressure() const;
    double injectionBhp() const;// {return injectionPressure() - dp_perf_;};
    void setPerfProps(double perfpressure,double perf_depth, double perfrate);//{perf_pressure_ = perfpressure;}
    void setWellProps(double wellrate, double WI, double wi_dp, double wi_respress, double ref_depth);//{well_rate_ = wellrate; total_wellindex_ = WI;}
    Dune::FieldVector<double, 6> stress(Dune::FieldVector<double, 3> obs) const;
    Dune::FieldVector<double, 6> strain(Dune::FieldVector<double, 3> obs) const;
    Dune::FieldVector<double, 3> disp(Dune::FieldVector<double, 3> obs) const;

    template <typename Scalar>
    void assignGeoMechWellState(PerfData<Scalar>& perfData) const;
    void setActive(bool active);
    bool isActive() const;
    const FractureSolveStats& lastSolveStats() const;
    std::array<double,2> hightAndWidth() const;
    double maxFlowTimeStep() const;
    double filterCakeVolume() const;
    void resetFracture();
    void moveForwardInTime();

private:
   struct CouplingMixState
   {
     bool initialized{false};
     bool has_previous_residual{false};
     double previous_target{0.0};
     double previous_residual{0.0};
     double omega{1.0};
   };

   double applyCouplingUpdate(double current,
                double target,
                double damping_factor,
                const std::string& mode,
                const std::string& channel,
                CouplingMixState& state);

   // Vector fixed-point state for the per-cell CTF coupling (opt-in
   // solver.wi_vector_acceleration). Keyed by reservoir cell index so it is
   // stable across fracture-grid changes; new cells take the target directly.
   struct VectorMixState
   {
     bool initialized{false};
     bool has_prev_residual{false};
     double omega{1.0};
     std::map<int, double> prev_mixed;    // accepted CTF of the previous solve
     std::map<int, double> prev_target;   // raw CTF target of the previous solve
     std::map<int, double> prev_residual; // target - prev_mixed of the previous solve
   };

   std::vector<RuntimePerforation>
   applyVectorCouplingUpdate(const std::vector<RuntimePerforation>& target);

   double reservoirTraction(int i) const;
   double fractureForce(int i) const;
   double leakofDp(int i) const;
   void summary_of_solve();
   std::vector<RuntimePerforation> wellIndices_() const;
   bool  expantionMax(const FractureProperties& fprop);
   bool removeNewZeroWithCells(RegularTrimesh& mesh,int cur_level,const RegularTrimesh& original_mesh) const;
   std::vector<double> redistribute_values(const std::vector<double>& values,
                                        const std::vector<std::vector<CellRef>>& map1,
                                        const std::vector<std::vector<CellRef>>& map2,
                                        const int level,
                                                     bool point_wise);

   void redistribute_values(Dune::BlockVector<Dune::FieldVector<double, 1>>& values,
                                        const std::vector<std::vector<CellRef>>& map1,
                                        const std::vector<std::vector<CellRef>>& map2,
                                        const int level,
                             bool point_wise);
  
    using ResVector = ::Opm::Fracture::Vector;
    using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

  // using Krull = Dune::MultiTypeBlockVector<Dune::BlockVector<double>>;
  // Dune::MultiTypeBlockVector<Dune::BlockVector<double, std::allocator<double>>> dill;
  // Krull tull(2);  
  
    using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
    using FMatrix = Dune::DynamicMatrix<double>;                       // full matrix

    size_t numFractureCells() const { return grid_->leafGridView().size(0); }
    std::vector<int> identify_closed(const FMatrix& A,
                                     const VectorHP& x,
                                     const ResVector& rhs,
                                     const int nwells);
    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator);

    void resetWriters();
    Dune::BlockVector<Dune::FieldVector<double, 3>> all_slips() const;
    // helpers for growing grid
    void insertLinear(const std::vector<unsigned int>& inner_indices);
    void insertExp(const std::vector<unsigned int>& inner_indices);
    void initFracture(); // create a new fracture grid from scratch
    Point3D surfaceMap(double x, double y);

    std::unique_ptr<GridStretcher> grid_stretcher_; //@@ experimental, for stretching grids
    std::unique_ptr<Opm::RegularTrimesh> trimesh_; // @@ experimental, implicitly defined grids

    std::vector<std::vector<CellRef>> grid_mesh_map_; // @@ index mapping from cells in grid_ to trimesh_.
                                                      // @@ NB: in general many-to-many
    std::unique_ptr<Grid> grid_;
    // grid definition in case of step fails
    std::unique_ptr<Grid> grid_prev_;
    std::unique_ptr<GridStretcher> grid_stretcher_prev_; //@@ experimental, for stretching grids
    std::unique_ptr<Opm::RegularTrimesh> trimesh_prev_; // @@ experimental, implicitly defined grids
    std::vector<std::vector<CellRef>> grid_mesh_map_prev_; // @@ index mapping from cells in grid_ to trimesh_.
                                                      // @@ NB: in general many-to-many
    std::vector<double> filtercake_thikness_prev_; // properties of filter cake, if any
    // well source checkpoints (needed to rebuild the dune grid on reset)
    std::vector<CellRef> well_source_cellref_prev_;
    std::vector<int> well_source_prev_;
    // copy of state variables
  
    

    Point3D origo_;
    std::array<Point3D, 3> axis_;
    std::array<Point3D, 3> naxis_;
    WellInfo wellinfo_;
    bool active_{false}; // is fracture active?
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    static constexpr auto VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<::Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_;
    std::vector<unsigned int> out_indices_;

    std::vector<double> filtercake_thikness_; // properties of filter cake, if any
    double filtercake_perm_{0.0}; // permeability of filter cake, if any
    double filtercake_poro_{0.0}; // permeability of filter cake, if any
    bool has_filtercake_{false}; // if true, filter cake is used in the model

    // should probably not be needed
    int layers_;
    int nlinear_;

    // help function for solving
    FracturePressureInput makePressureAssemblyInput() const;
    void assemblePressure();
    void assemblePressureAndCouplingAD(const std::vector<int>& closed_cells);
    void addSource();
    void initPressureMatrix();
    void setupPressureSolver();
    void updateFractureRHS();
    void updateLeakoff();
    void updateCellNormals();
    void normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction, bool resize=true) const;
    double normalFractureTraction(size_t ix) const;

    // one nonlinear iteration of fully coupled system.  Returns 'true' if converged
  bool fullSystemIteration(const double tol,const int iteration);

    void assembleFractureMatrix() const;
    std::vector<double> stressIntensityK1() const;
    int numWellEquations() const;
  
    //double well_pressure_;// for now using prm object for definition
    std::vector<CellRef> well_source_cellref_; // references to well cells in the fully resolved TriMesh
    std::vector<int> well_source_; // indices to well cells
    // for reservoir
    std::vector<int> reservoir_cells_;
    std::vector<std::vector<int>> all_reservoir_cells_;
    std::vector<std::vector<double>> all_reservoir_areas_;
    std::vector<std::vector<Dune::FieldVector<double,3>>> all_reservoir_centers_;
    // std::vector< Dune::FieldMatrix<double, 3, 3> > reservoir_perm_;
    std::vector<double> reservoir_perm_;
    std::vector<double> reservoir_cstress_;
    std::vector<double> reservoir_mobility_;
    std::vector<double> reservoir_density_;
    std::vector<double> reservoir_cell_z_;
    using WaterPropertyEvaluator = std::function<std::pair<CellFluidProperty, CellFluidProperty>(size_t, double)>;
    WaterPropertyEvaluator fracture_water_property_evaluator_;
    //
    //std::vector<int> unique_reservoir_cells_;
    std::map<int,double> map_reservoir_mobility_;
    std::map<int,double> map_reservoir_density_;
    std::map<int,double> map_reservoir_cell_z_;
    std::map<int,double> map_reservoir_pressure_;
    //

    std::vector<double> reservoir_dist_;
    std::vector<double> reservoir_pressure_;
    std::vector<Dune::FieldVector<double, 6>> reservoir_stress_;

    // only for radom access need to be updater after trid change
    Dune::BlockVector<Dune::FieldVector<double, 3>> cell_normals_;

    // solution variables (only to avoid memory allocation, do not trust their state)

    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_pressure_;

// State variables
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_pressure_;


    // copy of state variables
    Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_width_prev_;
    Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_pressure_prev_;

    // transmissibilities
    using Htrans = std::tuple<size_t, size_t, double, double>;
    std::vector<Htrans> htrans_;
    std::vector<std::tuple<int,double>> perfinj_;
    double perf_pressure_;
    std::vector<double> leakof_;
    
    PropertyTree prmpressure_;
    using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
    mutable std::unique_ptr<Matrix> pressure_matrix_;
    mutable std::unique_ptr<PressureOperatorType> pressure_operator_;
    mutable std::unique_ptr<FlexibleSolverType> pressure_solver_;
    mutable std::unique_ptr<Matrix> coupling_matrix_; // will be updated by `fullSystemIteration`
  
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DynamicMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    mutable std::unique_ptr<DynamicMatrix> fracture_matrix_; 

    // function ensuring that the fracture matrix exists, and returning a reference to it
    DynamicMatrix& fractureMatrix() const;
  
    double E_;
    double nu_;
    double min_width_; // minimum width of fracture, used for convergence criterion
    double gravity_{Opm::unit::gravity};//{9.81}; // gravity acceleration, used for leakoff calculations
    std::vector<double> fracture_dgh_; // gravity contribution to fracture pressure, used for leakoff calculations  
    PropertyTree prm_;
    double total_WI_well_{0.0}; // total well index for the well, used for leakoff calculations
    // for coupled solve
  // using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
  //using FMatrix = Dune::DynamicMatrix<double>; 
    using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                     Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
    using SystemMatrixRef = Dune::MultiTypeBlockMatrixRef<Dune::MultiTypeBlockVectorRef<FMatrix, SMatrix>,
                                                     Dune::MultiTypeBlockVectorRef<SMatrix, SMatrix>>;
    //using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<std::reference_wrapper<FMatrix>,std::reference_wrapper<SMatrix>>,
    //                                                Dune::MultiTypeBlockVector<std::reference_wrapper<SMatrix>, std::reference_wrapper<SMatrix>>>;                                                 
    std::unique_ptr<SystemMatrix> S_;
    std::unique_ptr<FMatrix> A_; // system matrix without coupling terms    
    //std::unique_ptr<SMatrix> C_; // system matrix without coupling terms    
    std::unique_ptr<SMatrix> I_; // system matrix without coupling terms    
  //std::unique_ptr<SMatrix> M_; // system matrix without coupling terms
    // using SystemOperator = Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>;
    // std::unique_ptr<SystemOperator> S_linop_;
    // using CoupledFractureSolver = Dune::InverseOperator<VectorHP, VectorHP>;
    // std::unique_ptr<CoupledFractureSolver> coupled_solver_;
    std::vector<std::vector<RuntimePerforation>> well_indices_; // for each perf, the corresponding well index info
    double max_flow_time_step_{-1.0};
    double well_rate_;
    double total_wellindex_;
    //double dz_perf_{0.0};
    double wi_dz_{0.0};
    double wi_respress_{0.0};
    double well_ref_depth_{0.0};
    double perf_ref_depth_{0.0};
  //double ref_depth_{0.0};
    double well_perf_rate_{0.0};
    CouplingMixState perf_pressure_mix_{};
    CouplingMixState well_rate_mix_{};
    CouplingMixState total_wellindex_mix_{};
    VectorMixState wi_vec_mix_{};
    // mixed CTF vector returned by wellIndices() when wi_vector_acceleration is on
    std::vector<RuntimePerforation> well_indices_accel_;
    double density_{1000.0};// default to water
    double density_perf_{1000.0};
    double mobility_water_perf_{1000.0};

    // handling preconditioner
    std::vector<int> closed_cells_; // indices of cells that are closed (i.e.
    // per-cell contact state changes within the current solve, used by the
    // opt-in anti-chatter guard (solver.toggle_guard_mode = "chatter")
    std::vector<int> cell_flip_counts_;
    //using AbstractPreconditioner = Dune::PreconditionerWithUpdate<VectorHP, VectorHP>;
    using AbstractPreconditioner = FractureMechanicsPreconditioner;//<VectorHP, VectorHP>;
    std::unique_ptr<AbstractPreconditioner> frac_flow_precond_; // have zero width) in the current iteration, used for modifying the system matrix and convergence criterion
    std::unique_ptr<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>> S_linop_;
    using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;
    std::unique_ptr<LinearSolverBase> psolver_;
    FractureSolveStats last_solve_stats_{};

};
} // namespace Opm

#include "Fracture_impl.hpp"
#endif
