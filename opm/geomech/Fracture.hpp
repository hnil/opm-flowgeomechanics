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



#include <opm/geomech/CentreStressInterpolation.hpp>
#include <opm/geomech/FractureConnectionId.hpp>
#include <limits>
#include <optional>
#include <unordered_map>

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
#include <opm/geomech/GeomechCheckpoint.hpp>


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
    int global_index; // Level-zero WSEED identity, used for activation.
    FractureConnectionId connection_id{};
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
    // set is oscillating and driving up the iteration count.
    int closed_cell_toggles = 0;
    // Robustness signals: coupled linear solves that did not converge,
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

    // The completion identity is required; the WSEED parent index cannot
    // safely provide a default, even when a caller expects grid zero.
    void init(const std::string& well,
              const int perf,
              const int well_cell,
              const int global_index,
              const int segment,
              const std::optional<std::pair<double, double>>& perf_range,
              const Point3D& origo,
              const Point3D& normal,
              const PropertyTree& prm,
              FractureConnectionId connection_id);

    std::string name() const;
    void write(int reportStep = -1) const;
    void writemulti(double time) const;
    // Opt-in (solver.write_iteration_snapshots): append the current fracture
    // state to <output>/<name>_iter.csv (per cell) and <name>_iter_summary.csv
    // (one line) tagged with timestep, outer coupling round and a label — lets
    // the open/closed branches the coupling visits be plotted round by round.
    void writeIterationSnapshot(int step, int round, const std::string& tag) const;
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

    void updateFilterCakeProps(const Opm::Connection& connection,
                               const Opm::SingleWellState<double,IndexTraits>& wellstate,double dt);
    void initFracturePressureFromReservoir();
    double hydrostaticReservoirPressure(size_t i) const;
    void initFractureStates();
    void initFractureWidth();

    // GrowOnly: one growth round on the current state, without the fracture's own
    // solve; new cells get no width (the joint well-fracture solve decides them).
    enum class SolveMode { Full, GrowOnly };
    template <class TypeTag, class Simulator>
    void solve(const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
               const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
               const Simulator& simulator,
               SolveMode mode = SolveMode::Full);
    FractureProperties calculateFractureProperties() const;
    void setFractureGrid(std::unique_ptr<Fracture::Grid> gptr);
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

    /*!
     * \brief Recompute the leak-off coefficients if they lag the grid.
     *
     * The propagation may rebuild the grid after the last leak-off update -- an attempt
     * that was rolled back leaves leakof_ sized for the grid that was tried -- and the
     * embedded flow representation reads it per current cell.  Returns whether leak-off
     * now matches the grid; false when the reservoir property arrays lag it as well, in
     * which case nothing sensible can be computed and the caller keeps what it has.
     */
    bool ensureLeakoffCurrent();

    //! Reservoir cell each fracture cell leaks into, by fracture cell index.
    const std::vector<int>& reservoirCells() const { return reservoir_map_.cells; }

    //! Leak-off coefficient per fracture cell.  Includes the reservoir mobility.
    const std::vector<double>& leakOf() const { return leakof_; }

    //! Reservoir mobility used to form leakOf(), so that it can be divided back out.
    const std::vector<double>& reservoirMobility() const { return reservoir_.mobility; }

    //! Distance from the fracture face to the reservoir cell's pressure point used by
    //! the leak-off, and the number of faces that leak.
    const std::vector<double>& reservoirDistance() const { return reservoir_map_.dist; }
    int leakingSides() const { return leakoffSides(); }

    //! The mobility actually folded into leakOf(), which is reservoirMobility() only
    //! under the legacy convention. Divide by this, not by reservoirMobility().
    const std::vector<double>& leakoffMobility() const { return leakoff_mobility_; }
    //! Mobility the flow applies on the embedded leak-off connection of a cell.
    double flowLeakoffMobility(int eIdx) const;

    //! Pressure per fracture cell as the fracture's own solver last left it.  In the
    //! embedded representation the same quantity is also a degree of freedom of the
    //! flow problem, and the two must agree; this is what makes that checkable.
    const Dune::BlockVector<Dune::FieldVector<double, 1>>& fracturePressure() const
    { return fracture_pressure_; }

    //! The width floor this fracture's own cubic-law assembly applies
    //! (config.min_width).  The embedded representation must apply the same floor to
    //! the same law, or the two solve differently conductive fractures.
    double cubicLawMinWidth() const { return min_width_; }

    //! Well-to-fracture connections: (fracture cell, well index).  The well injects
    //! into these cells directly, which in the embedded representation is a perforation
    //! of the fracture's own degrees of freedom rather than of the reservoir.
    const std::vector<std::tuple<int, double>>& wellPerforations() const { return perfinj_; }

    //! Half transmissibilities (i, j, t_i, t_j) between neighbouring fracture cells.
    const std::vector<Htrans>& halfTrans() const { return htrans_; }

    /*!
     * \brief Half transmissibilities computed fresh from the current grid.
     *
     * The cached htrans_ is rebuilt at particular points of the pressure solve and can
     * lag a re-gridding with indices that are still in range -- a wrong answer that no
     * size check catches.  A consumer describing the fracture to another system reads
     * this instead: pure geometry of the grid as it stands, nothing cached.
     */
    std::vector<Htrans> currentHalfTrans() const;

    //! Aperture per fracture cell.
    const Dune::BlockVector<Dune::FieldVector<double, 1>>& fractureWidth() const
    { return fracture_width_; }

    //! Area of each fracture cell.
    std::vector<double> cellAreas() const;
    //! Indices of the cells contact currently holds shut; empty before the first
    //! solve. Used by the fracture <-> mechanics coupling blocks for their mask.
    const std::vector<int>& closedCells() const { return closed_cells_; }

    //! Depth of each fracture cell's centre, positive downwards.
    std::vector<double> cellDepths() const;
    double injectionPressure() const;
    double injectionBhp() const;
    //! true once fracture_pressure_ holds a solved state (injectionPressure/
    //! injectionBhp index into it and must not be called before)
    bool hasPressureState() const { return fracture_pressure_.size() > 0; }
    //! Warn when the fracture is mechanically open yet carries no flow - the
    //! signature of a lost well connection (e.g. CTFs zeroed by the legacy WI
    //! upscaling), which otherwise passes silently as a pressure transient.
    void warnIfOpenButNotConducting() const;
    void setPerfProps(double perfpressure,double perf_depth, double perfrate);
    //! resolved fracture-BC control: "control.type" verbatim, or picked from the
    //! well's active constraint when control.type == "well"
    std::string effectiveControlType() const;
    void setWellControl(bool is_rate, double bhp)
    {
        well_control_is_rate_ = is_rate;
        well_target_bhp_ = bhp;
    }
    void setWellProps(double wellrate, double WI, double wi_dp, double wi_respress, double ref_depth);
    // Reservoir cells the well perforates (all perforations, not just this
    // fracture's), used by solver.well_source_all_perfs to feed the fracture
    // wherever the well crosses it.
    void setWellPerfCells(std::vector<int> cells) { well_perf_cells_ = std::move(cells); }
    //! current flow timestep (s); needed by the opt-in fracture storage term
    void setTimeStep(double dt) { current_dt_ = dt; }
    //! Back to solving the fracture's own pressure (see setExternalPressure).
    void clearExternalPressure() { external_pressure_ = false; }
    //! Whether the last SolveMode::GrowOnly solve added cells.
    bool grewInLastSolve() const { return joint_grew_; }

    // State and blocks for WellFractureJointSolver.
    const DynamicMatrix& jointElasticMatrix() const { return fractureMatrix(); }
    Dune::BlockVector<Dune::FieldVector<double, 1>> jointNormalTraction()
    {
        Dune::BlockVector<Dune::FieldVector<double, 1>> traction(fracture_width_);
        normalFractureTraction(traction, false);
        return traction;
    }
    void setJointWidth(const Dune::BlockVector<Dune::FieldVector<double, 1>>& width)
    { fracture_width_ = width; }
    void setJointContactState(const std::vector<double>& pressure, const std::vector<int>& closed)
    {
        for (std::size_t i = 0; i < pressure.size() && i < fracture_pressure_.size(); ++i) {
            fracture_pressure_[i] = pressure[i];
        }
        closed_cells_ = closed;
    }
    //! Per cell: the cell of the layout at the last markFlowStateBound() whose flow
    //! state it inherits across regrids (-1: none). Only valid if flowStateRemapped().
    const std::vector<int>& flowStateDonor() const { return flow_state_donor_; }
    bool flowStateRemapped() const { return flow_state_remapped_; }
    void markFlowStateBound() { flow_state_donor_.clear(); flow_state_remapped_ = false; }
    //! Re-solve on the current mesh without growing it (one call, set by the caller).
    void suppressPropagation(bool on) { suppress_propagation_ = on; }
    //! Set the fracture pressure per cell (and the well DOF, if any) from outside
    //! and switch to external-pressure mode: the coupled iteration then keeps the
    //! pressure fixed and solves mechanics, contact and propagation at it.
    //! False (mode off) if the sizes do not match the current grid or no cell is pinned.
    //! A NaN entry marks a cell with no flow DOF (no reservoir partner, or born
    //! since the last bind): it keeps the fracture's own pressure equation.
    bool setExternalPressure(const std::vector<double>& cellPressures, double wellPressure)
    {
        const std::size_t nc = numFractureCells();
        if (cellPressures.size() != nc || fracture_pressure_.size() < nc) {
            external_pressure_ = false;
            return false;
        }
        external_cell_mask_.assign(nc, 0);
        std::size_t pinned = 0;
        for (std::size_t i = 0; i < nc; ++i) {
            if (std::isfinite(cellPressures[i])) {
                fracture_pressure_[i][0] = cellPressures[i];
                external_cell_mask_[i] = 1;
                ++pinned;
            }
        }
        if (numWellEquations() > 0 && fracture_pressure_.size() > nc) {
            fracture_pressure_[nc][0] = wellPressure;
        }
        external_pressure_ = pinned > 0;
        return external_pressure_;
    }
    Dune::FieldVector<double, 6> stress(Dune::FieldVector<double, 3> obs) const;
    Dune::FieldVector<double, 6> strain(Dune::FieldVector<double, 3> obs) const;
    Dune::FieldVector<double, 3> disp(Dune::FieldVector<double, 3> obs) const;

    template <typename Scalar>
    void assignGeoMechWellState(PerfData<Scalar>& perfData) const;
    void setActive(bool active);
    bool isActive() const;
    const FractureSolveStats& lastSolveStats() const;
    double maxFlowTimeStep() const;
    // Per-step growth guard: area now vs at the step checkpoint. Returns a
    // non-empty message if growth exceeded solver.max_area_growth_per_step
    // (relative factor; the absolute floor solver.min_area_growth_guard lets a
    // small seed open freely). Empty string = fine / guard off.
    std::string growthGuardViolation() const;
    // Onset dt-hold: true while the fracture is still an unopened seed and the
    // perf pressure is rising (solver.onset_dt_limit > 0 arms it).
    bool onsetHoldActive() const;
    double filterCakeVolume() const;
    void resetFracture();
    void moveForwardInTime(double dt_last = -1.0);

    //! Checkpoint of everything a solve leaves behind. Loading rebuilds the grid from
    //! the mesh; matrices and solvers are rebuilt by the next solve.
    template <class Serializer>
    void serializeState(Serializer& serializer);

private:
    void rebuildGridForLoad(const std::vector<std::vector<CellRef>>& saved_map);

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
   enum class WiMode { Legacy, Conductivity, FracturePressure };
   struct WiOptions
   {
       WiMode mode{WiMode::Legacy};
       bool divideWellIdx{false};
       bool fluxNorm{true};
       bool signGate{false};
       double dpFloor{1.0e2};  // Pa
       double alphaMax{2.0};
       static WiOptions read(const PropertyTree& prm);
   };
   // Fracture-to-reservoir exchange gathered per reservoir cell.
   struct CellExchange
   {
       int cell{};
       double q{}, p{}, mob{}, dens{}, z{}, trans{}, area{};
       double pfracLeak{}, leak{}, pfracArea{}; // WiMode::FracturePressure only
   };
   std::vector<CellExchange> gatherCellExchange_(const WiOptions& opt) const;
   RuntimePerforation perforation_(const CellExchange& c, double inj_press) const;
   std::vector<RuntimePerforation> wiLegacy_(const std::vector<CellExchange>& cells,
                                             double inj_press) const;
   std::vector<RuntimePerforation> wiConductivity_(const std::vector<CellExchange>& cells,
                                                   double inj_press, const WiOptions& opt) const;
   std::vector<RuntimePerforation> wiFracturePressure_(const std::vector<CellExchange>& cells,
                                                       double inj_press) const;
   std::vector<RuntimePerforation> wellIndices_() const;
   bool  expantionMax(const FractureProperties& fprop);

   // ---- trimesh propagation solve (solve()) ----
   struct PropagationOptions
   {
       int maxIter{};
       double tol{};
       int nlinVerbosity{};
       int maxNumCoarsening{};
       int numcellThreshold{};
       bool smoothBoundary{false};
       bool remapSolution{true};
       bool conservativeCake{false};
       double forceLimit{0.0};
       bool conservative{false};
       bool convergedCellsOnly{false};
       double fbTol{1e-3};
       int stableIters{2};
       int verbosity{0};
       bool disable{false};
       bool openFrontOnly{false};
       static PropagationOptions read(const PropertyTree& prm);
   };
   // Stops expansion once new fracture volume exceeds what the well can deliver
   // this step; optionally draws the BC pressure down by the over-draw.
   struct VolumePacing
   {
       bool on{false};
       bool pressure{false};
       bool wholeWell{false};
       double fac{1.0};
       double dt{0.0};
       double v0{0.0};
       double perfP0{0.0};
       int verbosity{0};
       bool budgetSpent{false};
       static VolumePacing start(const Fracture& f, double dt);
       void afterSolve(Fracture& f);
       void restorePerfPressure(Fracture& f) const;
   };
   struct FilterCakeSnapshot
   {
       std::vector<double> thickness;
       std::vector<std::vector<CellRef>> gridMeshMap;
       double total{0.0};
       RegularTrimesh::CellVolumes volumes;
   };
   struct RegridMaps
   {
       std::vector<std::vector<CellRef>> before; // grid_mesh_map_ before the regrid
       std::vector<std::vector<CellRef>> after;
       std::map<CellRef, int> boundary;          // trimesh boundary cell -> grid cell
   };
   struct GridSolve
   {
       bool converged{true};
       int iterations{0};
   };
   void reinitializeTrimeshState_();
   FilterCakeSnapshot snapshotFilterCake_(const PropagationOptions& opt) const;
   template <class TypeTag, class Simulator>
   RegridMaps regridTo_(const RegularTrimesh& trimesh, int level, const PropagationOptions& opt,
                        const external::cvf::ref<external::cvf::BoundingBoxTree>& cell_search_tree,
                        const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds,
                        const Simulator& simulator);
   void remapState_(const RegridMaps& maps, int level, const PropagationOptions& opt,
                    const FilterCakeSnapshot& cake0);
   GridSolve solveOnCurrentGrid_(const PropagationOptions& opt);
   void reportGridSolve_(const GridSolve& solve, const PropagationOptions& opt);
   // K1/K1c per trimesh front cell; -1 where a veto forbids growth
   std::vector<double> propagationScores_(std::map<CellRef, int>& boundary, const GridSolve& solve,
                                          const PropagationOptions& opt) const;
   void vetoForceLimit_(std::vector<double>& score, const std::vector<int>& front,
                        const PropagationOptions& opt) const;
   void vetoUnconvergedCells_(std::vector<double>& score, const std::vector<int>& front,
                              const PropagationOptions& opt) const;
   void vetoClosedFront_(std::vector<double>& score, const std::vector<int>& front) const;
   bool removeNewZeroWithCells(RegularTrimesh& mesh,int cur_level,const RegularTrimesh& original_mesh) const;
   bool suppress_propagation_ = false;
   void composeFlowStateDonor(const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level);
   RegularTrimesh::CellVolumes filterCakeVolumesByCell() const;
   std::vector<double> remapFilterCakeConservative(const RegularTrimesh::CellVolumes& source,
                                                   const std::vector<std::vector<CellRef>>& source_map,
                                                   const int level) const;

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

    std::unique_ptr<Opm::RegularTrimesh> trimesh_; // experimental, implicitly defined grids

    std::vector<std::vector<CellRef>> grid_mesh_map_; // index mapping from cells in grid_ to trimesh_.
                                                      // NB: in general many-to-many
    std::vector<CellRef> grid_fixed_cells_; // cells kept fine when grid_ was built from trimesh_
    std::unique_ptr<Grid> grid_;
    // grid definition in case of step fails
    std::unique_ptr<Opm::RegularTrimesh> trimesh_prev_; // experimental, implicitly defined grids
    std::vector<std::vector<CellRef>> grid_mesh_map_prev_; // index mapping from cells in grid_ to trimesh_.
                                                      // NB: in general many-to-many
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
    bool joint_growth_round_{false};
    int joint_scores_{0};
    bool joint_grew_{false};
    bool joint_level_held_{false}; // the seed was coarsened once; growth keeps that level
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    static constexpr auto VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<::Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_;

    std::vector<double> filtercake_thikness_; // properties of filter cake, if any
    double filtercake_perm_{0.0}; // permeability of filter cake, if any
    double filtercake_poro_{0.0}; // permeability of filter cake, if any
    bool has_filtercake_{false}; // if true, filter cake is used in the model

    // help function for solving
    FracturePressureInput makePressureAssemblyInput() const;
    std::vector<double> faceGravityHeads() const;
    void assemblePressure();
    void assemblePressureAndCouplingAD(const std::vector<int>& closed_cells);
    void addSource();
    void initPressureMatrix();
    void updateLeakoff();
    double leakoffMobilityFor(int eIdx, bool upwind) const;
    double leakoffPotentialDifference(int eIdx) const;
    // leak-off convention (solver.leakoff_model): number of leaking faces and the
    // reservoir distance as a fraction of the cell extent normal to the fracture
    int leakoffSides() const;
    double leakoffDistanceFactor() const;
    void updateCellNormals();
    void normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction, bool resize=true) const;
    double normalFractureTraction(size_t ix) const;

    // one nonlinear iteration of fully coupled system.  Returns 'true' if converged
  bool fullSystemIteration(const double tol,const int iteration);

    void assembleFractureMatrix() const;
    std::vector<double> stressIntensityK1() const;
    int numWellEquations() const;

    std::vector<CellRef> well_source_cellref_; // references to well cells in the fully resolved TriMesh
    std::vector<int> well_source_; // indices to well cells
    // Reservoir cells each fracture cell intersects; rebuilt when the fracture grid changes.
    struct ReservoirMap
    {
        std::vector<int> cells; // the cell holding the fracture cell's centroid
        std::vector<std::vector<int>> allCells;
        std::vector<std::vector<double>> allAreas;
        std::vector<std::vector<Dune::FieldVector<double, 3>>> allCenters;
        std::vector<double> dist;
    };
    // Reservoir state per fracture cell; refreshed from flow and mechanics.
    struct ReservoirSamples
    {
        std::vector<double> perm, cstress, mobility, waterMobility, density, cellZ, pressure;
        std::vector<Dune::FieldVector<double, 6>> stress;
    };
    // Flow state of an intersected reservoir cell that holds no fracture-cell centroid.
    struct ReservoirCellValues
    {
        double mobility{}, density{}, cellZ{}, pressure{};
    };
    ReservoirMap reservoir_map_;
    ReservoirSamples reservoir_;
    std::map<int, ReservoirCellValues> off_centroid_cells_;
    using WaterPropertyEvaluator = std::function<std::pair<CellFluidProperty, CellFluidProperty>(size_t, double)>;
    WaterPropertyEvaluator fracture_water_property_evaluator_;
    std::optional<centre_stress::DeclaredFaults> stress_faults_;
    std::unordered_map<int, centre_stress::CellSides> stress_sides_;

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
    // width exposed to the coupling last outer round (solver.coupling_state_relaxation);
    // valid only within one timestep's outer iterations
    Dune::BlockVector<Dune::FieldVector<double, 1>> width_round_prev_;
    bool width_round_valid_{false};
    // contact set at the end of the previous coupling round in this timestep
    // (solver.propagate_on_stable_contact): growth is only allowed once the
    // open/closed set has stopped changing between rounds
    // fracture area and perf pressure at the timestep checkpoint (moveForwardInTime);
    // used by the per-step growth guard and the onset dt-hold
    double area_step_start_{-1.0};
    double perf_pressure_step_start_{-1.0};
    double volume_step_start_{-1.0}; // V0 for volume-paced propagation
    // coupling-rate dt controller (opt-in): cap on the next flow step, set at the
    // step checkpoint from how fast perf pressure / area moved in the last step
    double coupling_dt_cap_{-1.0};
    // stress-rate dt controller (opt-in, solver.stress_dt_target bar/step): the
    // thermal closure-stress transient is smooth, so a proportional cap on the
    // per-step traction change paces dt through cooling without the episode
    // problem of the coupling-rate controller
    double traction_step_start_{std::numeric_limits<double>::quiet_NaN()};
    double stress_dt_cap_{-1.0};
    bool well_control_is_rate_{true}; // active well constraint (from the well state)
    double well_target_bhp_{-1.0};    // active BHP (bhp_well constraint target)
    // control.type=="well": the BC resolved at the start of the current solve.
    // Latched so the system size (well DOF) is stable within a solve; transitions
    // between solves resize the persistent state.
    std::string effective_control_latched_;
    int coupling_quiet_steps_{0};

    // transmissibilities
    using Htrans = std::tuple<size_t, size_t, double, double>;
    std::vector<Htrans> htrans_;
    std::vector<std::tuple<int,double>> perfinj_;
    std::vector<int> well_perf_cells_; // see setWellPerfCells
    //! Per-cell scaled FB complementarity residual |phi|/(1+|a|+|b|) from the last
    //! assembled iterate; empty unless closed_cell_policy=fischer_burmeister.
    //! Lets propagation veto individual untrustworthy cells instead of all growth.
    std::vector<double> fb_cell_residual_;
    double current_dt_{-1.0}; // flow timestep (s), see setTimeStep
    bool external_pressure_{false}; // see setExternalPressure
    std::vector<char> external_cell_mask_; // per cell: pressure pinned from outside
    std::vector<int> flow_state_donor_; // see flowStateDonor()
    bool flow_state_remapped_{false};
    //! Nonlinear iteration at which each cell last flipped open/closed in the
    //! current solve; the binary active set's analogue of fb_cell_residual_ for
    //! the propagation veto (recent stability, not lifetime stability).
    std::vector<int> cell_last_toggle_iter_;
    int current_solve_iter_{0}; //!< last nonlinear iteration index of the solve
    double perf_pressure_;
    std::vector<double> leakof_;
    std::vector<double> leakoff_mobility_;

    mutable std::unique_ptr<Matrix> pressure_matrix_;
    mutable std::unique_ptr<Matrix> coupling_matrix_; // will be updated by `fullSystemIteration`

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
    mutable std::unique_ptr<DynamicMatrix> fracture_matrix_;

    // function ensuring that the fracture matrix exists, and returning a reference to it
    DynamicMatrix& fractureMatrix() const;

    double E_;
    double nu_;
    double min_width_; // minimum width of fracture, used for convergence criterion
    double gravity_{Opm::unit::gravity}; // gravity acceleration, used for leakoff calculations
    std::vector<double> fracture_dgh_; // gravity contribution to fracture pressure, used for leakoff calculations
    PropertyTree prm_;
    double total_WI_well_{0.0}; // total well index for the well, used for leakoff calculations
    // for coupled solve
    using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                     Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
    using SystemMatrixRef = Dune::MultiTypeBlockMatrixRef<Dune::MultiTypeBlockVectorRef<FMatrix, SMatrix>,
                                                     Dune::MultiTypeBlockVectorRef<SMatrix, SMatrix>>;
    std::unique_ptr<SystemMatrix> S_;
    std::unique_ptr<FMatrix> A_; // system matrix without coupling terms
    std::unique_ptr<SMatrix> I_; // system matrix without coupling terms
    std::vector<std::vector<RuntimePerforation>> well_indices_; // for each perf, the corresponding well index info
    double max_flow_time_step_{-1.0};
    double well_rate_;
    double total_wellindex_;
    double wi_dz_{0.0};
    double wi_respress_{0.0};
    double well_ref_depth_{0.0};
    double perf_ref_depth_{0.0};
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
    std::vector<int> closed_cells_; // indices of cells that are closed (zero width) in the current iteration
    // per-cell contact state changes within the current solve, used by the
    // opt-in anti-chatter guard (solver.toggle_guard_mode = "chatter")
    std::vector<int> cell_flip_counts_;
    //using AbstractPreconditioner = Dune::PreconditionerWithUpdate<VectorHP, VectorHP>;
    using AbstractPreconditioner = FractureMechanicsPreconditioner;
    std::unique_ptr<AbstractPreconditioner> frac_flow_precond_;
    std::unique_ptr<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>> S_linop_;
    using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;
    std::unique_ptr<LinearSolverBase> psolver_;
    FractureSolveStats last_solve_stats_{};

};
} // namespace Opm

#include "Fracture_impl.hpp"
#endif
