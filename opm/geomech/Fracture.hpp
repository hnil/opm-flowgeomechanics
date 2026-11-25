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
#include <algorithm>
#include <limits>
#include <optional>
#include <array>
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
#include <set>
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
#include <opm/geomech/FractureCoupledLinearSolver.hpp>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>
// for linear solve
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/InjectedPhase.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <dune/istl/matrixmarket.hh>

#include <opm/geomech/DiscreteDisplacement.hpp>
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
    // Dense mechanics factorisations, Woodbury-updated iterations and their time split.
    int mech_factorizations = 0;
    double mech_factor_seconds = 0.0;
    double linear_solve_seconds = 0.0;
    double precond_apply_seconds = 0.0;

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
      mech_factorizations += other.mech_factorizations;
      mech_factor_seconds += other.mech_factor_seconds;
      linear_solve_seconds += other.linear_solve_seconds;
      precond_apply_seconds += other.precond_apply_seconds;
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
    //! Cell search on a separate mechanics grid, for fracture cells outside the flow grid.
    struct MechLookup
    {
        const external::cvf::ref<external::cvf::BoundingBoxTree>* tree = nullptr;
        const Dune::CpGrid* grid = nullptr;
        const std::vector<EntitySeed>* seeds = nullptr;
    };
    void setMechLookup(const MechLookup& lookup) { mech_lookup_ = lookup; }
    template <class TypeTag, class Simulator>
    void updateFilterCakePropertiesPost(const Simulator& simulator);
    // solver related
    void updateReservoirProperties();
    bool removeCells();

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator, bool init_constant_vals = false);

    // Max |change| of water pressure and temperature in the fracture's reservoir cells
    // since the last call with record=true; infinite if nothing was recorded on this mesh.
    template <class TypeTag, class Simulator>
    std::pair<double, double> drivingStateChange(const Simulator& simulator, bool record);
    // > 0 replaces solver.area_change_fac; <= 0 restores it.
    void setAreaChangeFactorOverride(double fac) { area_change_fac_override_ = fac; }

    //! Deposit filtrate on the fracture faces over a step of length dt. The
    //! concentration is per surface volume of injected fluid if surface_conc_basis,
    //! else per reservoir volume (WINJFCNC item 3).
    void updateFilterCakeProps(const Opm::Connection& connection,
                               const Opm::SingleWellState<double,IndexTraits>& wellstate,
                               double dt,
                               bool surface_conc_basis = true);
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
    //! well indices as of the last converged solve; empty before the first one
    const std::vector<RuntimePerforation>& convergedWellIndices() const { return converged_well_indices_; }
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
    //! Elements outside this rank's cells but inside another rank's (last updateReservoirCells).
    int numPartitionClipped() const { return partition_clipped_; }

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

    //! Canonical phase index (water, oil or gas) of the fluid the well injects into,
    //! and that fills, the fracture.  Water until the well is known as an injector of
    //! another single phase.
    unsigned injectedPhaseIdx() const { return injected_phase_idx_; }
    bool injectsWater() const { return injected_phase_idx_ == IndexTraits::waterPhaseIdx; }

    //! Mobility per fracture cell used for flow inside the fracture (cubic law faces).
    std::vector<double> fractureFlowMobilities() const;

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
    //! Well path as a polyline in connection order (solver.well_source_mode=trajectory).
    void setWellPath(std::vector<Point3D> points) { well_path_ = std::move(points); }
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
    // Last solve stopped on the front criterion (not the area cap or volume
    // budget) and changed the area by less than rel_tol.
    bool growthSettled(double rel_tol) const;
    // (area now - area at the step checkpoint) / checkpoint area; 0 before a checkpoint.
    double relativeStepGrowth() const;
    // The step ending now was accepted as a burst: the area dt controller holds dt.
    void markBurstStep() { burst_step_ = true; }
    // Largest K1/K1c on the front the last solve left (vetoed cells score -1).
    double lastFrontScore() const { return last_front_score_; }
    // |well rate - rate at the step checkpoint| / |rate at the checkpoint|; 0 before a checkpoint.
    double relativeStepRateChange() const;
    // Onset dt-hold: true while the fracture is still an unopened seed and the
    // perf pressure is rising (solver.onset_dt_limit > 0 arms it).
    bool onsetHoldActive() const;
    double filterCakeVolume() const;
    void resetFracture();
    // after resetFracture and the reservoir-cell update: rebuild well sources and pressure matrix if the grid changed
    void finishReset();
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
   RuntimePerforation perforation_(const CellExchange& c) const;
   std::vector<RuntimePerforation> wiLegacy_(const std::vector<CellExchange>& cells,
                                             double inj_press) const;
   std::vector<RuntimePerforation> wiConductivity_(const std::vector<CellExchange>& cells,
                                                   double inj_press, const WiOptions& opt) const;
   std::vector<RuntimePerforation> wiFracturePressure_(const std::vector<CellExchange>& cells,
                                                       double inj_press) const;
   std::vector<RuntimePerforation> wellIndices_() const;
   bool  expantionMax(const FractureProperties& fprop);
   bool capFlowTimeStep_(double area, double ref_area); // sets max_flow_time_step_; true if over area_change_fac

   // ---- trimesh propagation solve (solve()) ----
   struct PropagationOptions
   {
       int maxIter{};
       double tol{};
       int nlinVerbosity{};
       int maxNumCoarsening{};
       int numcellThreshold{};
       int bandRings{1};
       bool smoothBoundary{false};
       bool remapSolution{true};
       bool conservativeCake{false};
       bool majorityCoarsen{false}; // width remapped by volume
       bool levelConsistent{false}; // expansion ends on the last solved grid and level
       double forceLimit{0.0};
       bool conservative{false};
       bool convergedCellsOnly{false};
       double fbTol{1e-3};
       int stableIters{2};
       int verbosity{0};
       bool disable{false};
       bool openFrontOnly{false};
       double frontAdvanceFraction{0.0};
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
       int before_level{0};                      // trimesh level of the refs in before
       std::vector<std::vector<CellRef>> after;
       std::map<CellRef, int> boundary;          // trimesh boundary cell -> grid cell
       RegularTrimesh::CellVolumes width_volumes; // on the grid before the regrid, majorityCoarsen only
       double open_area_before{0.0};
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
   std::vector<char> initNewCellsFromNeighbours_(size_t num_press_cells);
   void initNewWidthFromElasticity_(const std::vector<char>& is_new);
   void remapState_(const RegridMaps& maps, int level, const PropagationOptions& opt,
                    const FilterCakeSnapshot& cake0);
   GridSolve solveOnCurrentGrid_(const PropagationOptions& opt);
   void reportGridSolve_(const GridSolve& solve, const PropagationOptions& opt);
   // K1/K1c per trimesh front cell; -1 where a veto forbids growth
   std::vector<double> propagationScores_(std::map<CellRef, int>& boundary, const GridSolve& solve,
                                          const PropagationOptions& opt) const;
   void logDtCaps_() const;
   void clearContactMemoryOnRegrid_();
   void vetoForceLimit_(std::vector<double>& score, const std::vector<int>& front,
                        const PropagationOptions& opt) const;
   void vetoUnconvergedCells_(std::vector<double>& score, const std::vector<int>& front,
                              const PropagationOptions& opt) const;
   void vetoClosedFront_(std::vector<double>& score, const std::vector<int>& front) const;
   bool removeNewZeroWithCells(RegularTrimesh& mesh,int cur_level,const RegularTrimesh& original_mesh) const;
   bool suppress_propagation_ = false;
   void composeFlowStateDonor(const std::vector<std::vector<CellRef>>& map1,
                              const std::vector<std::vector<CellRef>>& map2,
                              const int level, const int map1_level = 0);
   RegularTrimesh::CellVolumes filterCakeVolumesByCell() const;
   RegularTrimesh::CellVolumes volumesByCell_(const std::vector<double>& thickness) const;
   std::vector<double> remapFilterCakeConservative(const RegularTrimesh::CellVolumes& source,
                                                   const std::vector<std::vector<CellRef>>& source_map,
                                                   const int level,
                                                   const int source_level = 0) const;

   std::vector<double> redistribute_values(const std::vector<double>& values,
                                        const std::vector<std::vector<CellRef>>& map1,
                                        const std::vector<std::vector<CellRef>>& map2,
                                        const int level,
                                                     bool point_wise, const int map1_level = 0);

   void redistribute_values(Dune::BlockVector<Dune::FieldVector<double, 1>>& values,
                                        const std::vector<std::vector<CellRef>>& map1,
                                        const std::vector<std::vector<CellRef>>& map2,
                                        const int level,
                             bool point_wise, const int map1_level = 0);

    using ResVector = ::Opm::Fracture::Vector;
    using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;


    using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
    using FMatrix = Dune::DynamicMatrix<double>;                       // full matrix

    size_t numFractureCells() const { return grid_->leafGridView().size(0); }
    std::vector<int> identify_closed(const FMatrix& A,
                                     const VectorHP& x,
                                     const ResVector& rhs);
    static ResVector contactForce(const FMatrix& A,
                                  const ResVector& w,
                                  const ResVector& p,
                                  const ResVector& traction);

    // Inner Newton settings, read once per grid solve.
    struct InnerSolveOptions
    {
        int verbosity{0};
        bool isFb{false}; // closed_cell_policy fischer_burmeister; otherwise sticky
        double closeForceTol{0.0};
        double reopenForceTol{0.0};
        bool fbScaledConvergence{false};
        bool fbLineSearch{false};
        double fbLineSearchMin{0.0625};
        bool fbAcceptConverged{false};
        bool fbNonmonotone{false};
        bool fbClampBeforeLineSearch{false};
        double fbPhiAbsTol{0.0}; // [Pa]; |phi| below it passes regardless of the cell's scale
        bool fbFlowMerit{false};
        // > 0: coupling_matrix_ kept only once the open set is this many iterations stable (and, with FB, phi < 10 tol)
        int keepCouplingAfterStableIters{0};
        bool fbTrace{false};
        double flowTolerance{-1.0}; // < 0: same as the mechanics tolerance
        bool dropCoupling{true};
        bool storageTerm{false};
        bool storageJacobian{false}; // storage d/dw in C also when the rest of C is dropped
        FractureCoupledLinearSolver::Options linear;
        int dumpMinNonlinIter{0};
        bool writeCoupledSystem{false};
        bool dumpOnFailure{true};
        double maxDwidth{0.0};
        static InnerSolveOptions read(const PropertyTree& prm);
    };
    struct ContactUpdate
    {
        std::vector<int> closed;
        bool changed{false};
    };
    // Fischer-Burmeister residual phi and its per-cell scale 1+|a|+|b|; empty otherwise.
    struct FbRows
    {
        std::vector<double> phi;
        std::vector<double> scale;
    };
    InnerSolveOptions inner_;
    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator);

    void resetWriters();
    Dune::BlockVector<Dune::FieldVector<double, 3>> all_slips() const;

    std::unique_ptr<Opm::RegularTrimesh> trimesh_; // experimental, implicitly defined grids

    std::vector<std::vector<CellRef>> grid_mesh_map_; // index mapping from cells in grid_ to trimesh_.
    int grid_mesh_level_ = 0; // level of grid_mesh_map_ refs; nonzero only inside a propagation
                                                      // NB: in general many-to-many
    std::vector<CellRef> grid_fixed_cells_; // cells kept fine when grid_ was built from trimesh_
    struct TriangulationArgs {
        int coarsen_levels = 1;
        bool smooth = true;
        int cellnum_threshold = 10;
        int band_rings = 1;
    };
    TriangulationArgs grid_tri_args_; // with grid_fixed_cells_: the createDuneGrid call grid_ came from
    std::optional<std::map<CellRef, int>> removal_boundary_map_; // boundary map of the last removeCells
    std::set<CellRef> removed_cells_; // trimesh_ cells dropped by removeCells since the last regrid
    bool reset_pending_ = false; // set by resetFracture, consumed by finishReset
    bool grid_changed_since_checkpoint_ = false; // grid_ changed since the checkpoint: reset must rebuild perfinj_, htrans_, pressure_matrix_
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
    void assemblePressureAndCouplingAD(const std::vector<int>& closed_cells);
    void addSource();
    void initPressureMatrix();
    bool halfTransOwnCell_() const;
    void restoreLeakoffOfIsolatedPieces_(const std::vector<double>& full_leakof);
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
    void updateFlowCoefficients_();
    ContactUpdate updateContactSet_(const VectorHP& x, const ResVector& traction, int nlin_iteration);
    FbRows assembleMechanicsRows_(const VectorHP& x,
                                  ResVector& mech_rhs,
                                  const std::vector<int>& closed_cells,
                                  bool first_iteration);
    void assembleFlowRows_(const std::vector<int>& closed_cells);
    void assembleSystemMatrix_(const std::vector<int>& closed_cells, bool first_iteration);
    VectorHP assembleResidual_(const VectorHP& x, const ResVector& mech_rhs, const FbRows& fb) const;
    bool isConverged_(const VectorHP& rhs, const FbRows& fb, double tol) const;
    void limitStep_(const VectorHP& x, VectorHP& dx, const ResVector& mech_rhs, const FbRows& fb,
                    double tol, double flow_res, int nlin_iteration,
                    const std::vector<int>& closed_cells);
    double flowMeritAt_(const VectorHP& xt, const std::vector<int>& closed_cells, double tol_flow);
    void applyUpdate_(VectorHP& x, const VectorHP& dx);

    void assembleFractureMatrix() const;
    std::vector<double> stressIntensityK1() const;
    // Output only: (flag, WI) per fracture cell; flag 2 = seed ring, 1 = other fed cell.
    std::pair<std::vector<double>, std::vector<double>> wellSourceFields() const;
    int numWellEquations() const;

    std::vector<CellRef> well_source_cellref_; // references to well cells in the fully resolved TriMesh
    std::vector<int> well_source_; // indices to well cells
    // Reservoir cells each fracture cell intersects; rebuilt when the fracture grid changes.
    struct ReservoirMap
    {
        std::vector<int> cells; // the cell holding the fracture cell's centroid
        std::vector<int> mechCells; // mechanics cell of a centroid outside the flow grid, else -1
        std::vector<std::vector<int>> allCells;
        std::vector<std::vector<double>> allAreas;
        std::vector<std::vector<Dune::FieldVector<double, 3>>> allCenters;
        std::vector<double> dist;
    };
    // Reservoir state per fracture cell; refreshed from flow and mechanics.
    struct ReservoirSamples
    {
        std::vector<double> perm, cstress, mobility, density, cellZ, pressure;
        // mobility of the injected phase in the reservoir cell (kr included)
        std::vector<double> injectedMobility;
        // 1/mu of the injected fluid filling the fracture (kr = 1)
        std::vector<double> fluidMobility;
        std::vector<Dune::FieldVector<double, 6>> stress;
    };
    // Flow state of an intersected reservoir cell that holds no fracture-cell centroid.
    struct ReservoirCellValues
    {
        double mobility{}, density{}, cellZ{}, pressure{};
    };
    ReservoirMap reservoir_map_;
    int partition_clipped_ = 0;
    MechLookup mech_lookup_;
    // Reservoir-cell mapping per fracture triangle, keyed on its exact corner coordinates
    // (solver.cache_cell_mapping); valid while the reservoir grid is fixed.
    struct TriCellMapping {
        int centre_cell;
        std::vector<int> cells;
        std::vector<double> areas;
        std::vector<Dune::FieldVector<double,3>> centers;
    };
    struct TriKeyHash {
        std::size_t operator()(const std::array<double, 9>& k) const noexcept;
    };
    std::unordered_map<std::array<double, 9>, TriCellMapping, TriKeyHash> tri_mapping_cache_;
    ReservoirSamples reservoir_;
    std::map<int, ReservoirCellValues> off_centroid_cells_;
    using FluidPropertyEvaluator = std::function<std::pair<CellFluidProperty, CellFluidProperty>(size_t, double)>;
    FluidPropertyEvaluator fracture_fluid_property_evaluator_;
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
    // area-growth dt controller (opt-in, solver.area_dt_target)
    double area_dt_cap_{-1.0};
    bool burst_step_{false};
    // outcome of the last solve, for solver.growth_stop=criterion
    bool last_growth_capped_{false};
    double last_solve_rel_growth_{0.0};
    double last_front_score_{-1.0};
    double rate_step_start_{std::numeric_limits<double>::quiet_NaN()};
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
    std::vector<Point3D> well_path_; // see setWellPath
    bool well_path_fallback_logged_ = false;
    std::vector<int> wellPathSourceCells() const;
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
    mutable double fb_prev_flow_res_ = 0.0; // flow residual at the previous FB iteration
    mutable double fb_last_step_ = 1.0; // line-search factor of the last FB step (trace only)
    std::vector<char> stage_open_; // w >= min_width per cell at the previous iteration
    int stage_stable_iters_ = 0;
    bool stage_coupling_on_ = false;
    long stage_iters_on_ = 0;
    long stage_iters_total_ = 0;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
    mutable std::unique_ptr<DynamicMatrix> fracture_matrix_;
    mutable ddm::KernelTable ddm_table_; // DDM kernel values on the trimesh lattice
    mutable std::uint64_t ddm_generation_{0}; // bumped at every assembly of fracture_matrix_
    mutable std::vector<double> mech_centres_; // xyz per cell, element-mapper order
    mutable std::uint64_t mech_centres_generation_{0};
    // mechanics block of the coupled system = diag(mech_alpha_) * fractureMatrix() + diag(mech_beta_)
    std::vector<double> mech_alpha_, mech_beta_;

    // function ensuring that the fracture matrix exists, and returning a reference to it
    DynamicMatrix& fractureMatrix() const;

    double E_;
    double nu_;
    double min_width_; // minimum width of fracture, used for convergence criterion
    double max_flow_width_{std::numeric_limits<double>::infinity()}; // cubic-law width cap
    double flowWidth(double w) const { return std::clamp(w, min_width_, max_flow_width_); }
    // cubic_sum: h^3 = max(w,0)^3 + floor^3, smooth across the floor
    bool flow_width_cubic_sum_{false};
    double flowWidthCubed(double w) const
    {
        if (!flow_width_cubic_sum_) {
            const double h = flowWidth(w);
            return h * h * h;
        }
        const double wp = std::max(w, 0.0);
        return std::min(wp * wp * wp + min_width_ * min_width_ * min_width_,
                        max_flow_width_ * max_flow_width_ * max_flow_width_);
    }
    bool closed_cell_no_leakoff_{false};
    bool closed_cell_leakoff_ramp_{false};
    double closed_cell_leakoff_width_{1e-6};
    // leakof_ with cells at or below closed_cell_leakoff_width_ zeroed (none) or scaled by w/width (ramp)
    double activeLeakof(std::size_t i) const;
    std::vector<double> activeLeakof() const;
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
    std::unique_ptr<SMatrix> I_; // system matrix without coupling terms
    std::vector<std::vector<RuntimePerforation>> well_indices_; // for each perf, the corresponding well index info
    double max_flow_time_step_{-1.0};
    double area_change_fac_override_{-1.0};
    std::vector<double> driving_pressure_;
    std::vector<double> driving_temperature_;
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
    std::vector<RuntimePerforation> converged_well_indices_;
    double density_{1000.0};// default to water
    double density_perf_{1000.0};
    double mobility_injected_perf_{1000.0};
    // b of the injected fluid in the well cell (surface volume per reservoir volume)
    double inv_b_perf_{1.0};
    // canonical index of the phase the well injects and the fracture is filled with
    unsigned injected_phase_idx_{IndexTraits::waterPhaseIdx};

    template <class FluidSystem, class Simulator>
    void updateInjectedPhase(const Simulator& simulator);
    // solver.fracture_fluid_mobility: whether fracture flow, well-to-fracture flow and
    // cake resistance use the injected fluid's 1/mu (kr = 1) instead of the legacy
    // reservoir mobilities
    bool useInjectedFluidMobility() const;
    double fractureFlowMobility(std::size_t cell) const;

    // handling preconditioner
    std::vector<int> closed_cells_; // indices of cells that are closed (zero width) in the current iteration
    // per-cell contact state changes within the current solve, used by the
    FractureCoupledLinearSolver coupledSolver_;
    FractureSolveStats last_solve_stats_{};

};
} // namespace Opm

#include "Fracture_impl.hpp"
#endif
