#pragma once

#include "Fracture.hpp"
#include "FractureWell.hpp"
#include "GeometryHelpers.hpp"

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFractureSeeds.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Opm {
    class RuntimePerforation;

    class ScheduleState;

  template <typename Scalar,typename IndexTraits>
    class WellState;
} // namespace Opm

namespace Opm{
    PropertyTree makeDefaultFractureParam(bool rate_control = false);
class FractureModel{
    //using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
public:
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned int,2>;
    using EntitySeed = Dune::CpGrid::Codim<0>::Entity::EntitySeed;
    template<class Grid>
    FractureModel(const Grid& grid,
                  const std::vector<Well>& wells,
                  const PropertyTree&);

    /// Initialise fracture objects.
    ///
    /// Initialisation method selected by the property tree's "type" node,
    /// which defaults to the "well_seed" method (WSEED keyword).
    ///
    /// \param[in] sched Dynamic objects in current run, especially
    /// including well fracturing seed points and fracturing plane normal
    /// vectors in addition to all current well objects.
    void addFractures(const ScheduleState& sched);

     void updateFractureReservoirCells(const Dune::CpGrid& cpgrid);
  
  //void updateFractureReservoirCells(const Dune::CpGrid& cpgrid);
  
  // template<class Grid>
  // void updateFractureReservoirCells(const Grid& cpgrid)
  //   {
  //       external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
  //       external::buildBoundingBoxTree(cellSearchTree, cpgrid);
  //       for (auto& well_fracture : well_fractures_) {
  //           for (auto& fracture : well_fracture) {
  //               fracture.updateReservoirCells(cellSearchTree, cpgrid);
  //           }
  //       }
  //   }

    void addWell(const std::string& name,
                 const std::vector<FractureWell::Connection>& conns,
                 const std::vector<Point3D>& points,
                 const std::vector<std::array<unsigned,2>>& segments);

    void write(int ReportStep = -1) const;
    void writemulti(double time) const;

    template<class wseed_collection>
    void updateActive(const wseed_collection& current_wseed);

    template <class TypeTag, class Simulator>
    void solve(const Simulator& simulator);

    template <class TypeTag, class Simulator>
    void updateFilterCakeProperties(const Simulator& simulator);

    template <class TypeTag, class Simulator>
    double maxFlowTimeStep(const Simulator& simulator) const;

    // Per-step growth guard: first violation message across all active
    // fractures (empty = none). Not MPI-reduced; caller decides.
    std::string growthGuardViolation() const;
    //! true when any active fracture's last solve ended unconverged
    bool anyLastSolveUnconverged() const;
    void writeIterationSnapshots(int step, int round, const std::string& tag) const;

        void updateReservoirProperties(); // for testing without simulator
        void initFractureStates();

        template <class TypeTag, class Simulator>
        void initReservoirProperties(const Simulator& simulator);
        template <class TypeTag, class Simulator>
        void updateReservoirAndWellProperties(const Simulator& simulator);
        template <class TypeTag, class Simulator>
        void resetFractures(const Simulator& simulator);
        void moveForwardInTime(double dt_last = -1.0);

        template <class TypeTag, class Simulator>
        void updateReservoirWellProperties(const Simulator& simulator);


        std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname) const;

        template <typename Scalar, typename IndexTraits>
        void assignGeoMechWellState(WellState<Scalar, IndexTraits> & wellState) const;

        bool addPertsToSchedule();
        // probably this should be collected in one loop since all do full loop over fracture ++ well
        Dune::FieldVector<double, 6> stress(Dune::FieldVector<double, 3> obs) const;
        Dune::FieldVector<double, 6> strain(Dune::FieldVector<double, 3> obs) const;
        Dune::FieldVector<double, 3> disp(Dune::FieldVector<double, 3> obs) const;
        Opm::PropertyTree& getParam();
        const FractureSolveStats& lastSolveStats() const;
        const FractureSolveStats& totalSolveStats() const;
        static Opm::DeferredLogger fractureLogger;

        //! The fractures of each well, in well order.
        const std::vector<std::vector<Fracture>>& wellFractures() const
        { return well_fractures_; }

        /*!
         * \brief Bring every fracture's leak-off in step with its grid.
         *
         * Called before the embedded flow representation reads the fractures, so that a
         * propagation attempt rolled back after the last leak-off update does not leave
         * a stale array behind.  Returns whether every fracture is now consistent.
         */
        bool ensureFlowDescriptionCurrent()
        {
            bool ok = true;
            for (auto& well_fracture : well_fractures_) {
                for (auto& fracture : well_fracture) {
                    ok = fracture.ensureLeakoffCurrent() && ok;
                }
            }
            return ok;
        }

    private:
        bool vtkwritewells_ = false; // write wells to VTK files
        template <class TypeTag, class Simulator>
        void updateReservoirProperties(const Simulator& simulator);

        template <class TypeTag, class Simulator>
        void updateWellProperties(const Simulator& simulator); // update all well related properties

        std::vector<FractureWell> wells_;
        std::vector<std::vector<Fracture>> well_fractures_;
        Opm::PropertyTree prm_;
        FractureSolveStats last_solve_stats_{};
        FractureSolveStats total_solve_stats_{};
        external::cvf::ref<external::cvf::BoundingBoxTree> cell_search_tree_;
        std::vector<EntitySeed> cell_seeds_;
        /// Initialise fractures perpendicularly to each reservoir connection.
        void addFracturesPerpWell();
        void addFracturesTensile();

        /// Initialise fractures in each seed identified in the WSEED keyword.
        ///
        /// \param[in] sched Dynamic objects in current run, especially
        /// including well fracturing seed points and fracturing plane normal
        /// vectors in addition to all current well objects.
        void addFracturesWellSeed(const ScheduleState& sched);

    };
}
#include "FractureModel_impl.hpp"
