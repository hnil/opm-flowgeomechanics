#pragma once
#include "Fracture.hpp"
#include "FractureWell.hpp"
#include "GeometryHelpers.hpp"

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

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

    template <typename Scalar>
    class WellState;
} // namespace Opm

namespace Opm{
  Opm::PropertyTree makeDefaultFractureParam();
class FractureModel{
    //using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
public:
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned int,2>;
    template<class Grid>
    FractureModel(const Grid& grid,
                  const std::vector<Opm::Well>& wells,
                   const Opm::PropertyTree&
        );

    /// Initialise fracture objects.
    ///
    /// Initialisation method selected by the property tree's "type" node,
    /// which defaults to the "well_seed" method (WSEED keyword).
    ///
    /// \param[in] sched Dynamic objects in current run, especially
    /// including well fracturing seed points and fracturing plane normal
    /// vectors in addition to all current well objects.
    void addFractures(const ScheduleState& sched);

    void updateFractureReservoirCells()
    {
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cell_search_tree_);
            }
        }
    }
  
  //void updateFractureReservoirCells(const Dune::CpGrid& cpgrid)
  template<class Grid>
  void updateFractureReservoirCells(const Grid& cpgrid)
    {
        external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
        external::buildBoundingBoxTree(cellSearchTree, cpgrid);
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cellSearchTree, cpgrid);
            }
        }
    }

  void addWell(std::string name,
	       const std::vector<Point3D>& points,
	       const std::vector<std::array<unsigned,2>>& segments,
	       const std::vector<int>& reservoir_cells );

    void write(int ReportStep = -1) const;
    void writemulti(double time) const;

  template <class TypeTag, class Simulator> void solve(const Simulator& simulator) {
      for(size_t i=0; i < wells_.size(); ++i){
        std::vector<Fracture>& fractures = well_fractures_[i];
        for(size_t j=0; j < fractures.size(); ++j){
          if(fractures[j].isActive()){  
            std::cout << "Solving fracture " << fractures[j].name() << std::endl;
            fractures[j].solve<TypeTag, Simulator>(cell_search_tree_, simulator);
          }
        }
      }
    }
    void updateReservoirProperties();
    void initFractureStates();
    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator)
    {
        for (size_t i=0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]) {
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator, true);
            }
        }
    }
    template <class TypeTag, class Simulator>   
    void updateReservoirAndWellProperties(const Simulator& simulator){
        this->updateReservoirProperties<TypeTag, Simulator>(simulator);
        this->updateWellProperties<TypeTag, Simulator>(simulator);
    }
    

    template<class TypeTag, class Simulator>
    void updateReservoirWellProperties(const Simulator& simulator) {
        for(auto& well : wells_){
             well.updateReservoirProperties<TypeTag,Simulator>(simulator);
        }
    }
    
    
    std::vector<RuntimePerforation>
    getExtraWellIndices(const std::string& wellname) const;

    template <typename Scalar>
    void assignGeomechWellState(WellState<Scalar>& wellState) const;

    bool addPertsToSchedule(){return prm_.get<bool>("addperfs_to_schedule");}
    // probably this should be collected in one loop since all do full loop over fracture ++ well
    Dune::FieldVector<double,6> stress(Dune::FieldVector<double,3> obs) const;
    Dune::FieldVector<double,6> strain(Dune::FieldVector<double,3> obs) const;
    Dune::FieldVector<double,3> disp(Dune::FieldVector<double,3> obs) const;
    Opm::PropertyTree& getParam(){return prm_;}
private:
    bool vtkwritewells_ = false; // write wells to VTK files
    template <class TypeTag, class Simulator>   
    void updateReservoirProperties(const Simulator& simulator){
        for (size_t i=0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]){
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator);
            }
        }
    }

    template <class TypeTag, class Simulator>
    void updateWellProperties(const Simulator& simulator);//update all well related properties 

    std::vector<FractureWell> wells_;
    std::vector<std::vector<Fracture>> well_fractures_;
    Opm::PropertyTree prm_;
    external::cvf::ref<external::cvf::BoundingBoxTree> cell_search_tree_;

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
