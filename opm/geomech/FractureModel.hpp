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
#include <iostream>
#include <sstream>
#include <string>

namespace Opm {
    class RuntimePerforation;
} // namespace Opm

namespace Opm{
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
    void addFractures();

    void updateFractureReservoirCells()
    {
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cell_search_tree_);
            }
        }
    }

    template <class Grid>
    void addFractures(const Grid& grid, const Opm::EclipseGrid& eclgrid);

    template <class Grid>
    void addDefaultsWells(const Grid& grid, const Opm::EclipseGrid& eclgrid);


  
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
          fractures[j].solve<TypeTag, Simulator>(cell_search_tree_, simulator);
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
    void updateReservoirProperties(const Simulator& simulator)
    {
        for (size_t i=0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]){
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator);
                // set well properties
                WellInfo wellinfo = fracture.wellInfo();
                // need to be double checked how to assosiate correct perforation/segment
                int perf_index_frac = wellinfo.perf;
                int cell_index_frac = wellinfo.well_cell;
                std::size_t well_index = simulator.problem().wellModel().wellState().index(wellinfo.name).value();
                const auto& wellstate = simulator.problem().wellModel().wellState().well(well_index);
                const auto& perf_data = wellstate.perf_data;
                auto it = std::find(perf_data.cell_index.begin(),
                                    perf_data.cell_index.end(),
                                    cell_index_frac);
                int perf_index = it- perf_data.cell_index.begin();
                std::cout << "Perf index flow " << perf_index << " fracture " << perf_index_frac << std::endl;
                double perf_pressure = perf_data.pressure[perf_index];
                fracture.setPerfPressure(perf_pressure);
                wells_[i].setPerfPressure(perf_index_frac, perf_pressure);
                //NB do we need some rates? need to be summed over "peforations of the fractures"
            }
        }
    }

    template<class TypeTag, class Simulator>
    void updateReservoirWellProperties(const Simulator& simulator) {
        for(auto& well : wells_){
             well.updateReservoirProperties<TypeTag,Simulator>(simulator);
        }
    }
    std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname) const;
    bool addPertsToSchedule(){return prm_.get<bool>("addperfs_to_schedule");}
    // probably this should be collected in one loop since all do full loop over fracture ++ well
    Dune::FieldVector<double,6> stress(Dune::FieldVector<double,3> obs) const;
    Dune::FieldVector<double,6> strain(Dune::FieldVector<double,3> obs) const;
    Dune::FieldVector<double,3> disp(Dune::FieldVector<double,3> obs) const;
private:
    std::vector<FractureWell> wells_;
    std::vector<std::vector<Fracture>> well_fractures_;
    Opm::PropertyTree prm_;
    external::cvf::ref<external::cvf::BoundingBoxTree> cell_search_tree_;
};
}
#include "FractureModel_impl.hpp"
