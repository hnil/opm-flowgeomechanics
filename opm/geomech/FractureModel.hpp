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
namespace Opm{
class FractureModel{
    //using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
public:
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned int,2>;
    template<class Grid>
    FractureModel(const Grid& grid,
                  const std::vector<Opm::Well>& wells,
                  const Opm::EclipseGrid& eclgrid,
                  const Opm::PropertyTree& param,
                  const bool default_fractures
        );
    void addFractures();

    template <class Grid>
    void updateFractureReservoirCells(const Grid& grid, const Opm::EclipseGrid& eclgrid)
    {
        external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
        external::buildBoundingBoxTree(cellSearchTree, eclgrid);
        for (auto& well_fracture : well_fractures_) {
            for (auto& fracture : well_fracture) {
                fracture.updateReservoirCells(cellSearchTree, grid);
            }
        }
    }
    void addWell(std::string name, const std::vector<Point3D>& points,const std::vector<std::array<unsigned,2>>& segments,const std::vector<int>& reservoir_cells );
    void write(int ReportStep = -1) const;
    void writemulti(double time) const;
    void solve();
    void updateReservoirProperties();

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator)
    {
        for (size_t i=0; i < wells_.size(); ++i) {
            for (auto& fracture : well_fractures_[i]){
                fracture.updateReservoirProperties<TypeTag, Simulator>(simulator);
            }
        }
    }

    template<class Problem>
    void updateReservoirWellProperties(const Problem& problem) {
        for(auto& well : wells_){
             well.updateReservoirStress(problem);
        }
    }
private:
    std::vector<FractureWell> wells_;
    std::vector<std::vector<Fracture>> well_fractures_;
    Opm::PropertyTree prm_;
};
}
#include "FractureModel_impl.hpp"
