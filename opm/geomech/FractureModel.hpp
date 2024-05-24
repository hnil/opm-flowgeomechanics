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
                  const Opm::PropertyTree& param
        ):
        prm_(param)
    {
        //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        //CartesianIndexMapper cartmapper(grid);
        //const std::array<int, dimension>
        //auto cartdims = cartmapper.cartesianDimensions();
        GeometryHelper geomhelp(grid);
        for(int wellIdx=0; wellIdx < wells.size(); ++wellIdx){
            std::vector<Point3D> vertices;
            std::vector<Segment> wellsegment;
            auto well = wells[wellIdx];
            if(true){//well.isInjector()){// only look at injectors
            std::vector<int> cells;
            // should be done property with welltraj
            for(const auto& connection : well.getConnections()){
                int global_index = connection.global_index();
                int cell_index = geomhelp.compressedIndex(global_index);
                Point3D vertex = geomhelp.centroid(cell_index);
                vertices.push_back(vertex);
                cells.push_back(cell_index);
            }
            Point3D refpoint(vertices[0]);
            if(well.hasRefDepth()){
                if(std::abs(well.getRefDepth()-refpoint[2])< 10){
                    refpoint[2] = well.getRefDepth()-10;//avoid zero
                }else{
                    refpoint[2] = well.getRefDepth();
                }
            }else{
                refpoint[2] -= 10;
            }
            vertices.insert(vertices.begin(),refpoint);
            std::vector<Segment> segments;
            for(size_t i=0; i < cells.size(); ++i){
                segments.push_back(Segment({i,i+1}));
            }
            // NB should add gri cells
            this->addWell(well.name(),vertices, segments);
            }
        }
        this->addFractures();
        external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
        external::buildBoundingBoxTree(cellSearchTree, eclgrid);
        for(auto& well_fracture : well_fractures_){
            for(auto& fracture : well_fracture){
                fracture.updateReservoirCells(cellSearchTree, grid);
            }
        }

    }
    void addFractures();
    void addWell(std::string name, const std::vector<Point3D>& points,const std::vector<std::array<unsigned,2>>& segments );
    void write(int ReportStep = -1) const;
    void solve();
    void updateReservoirProperties();
private:
    std::vector<FractureWell> wells_;
    std::vector<std::vector<Fracture>> well_fractures_;
    Opm::PropertyTree prm_;
};
}
