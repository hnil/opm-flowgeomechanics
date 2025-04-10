#pragma once
#include <opm/simulators/linalg/PropertyTree.hpp>
namespace Opm{
    template<class Grid>
    FractureModel::FractureModel(const Grid& grid,
                  const std::vector<Opm::Well>& wells,
                  const Opm::PropertyTree& param)
                  :
        prm_(param)
    {
        //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        //CartesianIndexMapper cartmapper(grid);
        //const std::array<int, dimension>
        //auto cartdims = cartmapper.cartesianDimensions();
        GeometryHelper geomhelp(grid);
	// NB: need to be carefull in parallel
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
            //NB unsinged so grid can not be to big
            for(unsigned i=0; i < unsigned(cells.size()); ++i){
                segments.push_back(Segment({i,i+1}));
            }
            // NB should add gri cells
            this->addWell(well.name(),vertices, segments, cells);
            }
        }

        external::buildBoundingBoxTree(cell_search_tree_, grid);
    }
}
