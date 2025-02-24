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
      using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
      CartesianIndexMapper cartmapper(grid);
      const auto& cartdim = cartmapper.cartesianDimensions();
      std::vector<int> cart(cartdim[0]*cartdim[1]*cartdim[2], -1);
      using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView>;
      ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
      for (const auto& elem : elements(grid.leafGridView())) {
	int eIdx = mapper.index(elem);
	cart[ cartmapper.cartesianIndex( eIdx ) ] = eIdx;
      }
      const auto& config = prm_.get_child("config");
      std::string type = config.get<std::string>("type");
      if( type  == "well_seed"){
	const auto cell_ijk = config.get_child_items_as_vector<int>("cell_ijk");
	std::array<int, 3> ijk;
	   assert(cell_ijk.has_value() && cell_ijk->size() == 3);
	for(int i=0; i < 3; ++i){
	  // map to 0 based
	  ijk[i] = (*cell_ijk)[i]-1;
	}
	int cartIdx =  ijk[0]+(cartdim[0])*ijk[1]+(cartdim[0]*cartdim[1])*ijk[2];//NB assumes ijk ordering
	int reservoir_cell = cart[cartIdx];
	assert(cartIdx>=0);
	prm_.put("config.cell",reservoir_cell);
      }
      
     
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
  
    template<class Grid>
    void FractureModel::addDefaultsWells(const Grid& grid, const Opm::EclipseGrid& eclgrid) {
      this->addFractures();
      this->updateFractureReservoirCells();
    }
}
