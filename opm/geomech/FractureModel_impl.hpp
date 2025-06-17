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
        vtkwritewells_= prm_.get<bool>("vtkwritewells", false);
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
            if(cells.size() == 0){
                std::cerr << "Warning: No connections found for well " << well.name() << std::endl;
                continue; // skip if no connections
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
    template <class TypeTag, class Simulator>
    void FractureModel::updateWellProperties(const Simulator& simulator)
    {
        for (size_t i=0; i < wells_.size(); ++i) {
            // TO DO set wells to active even without fractures
            {
                
            }  
            for (auto& fracture : well_fractures_[i]){
                // do update wells
                // set well properties
                WellInfo wellinfo = fracture.wellInfo();
                // need to be double checked how to assosiate correct perforation/segment
                int perf_index_frac = wellinfo.perf;
                int cell_index_frac = wellinfo.well_cell;
                // see if well exist
                fracture.setActive(false);
                auto well_index = simulator.problem().wellModel().wellState().index(wellinfo.name);
                if(!well_index.has_value()){
                    wells_[i].setActive(false);
                    fracture.setActive(false);
                    continue;
                }
                const auto& wellstate = simulator.problem().wellModel().wellState().well(*well_index);
                // check if well is open
                if(wellstate.status != Opm::WellStatus::OPEN) {
                    wells_[i].setActive(false);
                    std::cerr << "Warning: Well " << wellinfo.name << " is not open, skipping update." << std::endl;
                    fracture.setActive(false);
                    continue; // skip if not open
                }
                // well is active
                wells_[i].setActive(true);
                // get well perforation
                const auto& perf_data = wellstate.perf_data;
                auto it = std::find(perf_data.cell_index.begin(),
                                    perf_data.cell_index.end(),
                                    cell_index_frac);
                // check if perforation exists
                if(it == perf_data.cell_index.end()) {
                    std::cerr << "Warning: Could not find perforation for well " << wellinfo.name
                              << " in cell index " << cell_index_frac << std::endl;
                    fracture.setActive(false);
                    wells_[i].setPerfActive(perf_index_frac,false);          
                    continue; // skip if not found
                }
                int perf_index = it- perf_data.cell_index.begin();
                double perf_pressure = perf_data.pressure[perf_index];
                std::cout << "Perf index flow " << perf_index << " fracture " << perf_index_frac << " pressure " << perf_pressure << std::endl;
                fracture.setPerfPressure(perf_pressure);
                wells_[i].setPerfPressure(perf_index_frac, perf_pressure);
                fracture.setActive(true);
                wells_[i].setPerfActive(perf_index_frac,true);
                //NB do we need some rates? need to be summed over "peforations of the fractures"
            }
        }
    }
    

}
