#pragma once

#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm {
    template<class Grid>
    FractureModel::FractureModel(const Grid&              grid,
                                 const std::vector<Well>& wells,
                                 const PropertyTree&      param)
        : prm_(param)
    {
        //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        //CartesianIndexMapper cartmapper(grid);
        //const std::array<int, dimension>
        //auto cartdims = cartmapper.cartesianDimensions();
        GeometryHelper geomhelp(grid);
	// NB: need to be carefull in parallel
        vtkwritewells_= prm_.get<bool>("vtkwritewells", false);

        for (const auto& well : wells) {
            std::vector<Point3D> vertices;

            if (true) {//well.isInjector()){// only look at injectors
                auto conns = std::vector<FractureWell::Connection>{};

                // should be done property with welltraj
                for (const auto& connection : well.getConnections()) {
                    const int cell_index = geomhelp
                        .compressedIndex(connection.global_index());

                    vertices.push_back(geomhelp.centroid(cell_index));

                    auto& conn = conns.emplace_back();
                    conn.cell = cell_index;
                    conn.segment = connection.segment();
                    conn.perf_range = connection.perf_range();
                }

                if (conns.empty()) {
                    std::cerr << "Warning: No connections found for well " << well.name() << std::endl;
                    continue;
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
                for (unsigned i=0; i < unsigned(conns.size()); ++i) {
                    segments.push_back(Segment({i,i+1}));
                }

                // NB should add gri cells
                this->addWell(well.name(), conns, vertices, segments);
            }
        }

        external::buildBoundingBoxTree(cell_search_tree_, cell_seeds_, grid);
    }

  template<class WellState>
    int findPerf(const WellState& wellstate, int cell_index_frac){
        const auto& perf_data = wellstate.perf_data;
        auto it = std::find(perf_data.cell_index.begin(),
                            perf_data.cell_index.end(),
                            cell_index_frac);
                            // check if perforation exists
        if(it == perf_data.cell_index.end()) {
            std::cerr << "Warning: Could not find perforation for well " << wellstate.name
                    << " in cell index " << cell_index_frac << std::endl;         
            return -1;
        }// skip if not found
                            
        int perf_index = it - perf_data.cell_index.begin();
        return perf_index;
    }

    template <class TypeTag, class Simulator>
    void FractureModel::updateWellProperties(const Simulator& simulator)
    {
       using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
       //const int water_index = 0;//FluidSystem::waterPhaseIdx
        for (size_t i=0; i < wells_.size(); ++i) {
            // TO DO set wells to active even without fractures
            double injection_rate = 0.0;
            double well_depth = 0.0;
            double total_wellindex = 0.0;
            double wi_dz = 0.0;
            double wi_respress = 0.0;
            std::vector<double> perf_depths;// well.perfDepth()[perf_index]; 
            {
              auto well_index = simulator.problem().wellModel().wellState().index(wells_[i].name());
                if(well_index.has_value()){
                    const auto& wellmodel = simulator.problem().wellModel();
                    const auto& wellstate = wellmodel.wellState().well(*well_index);
                    //const auto& reportStepIndex = simulator.episodeIndex();
                    const auto& wellstate_nupcol = wellmodel.nupcolWellState().well(*well_index);
                    const auto& wells = wellmodel.localNonshutWells();
                    const auto& well = wells[*well_index];
                    assert(well->name() == wells_[i].name());
                    // get cell index from well perforation data
                    // check if well is open
                    if(wellstate.status == ::Opm::WellStatus::OPEN) {
                      //for (const auto& perf : wellstate.perf_data) {
                      for (int perf_index=0; perf_index < wellstate.perf_data.cell_index.size(); ++perf_index) {    
                            //water_rate += perf.flux[TypeTag::FluidSystem    ::waterPhaseIdx];
                            //total_wellindex += perf.well_index;// nead to multiply with formation damange
                            //int  perf_index =  findPerf(wellstate, cell_idx);
                            ///if(perf_index == -1) continue; // skip if not found
                            const int cell_idx = wellstate.perf_data.cell_index[perf_index];
                            const auto& intQuants = simulator.model()
                                .intensiveQuantities(cell_idx, /*timeIdx=*/0);
                            using Scalar = double;
                            const auto trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(intQuants,cell_idx);
                            const auto effective_well_indexs = well->wellIndex(perf_index, 
                                                                              intQuants, 
                                                                              trans_mult, 
                                                                              wellstate_nupcol, 
                                                                              /*with_fracture*/ false);
                            
                             const auto effective_well_index   =  effective_well_indexs[FluidSystem::waterPhaseIdx];
                            double density = intQuants.fluidState().density(FluidSystem::waterPhaseIdx).value();
                            total_wellindex += effective_well_index;
                            double dzwell = (well->perfDepth()[perf_index]- well->refDepth());  
                            wi_dz += effective_well_index*dzwell;
                            double pressure = intQuants.fluidState().pressure(FluidSystem::waterPhaseIdx).value();
                            wi_respress += effective_well_index*pressure;

                        }
                        injection_rate = wellstate.reservoir_rates[FluidSystem::waterPhaseIdx];
                        well_depth = well->refDepth();
                        perf_depths = well->perfDepth();
                    }   
                 }  
            }
            assert(well_fractures_[i].size() < 2);// for now asser this if not better "well model is need"
            for (auto& fracture : well_fractures_[i]){
              WellInfo wellinfo = fracture.wellInfo();
              std::cout << " Well " << wellinfo.name << " injection "
                        << injection_rate << " WI " << total_wellindex << std::endl;  
              fracture.setWellProps(injection_rate,  total_wellindex,  wi_dz,  wi_respress,  well_depth);
                // do update wells
                // set well properties

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
                if(wellstate.status != ::Opm::WellStatus::OPEN) {
                    wells_[i].setActive(false);
                    std::cerr << "Warning: Well " << wellinfo.name << " is not open, skipping update." << std::endl;
                    fracture.setActive(false);
                    continue; // skip if not open
                }
                // well is active
                wells_[i].setActive(true);
                // get well perforation
                const auto& perf_data = wellstate.perf_data;
                int  perf_index =  findPerf(wellstate,cell_index_frac);
                if(perf_index == -1) {
                    fracture.setActive(false);
                    wells_[i].setPerfActive(perf_index_frac,false);          
                    continue; // skip if not found
                }
                double perf_pressure = perf_data.pressure[perf_index];
                double perf_rate = perf_data.rates[perf_index];//[FluidSystem::waterPhaseIdx];
                std::cout << "Perf index perf " << perf_index << " fracture " << perf_index_frac << " pressure " << perf_pressure << " rate " << perf_rate << std::endl;
                // std::cout << 
                //std::cout << " Well " << perf_pressure << std::endl;
                //std::cout << "Perf index åerf " << perf_index << " fracture " << perf_index_frac << " pressure " << perf_pressure << std::endl;
                double perf_depth = perf_depths[perf_index];
                fracture.setPerfProps(perf_pressure, perf_depth, perf_rate);
                wells_[i].setPerfPressure(perf_index_frac, perf_pressure);
                fracture.setActive(true);
                wells_[i].setPerfActive(perf_index_frac,true);
                //NB do we need some rates? need to be summed over "peforations of the fractures"
            }
        }
    }
    

}
