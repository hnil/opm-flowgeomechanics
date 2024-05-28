#include "config.h"
#include <opm/geomech/FractureModel.hpp>
namespace Opm{
    void FractureModel::addWell(std::string name,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments){
        wells_.push_back(FractureWell(name,points,segments));
        // add witth no fractures
        well_fractures_.push_back(std::vector<Fracture>());
    }

    void FractureModel::addFractures(){
        for(size_t i=0; i < wells_.size(); ++i){
            const FractureWell& well = wells_[i];
             auto& grid = well.grid();
            for(const auto& elem : elements(grid.leafGridView())){
                auto geo = elem.geometry();
                assert(geo.corners() == 2);
                //auto origo = geo.center();
                auto origo = geo.corner(1);// assume this is cell center
                auto normal = geo.corner(1) - geo.corner(0);
                int perf = well_fractures_[i].size();
                int well_cell = -1;
                Fracture fracture;
                fracture.init(well.name(),
                              perf,
                              well_cell,
                              origo,
                              normal,
                              prm_
                    );
                well_fractures_[i].push_back(std::move(fracture));
            }
        }
    }
    void FractureModel::write(int reportStep) const{
        for(size_t i=0; i < wells_.size(); ++i){
            wells_[i].write();
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].write(reportStep);
            }
        }
    }
    void FractureModel::writemulti(double time) const{
        for(size_t i=0; i < wells_.size(); ++i){
            //wells_[i].write();
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].writemulti(time);
            }
        }
    }
    void FractureModel::solve() {
        for(size_t i=0; i < wells_.size(); ++i){
            std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].solve();
            }
        }
    }
    void FractureModel::updateReservoirProperties() {
        for(size_t i=0; i < wells_.size(); ++i){
            std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].updateReservoirProperties();
            }
        }
    }
}
