#include "config.h"
#include <opm/geomech/FractureModel.hpp>
#include <dune/common/fmatrixev.hh>
#include <opm/geomech/DiscreteDisplacement.hpp>
namespace Opm{
    void FractureModel::addWell(std::string name,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments,
                           const std::vector<int>& res_cells){
        std::string outputdir = prm_.get<std::string>("outputdir");                   
        wells_.push_back(FractureWell(outputdir,name,points,segments, res_cells));
        // add witth no fractures
        well_fractures_.push_back(std::vector<Fracture>());
    }

    void
    FractureModel::addFractures()
    {
        auto config = prm_.get_child("config");
        std::string type = config.get<std::string>("type");
        for (size_t i = 0; i < wells_.size(); ++i) {
            const FractureWell& well = wells_[i];
            auto& grid = well.grid();
            using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;
            ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
            for (const auto& elem : elements(grid.leafGridView())) {
                int eIdx = mapper.index(elem);
                auto geo = elem.geometry();
                assert(geo.corners() == 2);
                // auto origo = geo.center();
                Dune::FieldVector<double, 3> origo, normal;
                int perf = eIdx;
                assert(eIdx == well_fractures_[i].size());
                int well_cell = well.reservoirCell(eIdx);

                if (type == "perp_well") {
                    origo = geo.corner(1); // assume this is cell center
                    normal = geo.corner(1) - geo.corner(0);
                    // perf = well_fractures_[i].size();
                } else if (type == "tensile_fracture") {
                    // https://link.springer.com/article/10.1007/s40948-023-00694-1
                    double fractureThoughness = 1.0e6; // reservoir_fractureThoughness_[eIdx]; // ca 1.0 MPa m^1/2
                    double tensilestrength = 5e6; // reservoir_tensilestrength_[eIdx]; //  5 MPa
                    double criticallength = (fractureThoughness / tensilestrength); // ca (1/5)^2 5 mm.
                    criticallength *= criticallength;
                    auto stressmat = ddm::symTensor2Matrix(well.reservoirStress(eIdx));
                    Dune::FieldMatrix<double, 3, 3> eigenVectors;
                    Dune::FieldVector<double, 3> eigenValues;
                    Dune::FMatrixHelp::eigenValuesVectors(stressmat, eigenValues, eigenVectors);
                    int min_dir = -1;
                    int max_dir = -1;
                    double min_eig = 1e99;
                    double max_eig = -1e99;
                    for (int i = 0; i < 3; ++i) {
                        if (eigenValues[i] < min_eig) {
                            min_dir = i;
                            min_eig = eigenValues[i];
                        }
                        if (eigenValues[i] > max_eig) {
                            max_dir = i;
                            max_eig = eigenValues[i];
                        }
                    }
                    normal = eigenVectors[min_dir];
                    // take midpoint
                    origo = geo.corner(1); //-geo.corner(0);
                    // origo /= 2.0;
                    //  expression for size;
                } else {
                    OPM_THROW(std::runtime_error, "Invalid fracture type");
                }


                Fracture fracture;
                fracture.init(well.name(), perf, well_cell, origo, normal, prm_);
                well_fractures_[i].push_back(std::move(fracture));
            }
        }
    }

    void FractureModel::initFractureStates(){
        for(size_t i=0; i < wells_.size(); ++i){
            std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].initFractureStates();
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
    std::vector<std::tuple<int,double,double>> 
    FractureModel::getExtraWellIndices(std::string wellname) const{
        // for now just do a search
        bool addconnections = prm_.get<bool>("addconnections");
        if (addconnections) {
            for (size_t i = 0; i < wells_.size(); ++i) {
                if (wells_[i].name() == wellname) {
                    // collect all from a well
                    std::vector<std::tuple<int, double, double>> wellindices;
                    for (const auto& frac : well_fractures_[i]) {
                        auto perfs = frac.wellIndices();
                        wellindices.insert(wellindices.begin(),perfs.begin(), perfs.end());
                    }
                    return wellindices;
                }
            }
            std::string message = "Now fractures on this well found";
            message += wellname;
            OPM_THROW(std::runtime_error, message.c_str());
        }
        std::vector<std::tuple<int, double, double>> empty; 
        return empty;
    }
}
