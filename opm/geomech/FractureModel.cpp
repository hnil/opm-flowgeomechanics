#include "config.h"

#include <opm/geomech/FractureModel.hpp>
#include <opm/geomech/DiscreteDisplacement.hpp>

#include <dune/common/fmatrixev.hh>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <stddef.h>

namespace Opm{
    void FractureModel::addWell(const std::string& name,
                                const std::vector<Point3D>& points,
                                const std::vector<Segment>& segments,
                                const std::vector<int>& res_cells)
    {
        const std::string outputdir = prm_.get<std::string>("outputdir");
        const std::string casename = prm_.get<std::string>("casename");

        wells_.emplace_back(outputdir, casename, name, points, segments, res_cells);

        // add with no fractures
        well_fractures_.emplace_back();
    }

    void FractureModel::addFractures()
    {
        const auto config = prm_.get_child("config");
        const std::string type = config.get<std::string>("type");

        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper
            <FractureWell::Grid::LeafGridView>;

        auto wellFrac = this->well_fractures_.begin();
        for (const auto& well : this->wells_) {
            const auto& grid = well.grid();

            ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
            for (const auto& elem : elements(grid.leafGridView())) {
                const int eIdx = mapper.index(elem);
                const auto geo = elem.geometry();
                assert(geo.corners() == 2);

                // auto origo = geo.center();
                Dune::FieldVector<double, 3> origo, normal;
                int perf = eIdx;
                int well_cell = well.reservoirCell(eIdx);

                if (type == "perp_well") {
                    origo = geo.corner(1); // assume this is cell center
                    normal = geo.corner(1) - geo.corner(0);
                    // perf = well_fractures_[i].size();
		}
                else if (type == "well_seed") {
                    if (config.get<std::string>("well") == well.name()) {
                        std::cout << "Fracure added for well " << well.name() << std::endl;
                        //std::vector<int> cell_ijk = config.get< std::vector<int> > ("cell_ijk");
                        int cell = config.get< int> ("cell");
                        if (well.reservoirCell(eIdx) == cell) {
                            origo = geo.corner(1);
                            auto config_bst = config.getBoostParamPtr();
                            //double tmp = config_bst->get<double > ("normal",1);
                            //for (auto i : as_vector<int>(pt, "a")) std::cout << i << ' ';
                            std::vector<double> tmp_normal = as_vector<double>(*config_bst,"normal");
                            assert(tmp_normal.size() == 3); // wrong use of tmpmal.
                            for (int i=0; i < 3; ++i) {
                                normal[i] = tmp_normal[i];
                            }
                        }
                        else {
                            continue;
                        }
                    }
                    else {
                        continue;
                    }
                }
                else if (type == "tensile_fracture") {
                    // https://link.springer.com/article/10.1007/s40948-023-00694-1
                    const double fractureThoughness = 1.0e6; // reservoir_fractureThoughness_[eIdx]; // ca 1.0 MPa m^1/2
                    const double tensilestrength = 5e6; // reservoir_tensilestrength_[eIdx]; //  5 MPa

                    double criticallength = (fractureThoughness / tensilestrength); // ca (1/5)^2 5 mm.
                    criticallength *= criticallength;

                    Dune::FieldMatrix<double, 3, 3> eigenVectors;
                    Dune::FieldVector<double, 3> eigenValues;
                    auto stressmat = ddm::symTensor2Matrix(well.reservoirStress(eIdx));
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
                }
                else {
                    OPM_THROW(std::runtime_error, "Invalid fracture type");
                }

                wellFrac->emplace_back()
                    .init(well.name(), perf, well_cell, origo, normal, prm_);
            }

            ++wellFrac;
        }
    }

    void FractureModel::initFractureStates()
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                fracture.initFractureStates();
            }
        }
    }

    // probably this should be collected in one loop
    Dune::FieldVector<double,6>
    FractureModel::stress(const Dune::FieldVector<double,3>& obs) const
    {
        Dune::FieldVector<double,6> stress{};

        for (const auto& fractures : this->well_fractures_) {
            for (const auto& fracture : fractures) {
                stress += fracture.stress(obs);
            }
        }

        return stress;
    }

    Dune::FieldVector<double,6>
    FractureModel::strain(const Dune::FieldVector<double,3>& obs) const
    {
        Dune::FieldVector<double,6> strain{};

        for (const auto& fractures : this->well_fractures_) {
            for (const auto& fracture : fractures) {
                strain += fracture.strain(obs);
            }
        }

        return strain;
    }

    Dune::FieldVector<double,3>
    FractureModel::disp(const Dune::FieldVector<double,3>& obs) const
    {
        Dune::FieldVector<double,3> disp{};

        for (const auto& fractures : this->well_fractures_) {
            for (const auto& fracture : fractures) {
                disp += fracture.disp(obs);
            }
        }

        return disp;
    }

    void FractureModel::write(int reportStep) const
    {
        for (auto i = 0*this->wells_.size(); i < this->wells_.size(); ++i) {
            this->wells_[i].write();

            for (const auto& fracture : this->well_fractures_[i]) {
                fracture.write(reportStep);
            }
        }
    }

    void FractureModel::writemulti(double time) const
    {
        for (auto i = 0*this->wells_.size(); i < this->wells_.size(); ++i) {
            this->wells_[i].writemulti(time);

            for (const auto& fracture : this->well_fractures_[i]) {
                fracture.writemulti(time);
            }
        }
    }

    void FractureModel::solve()
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                fracture.solve();
            }
        }
    }

    void FractureModel::updateReservoirProperties()
    {
        for (auto& fractures : this->well_fractures_) {
            for (auto& fracture : fractures) {
                fracture.updateReservoirProperties();
            }
        }
    }

    std::vector<std::tuple<int,double,double>> 
    FractureModel::getExtraWellIndices(const std::string& wellname) const
    {
        auto wellindices = std::vector<std::tuple<int, double, double>>{};

        if (! prm_.get<bool>("addconnections")) {
            return wellindices;
        }

        // For now just do a search
        auto pos = std::find_if(this->wells_.begin(), this->wells_.end(),
                                [&wellname](const auto& well)
                                { return well.name() == wellname; });

        if (pos == this->wells_.end()) {
            OPM_THROW(std::runtime_error, "No fractures found for well " + wellname);
        }

        const auto i = std::distance(this->wells_.begin(), pos);
        for (const auto& frac : well_fractures_[i]) {
            auto perfs = frac.wellIndices();

            wellindices.insert(wellindices.begin(),
                               std::make_move_iterator(perfs.begin()),
                               std::make_move_iterator(perfs.end()));
        }

        return wellindices;
    }
}
