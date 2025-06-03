#include "config.h"
#include <opm/geomech/FractureModel.hpp>
#include <dune/common/fmatrixev.hh>
#include <opm/geomech/DiscreteDisplacement.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFractureSeeds.hpp>

#include <algorithm>
#include <cstddef>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <stddef.h>

namespace Opm{
    Opm::PropertyTree makeDefaultFractureParam(){
        using namespace std::string_literals;
        Opm::PropertyTree fracture_param;
        fracture_param.put("hasfractures",false);
        fracture_param.put("add_perfs_to_schedule", true);
        // solution method    
        fracture_param.put("solver.method", "PostSolve"s);
        fracture_param.put("solver.implicit_flow", false);
        fracture_param.put("solver.max_mech_it", 2);

        // solve method for fracture
        fracture_param.put("fractureparam.method.iterate", true);
        fracture_param.put("fractureparam.method.max_it", 3);

        //
        fracture_param.put("fractureparam.reduce_boundary", false);
        fracture_param.put("fractureparam.addconnections", true);
        // very experimental to calculate stress contributions from fracture to cell values
        fracture_param.put("fractureparam.include_fracture_contributions", false);

        // seed to be in input file   
        // fracture_param.put("fractureparam.config.type", "well_seed");
        // fracture_param.put("fractureparam.config.well", "P1");
        // fracture_param.put("fractureparam.config.cell_ijk", std::vector<int>{8, 8, 10});
        // fracture_param.put("fractureparam.config.normal", std::vector<double>{0, 1, 0});
        // fracture_param.put("fractureparam.config.initial_fracture_width", 0.0);
        // fracture_param.put("fractureparam.config.num_exp", 3);
        // fracture_param.put("fractureparam.config.num_lin", 2);
        // fracture_param.put("fractureparam.config.axis_scale", 4.50);

        // propagation properties
        //fracture_param.put("fractureparam.solver.method", "if_propagate_trimesh"s);
        fracture_param.put("fractureparam.solver.method", "if_propagate"s);
        fracture_param.put("fractureparam.solver.efac", 0.5);
        fracture_param.put("fractureparam.solver.rfac", 0.1);
        fracture_param.put("fractureparam.solver.max_iter", 100);
        fracture_param.put("fractureparam.solver.damping", 1e0);
        fracture_param.put("fractureparam.solver.min_width", 0.0);
        fracture_param.put("fractureparam.solver.max_width", 1e10);
        fracture_param.put("fractureparam.solver.max_dwidth", 1e-3);
        fracture_param.put("fractureparam.solver.max_dp", 1e6);
        fracture_param.put("fractureparam.solver.max_change", 1e-3);
        fracture_param.put("fractureparam.solver.verbosity", 0);

        // fracture linear solve
        fracture_param.put("fractureparam.solver.linsolver.tol", 1e-10);
        fracture_param.put("fractureparam.solver.linsolver.max_iter", 1000);
        fracture_param.put("fractureparam.solver.linsolver.verbosity", 1);

        // reservoir fracture coupling
        fracture_param.put("fractureparam.reservoir.dist", 1e1);
        fracture_param.put("fractureparam.reservoir.calculate_dist", true);
        //fracture_param.put("fractureparam.reservoir.mobility", 1.3e-3);
        //fracture_param.put("fractureparam.reservoir.perm", 1e-13);

        // well fracture coupling
        fracture_param.put("fractureparam.control.type", "perf_pressure"s);
        //fracture_param.put("fractureparam.control.rate", 2.9e-2);
        fracture_param.put("fractureparam.control.WI", 1.0e-11);


        //fracture_param.put("fractureparam.KMax", 1e6);// in input file
        fracture_param.put("fractureparam.extended_fractures", true);
        fracture_param.put("fractureparam.fractureWI", 0.1);
        fracture_param.put("fractureparam.write_pressure_system", false);
        fracture_param.put("fractureparam.write_fracture_system", false);
        fracture_param.put("fractureparam.pressuresolver", "umfpack"s);
        fracture_param.put("fractureparam.fracturesolver", "notused"s);
        return fracture_param;
   }

    void FractureModel::addWell(std::string name,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments,
                           const std::vector<int>& res_cells){
        std::string outputdir = prm_.get<std::string>("outputdir");
        std::string casename = prm_.get<std::string>("casename"); 
        wells_.push_back(FractureWell(outputdir,casename, name,points,segments, res_cells));
        // add witth no fractures
        well_fractures_.push_back(std::vector<Fracture>());
    }

    void FractureModel::addFractures(const ScheduleState& sched)
    {
        const auto fracture_type = this->prm_
            .get<std::string>("config.type", "well_seed");

        if (fracture_type == "perp_well") {
            this->addFracturesPerpWell();
        }
        else if (fracture_type == "tensile_fracture") {
            this->addFracturesTensile();
        //}
        // else if (fracture_type == "well_seed_json") {
        //     // old method
        //     Fracture fracture;
        //     fracture.init(well.name(), perf, well_cell, origo, normal, prm_);
        //     well_fractures_[i].push_back(std::move(fracture));
        }else if (fracture_type == "well_seed") {
            this->addFracturesWellSeed(sched);
        }
        else {
            OPM_THROW(std::runtime_error, "Fracture type '" + fracture_type + "' is not supported");
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
    // probably this should be collected in one loop
    Dune::FieldVector<double,6> FractureModel::stress(Dune::FieldVector<double,3> obs) const{
        Dune::FieldVector<double,6> stress;
        for(size_t i=0; i < wells_.size(); ++i){
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                Dune::FieldVector<double,6> tmp = fractures[j].stress(obs);
                stress += tmp;
            }
        }
        return stress;
    }
    Dune::FieldVector<double,6> FractureModel::strain(Dune::FieldVector<double,3> obs) const{
        Dune::FieldVector<double,6> strain;
        for(size_t i=0; i < wells_.size(); ++i){
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                Dune::FieldVector<double,6> tmp = fractures[j].strain(obs);
                strain += tmp;
            }
        }
        return strain;
    }
    Dune::FieldVector<double,3> FractureModel::disp(Dune::FieldVector<double,3> obs) const{
        Dune::FieldVector<double,3> disp;
        for(size_t i=0; i < wells_.size(); ++i){
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                Dune::FieldVector<double,3> tmp = fractures[j].disp(obs);
                disp += tmp;
            }
        }
        return disp;

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
            wells_[i].writemulti(time);
            const std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].writemulti(time);
            }
        }
    }
    // void FractureModel::solve() {
    //     for(size_t i=0; i < wells_.size(); ++i){
    //         std::vector<Fracture>& fractures = well_fractures_[i];
    //         for(size_t j=0; j < fractures.size(); ++j){
    //             fractures[j].solve(cell_search_tree_);
    //         }
    //     }
    // }
    void FractureModel::updateReservoirProperties() {
        for(size_t i=0; i < wells_.size(); ++i){
            std::vector<Fracture>& fractures = well_fractures_[i];
            for(size_t j=0; j < fractures.size(); ++j){
                fractures[j].updateReservoirProperties();
            }
        }
    }

    std::vector<RuntimePerforation>
    FractureModel::getExtraWellIndices(const std::string& wellname) const
    {
        // for now just do a search
        bool addconnections = prm_.get<bool>("addconnections");
        if (addconnections) {
            for (size_t i = 0; i < wells_.size(); ++i) {
                if (wells_[i].name() == wellname) {
                    // collect all from a well
                    std::vector<RuntimePerforation> wellindices;
                    for (const auto& frac : well_fractures_[i]) {
                        auto perfs = frac.wellIndices();
                        wellindices.insert(wellindices.end(),perfs.begin(), perfs.end());
                    }
                    return wellindices;
                }
            }
            std::string message = "Now fractures on this well found";
            message += wellname;
            OPM_THROW(std::runtime_error, message.c_str());
        }
        return {};
    }

    template <typename Scalar>
    void FractureModel::assignGeomechWellState(WellState<Scalar>& wellState) const
    {
        const auto nWells = this->wells_.size();
        for (auto i = 0*nWells; i < nWells; ++i) {
            const auto wsIx = wellState.index(this->wells_[i].name());
            if (! wsIx.has_value()) { continue; }

            auto& perfData = wellState[*wsIx].perf_data;

            if (perfData.connFracStatistics.size() != perfData.cell_index.size()) {
                perfData.connFracStatistics.resize(perfData.cell_index.size());
            }

            for (const auto& fracture : this->well_fractures_[i]) {
                auto perfPos = std::find(perfData.cell_index.begin(),
                                         perfData.cell_index.end(),
                                         fracture.wellInfo().well_cell);
                if (perfPos == perfData.cell_index.end()) { continue; }

                // Possibly just "fracture.wellInfo().perf" instead.
                const auto perfIx = std::distance(perfData.cell_index.begin(), perfPos);

                fracture.assignGeomechWellState(perfData.connFracStatistics[perfIx]);
            }
        }
    }
} // namespace Opm

// ===========================================================================
// Private member functions
// ===========================================================================

void Opm::FractureModel::addFracturesPerpWell()
{
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper
        <FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0*this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {
            fracWell.grid().leafGridView(), Dune::mcmgElementLayout()
        };

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto& geom = elem.geometry();
            const auto& origin = geom.corner(1);
            const auto  normal = origin - geom.corner(0);
            const auto  elemIdx = emap.index(elem);

            this->well_fractures_[wellIx].emplace_back()
                .init(fracWell.name(),
                      /* connection index */ elemIdx,
                      /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                      /* seed origin */ origin,
                      /* fracturing plane's normal vector */ normal,
                      this->prm_);
        }
    }
}

void Opm::FractureModel::addFracturesTensile()
{
    // https://link.springer.com/article/10.1007/s40948-023-00694-1

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper
        <FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0*this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {
            fracWell.grid().leafGridView(), Dune::mcmgElementLayout()
        };

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            auto eigenVectors = Dune::FieldMatrix<double, 3, 3>{};
            auto eigenValues  = Dune::FieldVector<double, 3>{};

            const auto elemIdx = emap.index(elem);

            const auto stressmat = ddm::symTensor2Matrix(fracWell.reservoirStress(elemIdx));
            Dune::FMatrixHelp::eigenValuesVectors(stressmat, eigenValues, eigenVectors);

            // Note: documentation for eigenValuesVectors() seems to imply
            // that minPos == eigenValues.begin() here.
            const auto  minPos = std::min_element(eigenValues.begin(), eigenValues.end());
            const auto& normal = eigenVectors[std::distance(eigenValues.begin(), minPos)];

            this->well_fractures_[wellIx].emplace_back()
                .init(fracWell.name(),
                      /* connection index */ elemIdx,
                      /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                      /* seed origin */ elem.geometry().corner(1),
                      /* fracturing plane's normal vector */ normal,
                      this->prm_);
        }
    }
}

namespace {
    std::unordered_map<int, std::size_t>
    localSeedCells(const Opm::FractureWell&      fracWell,
                   const Opm::WellConnections&   conns,
                   const Opm::WellFractureSeeds& seeds)
    {
        auto localSeedIxMap = std::unordered_map<int, std::size_t>{};

        auto connIx = [&fracWell, &conns](const std::size_t seedCellGlobal)
        {
            auto connPos = std::find_if(conns.begin(), conns.end(),
                                        [seedCellGlobal](const auto& conn)
                                        { return conn.global_index() == seedCellGlobal; });

            if (connPos == conns.end()) {
                return -1;
            }

            return fracWell.reservoirCell(std::distance(conns.begin(), connPos));
        };

        const auto& cells = seeds.seedCells();

        for (auto seedIx = 0*cells.size(); seedIx < cells.size(); ++seedIx) {
            if (const auto ix = connIx(cells[seedIx]); ix >= 0) {
                localSeedIxMap.insert_or_assign(ix, seedIx);
            }
        }
        return localSeedIxMap;
    }
}

void Opm::FractureModel::addFracturesWellSeed(const ScheduleState& sched)
{
    if (sched.wseed().empty()) {
        return;
    }

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper
        <FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0*this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        if (! sched.wseed.has(fracWell.name())) {
            continue;
        }

        const auto& wseed = sched.wseed(fracWell.name());
        const auto localSeeds = localSeedCells
            (fracWell, sched.wells(fracWell.name()).getConnections(), wseed);

        const auto emap = ElementMapper {
            fracWell.grid().leafGridView(), Dune::mcmgElementLayout()
        };

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto elemIdx = emap.index(elem);
            size_t res_cell = fracWell.reservoirCell(elemIdx);
            const auto seedPos = localSeeds.find(res_cell);
            if (seedPos == localSeeds.end()) {
                continue;
            }

            const auto& seedNormal = wseed.getNormal
                (WellFractureSeeds::SeedIndex { seedPos->second });

            const auto normal = Dune::FieldVector<double, 3> {
                seedNormal[0], seedNormal[1], seedNormal[2],
            };

            this->well_fractures_[wellIx].emplace_back()
                .init(fracWell.name(),
                      /* connection index */ elemIdx,
                      /* reservoir cell */ seedPos->first,
                      /* seed origin */ elem.geometry().corner(1),
                      /* fracturing plane's normal vector */ normal,
                      this->prm_);
        }
    }
}

// ===========================================================================
// Explicit specialisations.  No other code below separator.
// ===========================================================================

template void Opm::FractureModel::assignGeomechWellState(WellState<float>&) const;
template void Opm::FractureModel::assignGeomechWellState(WellState<double>&) const;
