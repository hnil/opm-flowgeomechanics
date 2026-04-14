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
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <stddef.h>

namespace Opm
{
PropertyTree
makeDefaultFractureParam(bool rate_control)
{
    using namespace std::string_literals;
    PropertyTree fracture_param;
    fracture_param.put("hasfractures", true); //used ??
    fracture_param.put("verbosity", 1); //used ??
    fracture_param.put("add_perfs_to_schedule", true);
    // solution method method for coupled system only postsolve tested so far
    fracture_param.put("solver.method", "PostSolve"s);
    fracture_param.put("solver.implicit_flow", false);
    fracture_param.put("solver.max_mech_it", 2);

    // solve method for fracture
    fracture_param.put("fractureparam.verbosity", 1);
    fracture_param.put("fractureparam.method.iterate", true);
    fracture_param.put("fractureparam.method.max_it", 3);
    fracture_param.put("fractureparam.method.tolerance", 1e-4);

    //
    fracture_param.put("fractureparam.reduce_boundary", false);
    fracture_param.put("fractureparam.vem_stability_choice", 3);
    fracture_param.put("fractureparam.smooth_force", false);
    fracture_param.put("fractureparam.addconnections", true);
    // very experimental to calculate stress contributions from fracture to cell values
    fracture_param.put("fractureparam.include_fracture_contributions", false);

    // seed to be in input file
    fracture_param.put("fractureparam.config.verbosity", 1);// well seed should not be used any more
    fracture_param.put("fractureparam.config.type", "well_seed"s);// well seed should not be used any more
    fracture_param.put("fractureparam.config.initial_fracture_width", 1e-4);
    fracture_param.put("fractureparam.config.min_width", 0.0);// normaly taken from deck
    fracture_param.put("fractureparam.config.trires", 5);
    fracture_param.put("fractureparam.config.gravity_off", false);
    fracture_param.put("fractureparam.config.scale_filtrate", true);

    fracture_param.put("fractureparam.solver.verbosity", 1);
    fracture_param.put("fractureparam.solver.method", "if_propagate_trimesh"s);
    fracture_param.put("fractureparam.solver.target_cellcount", 300);
    fracture_param.put("fractureparam.solver.cellcount_threshold", 300);
    fracture_param.put("fractureparam.solver.numcell_threshold", 100000);// probably 1200 should be ok
    fracture_param.put("fractureparam.solver.max_num_coarsening", 0);
    fracture_param.put("fractureparam.solver.max_iter_on_same_level", 100000000);
    
    // used ??
    fracture_param.put("fractureparam.solver.efac", 0.5);
    fracture_param.put("fractureparam.solver.rfac", 0.1);
    //fracture_param.put("fractureparam.solver.damping_w", 0.5);
    //fracture_param.put("fractureparam.solver.damping_p", 0.5);
    fracture_param.put("fractureparam.solver.max_expand_iter", 20);
    fracture_param.put("fractureparam.solver.max_iter", 100);
    fracture_param.put("fractureparam.solver.tolerance", 1e-6);
    fracture_param.put("fractureparam.solver.damping", 1e0);
    fracture_param.put("fractureparam.solver.min_width", 0.001);// NB should not be used
    fracture_param.put("fractureparam.solver.max_width", 0.5);
    fracture_param.put("fractureparam.solver.max_change", 1e5);
    // end used ??

    fracture_param.put("fractureparam.solver.max_dwidth", 5e-3);
    fracture_param.put("fractureparam.solver.max_dp", 1000000000);// NB probably better value 5e6
    fracture_param.put("fractureparam.solver.verbosity", 2);

    // fracture linear solve
    fracture_param.put("fractureparam.solver.remap_solution", false);
    fracture_param.put("fractureparam.solver.linsolver.tol", 1e-8);
    fracture_param.put("fractureparam.solver.area_change_fac", 3.0);
    fracture_param.put("fractureparam.solver.dt_limit", 0.1);// in days
    if(rate_control){
        fracture_param.put("fractureparam.solver.damping_factor_perf", 0.0);// doesnot matter
        fracture_param.put("fractureparam.solver.damping_factor_wi", 0.0);
    }else{
        fracture_param.put("fractureparam.solver.damping_factor_perf", 2.0);
        fracture_param.put("fractureparam.solver.damping_factor_wi", 2.0);
    }
    fracture_param.put("fractureparam.solver.failure_on_nonconvergence", false);
    fracture_param.put("fractureparam.solver.force_limit", 0.0);
    fracture_param.put("fractureparam.solver.smooth_boundary", false);
    fracture_param.put("fractureparam.solver.full_intersections", false);
    fracture_param.put("fractureparam.solver.divide_wellidx", false);
    fracture_param.put("fractureparam.solver.no_leakof_outercells",false);
    // linear solver
    fracture_param.put("fractureparam.solver.linsolver.atol", 1e-20);
    fracture_param.put("fractureparam.solver.linsolver.tol", 1e-10);
    fracture_param.put("fractureparam.solver.linsolver.max_iter", 1000);
    fracture_param.put("fractureparam.solver.linsolver.verbosity", 2);
    fracture_param.put("fractureparam.solver.linsolver.solver", "bicgstab"s);
    // preconditioner
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.diag_mech", false);// for large systems seem sto better with better preconditioner
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.diag_flow", false);
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.mech_first", false);
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.mech_press_coupling", true);
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.verbosity", 0);  
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.solver","umfpack"s);
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.verbosity",0 );
    // not used if not solver is iterative
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.tol",1 );
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.maxiter",1 );
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.preconditioner.type","DILU"s);
    fracture_param.put("fractureparam.solver.linsolver.preconditioner.flow_solver.preconditioner.verbosity",0);
    
    // reservoir fracture coupling
    fracture_param.put("fractureparam.reservoir.calculate_dist", true);
    // not used if calculate_dist is true
    fracture_param.put("fractureparam.reservoir.dist", 1e0);
    
    // normally taken from resrevoir
    fracture_param.put("fractureparam.reservoir.mobility", 1.3e-3);
    fracture_param.put("fractureparam.reservoir.perm", 1e-13);

    // well fracture coupling
    if(rate_control){
        fracture_param.put("fractureparam.control.type", "perf_rate"s);
    }else{
        fracture_param.put("fractureparam.control.type", "perf_pressure"s);
    }
    // only used for some of the options
    fracture_param.put("fractureparam.control.rate", 2.9e-2);
    fracture_param.put("fractureparam.control.WI", 1.0e-6); // calculated from well

    fracture_param.put("fractureparam.fractureWI", 1e-6);//NB to do need to be caluclated based on fracture properties
    
    

    fracture_param.put("fractureparam.extended_fractures", true);

    // nomrmally taken from reservoir
    fracture_param.put("fractureparam.KMax", 1e6); // in input file
    

    fracture_param.put("fractureparam.write_pressure_system", false);
    fracture_param.put("fractureparam.write_fracture_system", false);
    fracture_param.put("fractureparam.pressuresolver", "umfpack"s);
    fracture_param.put("fractureparam.fracturesolver", "notused"s);
    return fracture_param;
}

void
FractureModel::addWell(const std::string& name,
                       const std::vector<FractureWell::Connection>& conns,
                       const std::vector<Point3D>& points,
                       const std::vector<Segment>& segments)
{
    const std::string outputdir = prm_.get<std::string>("outputdir");
    const std::string casename = prm_.get<std::string>("casename");

    // add with no fractures
    wells_.emplace_back(outputdir, casename, name, conns, points, segments);
    well_fractures_.emplace_back();
}

void
FractureModel::addFractures(const ScheduleState& sched)
{
    const auto fracture_type = this->prm_.get<std::string>("config.type", "well_seed");

    if (fracture_type == "perp_well") {
        this->addFracturesPerpWell();
    } else if (fracture_type == "tensile_fracture") {
        this->addFracturesTensile();
    }
    // else if (fracture_type == "well_seed_json") {
    //     // old method
    //     Fracture fracture;
    //     fracture.init(well.name(), perf, well_cell, origo, normal, prm_);
    //     well_fractures_[i].push_back(std::move(fracture));
    //}
    else if (fracture_type == "well_seed") {
        this->addFracturesWellSeed(sched);
    } else {
        OPM_THROW(std::runtime_error, "Fracture type '" + fracture_type + "' is not supported");
    }

    std::cout << "Added fractures to " << wells_.size() << " wells\n"
              << "Total number of fractures_wells: " << well_fractures_.size() << std::endl;

    int count_frac = 0;
    for (auto i = 0 * well_fractures_.size(); i < well_fractures_.size(); ++i) {
        const auto nfrac = well_fractures_[i].size();
        const auto* pl = (nfrac != 1) ? "s" : "";

        count_frac += static_cast<int>(nfrac);

        std::cout << "Well " << wells_[i].name() << " has " << nfrac << " fracture" << pl << '\n';
    }

    std::cout << "Total number of fractures: " << count_frac << std::endl;
}

void
FractureModel::initFractureStates()
{
    for (size_t i = 0; i < wells_.size(); ++i) {
        std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
            fractures[j].initFractureStates();
        }
    }
}
// probably this should be collected in one loop
Dune::FieldVector<double, 6>
FractureModel::stress(Dune::FieldVector<double, 3> obs) const
{
    Dune::FieldVector<double, 6> stress;
    for (size_t i = 0; i < wells_.size(); ++i) {
        const std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
            Dune::FieldVector<double, 6> tmp = fractures[j].stress(obs);
            stress += tmp;
        }
    }
    return stress;
}
Dune::FieldVector<double, 6>
FractureModel::strain(Dune::FieldVector<double, 3> obs) const
{
    Dune::FieldVector<double, 6> strain;
    for (size_t i = 0; i < wells_.size(); ++i) {
        const std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
            Dune::FieldVector<double, 6> tmp = fractures[j].strain(obs);
            strain += tmp;
        }
    }
    return strain;
}
Dune::FieldVector<double, 3>
FractureModel::disp(Dune::FieldVector<double, 3> obs) const
{
    Dune::FieldVector<double, 3> disp;
    for (size_t i = 0; i < wells_.size(); ++i) {
        const std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
            Dune::FieldVector<double, 3> tmp = fractures[j].disp(obs);
            disp += tmp;
        }
    }
    return disp;
}

void
FractureModel::write(int reportStep) const
{
    for (size_t i = 0; i < wells_.size(); ++i) {
        wells_[i].write();
        const std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
            fractures[j].write(reportStep);
        }
    }
}
void
FractureModel::writemulti(double time) const
{
    for (size_t i = 0; i < wells_.size(); ++i) {
        if (vtkwritewells_) {
            wells_[i].writemulti(time);
        }
        const std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
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
void
FractureModel::updateReservoirProperties()
{
    for (size_t i = 0; i < wells_.size(); ++i) {
        std::vector<Fracture>& fractures = well_fractures_[i];
        for (size_t j = 0; j < fractures.size(); ++j) {
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
        int well_idx = -1;
        for (size_t i = 0; i < wells_.size(); ++i) {
            if (!wells_[i].isActive()) {
                if (wells_[i].name() == wellname) {
                    well_idx = i;
                }
                continue;
            }
            if (wells_[i].name() == wellname) {
                well_idx = i;
                // collect all from a well
                std::vector<RuntimePerforation> wellindices;
                for (const auto& frac : well_fractures_[i]) {
                    // assert(frac.isActive());
                    if (!frac.isActive()) {
                        continue;
                    }
                    auto perfs = frac.wellIndices();
                    wellindices.insert(wellindices.end(), perfs.begin(), perfs.end());
                }
                return wellindices;
            }
        }
        if (false) {
            std::cout << "Well " << wellname << " not connections added" << std::endl;
            if (well_idx < -1) {
                std::cout << "Well not found " << wellname << std::endl;
            } else {
                std::cout << "Well " << wellname << " active " << wells_[well_idx].isActive()
                          << " fractures " << well_fractures_[well_idx].size() << std::endl;
                if (well_fractures_[well_idx].size() > 0) {
                    for (const auto& frac : well_fractures_[well_idx]) {
                        std::cout << "Fracture " << frac.name() << " active " << frac.isActive()
                                  << std::endl;
                    }
                }
            }
        }
        // message += wellname;
        // OPM_THROW(std::runtime_error, message.c_str());
    }
    return {};
}

template <typename Scalar,typename IndexTraits>
void
FractureModel::assignGeomechWellState(WellState<Scalar,IndexTraits>& wellState) const
{
    const auto nWells = this->wells_.size();
    for (auto i = 0 * nWells; i < nWells; ++i) {
        if (!wells_[i].isActive()) {
            continue;
        }
        const auto wsIx = wellState.index(this->wells_[i].name());
        if (!wsIx.has_value()) {
            continue;
        }

        auto& perfData = wellState[*wsIx].perf_data;

        if (perfData.connFracStatistics.size() != perfData.cell_index.size()) {
            perfData.connFracStatistics.resize(perfData.cell_index.size());
        }

        for (const auto& fracture : this->well_fractures_[i]) {
            if (!fracture.isActive()) {
                continue;
            }
            auto perfPos = std::find(
                perfData.cell_index.begin(), perfData.cell_index.end(), fracture.wellInfo().well_cell);
            if (perfPos == perfData.cell_index.end()) {
                continue;
            }

            // Possibly just "fracture.wellInfo().perf" instead.
            //const auto perfIx = std::distance(perfData.cell_index.begin(), perfPos);
            fracture.assignGeomechWellState(perfData);//.connFracStatistics[perfIx]);
        }
    }
}
} // namespace Opm

// ===========================================================================
// Private member functions
// ===========================================================================

void
Opm::FractureModel::addFracturesPerpWell()
{
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto& geom = elem.geometry();
            const auto& origin = geom.corner(1);
            const auto normal = origin - geom.corner(0);
            const auto elemIdx = emap.index(elem);

            this->well_fractures_[wellIx].emplace_back().init(
                fracWell.name(),
                /* connection index */ elemIdx,
                /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                /*global_index */ -1, // dummy for now
                /* well segment */ fracWell.segment(elemIdx),
                /* distance range */ fracWell.perfRange(elemIdx),
                /* seed origin */ origin,
                /* fracturing plane's normal vector */ normal,
                this->prm_);
        }
    }
}

void
Opm::FractureModel::addFracturesTensile()
{
    // https://link.springer.com/article/10.1007/s40948-023-00694-1

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            auto eigenVectors = Dune::FieldMatrix<double, 3, 3> {};
            auto eigenValues = Dune::FieldVector<double, 3> {};

            const auto elemIdx = emap.index(elem);

            const auto stressmat = ddm::symTensor2Matrix(fracWell.reservoirStress(elemIdx));
            Dune::FMatrixHelp::eigenValuesVectors(stressmat, eigenValues, eigenVectors);

            // Note: documentation for eigenValuesVectors() seems to imply
            // that minPos == eigenValues.begin() here.
            const auto minPos = std::min_element(eigenValues.begin(), eigenValues.end());
            const auto& normal = eigenVectors[std::distance(eigenValues.begin(), minPos)];

            this->well_fractures_[wellIx].emplace_back().init(
                fracWell.name(),
                /* connection index */ elemIdx,
                /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                /*global_index */ -1, // dummy for now
                /* well segment */ fracWell.segment(elemIdx),
                /* distance range */ fracWell.perfRange(elemIdx),
                /* seed origin */ elem.geometry().corner(1),
                /* fracturing plane's normal vector */ normal,
                this->prm_);
        }
    }
}

namespace
{
std::unordered_map<int, std::size_t>
localSeedCells(const Opm::FractureWell& fracWell,
               const Opm::WellConnections& conns,
               const Opm::WellFractureSeeds& seeds)
{
    auto localSeedIxMap = std::unordered_map<int, std::size_t> {};

    auto connIx = [&fracWell, &conns](const std::size_t seedCellGlobal) {
        auto connPos = std::find_if(conns.begin(), conns.end(), [seedCellGlobal](const auto& conn) {
            return conn.global_index() == seedCellGlobal;
        });

        if (connPos == conns.end()) {
            return -1;
        }

        return fracWell.reservoirCell(std::distance(conns.begin(), connPos));
    };

    const auto& cells = seeds.seedCells();

    for (auto seedIx = 0 * cells.size(); seedIx < cells.size(); ++seedIx) {
        if (const auto ix = connIx(cells[seedIx]); ix >= 0) {
            localSeedIxMap.insert_or_assign(ix, seedIx);
        }
    }
    return localSeedIxMap;
}
} // namespace

void
Opm::FractureModel::addFracturesWellSeed(const ScheduleState& sched)
{
    if (sched.wseed().empty()) {
        return;
    }

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        if (!sched.wseed.has(fracWell.name())) {
            continue;
        }

        const auto& wseed = sched.wseed(fracWell.name());
        const auto localSeeds
            = localSeedCells(fracWell, sched.wells(fracWell.name()).getConnections(), wseed);

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto elemIdx = emap.index(elem);
            const auto seedPos = localSeeds.find(fracWell.reservoirCell(elemIdx));
            if (seedPos == localSeeds.end()) {
                continue;
            }

            const auto& seedNormal = wseed.getNormal(WellFractureSeeds::SeedIndex {seedPos->second});

            const auto& seedSize = wseed.getSize(WellFractureSeeds::SeedIndex {seedPos->second});

            const auto normal = Dune::FieldVector<double, 3> {
                seedNormal[0],
                seedNormal[1],
                seedNormal[2],
            };
            //const auto frac_size = Dune::FieldVector<double, 3> {seedSize[0], seedSize[1], seedSize[2]};
            const auto frac_size = Dune::FieldVector<double, 3> {seedSize.verticalExtent(),
                                                                 seedSize.horizontalExtent(),
                                                                 seedSize.width()};
            // hack
            prm_.put("config.axis_scale", frac_size[0]); // vertical scale
            // prm_.put("fractureparam.config.axis_scale", frac_size[1]);// horizontal scale
            prm_.put("config.min_width", frac_size[2]);

            assert(normal.two_norm() > 0.0);
            const auto& conn = sched.wells(fracWell.name()).getConnections();
            int globalIndex = conn[elemIdx].global_index();
            this->well_fractures_[wellIx].emplace_back().init(
                fracWell.name(),
                /* connection index */ elemIdx,
                /* reservoir cell */ seedPos->first,
                /*global_index */ globalIndex,
                /* well segment */ fracWell.segment(elemIdx),
                /* distance range */ fracWell.perfRange(elemIdx),
                /* seed origin */ elem.geometry().corner(1),
                /* fracturing plane's normal vector */ normal,
                this->prm_);
        }
    }
}

// ===========================================================================
// Explicit specialisations.  No other code below separator.
// ===========================================================================

template void Opm::FractureModel::assignGeomechWellState(WellState<float,Fracture::IndexTraits>&) const;
template void Opm::FractureModel::assignGeomechWellState(WellState<double,Fracture::IndexTraits>&) const;
Opm::DeferredLogger Opm::FractureModel::fractureLogger = Opm::DeferredLogger();
