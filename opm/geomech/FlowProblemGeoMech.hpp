#ifndef OPM_FLOW_PROBLEM_GEOMECH_HPP
#define OPM_FLOW_PROBLEM_GEOMECH_HPP

#include <opm/common/ErrorMacros.hpp>

#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/geomech/FlowGeoMechLinearSolverParameters.hpp>
#include <opm/geomech/FlowProblemMech.hpp>
#include <opm/geomech/BoundaryUtils.hpp>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/VtkGeoMechModule.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>


#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/elasticity/material.hh>
#include <opm/elasticity/materials.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm{

    template<typename TypeTag>
    class FlowProblemGeoMech: public FlowProblemMech<TypeTag, FlowProblemBlackoil<TypeTag>>{
    public:
        using MechParent = FlowProblemMech<TypeTag, FlowProblemBlackoil<TypeTag>>;
        using Parent = FlowProblemBlackoil<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper = AdaptiveTimeStepping<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
        enum { dim = GridView::dimension };
        enum { dimWorld = GridView::dimensionworld };
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;

      //using CellSeedType = typename GridView::template Codim<0>::EntitySeed;
        FlowProblemGeoMech(Simulator& simulator):
            MechParent(simulator),
            geoMechModel_(simulator)
        {
            if(this->simulator().vanguard().eclState().runspec().mech()){
              this->model().addOutputModule(std::make_unique<VtkGeoMechModule<TypeTag>>(simulator));
            }
        }

        static void registerParameters(){
            MechParent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
            FlowLinearSolverParametersGeoMech::registerParameters<TypeTag>();
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
        }

        void finishInit(){
            OPM_TIMEBLOCK(finishInit);
            MechParent::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
            if(eclState.runspec().mech()){
                const auto& initconfig = eclState.getInitConfig();
                geoMechModel_.init(initconfig.restartRequested());
                for(size_t i=0; i < this->ymodule_.size(); ++i){
                    using IsoMat = Opm::Elasticity::Isotropic;
                    elasticparams_.push_back(std::make_shared<IsoMat>(i,this->ymodule_[i],this->pratio_[i]));
                }
                // read mechanical boundary conditions
                const auto& vanguard = simulator.vanguard();
                const auto& bcconfigs = vanguard.eclState().getSimulationConfig().bcconfig();
                const auto& bcprops = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
                const auto& gv = this->gridView();
                const auto& cartesianIndexMapper = vanguard.cartesianIndexMapper();
                Opm::Elasticity::nodesAtBoundary(bc_nodes_,
                                                 bcconfigs,
                                                 bcprops,
                                                 gv,
                                                 cartesianIndexMapper);

                bool is_ok = checkBcConfig(bc_nodes_);
                if(!is_ok){
                  // this need to be fixed for parallel runs
                  std::cout << "Error in boundary condition specification not proper for mechanical problem" << std::endl;
                }
            }
        }

        // ///
        // Backend hooks for the shared initial-stress handling
        // ///
        void prepareMechForInit() override{
            this->geoMechModel_.setMaterial(this->ymodule_, this->pratio_);
            this->geoMechModel_.updatePotentialForces();//Neede only of output
        }

        void initializeStressFromMechSolve() override{
            geoMechModel_.solveGeomechanics(/*use_body_force*/ true, /*relative_solve*/ false);
            for (size_t i = 0; i < this->initstress_.size(); ++i) {
                this->initstress_[i] = geoMechModel_.stress(i);
            }
            // stress in output on first step maybe wrong i.e. 2*stress;
            this->geoMechModel_.setFirstSolveTrue();// to do full rebuild next time step
        }

        void applyInitialOutputStress() override{
            this->geoMechModel_.setOutputPutStress(this->initstress_);
        }

        void timeIntegration()
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start TimeIntegration-------------------\n"
                << std::flush;
            }
            Parent::timeIntegration();
        }
        void emptyFractureLogger(){
            const auto& comm = this->simulator().vanguard().grid().comm();
            Opm::DeferredLogger global_logger = gatherDeferredLogger(FractureModel::fractureLogger, comm);
            FractureModel::fractureLogger.clearMessages();
            global_logger.logMessages();
        }
        void beginTimeStep() override{
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start beginTimeStep-------------------\n"
                << std::flush;
            }
            Parent::beginTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures()){
                  if(!(int(this->cstress_.size()) == int(this->gridView().size(0)))){
                        OPM_THROW(std::runtime_error,"CSTRESS not set but fractures exists");
                    }
                }
                geoMechModel_.beginTimeStep();
                if(this->hasFractures()){
                    this->wellModel().beginTimeStep();// just to be sure well conteiner is reinitialized
                    this->addConnectionsToWell(); // modify wells WI wiht fracture well 
                }
                
            }
            OPM_END_PARALLEL_TRY_CATCH("Begin time step geomech failed:",this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
        }
        void endTimeStep() override{
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;
            }
            //Parent::FlowProblemType::endTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                geoMechModel_.endTimeStep();
                if(this->hasFractures() && this->geoMechModel().fractureModelActive()){
                    // method for handling extra connections from fractures
                    // it is options for not including them in fractures i.e. addconnections
                    //if(addPerfsToSchedule_){
                        this->addConnectionsToSchedual();
                        this->gridView().comm().barrier();
                    //}else{
                    // not not working ... more work...
                    // will only work if structure is ok
                    //    assert(false);
                    //    this->addConnectionsToWell();
                    //}
                
                    this->geoMechModel_.fractureModel().
                        assignGeoMechWellState(this->wellModel_.wellState());
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("End time step geomech failed: ", this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
            Parent::endTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures() ){
                    // need to be here ?? to have updated values
                    this->geoMechModel_.updateFilterCakePropertiesOnFractures();
                }
            }
            //
            // if(first_fracture_solve_){
            //     first_fracture_solve_ = false;
            //     geoMechModel_.writeFractureSolutionFirst();
            // }
            if(Parameters::Get<Parameters::EnableWriteAllSolutions>()){
              //OPM_BEGIN_PARALLEL_TRY_CATCH();
                geoMechModel_.writeFractureSolution();
                //  OPM_END_PARALLEL_TRY_CATCH("Writing fracture solution failed: ", this->simulator().vanguard().grid().comm());
            }
            //
            //Parent::FlowProblemType::endTimeStep();
            //Parent::endStepApplyAction();
        }

      double maxNextTimeStepSize() const override{
        double dt_max = Parent::maxNextTimeStepSize();
        if(geoMechModel_.fractureModelActive()){
            const auto& simulator = this->simulator();
            double dt_max_frac = geoMechModel_.fractureModel().template maxFlowTimeStep<TypeTag,Simulator>(simulator);
            dt_max = std::min(dt_max,dt_max_frac);
        }
        return dt_max;
      }

        std::vector<std::vector<RuntimePerforation> >getAllExtraWellIndices()
        {
            std::vector<std::vector<RuntimePerforation> > all_indices;
            auto& wellcontainer = this->wellModel().localNonshutWells();
            for (auto& wellPtr : wellcontainer) {
                auto wellName = wellPtr->name();
                const auto& wellcons = geoMechModel_.getExtraWellIndices(wellName);
                all_indices.push_back(wellcons);
            }
            return all_indices;
        }
      
        void addConnectionsToWell(){
            auto& wellcontainer = this->wellModel().localNonshutWells();
            for (auto& wellPtr : wellcontainer) {
                auto wellName = wellPtr->name();
                const auto& wellcons = geoMechModel_.getExtraWellIndices(wellName);

                wellPtr->addFracturePerforations(wellcons);
            }
        }

        void addConnectionsToSchedual()
        {
            const auto& prm = this->getFractureParam();
            int verbosity = prm.get("fractureparam.verbosity", 1);
            if(verbosity >0){
                std::stringstream os;
                os << "Adding connections to schedule" << std::endl;
                FractureModel::fractureLogger.info(os.str());
            }
            //return;
        // add extra connections to schedule and use the action framework to handle it
            //const auto& problem = simulator_.problem();
            auto& simulator = this->simulator();
            auto& schedule = simulator.vanguard().schedule();
            const int reportStep = this->episodeIndex();

            // const auto sim_time = simulator_.time() + simulator_.timeStepSize();
            //  const auto now = TimeStampUTC {schedule_.getStartTime()} + std::chrono::duration<double>(sim_time);
            // const auto ts = formatActionDate(now, reportStep);

            std::map<std::string, std::vector<Connection>> extra_perfs;

            //auto mapper = simulator.vanguard().cartesianMapper();
            //auto& wellcontainer = this->wellModel().localNonshutWells();
            //for (auto& wellPtr : wellcontainer) {
            int new_conns = 0;
            int old_conns = 0;
            for (const auto& wellName : schedule.wellNames(reportStep)) {
                const auto wellcons = geoMechModel_.getExtraWellIndices(wellName);

                if (wellcons.empty()) {
                    // No extra connections for this well.
                    continue;
                }

                const auto& origConns = this->schedule_[reportStep]
                    .wells(wellName).getConnections();

                auto extra = std::vector<Connection>{};

                for (const auto& wellconn : wellcons) {
                    // simple calculated with upscaling
                    // map to cartesian
                    const auto cartesianIdx = simulator.vanguard()
                        .cartesianIndex(wellconn.cell);

                    if (origConns.hasGlobalIndex(cartesianIdx)) {
                        
                        if (verbosity > 2) {
                            std::stringstream os;
                            os << "Connection already exists for cell: "
                               << wellconn.cell; // << std::endl;
                            FractureModel::fractureLogger.info(os.str());
                        }
                        old_conns += 1;
                        continue;
                    }

                    // get ijk
                    std::array<int, 3> ijk{};
                    simulator.vanguard().cartesianCoordinate(wellconn.cell, ijk);
                    new_conns += 1;
                    if(verbosity > 1) {
                        std::stringstream os;

                        os << "New connection for cell: "
                           << wellconn.cell << " ("
                           << (ijk[0] + 1)  << ", "
                           << (ijk[1] + 1)  << ", "
                           << (ijk[2] + 1)  << ')';

                        FractureModel::fractureLogger.info(os.str());
                    }

                    // Making preliminary connection to be added in schedule
                    // with correct numbering
                    // Dynamic fracturing completions enter the schedule as open
                    // zero-CF well connections; the fracture WI/CTF is added later
                    // through addFracturePerforations(). Even a zero-CF open
                    // completion still changes the rebuilt well topology, perf
                    // ordering, and first-perforation state initialization.
                    auto& connection = extra
                        .emplace_back(ijk[0], ijk[1], ijk[2], cartesianIdx,
                                      /*complnum*/ -1,
                                      Connection::State::OPEN,
                                      Connection::Direction::Z,
                                      Connection::CTFKind::DynamicFracturing,
                                      /* satTableId */ -1,
                                      wellconn.depth,
                                      Connection::CTFProperties{},
                                      /* sort_value */ -1,
                                      /* defaut sattable */ true);

                    //only add zero value 
                    // connection need to be modified later.
                    connection.setCF(wellconn.ctf * 0.0);

                    if (wellconn.perf_range.has_value()) {
                        const auto compseg_insert_index =
                            std::numeric_limits<std::size_t>::max();

                        connection.updateSegment(wellconn.segment,
                                                 wellconn.depth,
                                                 connection.thermalLength(),
                                                 compseg_insert_index,
                                                 wellconn.perf_range);
                    }
                }

                if (! extra.empty()) {
                    const auto* pl = (extra.size() != 1) ? "s" : "";

                    std::stringstream os;
                    os << "Adding " << extra.size()
                       << " extra connection" << pl << " for well: "
                       << wellName;
                    FractureModel::fractureLogger.info(os.str());

                    extra_perfs.insert_or_assign(wellName, std::move(extra));
                }
            }

            if (this->gridView().comm().sum(static_cast<int>(extra_perfs.size())) == 0) {
                return;
            }
            else {
                // add to schedule
                // structure will be changed erase matrix, maybe only rebuilding of linear solver is neede
                //std::cout << "Rebuilding linear solver for extra connections" << std::endl;
                this->simulator().model().linearizer().eraseMatrix();
                this->simulator().model().newtonMethod().linearSolver().eraseMatrix();
                //if (this->gridView().comm().rank() == 0) {
                if(verbosity > 0) {
                    std::stringstream os;
                    os << "Adding extra connections to schedule for report step: "
                   << reportStep << " old connections: " << old_conns << " new connections: " << new_conns;
                    FractureModel::fractureLogger.info(os.str());
                }
                //}
            }

            auto sim_update = schedule.modifyCompletions(reportStep, extra_perfs);

            {
                // Some, or all, of the 'extra_perfs' are entirely new
                // connections created by the fracturing process.  Inform the
                // summary vector calculation engine of these new connections so
                // that they may be included in the summary output if needed.

                const auto root = 0;
                const auto newConns = CollectDynamicConns {
                    Mpi::Packer { this->gridView().comm() }
                }(this->gridView().comm(), root, sim_update.new_frac_wconns);

                if (this->gridView().comm().rank() == root) {
                    this->eclWriter_->recordNewDynamicWellConns(newConns);
                }
            }

            // alwas rebuild wells
            sim_update.well_structure_changed = true;

            bool commit_wellstate = false;
            {
                // The well/group updates triggered by the simulator update
                // log through the group-state helper's deferred logger,
                // which must be established by the caller.
                auto logger_guard = this->wellModel().groupStateHelper().pushLogger();
                this->actionHandler_.applySimulatorUpdate(reportStep,
                                                          sim_update,
                                                          /* updateTrans = */ [](const bool) {},
                                                          commit_wellstate);
            }
            if (commit_wellstate) {
                this->wellModel().commitWGState();
            }
        }


        void endEpisode() override{
            Parent::endEpisode();
            if(!Parameters::Get<Parameters::EnableWriteAllSolutions>()){
                geoMechModel_.writeFractureSolution();
            }
        }

        const GeoMechModel<TypeTag>& geoMechModel() const
        { return geoMechModel_; }

        GeoMechModel<TypeTag>& geoMechModel()
        { return geoMechModel_; }

        //! Mechanical pore-volume change used by the TPSA coupling in the
        //! intensive quantities.  Only evaluated at runtime for the TPSA
        //! mechanics solver; this problem couples through the elasticity/
        //! fracture model instead, so no porosity feedback is reported.
        Scalar rockMechPoroChange(unsigned /*elementIdx*/, unsigned /*timeIdx*/) const
        { return 0.0; }

        const std::vector<std::tuple<size_t,MechBCValue>>& bcNodes() const{
            return bc_nodes_;
        }

        Dune::FieldVector<double,6> stress(size_t globalIdx) const{
            return geoMechModel_.stress(globalIdx);
        }
        // double getFieldProps(const std::string& field, unsigned globalIdx) const{
        //     const auto& eclState = this->simulator().vanguard().eclState();
        //     const auto& fp = eclState.fieldProps();
        //     const auto& myvec = fp.get_double(field);
        //     return myvec[globalIdx];
        // }
      //const std::vector< GridView::Codim<0>::EntitySeed >& elementEntitySeed(){return entitity_seed_;}
      //const std::vector< CellSeedType >& elementEntitySeed(){return entity_seed_;}
        void updateFailed(){
            Parent::updateFailed();
            geoMechModel_.resetFractureModel();
        }

        void advanceTimeLevel(){
            Parent::advanceTimeLevel();
            //this->simulator_.problem().geoMechModel().advanceTimeLevel();// this is done in begin timestep
            //wellModel_.serialize(res);
            //aquiferModel_.serialize(res);
        }
 
    private:
        class CollectDynamicConns : private Serializer<Mpi::Packer>
        {
        public:
            using DynamicConns =
                std::vector<std::pair<std::string, std::vector<std::size_t>>>;

            explicit CollectDynamicConns(const Mpi::Packer& pack)
                : Serializer<Mpi::Packer> { pack }
            {}

            DynamicConns
            operator()(const Parallel::Communication comm,
                       const int                     root,
                       const DynamicConns&           newConns)
            {
                if (comm.size() == 1) {
                    return newConns;
                }

                this->pack(newConns);

                std::tie(this->rankBuffers_, this->rankStart_) =
                    gatherv(this->m_buffer, comm, root);

                if (comm.rank() != root) {
                    // Non-root processes don't need any new connection objects.
                    return {};
                }

                auto allNewConns = DynamicConns{};
                for (auto rank = 0*comm.size(); rank < comm.size(); ++rank) {
                    auto rankNewConns = this->deserialise(rank);

                    allNewConns.insert(allNewConns.end(),
                                       std::make_move_iterator(rankNewConns.begin()),
                                       std::make_move_iterator(rankNewConns.end()));
                }

                return allNewConns;
            }

        private:
            std::vector<char> rankBuffers_{};
            std::vector<int> rankStart_{};

            DynamicConns deserialise(const std::size_t rank)
            {
                auto newConns = DynamicConns{};

                auto begin = this->rankBuffers_.begin() + this->rankStart_[rank + 0];
                auto end   = this->rankBuffers_.begin() + this->rankStart_[rank + 1];

                this->m_buffer.assign(begin, end);
                this->m_packSize = std::distance(begin, end);

                this->unpack(newConns);

                return newConns;
            }
        };

        // ----------------------------------------------------------------------------
        // Heuristic rotational constraint check based only on fixed-direction
        // masks in bc_nodes.
        //
        // Since coordinates are not available here, this is a structural check:
        // to constrain rotation around an axis, we require fixed DOFs in the two
        // transverse directions on at least two distinct nodes.
        // ----------------------------------------------------------------------------
        static bool checkRotationsConstrained(
            const std::vector<std::tuple<size_t, MechBCValue>>& bc_nodes)
        {
            auto constrainedAroundAxis = [&](int d1,
                                             int d2,
                                             const std::string& axisName) -> bool
            {
                std::set<size_t> nodes_d1;
                std::set<size_t> nodes_d2;

                for (const auto& [node_idx, bc] : bc_nodes) {
                    if (bc.fixeddir[d1]) {
                        nodes_d1.insert(node_idx);
                    }
                    if (bc.fixeddir[d2]) {
                        nodes_d2.insert(node_idx);
                    }
                }

                if (nodes_d1.empty() || nodes_d2.empty()) {
                    OpmLog::warning(
                        "Mechanical BC: rotation around " + axisName
                        + " may be unconstrained (missing fixed DOFs in one or "
                          "both transverse directions).");
                    return false;
                }

                std::set<size_t> union_nodes = nodes_d1;
                union_nodes.insert(nodes_d2.begin(), nodes_d2.end());

                if (union_nodes.size() < 2) {
                    OpmLog::warning(
                        "Mechanical BC: rotation around " + axisName
                        + " may be unconstrained (transverse constraints are "
                          "applied at a single node only).");
                    return false;
                }

                return true;
            };

            const bool rx = constrainedAroundAxis(1, 2, "X");
            const bool ry = constrainedAroundAxis(0, 2, "Y");
            const bool rz = constrainedAroundAxis(0, 1, "Z");

            return rx && ry && rz;
        }

        // ----------------------------------------------------------------------------
        // Check that the mechanical boundary conditions are self-consistent and
        // cover enough DOFs to prevent rigid-body translations.
        //
        // Returns true if the configuration is valid, false otherwise.
        // An empty bc_nodes list is accepted (e.g. stress-only BCs on all faces).
        // A warning is emitted for any spatial direction that has no fixed node,
        // since the resulting system may be singular.
        // ----------------------------------------------------------------------------
        static bool checkBcConfig(
            const std::vector<std::tuple<size_t, MechBCValue>>& bc_nodes)
        {
            // Count how many nodes fix each direction, and check that every
            // listed node actually constrains something.
            std::array<int, 3> fixed_count = {0, 0, 0};
            const std::array<std::string, 3> dir_name = {"X", "Y", "Z"};

            for (const auto& [node_idx, bc] : bc_nodes) {
                bool node_fixes_anything = false;
                for (int d = 0; d < 3; ++d) {
                    if (bc.fixeddir[d]) {
                        ++fixed_count[d];
                        node_fixes_anything = true;
                    }
                }
                if (!node_fixes_anything) {
                    OpmLog::warning("BC node " + std::to_string(node_idx)
                        + " is listed in bc_nodes but does not fix any"
                          " displacement direction – check BCMECH input.");
                    //return false;
                }
            }

            // Warn (but do not fail) when a direction has no fixed nodes.
            // The system may still be well-posed via stress BCs or symmetry.
            if(bc_nodes.empty()){
                OpmLog::warning(
                    "Mechanical BC configuration: no nodes are fixed.  The system may be singular");
                return false;
            }
            if (!bc_nodes.empty()) {
                for (int d = 0; d < 3; ++d) {
                    if (fixed_count[d] == 0) {
                        OpmLog::warning(
                            "Mechanical BC configuration: no nodes are fixed in "
                            + dir_name[d]
                            + "-direction.  The system may is singular");
                        return false;      
                    };
                    //  else {
                    //     OpmLog::info(
                    //         "Mechanical BC: " + std::to_string(fixed_count[d])
                    //         + " node(s) fixed in " + dir_name[d] + "-direction.");
                    // }
                }
            }

            if (!checkRotationsConstrained(bc_nodes)) {
                OpmLog::warning(
                    "Mechanical BC configuration: rotational modes are not "
                    "sufficiently constrained based on bc_nodes.");
                return false;
            }
            return true;
        }

        GeoMechModel<TypeTag> geoMechModel_;

        std::vector<std::tuple<size_t,MechBCValue>> bc_nodes_;
        //std::vector<Opm::Elasticity::Material> elasticparams_;
        std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;
      //std::vector< CellSeedType > entity_seed_;
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };
}
#endif
