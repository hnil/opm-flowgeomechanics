#ifndef OPM_FLOW_PROBLEM_MECH_HPP
#define OPM_FLOW_PROBLEM_MECH_HPP

// Shared problem layer for mechanics-enabled flow problems, extracted from
// FlowProblemGeoMech so that both the VEM-based problem (over
// FlowProblemBlackoil) and a TPSA-based problem (over FlowProblemTPSA) host
// the same geomechanical field properties, initial-stress machinery
// (STRESSEQUIL or backend equilibration solve) and fracture configuration.
//
// The layer deliberately stays free of opm-flowgeomechanics solver includes:
// the STRESSEQUIL/initial-state block and the field-property handling only
// need opm-common/opm-simulators facilities and are candidates for a future
// home in FlowProblemBlackoil upstream.
//
// Backend hooks (implemented by the derived problem):
//   prepareMechForInit()            - backend setup before equilibration
//   initializeStressFromMechSolve() - no-STRESSEQUIL path: initial mech solve
//                                     with body force, capture initstress_
//   applyInitialOutputStress()      - STRESSEQUIL path: seed the backend's
//                                     output stress snapshot

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/common/utility/Serializer.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/geomech/FractureModel.hpp>

#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Parameters {
    struct FractureParamFile {
        inline static std::string value{"notafile"};
    };
}

namespace Opm{
    namespace Detail {
        inline bool fractureParamFileExists(const std::string& filename)
        {
            std::ifstream stream(filename);
            return stream.good();
        }

        bool isBuiltinFractureParamAlias(const std::string& filename);
        Opm::PropertyTree builtinFractureParam(const std::string& filename);
    } // namespace Detail

    // Defined in FractureModel.cpp
    PropertyTree makeDefaultFractureParam(bool rate_control);

    template<typename TypeTag, class Base, class Derived>
    class FlowProblemMech: public Base {
    public:
        using MechBase = Base;

        Derived& derived()
        { return static_cast<Derived&>(*this); }

        const Derived& derived() const
        { return static_cast<const Derived&>(*this); }

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;

        template <class FluidState>
        static int referencePhaseIdx(const FluidState& fs)
        {
            if (fs.phaseIsActive(FluidSystem::waterPhaseIdx)) {
                return FluidSystem::waterPhaseIdx;
            }
            if (fs.phaseIsActive(FluidSystem::oilPhaseIdx)) {
                return FluidSystem::oilPhaseIdx;
            }
            return FluidSystem::gasPhaseIdx;
        }

        explicit FlowProblemMech(Simulator& simulator):
            Base(simulator)
        {
            std::string filename = Parameters::Get<Parameters::FractureParamFile>();
            if (filename == "notafile") {
                Opm::PropertyTree fracture_param = makeDefaultFractureParam(false);
                // fracture_param.put("fractureparams.numfractures",1);
                fracture_param_ = fracture_param;
            } else {
                std::string fractureParamSource = filename;
                try {
                    if (Detail::fractureParamFileExists(filename)) {
                        Opm::PropertyTree fracture_param(filename);
                        fracture_param_ = fracture_param;
                    } else if (Detail::isBuiltinFractureParamAlias(filename)) {
                        fractureParamSource = "built-in alias '" + filename + "'";
                        fracture_param_ = Detail::builtinFractureParam(filename);
                    } else {
                        Opm::PropertyTree fracture_param(filename);
                        fracture_param_ = fracture_param;
                    }
                } catch (...) {
                    std::stringstream ss;
                    ss << "No fracture parameter file or error reading it: " << fractureParamSource
                       << " : Simulation stopped: correct file or use default";
                    // ss << e.what();
                    OpmLog::warning(ss.str());
                    // stop simulation
                    OPM_THROW(std::runtime_error, ss.str());
                }
            }
            std::stringstream os;
            fracture_param_.write_json(os, true);
            OpmLog::info(os.str());

            hasFractures_ = this->simulator().vanguard().eclState().runspec().frac();//fracture_param_.get<bool>("hasfractures");
                        const bool waterActive = this->simulator().vanguard().eclState().runspec().phases().active(Opm::Phase::WATER);
                        if (hasFractures_ && !waterActive) {
                                OPM_THROW(std::runtime_error,
                                                    "Fracture code requires an active WATER phase; decks with MECH+FRAC but without WATER are not supported");
                        }
        }

        static void registerParameters(){
            Base::registerParameters();
	    Parameters::Register<Parameters::FractureParamFile>("json file defining fracture setting or alias: standard, sequential_implicit");
        }

        void finishInit(){
            Base::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
            if(eclState.runspec().mech()){
                const auto& fp = eclState.fieldProps();
                std::vector<std::string> needkeys = {"YMODULE","PRATIO"};
                for(size_t i=0; i < needkeys.size(); ++i){
                    bool ok = fp.has_double(needkeys[i]);
                std::stringstream ss;
                if(!ok){
                    ss << "Missing keyword " << needkeys[i];
                    OPM_THROW(std::runtime_error, ss.str());
                }
                }
                ymodule_ = fp.get_double("YMODULE");
                pratio_ = fp.get_double("PRATIO");
                if(fp.has_double("BIOTCOEF")){
                    biotcoef_ = fp.get_double("BIOTCOEF");
                    poelcoef_.resize(ymodule_.size());
                    for(std::size_t i=0; i < ymodule_.size(); ++i){
                        poelcoef_[i] = (1-2*pratio_[i])/(1-pratio_[i])*biotcoef_[i];
                    }
                }else{
                    if(!fp.has_double("POELCOEF")){
                        OPM_THROW(std::runtime_error,"Missing keyword BIOTCOEF or POELCOEF");
                    }
                    poelcoef_ = fp.get_double("POELCOEF");
                    biotcoef_.resize(ymodule_.size());
                    for(std::size_t i=0; i < ymodule_.size(); ++i){
                        biotcoef_[i] = poelcoef_[i]*(1-pratio_[i])/(1-2*pratio_[i]);
                    }
                }

                // thermal related
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();
                if(thermal_expansion){
                    if(fp.has_double("THELCOEF")){
                        thelcoef_ = fp.get_double("THELCOEF");
                        thermexr_.resize(ymodule_.size());
                        for(std::size_t i=0; i < ymodule_.size(); ++i){
                            thermexr_[i] = thelcoef_[i]*(1-pratio_[i])/ymodule_[i];
                        }
                    }else{
                        if(!fp.has_double("THERMEXR")){
                            OPM_THROW(std::runtime_error,"Missing keyword THELCOEF or THERMEXR");
                        }
                        thermexr_ = fp.get_double("THERMEXR");
                        thelcoef_.resize(ymodule_.size());
                        for(std::size_t i=0; i < ymodule_.size(); ++i){
                            thelcoef_[i] = thermexr_[i]*ymodule_[i]/(1-pratio_[i]);
                        }
                    }
                }
                for(size_t i=0; i < ymodule_.size(); ++i){
                    if(pratio_[i]>0.5 || pratio_[i] < 0.0){
                        OPM_THROW(std::runtime_error,"Pratio not valid");
                    }
                    if(biotcoef_[i]>1.0 || biotcoef_[i] < 0.0){
                        OPM_THROW(std::runtime_error,"BIOTCOEF not valid");
                    }
                }
                if(fp.has_double("CSTRESS")){
                    cstress_ = fp.get_double("CSTRESS");
                }

                const auto& initconfig = eclState.getInitConfig();
                const auto& gv = this->gridView();
                const auto& cartesianIndexMapper = simulator.vanguard().cartesianIndexMapper();
                if( initconfig.hasStressEquil()) {
                    size_t numCartDof = cartesianIndexMapper.cartesianSize();
                    unsigned numElems = gv.size(/*codim=*/0);
                    std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);
                    for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx){
                        cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
                    }
                    const auto& stressequil = initconfig.getStressEquil();
                    const auto& equilRegionData = fp.get_int("STRESSEQUILNUM");
                    //make lambda functions for each regaion
                    std::vector<std::function<std::array<double,6>()>> functors;
                    int recnum=1;
                    initstress_.resize(gv.size(0));
                    for (const auto& record : stressequil) {
                        const auto datum_depth = record.datumDepth();
                        const auto STRESSXX= record.stressXX();
                        const auto STRESSXXGRAD = record.stressXX_grad();
                        const auto STRESSYY= record.stressYY();
                        const auto STRESSYYGRAD = record.stressYY_grad();
                        const auto STRESSZZ= record.stressZZ();
                        const auto STRESSZZGRAD = record.stressZZ_grad();
                        const auto STRESSXY= record.stressXY();
                        const auto STRESSXZ= record.stressXZ();
                        const auto STRESSYZ= record.stressYZ();

                        const auto stressXYGRAD = record.stressXY_grad();
                        const auto stressXZGRAD = record.stressXZ_grad();
                        const auto stressYZGRAD = record.stressYZ_grad();
                        for(const auto& cell : elements(gv)){
                            const auto& center = cell.geometry().center();
                            const auto& cellIdx = gv.indexSet().index(cell);
                            assert(cellIdx < equilRegionData.size());
                            const auto& region = equilRegionData[cellIdx];//cartesianIndexMapper.cartesianIndex(cellIdx)];
                            assert(region <= stressequil.size());
                            if(region == recnum){
                                Dune::FieldVector<double, 6> initstress;
                                initstress[0] = STRESSXX +  STRESSXXGRAD*(center[2] - datum_depth);
                                initstress[1] = STRESSYY +  STRESSYYGRAD*(center[2] - datum_depth);
                                initstress[2] = STRESSZZ +  STRESSZZGRAD*(center[2] - datum_depth);
                                initstress[3] = STRESSYZ +  stressYZGRAD*(center[2] - datum_depth);
                                initstress[4] = STRESSXZ +  stressXZGRAD*(center[2] - datum_depth);
                                initstress[5] = STRESSXY +  stressXYGRAD*(center[2] - datum_depth);
                                // NB share stress not set to zero
                                // we operate with stress = C \grad d + \grad d^T in the matematics
                                initstress_[cellIdx] = initstress;
                            }
                        }
                        recnum +=1;
                    }
                }
            }
        }

        void initialSolutionApplied() override{
            Base::initialSolutionApplied();
            const auto& simulator = this->simulator();
            size_t numDof = simulator.model().numGridDof();
            initpressure_.resize(numDof);
            inittemperature_.resize(numDof);

            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = this->model().intensiveQuantities(dofIdx,0);
                const auto& fs = iq.fluidState();
                const int phaseIdx = referencePhaseIdx(fs);
                initpressure_[dofIdx] = Toolbox::value(fs.pressure(phaseIdx));
                inittemperature_[dofIdx] = Toolbox::value(fs.temperature(phaseIdx));
            }

            initstress_.resize(numDof);

            if (simulator.vanguard().eclState().runspec().mech()) {
                this->prepareMechForInit();
                const auto& eclState = simulator.vanguard().eclState();
                const auto& initconfig = eclState.getInitConfig();
                if (!initconfig.hasStressEquil()) {
                    this->model().invalidateAndUpdateIntensiveQuantities(0);
                    std::stringstream os;
                    os << "No stress equilibration specified .. try to equilibrate";// << std::endl;
                    OpmLog::info(os.str());
                    initstress_.resize(numDof);
                    this->initializeStressFromMechSolve();
                }else{
                    // used set this for initial output of stress all other initial values is zero
                    // which is always calculated relative to initial configuration
                    this->applyInitialOutputStress();
                }
            }
        }

        // ///
        // Backend hooks (overridden by the derived problem)
        // ///
        //! Backend setup before the initial-stress handling (material setup,
        //! potential forces for output, ...).
        virtual void prepareMechForInit()
        {}

        //! No-STRESSEQUIL path: run an initial mechanics solve with body
        //! force and capture the resulting stress in initstress_.
        virtual void initializeStressFromMechSolve()
        {
            OPM_THROW(std::runtime_error,
                      "Stress equilibration without STRESSEQUIL is not "
                      "implemented for this mechanics backend");
        }

        //! STRESSEQUIL path: seed the backend's output stress snapshot with
        //! the initial stress.
        virtual void applyInitialOutputStress()
        {}

        // ///
        // Fracture connection plumbing (shared between mechanics backends;
        // reaches the backend through derived().fractureHost())
        // ///
      double maxNextTimeStepSize() const override{
        double dt_max = Base::maxNextTimeStepSize();
        if(this->derived().fractureHost().fractureModelActive()){
            const auto& simulator = this->simulator();
            double dt_max_frac = this->derived().fractureHost().fractureModel().template maxFlowTimeStep<TypeTag,Simulator>(simulator);
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
                const auto& wellcons = this->derived().fractureHost().getExtraWellIndices(wellName);
                all_indices.push_back(wellcons);
            }
            return all_indices;
        }
      
        void addConnectionsToWell(){
            auto& wellcontainer = this->wellModel().localNonshutWells();
            for (auto& wellPtr : wellcontainer) {
                auto wellName = wellPtr->name();
                const auto& wellcons = this->derived().fractureHost().getExtraWellIndices(wellName);

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
                const auto wellcons = this->derived().fractureHost().getExtraWellIndices(wellName);

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
            Base::endEpisode();
            if(!Parameters::Get<Parameters::EnableWriteAllSolutions>()){
                this->derived().fractureHost().writeFractureSolution();
            }
        }

        void updateFailed(){
            Base::updateFailed();
            this->derived().fractureHost().resetFractureModel();
        }


        // ///
        // Accessors
        // ///
        double initPressure(unsigned dofIdx) const{
            return initpressure_[dofIdx];
        }

        double initTemperature(unsigned dofIdx) const{
            return inittemperature_[dofIdx];
        }

        double initStress(unsigned dofIdx,int comp) const{
            return initstress_[dofIdx][comp];
        }

        const SymTensor& initStress(unsigned dofIdx) const{
            return initstress_[dofIdx];
        }

        double biotCoef(unsigned globalIdx) const{
            return biotcoef_[globalIdx];
        }
        double thelCoef(unsigned globalIdx) const{
            return thelcoef_[globalIdx];
        }
        double thermExr(unsigned globalIdx) const{
            return thermexr_[globalIdx];
        }
        double poelCoef(unsigned globalIdx) const{
            return poelcoef_[globalIdx];
        }

        bool hasFractures() const{ return hasFractures_;}
        Opm::PropertyTree getFractureParam() const{return fracture_param_.get_child("fractureparam");};
        Opm::PropertyTree getGeoMechParam() const{return fracture_param_;};
        // used for fracture model
        double yModule(size_t idx) const {return ymodule_[idx];}
        double pRatio(size_t idx) const {return pratio_[idx];}
        double cStress(size_t idx) const {return cstress_[idx];}

    protected:
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


        std::vector<double> ymodule_;
        std::vector<double> pratio_;
        std::vector<double> biotcoef_;
        std::vector<double> poelcoef_;
        std::vector<double> thermexr_;
        std::vector<double> thelcoef_;
        std::vector<double> cstress_;

        std::vector<double> initpressure_;
        std::vector<double> inittemperature_;
        Dune::BlockVector< SymTensor > initstress_;

        // for fracture calculation
        bool hasFractures_;
        Opm::PropertyTree fracture_param_;
        bool first_fracture_solve_ {true};
    };
}
#endif
