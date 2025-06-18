#ifndef OPM_ECLPROBLEM_GEOMECH_HH
#define OPM_ECLPROBLEM_GEOMECH_HH

#include <opm/common/ErrorMacros.hpp>

#include <opm/elasticity/material.hh>
#include <opm/elasticity/materials.hh>

#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/vtkgeomechmodule.hh>
#include <opm/geomech/boundaryutils.hh>
#include <opm/geomech/FlowGeomechLinearSolverParameters.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/elasticity/material.hh>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <array>
#include <functional>
#include <iostream>
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
    template<typename TypeTag>
    class EclProblemGeoMech: public FlowProblemBlackoil<TypeTag>{
    public:
        using Parent = FlowProblemBlackoil<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper = AdaptiveTimeStepping<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
        enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
        enum { dim = GridView::dimension };
        enum { dimWorld = GridView::dimensionworld };
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;
        using GeomechModel = EclGeoMechModel<TypeTag>;
      //using CellSeedType = typename GridView::template Codim<0>::EntitySeed;
        EclProblemGeoMech(Simulator& simulator):
            FlowProblemBlackoil<TypeTag>(simulator),
            geomechModel_(simulator)
        {
            std::string filename = Parameters::Get<Parameters::FractureParamFile>();
            try{
                Opm::PropertyTree fracture_param(filename);
                // set seed values
                fracture_param_ = fracture_param;
            }
            catch(...){
                std::stringstream ss;
                ss << "No fracture parameter file: " << filename << " : no fractures added ";
                //ss << e.what();
                OpmLog::warning(ss.str());
                Opm::PropertyTree fracture_param = makeDefaultFractureParam();
                fracture_param.put("hasfractures",false);
                //fracture_param.put("fractureparams.numfractures",1);
                fracture_param_ = fracture_param;
                
            }
            fracture_param_.write_json(std::cout, true);

            hasFractures_ = this->simulator().vanguard().eclState().runspec().frac();//fracture_param_.get<bool>("hasfractures");
            addPerfsToSchedule_ = fracture_param_.get<bool>("add_perfs_to_schedule");
            if(this->simulator().vanguard().eclState().runspec().mech()){
              this->model().addOutputModule(std::make_unique<VtkGeoMechModule<TypeTag>>(simulator));
            }
        }

        static void registerParameters(){
            Parent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
            FlowLinearSolverParametersGeoMech::registerParameters<TypeTag>();
            Parameters::Register<Parameters::FractureParamFile>("json file defining fracture setting");
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
        }

        void finishInit(){
            OPM_TIMEBLOCK(finishInit);
            Parent::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
            if(eclState.runspec().mech()){
                const auto& initconfig = eclState.getInitConfig();
                geomechModel_.init(initconfig.restartRequested());
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
                    for(int i=0; i < ymodule_.size(); ++i){
                        poelcoef_[i] = (1-2*pratio_[i])/(1-pratio_[i])*biotcoef_[i];
                    }
                }else{
                    if(!fp.has_double("POELCOEF")){
                        OPM_THROW(std::runtime_error,"Missing keyword BIOTCOEF or POELCOEF");
                    }
                    poelcoef_ = fp.get_double("POELCOEF");
                    biotcoef_.resize(ymodule_.size());
                    for(int i=0; i < ymodule_.size(); ++i){
                        biotcoef_[i] = poelcoef_[i]*(1-pratio_[i])/(1-2*pratio_[i]); 
                    }
                }
            
                // thermal related
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();
                if(thermal_expansion){
                    if(fp.has_double("THELCOEF")){
                        thelcoef_ = fp.get_double("THELCOEF");
                        thermexr_.resize(ymodule_.size());
                        for(int i=0; i < ymodule_.size(); ++i){
                            thermexr_[i] = thelcoef_[i]*(1-pratio_[i])/ymodule_[i];
                        }
                    }else{
                        if(!fp.has_double("THERMEXR")){
                            OPM_THROW(std::runtime_error,"Missing keyword THELCOEF or THERMEXR");
                        }
                        thermexr_ = fp.get_double("THERMEXR");
                        thelcoef_.resize(ymodule_.size());
                        for(int i=0; i < ymodule_.size(); ++i){
                            thelcoef_[i] = thermexr_[i]*ymodule_[i]/(1-pratio_[i]);
                        }
                    }
                }
                for(size_t i=0; i < ymodule_.size(); ++i){
                    using IsoMat = Opm::Elasticity::Isotropic;
                    if(pratio_[i]>0.5 || pratio_[i] < 0.0){
                        OPM_THROW(std::runtime_error,"Pratio not valid");
                    }
                    if(biotcoef_[i]>1.0 || biotcoef_[i] < 0.0){
                        OPM_THROW(std::runtime_error,"BIOTCOEF not valid");
                    }
                    elasticparams_.push_back(std::make_shared<IsoMat>(i,ymodule_[i],pratio_[i]));
                }
                if(fp.has_double("CSTRESS")){
                    cstress_ = fp.get_double("CSTRESS");
                }
                // read mechanical boundary conditions
                //const auto& simulator = this->simulator();
                const auto& vanguard = simulator.vanguard();
                const auto& bcconfigs = vanguard.eclState().getSimulationConfig().bcconfig();
                const auto& bcprops = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
                //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
                const auto& gv = this->gridView();
                //const auto& grid = simulator.grid();
                const auto& cartesianIndexMapper = vanguard.cartesianIndexMapper();
                //CartesianIndexMapper cartesianIndexMapper(grid);
                Opm::Elasticity::nodesAtBoundary(bc_nodes_,
                                                 bcconfigs,
                                                 bcprops,
                                                 gv,
                                                 cartesianIndexMapper);



                //using Opm::ParserKeywords::;
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
                                // functors.push_back([&]{

                                //     return center;
                                // }
                                //     );
                            }
                        }
                        recnum +=1;
                    }
                    // NB setting initial stress
                    //this->geomechModel_.setStress(initstress_);
                }else{
                    OPM_THROW(std::runtime_error, "Missing stress initialization keywords");
                }
                // entity_seed_.resize(gv.size(0));
                // for(const auto& elem : elements(gv)){
                //   const auto& cellIdx = gv.indexSet().index(elem);
                //   entity_seed_[cellIdx] = elem.seed();
                // }
            }
        }
        void initialSolutionApplied(){
            OPM_TIMEBLOCK(initialSolutionApplied);
            Parent::initialSolutionApplied();
            const auto& simulator = this->simulator();
            size_t numDof = simulator.model().numGridDof();
            initpressure_.resize(numDof);
            inittemperature_.resize(numDof);

            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = this->model().intensiveQuantities(dofIdx,0);
                const auto& fs = iq.fluidState();
                initpressure_[dofIdx] = Toolbox::value(fs.pressure(waterPhaseIdx));
                inittemperature_[dofIdx] = Toolbox::value(fs.temperature(waterPhaseIdx));
            }
            initstress_.resize(numDof);

            // for now make a copy
            if(simulator.vanguard().eclState().runspec().mech()){
                //this->geomechModel_.setMaterial(elasticparams_);
                this->geomechModel_.setMaterial(ymodule_,pratio_);
                this->geomechModel_.updatePotentialForces();
            }

        }
        void timeIntegration()
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start TimeIntegration-------------------\n"
                << std::flush;
            }
            Parent::timeIntegration();
        }
        void beginTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start beginTimeStep-------------------\n"
                << std::flush;
            }
            Parent::beginTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures()){
                  if(!(cstress_.size() == this->gridView().size(0))){
                        OPM_THROW(std::runtime_error,"CSTRESS not set but fractures exists");
                    }
                }
                geomechModel_.beginTimeStep();
                if(this->hasFractures()){
                    this->wellModel().beginTimeStep();// just to be sure well conteiner is reinitialized
                    this->addConnectionsToWell(); // modify wells WI wiht fracture well 
                }
                
            }
        }
        void endTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;
            }
            //Parent::FlowProblemType::endTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                geomechModel_.endTimeStep();
                if(this->hasFractures() && this->geomechModel().fractureModelActive()){
                    // method for handling extra connections from fractures
                    // it is options for not including them in fractures i.e. addconnections
                    if(addPerfsToSchedule_){
                        this->addConnectionsToSchedual();
                    }else{
                    // not not working ... more work...
                    // will only work if structure is ok
                        assert(false);
                        this->addConnectionsToWell();
                    }
                
                    this->geomechModel_.fractureModel().
                        assignGeomechWellState(this->wellModel_.wellState());
                }
            }

            Parent::endTimeStep();
            //Parent::FlowProblemType::endTimeStep();
            //Parent::endStepApplyAction();
 
        }

        void addConnectionsToWell(){
            //return;
        // add extra connections from fractures direclty to the well structure
            auto& wellcontainer = this->wellModel().localNonshutWells();
            for (auto& wellPtr : wellcontainer) {
                auto wellName = wellPtr->name();
                const auto& wellcons = geomechModel_.getExtraWellIndices(wellName);
                wellPtr->addPerforations(wellcons);
            }
        }

        void addConnectionsToSchedual()
        {
            //return;
        // add extra connections to schedule and use the action framework to handle it
            //const auto& problem = simulator_.problem();
            auto& simulator = this->simulator();
            auto& schedule = simulator.vanguard().schedule();
            const int reportStep = this->episodeIndex();
            // const auto sim_time = simulator_.time() + simulator_.timeStepSize();
            //  const auto now = TimeStampUTC {schedule_.getStartTime()} + std::chrono::duration<double>(sim_time);
            // const auto ts = formatActionDate(now, reportStep);
            std::map<std::string, std::vector<Opm::Connection>> extra_perfs;
            //auto mapper = simulator.vanguard().cartesianMapper();
            //auto& wellcontainer = this->wellModel().localNonshutWells();
            //for (auto& wellPtr : wellcontainer) {
            bool hasExtraPerfs = false;
            //const auto& wells = schedule.getWells(reportStep);
            //for (const auto& wellName : schedule.wellNames(reportStep) ){   
            for (const auto& well : schedule.getWells(reportStep) ){   
                //auto wellName = wellPtr->name();
                auto wellName = well.name();
                const auto& wellcons = geomechModel_.getExtraWellIndices(wellName);
                std::cout << "Adding extra connections for well: " << wellName << " number " << wellcons.size() << std::endl;
                // check if well connections already exist
                if (wellcons.empty()) {
                    continue; // no extra connections for this well
                }
                //hasExtraPerfs = false;
                for (const auto& cons : wellcons) {
                    // simple calculated with upscaling
                    const auto [cell, WI, depth] = cons;
                    // map to cartesian
                    const auto cartesianIdx = simulator.vanguard().cartesianIndex(cell);
                    // get ijk
                    if (well.getConnections().hasGlobalIndex(cell)) {
                        // already have connection for this cell
                        std::cout << "Connection already exists for cell: " << cell << std::endl;
                        continue;
                    }
                    hasExtraPerfs = true;
                    std::array<int, 3> ijk;
                    simulator.vanguard().cartesianCoordinate(cell, ijk);
                    // makeing preliminary connection to be added in schedual with correct numbering

                    Opm::Connection::CTFProperties ctfprop;
                    Opm::Connection connection(ijk[0],
                                               ijk[1],
                                               ijk[2],
                                               cartesianIdx,
                                               /*complnum*/ -1,
                                               Opm::Connection::State::OPEN,
                                               Opm::Connection::Direction::Z,
                                               Opm::Connection::CTFKind::Defaulted,
                                               /*sort_value*/ -1,
                                               depth,
                                               ctfprop,
                                               /*sort_value*/ -1,
                                               /*defaut sattable*/ true);
                    //only add zero value 
                    // connection need to be modified lager
                    connection.setCF(WI*0.0);
                    extra_perfs[wellName].push_back(connection);
                }
            }
            if(!hasExtraPerfs){
                // add to schedule
                return;
            }else{
            // add to schedule
                // structure will be changed erase matrix, maybe only rebuilding of linear solver is neede
                //std::cout << "Rebuilding linear solver for extra connections" << std::endl;
                //this->simulator().model().linearizer().eraseMatrix();
                //this->simulator().model().newtonMethod().linearSolver().setForceRecreate();
                if (this->gridView().comm().rank() == 0) {
                    std::cout << "Adding extra connections to schedule for report step: " << reportStep << std::endl;
                }
            }
            bool commit_wellstate = false;
            auto sim_update = schedule.modifyCompletions(reportStep, extra_perfs);
            // shouldnot be used
            auto updateTrans = [this](const bool global)
            {
                using TransUpdateQuantities = typename Vanguard::TransmissibilityType::TransUpdateQuantities;
                this->transmissibilities_
                    .update(global, TransUpdateQuantities::All, [&vg = this->simulator().vanguard()]
                            (const unsigned int i)
                    {
                        return vg.gridIdxToEquilGridIdx(i);
                    });
            };
            // alwas rebuild wells
            sim_update.well_structure_changed = true;
            this->actionHandler_.applySimulatorUpdate(reportStep,
                                                      sim_update,
                                                      updateTrans,
                                                      commit_wellstate);
            if (commit_wellstate) {
                this->wellModel().commitWGState();
            }
        }


        void endEpisode(){
            Parent::endEpisode();
            geomechModel_.writeFractureSolution();
        }

        const EclGeoMechModel<TypeTag>& geoMechModel() const
        { return geomechModel_; }

        EclGeoMechModel<TypeTag>& geoMechModel()
        { return geomechModel_; }

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

        // double pRatio(unsigned globalIdx) const{
        //     return pratio_[globalIdx];
        // }


        const std::vector<std::tuple<size_t,MechBCValue>>& bcNodes() const{
            return bc_nodes_;
        }

        Dune::FieldVector<double,6> stress(size_t globalIdx) const{
            return geomechModel_.stress(globalIdx);
        }
        // double getFieldProps(const std::string& field, unsigned globalIdx) const{
        //     const auto& eclState = this->simulator().vanguard().eclState();
        //     const auto& fp = eclState.fieldProps();
        //     const auto& myvec = fp.get_double(field);
        //     return myvec[globalIdx];
        // }
        bool hasFractures() const{ return hasFractures_;}
        Opm::PropertyTree getFractureParam() const{return fracture_param_.get_child("fractureparam");};
        Opm::PropertyTree getGeomechParam() const{return fracture_param_;};
        GeomechModel& geomechModel() {return geomechModel_;}
        const GeomechModel& geomechModel()const {return geomechModel_;}
        // used for fracture model
        double yModule(size_t idx) const {return ymodule_[idx];}
        double pRatio(size_t idx) const {return pratio_[idx];}
        double cStress(size_t idx) const {return cstress_[idx];}
        
      //const std::vector< GridView::Codim<0>::EntitySeed >& elementEntitySeed(){return entitity_seed_;}
      //const std::vector< CellSeedType >& elementEntitySeed(){return entity_seed_;}
    private:
        GeomechModel geomechModel_;

        std::vector<double> ymodule_;
        std::vector<double> pratio_;
        std::vector<double> biotcoef_;
        std::vector<double> poelcoef_;
        std::vector<double> thermexr_;
        std::vector<double> thelcoef_;
        std::vector<double> cstress_;
        
        std::vector<double> initpressure_;
        std::vector<double> inittemperature_;
        std::vector<std::tuple<size_t,MechBCValue>> bc_nodes_;
        Dune::BlockVector< SymTensor > initstress_;
        //std::vector<Opm::Elasticity::Material> elasticparams_;
        std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;

        // for fracture calculation
        bool hasFractures_;
        bool addPerfsToSchedule_;
        Opm::PropertyTree fracture_param_;
      //std::vector< CellSeedType > entity_seed_;
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };
}
#endif
