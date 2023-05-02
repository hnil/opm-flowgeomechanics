#ifndef OPM_ECLPROBLEM_GEOMECH_HH
#define OPM_ECLPROBLEM_GEOMECH_HH
#include "ebos/eclproblem.hh"
#include "opm/geomech/vtkgeomechmodule.hh"
#include "opm/geomech/boundaryutils.hh"
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/elasticity/material.hh>
namespace Opm{
    template<typename TypeTag>
    class EclProblemGeoMech: public EclProblem<TypeTag>{
    public:
        using Parent = EclProblem<TypeTag>;
        
        
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper =  AdaptiveTimeSteppingEbos<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
        enum { dim = GridView::dimension };
        enum { dimWorld = GridView::dimensionworld };
        using Toolbox = MathToolbox<Evaluation>;
        EclProblemGeoMech(Simulator& simulator):
            EclProblem<TypeTag>(simulator),
            geomechModel_(simulator)
        {
            if(this->simulator().vanguard().eclState().runspec().mech()){
                this->model().addOutputModule(new VtkGeoMechModule<TypeTag>(simulator));
            }
        }

        static void registerParameters(){
            Parent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
        }
        
        void finishInit(){
            Parent::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
            if(eclState.runspec().mech()){
                const auto& initconfig = eclState.getInitConfig();
                geomechModel_.init(initconfig.restartRequested());
                const auto& fp = eclState.fieldProps();
                std::vector<std::string> needkeys = {"YMODULE","PRATIO","BIOTCOEF"};
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
                biotcoef_ = fp.get_double("BIOTCOEF");
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
                // read mechanical boundary conditions
                //const auto& simulator = this->simulator();
                const auto& vanguard = simulator.vanguard();
                const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
                //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
                const auto& gv = this->gridView();
                //const auto& grid = simulator.grid();
                const auto& cartesianIndexMapper = vanguard.cartesianIndexMapper();
                //CartesianIndexMapper cartesianIndexMapper(grid);
                Opm::Elasticity::fixNodesAtBoundary(fixed_nodes_,
                                                    bcconfig,
                                                    gv,
                                                    cartesianIndexMapper);
            }
        }
        void initialSolutionApplied(){
            Parent::initialSolutionApplied();
            const auto& simulator = this->simulator();
            size_t numDof = simulator.model().numGridDof();
            initpressure_.resize(numDof);
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = this->model().intensiveQuantities(dofIdx,0);
                const auto& fs = iq.fluidState();
                initpressure_[dofIdx] = Toolbox::value(fs.pressure(waterPhaseIdx));
            }
            // for now make a copy
            if(simulator.vanguard().eclState().runspec().mech()){
                //this->geomechModel_.setMaterial(elasticparams_);
                this->geomechModel_.setMaterial(ymodule_,pratio_);
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
                geomechModel_.beginTimeStep();
            }
        }
        void endTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;                                                     
            }
            Parent::endTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                geomechModel_.endTimeStep();
            }
        }
        
        const EclGeoMechModel<TypeTag>& geoMechModel() const
        { return geomechModel_; }
        
        EclGeoMechModel<TypeTag>& geoMechModel()
        { return geomechModel_; }

        double initPressure(unsigned dofIdx) const{
            return initpressure_[dofIdx];
        }

        double biotCoef(unsigned globalIdx) const{
            return biotcoef_[globalIdx];
        }
        
        const std::vector<size_t>& fixedNodes() const{
            return fixed_nodes_;
        }
    private:
        using GeomechModel = EclGeoMechModel<TypeTag>;
        GeomechModel geomechModel_;
        std::vector<double> ymodule_;
        std::vector<double> pratio_;
        std::vector<double> biotcoef_;
        std::vector<double> poelcoef_;
        std::vector<double> initpressure_;
        std::vector<size_t> fixed_nodes_;
        Dune::BlockVector<Dune::FieldVector<double,6>> initstress_;
        //std::vector<Opm::Elasticity::Material> elasticparams_;
        std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };
}
#endif
