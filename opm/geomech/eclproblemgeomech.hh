#ifndef OPM_ECLPROBLEM_GEOMECH_HH
#define OPM_ECLPROBLEM_GEOMECH_HH
#include "ebos/eclproblem.hh"
#include "opm/geomech/vtkgeomechmodule.hh"
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
        enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
        enum { dim = GridView::dimension };
        enum { dimWorld = GridView::dimensionworld };
        using Toolbox = MathToolbox<Evaluation>;
        EclProblemGeoMech(Simulator& simulator):
            EclProblem<TypeTag>(simulator),
            geomechModel_(simulator)
        {
            this->model().addOutputModule(new VtkGeoMechModule<TypeTag>(simulator));
        }

        static void registerParameters(){
            Parent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
        }
        
        void finishInit(){
            Parent::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
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
                if(pratio_[i]>0.5 || pratio_[i] < 0){
                    OPM_THROW(std::runtime_error,"Pratio not valid");
                }
                elasticparams_.push_back(std::make_shared<IsoMat>(i,ymodule_[i],pratio_[i]));
            }
            // read mechanical boundary conditions
            //const auto& simulator = this->simulator();
            const auto& vanguard = simulator.vanguard();
            const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
            if (bcconfig.size() > 0) {
                //nonTrivialBoundaryConditions_ = true;

                size_t numCartDof = vanguard.cartesianSize();
                unsigned numElems = vanguard.gridView().size(/*codim=*/0);
                std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

                for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx){
                    cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;
                }
                std::array<int, 3> cartdim = cartesianIndexMapper.cartesianDimensions();
                if((bcface.i1 < 0) || (bcface.j1<0) || (bcface.k1<0)){
                 throw std::logic_error("Lower range of BC wrong");
                }
                if( (bcface.i2 > cartdim[0]) || (bcface.j1> cartdim[1]) || (bcface.k1 > cartdim[2])){
                    throw std::logic_error("Upper range of BC wrong");
                }
                for (const auto& bcface : bcconfig) {
                    const auto& type = bcface.bcmechtype;
                    if (type == BCMECHType::FREE) {
                        // do nothing
                    }else if (type == BCMECHType::FIXED) {
                        std::set<size_t> effected_cells;
                        for (int i = bcface.i1; i <= bcface.i2; ++i) {
                            for (int j = bcface.j1; j <= bcface.j2; ++j) {
                                for (int k = bcface.k1; k <= bcface.k2; ++k) {
                                    std::array<int, 3> tmp = {i,j,k};
                                    auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                                    if (elemIdx>-1){
                                        effected_cells.insert(elemIdx);
                                    }
                                }
                            }
                        }
                        const auto& gv = this->gridView();
                        for(const auto& cell:elements(gv)){
                            auto index = gv.indexSet().index(cell);
                            auto it = effected_cells.find(index);
                            if(!(it == effected_cells.end())){
                                // fix all noted for now
                                 for (const auto& vertex : Dune::subEntities(cell, Dune::Codim<dim>{})){
                                     fixed_nodes_.push_back(gv.indexSet().index(vertex));
                                 }
                            }
                        }                    
                    } else {    
                        throw std::logic_error("invalid type for BC. Use FREE or RATE");
                    }
                }
            }
            std::sort(fixed_nodes_.begin(), fixed_nodes_.end()); // {1 1 2 3 4 4 5}
            auto last = std::unique(fixed_nodes_.begin(), fixed_nodes_.end());
            // v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
            fixed_nodes_.erase(last, fixed_nodes_.end());
        }
        void initialSolutionApplied(){
            Parent::initialSolutionApplied();
            auto& simulator = this->simulator();
            size_t numDof = simulator.model().numGridDof();
            initpressure_.resize(numDof);
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = this->model().intensiveQuantities(dofIdx,0);
                const auto& fs = iq.fluidState();
                initpressure_[dofIdx] = Toolbox::value(fs.pressure(waterPhaseIdx));
            }
            // for now make a copy
            this->geomechModel_.setMaterial(elasticparams_);
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
            geomechModel_.beginTimeStep();
        }
        void endTimeStep(){
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;                                                     
            }
            Parent::endTimeStep();
            geomechModel_.endTimeStep();
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
