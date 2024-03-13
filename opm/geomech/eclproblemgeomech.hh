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

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm{
    template<typename TypeTag>
    class EclProblemGeoMech: public FlowProblem<TypeTag>{
    public:
        using Parent = FlowProblem<TypeTag>;


        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper = AdaptiveTimeStepping<TypeTag>;
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
            FlowProblem<TypeTag>(simulator),
            geomechModel_(simulator)
        {
            if(this->simulator().vanguard().eclState().runspec().mech()){
                this->model().addOutputModule(new VtkGeoMechModule<TypeTag>(simulator));
            }
        }

        static void registerParameters(){
            Parent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
            FlowLinearSolverParametersGeoMech::registerParameters<TypeTag>();
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
                poelcoef_ = fp.get_double("POELCOEF");

                // thermal related
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();
                if(thermal_expansion){
                    thelcoef_ = fp.get_double("THELCOEF");
                    thermexr_ = fp.get_double("THERMEXR");
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
                        for(const auto& cell : elements(gv)){
                            const auto& center = cell.geometry().center();
                            const auto& cellIdx = gv.indexSet().index(cell);
                            const auto& region = equilRegionData[cartesianIndexMapper.cartesianIndex(cellIdx)];
                            if(region == recnum){
                                Dune::FieldVector<double, 6> initstress;
                                initstress[0] = STRESSXX +  STRESSXXGRAD*(center[0] - datum_depth);
                                initstress[1] = STRESSYY +  STRESSYYGRAD*(center[1] - datum_depth);
                                initstress[2] = STRESSZZ +  STRESSZZGRAD*(center[2] - datum_depth);
                                initstress[3] = 0.0;
                                initstress[4] = 0.0;
                                initstress[5] = 0.0;
                                // NB share stress not set to zero
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
                    this->geomechModel_.setStress(initstress_);
                }else{
                    OPM_THROW(std::runtime_error, "Missing stress initialization keywords");
                }
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

        double initTemperature(unsigned dofIdx) const{
            return inittemperature_[dofIdx];
        }

        double initStress(unsigned dofIdx,int comp) const{
            return initstress_[dofIdx][comp];
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
        double pRatio(unsigned globalIdx) const{
            return pratio_[globalIdx];
        }
        const std::vector<std::tuple<size_t,MechBCValue>>& bcNodes() const{
            return bc_nodes_;
        }

        // double getFieldProps(const std::string& field, unsigned globalIdx) const{
        //     const auto& eclState = this->simulator().vanguard().eclState();
        //     const auto& fp = eclState.fieldProps();
        //     const auto& myvec = fp.get_double(field);
        //     return myvec[globalIdx];
        // }
    private:
        using GeomechModel = EclGeoMechModel<TypeTag>;
        GeomechModel geomechModel_;

        std::vector<double> ymodule_;
        std::vector<double> pratio_;
        std::vector<double> biotcoef_;
        std::vector<double> poelcoef_;
        std::vector<double> thermexr_;
        std::vector<double> thelcoef_;

        std::vector<double> initpressure_;
        std::vector<double> inittemperature_;
        std::vector<std::tuple<size_t,MechBCValue>> bc_nodes_;
        Dune::BlockVector<Dune::FieldVector<double,6>> initstress_;
        //std::vector<Opm::Elasticity::Material> elasticparams_;
        std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };
}
#endif
