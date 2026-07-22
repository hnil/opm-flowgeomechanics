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

#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

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

    template<typename TypeTag, class Base>
    class FlowProblemMech: public Base {
    public:
        using MechBase = Base;
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
