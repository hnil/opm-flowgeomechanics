#ifndef OPM_GEOMECH_MODEL_HPP
#define OPM_GEOMECH_MODEL_HPP
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/geomech/ElasticitySolver.hpp>
#include <opm/geomech/VemElasticitySolver.hpp>

#include <opm/geomech/FlowGeoMechLinearSolverParameters.hpp>
#include <opm/geomech/FractureMechHost.hpp>
#include <opm/geomech/MechFlowMap.hpp>
#include <opm/geomech/MechGridContext.hpp>
#include <opm/geomech/FractureModel.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <dune/common/timer.hh>

#include <iomanip>
#include <memory>
#include <sstream>

namespace Opm{
    template<typename TypeTag>
    class GeoMechModel : public BaseAuxiliaryModule<TypeTag>
    {
        using Parent = BaseAuxiliaryModule<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using Stencil = GetPropType<TypeTag, Properties::Stencil>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
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
    public:
        GeoMechModel(Simulator& simulator):
            //           Parent(simulator)
            first_solve_(true),
            write_system_(false),
            reduce_boundary_(false),
            simulator_(simulator)
        {

            //const auto& eclstate = simulator_.vanguard().eclState();
        };
        //ax model things
        void postSolve(GlobalEqVector&) {
            std::cout << "Geomech dummy PostSolve Aux" << std::endl;
        }
        void addNeighbors(std::vector<NeighborSet>&) const {
            std::cout << "Geomech add neigbors" << std::endl;
        };
        void applyInitial(){
            std::cout << "Geomech applyInitial" << std::endl;
         };
        unsigned numDofs() const{return 0;};
        void linearize(SparseMatrixAdapter&, GlobalEqVector&) {
            std::cout << "Geomech Dummy Linearize" << std::endl;
        };

        // model things
        void beginIteration(){
            //Parent::beginIteration();
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech begin iteration\n");
            }
        }
        void endIteration(){
            //Parent::endIteration();
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech end iteration\n");
            }
        }
        void beginTimeStep(){
            //Parent::beginIteration();
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech begin time step\n");
            }
            fracHost_.beginTimeStep();
        }
        void endTimeStep(){
            // always do post solve
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech model endstimeStep\n");
            }
            this->solveGeoMechAndFracture();
        }
        void solveGeoMechAndFracture(){
            //Parent::endIteration();
            const auto& problem = simulator_.problem();
            this->solveGeomechanics();
            if(problem.hasFractures()){
                this->solveFractures();
            }   
        }

       //! Fracture driving is hosted by FractureMechHost (shared between
       //! mechanics backends); the public surface is forwarded unchanged.
       void solveFractures(){ fracHost_.solveFractures(); }

       void writeFractureSolutionFirst(){ fracHost_.writeFractureSolutionFirst(); }

       void writeFractureSolution(){ fracHost_.writeFractureSolution(); }

       std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname){
           return fracHost_.getExtraWellIndices(wellname);
       }

       bool fractureModelActive() const{ return fracHost_.fractureModelActive(); }

       void updateFilterCakePropertiesOnFractures(){ fracHost_.updateFilterCakePropertiesOnFractures(); }

       const FractureModel& fractureModel() const{ return fracHost_.fractureModel(); }

       FractureModel& fractureModel() { return fracHost_.fractureModel(); }

       void resetFractureModel(){ fracHost_.resetFractureModel(); }

        //! Mechanics runs on its own grid when one was built, otherwise on
        //! the flow grid.
        const Grid& mechGrid() const {
            return mechCtx_ ? mechCtx_->grid() : simulator_.vanguard().grid();
        }

        //! Must be called before init().
        void setMechGrid(const MechGridContext& ctx){
            mechCtx_ = &ctx;
            mechMap_ = ctx.map();
        }

        void updatePotentialForces(bool relative_solve = true){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Update Forces for Geomechanics");
            }
            size_t numDof = simulator_.model().numGridDof();
            const auto& problem = simulator_.problem();
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = simulator_.model().intensiveQuantities(dofIdx,0);
                // pressure part
                const auto& fs = iq.fluidState();
                const int phaseIdx = referencePhaseIdx(fs);
                const auto& press = fs.pressure(phaseIdx);
                // const auto& biotcoef = problem.biotCoef(dofIdx); //NB not used
                //thermal part
                //Properties::EnableTemperature
                const auto& poelCoef = problem.poelCoef(dofIdx);
                double pressure = Toolbox::value(press);
                double diffpress = pressure;
                if(relative_solve){
                    diffpress -=  problem.initPressure(dofIdx);
                }
                
		        const auto& pratio = problem.pRatio(dofIdx);
                double fac = (1-pratio)/(1-2*pratio);
                double pcoeff = poelCoef*fac;
		//assert pcoeff == biot
                pressure_[dofIdx] = pressure;
                mechPotentialForce_[dofIdx] = diffpress*pcoeff;
                mechPotentialPressForce_[dofIdx] = diffpress*pcoeff;
                assert(pcoeff<=(1.0+1e-12));
                mechPotentialPressForceFracture_[dofIdx] = diffpress*(1.0-pcoeff);
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();

                if(thermal_expansion){
                    OPM_TIMEBLOCK(addTermalParametersToMech);
                    const auto& temp = fs.temperature(phaseIdx);//NB all phases have equal temperature
                    const auto& thelcoef = problem.thelCoef(dofIdx);
                    // const auto& termExpr = problem.termExpr(dofIdx); //NB not used
		    // tcoeff = (youngs*tempExp/(1-pratio))*fac;
                    double tcoeff = thelcoef*fac;
                    // assume difftemp = 0 for non termal runs
                    double difftemp = (Toolbox::value(temp));
                    if(relative_solve){
                        difftemp -= problem.initTemperature(dofIdx);
                    }
                    
                    mechPotentialForce_[dofIdx] += difftemp*tcoeff;
                    mechPotentialTempForce_[dofIdx] = difftemp*tcoeff;
                }
                //NB check sign !!
                mechPotentialForce_[dofIdx] *= 1.0;

                // Flow-only cells (numerical-aquifer cells etc.) carry no
                // mechanical loads: they stay in the stiffness as inert rock
                // but must not push on the mechanics (their pressure/PV/depth
                // are flow bookkeeping, not rock state).
                if (problem.mechFlowOnlyCell(dofIdx)) {
                    mechPotentialForce_[dofIdx] = 0.0;
                    mechPotentialPressForce_[dofIdx] = 0.0;
                    mechPotentialPressForceFracture_[dofIdx] = 0.0;
                    mechPotentialTempForce_[dofIdx] = 0.0;
                }
            }
            Opm::PropertyTree param = problem.getFractureParam();
            bool smooth = param.get<bool>("smooth_force",false);
            // smoothing length in meters: regularizes under-resolved (cell-wise sharp)
            // temperature/pressure fronts over a fixed *physical* scale, so the
            // resulting stress extrema are grid independent.  A natural choice is the
            // physical front width, e.g. sqrt(thermal diffusivity * time) for a
            // conductive thermal front.  Takes precedence over the legacy one-ring
            // 'smooth_force' smoothing, whose radius is grid dependent.
            double smooth_length = param.get<double>("smooth_force_length",0.0);
            // Safeguard: force smoothing averages neighbor values, so the
            // zeroed loads of flow-only cells would bleed into the adjacent
            // rock.  Warn loudly rather than silently distort the loads.
            if ((smooth || smooth_length > 0.0) && problem.hasMechFlowOnlyCells()
                && simulator_.gridView().comm().rank() == 0) {
                OpmLog::warning("smooth_force(_length) is active while flow-only "
                                "(mech-masked) cells exist: smoothing will mix their "
                                "zero loads into neighboring rock cells near the "
                                "aquifer attachment.");
            }
            if(smooth_length > 0.0){
              mechPotentialForce_ = vem::diffuseCellVector(simulator_.vanguard().grid(),mechPotentialForce_,smooth_length);
            } else if(smooth){
              mechPotentialForce_ = vem::smoothCellVector(simulator_.vanguard().grid(),mechPotentialForce_);
            }
            // Smoothing belongs on the flow grid; the mechanics grid sees the
            // volume-weighted mean, which keeps the load integral.
            if(!mechMap_.isIdentity()){
                mechMap_.restrict(mechPotentialForce_, mechPotentialForceMech_);
            }
        }

        //! The load the mechanics solver assembles: flow cells when the grids
        //! are the same, mechanics cells otherwise.
        const Dune::BlockVector<Dune::FieldVector<double,1>>& mechLoad() const {
            return mechMap_.isIdentity() ? mechPotentialForce_ : mechPotentialForceMech_;
        }

        //! Mechanics index of a flow cell; the identity while the grids agree.
        std::size_t mechIdx(std::size_t globalIdx) const {
            return static_cast<std::size_t>(mechMap_.mechCell(globalIdx));
        }
        void setupMechSolver(bool use_body_force = false){
                const auto& problem = simulator_.problem();
                //NB hack use fracture params to control mechancis
                Opm::PropertyTree param = problem.getFractureParam();
                reduce_boundary_ = param.get<bool>("reduce_boundary");
                int stability_choice_int = param.get<int>("vem_stability_choice",3);
                bool vem_source = param.get<bool>("vem_source",true);
                bool vem_force = param.get<bool>("vem_force",true);
                bool vem_stress = param.get<bool>("vem_stress",true);
                bool stab_on_stress = param.get<bool>("stab_on_stress",false);
                bool mech_diagonal_scaling = param.get<bool>("mech_diagonal_scaling",false);
                elacticitysolver_->setOptions(stability_choice_int,
                                             vem_source, vem_force,
                                             stab_on_stress,
                                             vem_stress,
                                             mech_diagonal_scaling);
                OPM_TIMEBLOCK(SetupMechSolver);
                bool do_matrix = true;//assemble matrix
                bool do_vector = true;//assemble matrix
                // set boundary
                double gravity = 0.0;
                if(use_body_force){
                    gravity = Opm::unit::gravity;
                    OpmLog::info("Using body force in geomechanics to 2000 kg/m^3 \n");
                }
                int nc = this->mechGrid().leafGridView().size(0);
                std::vector<double> density(nc, 2000.0);
                elacticitysolver_->setBodyForce(gravity, density);
                
                elacticitysolver_->fixNodes(problem.bcNodes());
                //
                elacticitysolver_->initForAssembly();
                elacticitysolver_->assemble(this->mechLoad(), do_matrix, do_vector,reduce_boundary_);
                //elacticitysolver_->assemble_fem(mechPotentialForce_, do_matrix, do_vector,reduce_boundary_);
                FlowLinearSolverParametersGeoMech p;
                p.init<TypeTag>();
                // Print parameters to PRT/DBG logs.

                PropertyTree prm = setupPropertyTree(p, true, true);
                if (p.linear_solver_print_json_definition_) {
                    if(simulator_.gridView().comm().rank() == 0){
                        std::ostringstream os;
                        os << "Property tree for mech linear solver:\n";
                        prm.write_json(os, true);
                        OpmLog::note(os.str());
                    }
                }
                elacticitysolver_->setupSolver(prm);
                elacticitysolver_->comm()->communicator().barrier();
                first_solve_ = false;
                write_system_ = prm.get<int>("verbosity", 0) > 10;
        }
        void writeMechSystem(){
            OPM_TIMEBLOCK(WriteMechSystem);
                    const auto& problem = simulator_.problem();
                    Opm::Helper::writeMechSystem(simulator_,
                    elacticitysolver_->getOperator(),
                    elacticitysolver_->getLoadVector(),
                    elacticitysolver_->comm());
                    {
                        int num_points = this->mechGrid().size(3);
                        Dune::BlockVector<Dune::FieldVector<double,1>> fixed(3*num_points);
                        fixed  = 0.0;
                        const auto& bcnodes = problem.bcNodes();
                        for(const auto& node: bcnodes){
                            int node_idx = std::get<0>(node);
                            auto bcnode = std::get<1>(node);
                            auto fixdir = bcnode.fixeddir;
                            for(int i = 0; i < 3; ++i){ 
                                fixed[3*node_idx+i][0] = fixdir[i];
                            }
                        }
                        Opm::Helper::writeVector(simulator_,
                                                fixed,
                                                "fixed_values_",
                                                elacticitysolver_->comm());
                    }
        }
      void setOutputPutStress(const Dune::BlockVector<Dune::FieldVector<double,6> >& initstress){
            // used at if initial stress not calculated
            outputstress_ = initstress;
      }
      void calculateOutputQuantitiesMech(){//bool relative_solve = true){
            OPM_TIMEBLOCK(CalculateOutputQuantitesMech);
            Opm::Elasticity::Vector field;
            const auto& grid = this->mechGrid();
            const auto& gv = grid.leafGridView();
            static constexpr int dim = Grid::dimension;
            field.resize(grid.size(dim)*dim);
            if(reduce_boundary_){
              elacticitysolver_->expandSolution(field,elacticitysolver_->getSolutionVector());
            }else{
              assert(field.size() == elacticitysolver_->getSolutionVector().size());
              field =  elacticitysolver_->getSolutionVector();
            }

            this->makeDisplacement(field);
            // update variables used for output to resinsight
            // NB TO DO
            {
            OPM_TIMEBLOCK(calculateStress);
            elacticitysolver_->calculateStressPrecomputed(field);
            elacticitysolver_->calculateStrainPrecomputed(field);
            }
            const auto& linstress = elacticitysolver_->stress();
            const auto& linstrain = elacticitysolver_->strain();

            for (const auto& cell: elements(gv)){
                auto cellindex = gv.indexSet().index(cell);
                // add initial stress
                //auto cellindex2 = gv.indexSet().index(cell);
                //stress_[cellindex] = linstress[cellindex];
                strain_[cellindex] = linstrain[cellindex];
                linstress_[cellindex] = linstress[cellindex];
            }
            // optional: replace the raw cell-wise stress by its patch recovery
            // (least-squares nodal fit averaged back to cells).  Reduces the
            // grid-dependent cell-to-cell jumps of the stress, in particular at
            // under-resolved temperature/pressure fronts, and thereby the
            // tractions sampled by the fracture model via problem.stress().
            {
                Opm::PropertyTree param = simulator_.problem().getFractureParam();
                const bool recover = param.get<bool>("patch_recovery_stress",false);
                if(recover){
                    OPM_TIMEBLOCK(patchRecoveryStress);
                    linstress_ = vem::patchRecoveryCells6(grid, linstress_);
                }
            }
            /// NB need to be updated after linstress to be correct
            // output stress is saved here in case init stress is changed before
            // output, and is reported per flow cell: it adds that cell's own
            // initial stress and load to the mechanics cell's solution.
            for(size_t dofIdx = 0; dofIdx < simulator_.model().numGridDof(); ++dofIdx){
                outputstress_[dofIdx] = this->stress(dofIdx);
            }
            //size_t lsdim = 6;
            //for(size_t i = 0; i < stress_.size(); ++i){
            //         stress_[i] += problem.initStress(i);
            //}
            //NB ssume initial strain is 0

            bool verbose = false;
            if(verbose){
                OPM_TIMEBLOCK(WriteMatrixMarket);
                // debug output to matrixmaket format
                Dune::storeMatrixMarket(elacticitysolver_->getOperator(), "A.mtx");
                Dune::storeMatrixMarket(elacticitysolver_->getLoadVector(), "b.mtx");
                Dune::storeMatrixMarket(elacticitysolver_->getSolutionVector(), "u.mtx");
                Dune::storeMatrixMarket(field, "field.mtx");
                Dune::storeMatrixMarket(mechPotentialForce_, "pressforce.mtx");
            }
        }


        void setupAndUpdateGemechanics(bool use_body_force = false, bool relative_solve = true){
            OPM_TIMEBLOCK(endTimeStepMech);
            this->updatePotentialForces(relative_solve);
            // for now assemble and set up solver her
            //const auto& problem = simulator_.problem();
            if(first_solve_){
                this->setupMechSolver(use_body_force);
            }else{
              // only need if first solve use different gravity
              double gravity = 0.0;// used for initialization
              if(use_body_force){
                gravity = 9.81;//normal case in forward simulation
              }
              int nc = this->mechGrid().leafGridView().size(0);
              std::vector<double> density(nc, 2000.0);
              elacticitysolver_->setBodyForce(gravity, density);
            }
            //else
            {
                // reset the rhs even in the first iteration maybe bug in rhs for reduce_boundary=false;
                OPM_TIMEBLOCK(AssembleRhs);
                // need "static boundary conditions is changing"
                //bool do_matrix = false;//assemble matrix
                //bool do_vector = true;//assemble matrix
                //elacticitysolver_->A.initForAssembly();
                //elacticitysolver_->assemble(mechPotentialForce_, do_matrix, do_vector);
                // need precomputed divgrad operator
                elacticitysolver_->updateRhsWithGrad(this->mechLoad());
            }    
        }
        void solveGeomechanics(bool use_body_force = false, bool relative_solve = true){
            setupAndUpdateGemechanics(use_body_force , relative_solve);
            {
                OPM_TIMEBLOCK(SolveMechanicalSystem);
                Dune::Timer solve_timer;
                elacticitysolver_->solve();
                last_mechanical_solve_time_seconds_ = solve_timer.stop();
                total_mechanical_solve_time_seconds_ += last_mechanical_solve_time_seconds_;
                if(simulator_.gridView().comm().rank() == 0){
                    std::ostringstream os;
                    os << "Mechanical solve stats: solves=" << elacticitysolver_->numSolves()
                       << ", linear_iterations=" << elacticitysolver_->lastLinearIterations()
                       << " (total " << elacticitysolver_->totalLinearIterations() << ")"
                       << ", solve_time_s=" << last_mechanical_solve_time_seconds_
                       << " (total " << total_mechanical_solve_time_seconds_ << ")"
                       << ", converged=" << (elacticitysolver_->lastLinearSolveConverged() ? "true" : "false");
                    OpmLog::info(os.str());
                }
                if(write_system_){
                    this->writeMechSystem();
                }
            }
            this->calculateOutputQuantitiesMech();
        }
        template<class Serializer>
        void serializeOp(Serializer&)
        {
            //serializer(tracerConcentration_);
            //serializer(wellTracerRate_);
         }

        // used in FlowProblemGeoMech
        void init(bool /*restart*/){
            if(simulator_.gridView().comm().rank() == 0){
                OpmLog::info("Geomech init\n");
            }
            size_t numDof = simulator_.model().numGridDof();
            pressure_.resize(numDof);
            mechPotentialForce_.resize(numDof);
            mechPotentialTempForce_.resize(numDof);
            mechPotentialPressForce_.resize(numDof);
            mechPotentialPressForceFracture_.resize(numDof);
            // hopefully temperature and pressure initilized
            if(!mechCtx_){
                mechMap_ = MechFlowMap::identity(numDof);
            }
            elacticitysolver_
                = std::make_unique<Opm::Elasticity::VemElasticitySolver<Grid>>(this->mechGrid());

            // Displacement, stress and strain live on the mechanics grid; the
            // loads above stay on the flow grid and are restricted to it.
            const size_t numMech = mechMap_.numMechCells();
            celldisplacement_.resize(numMech);
            std::fill(celldisplacement_.begin(),celldisplacement_.end(),0.0);
            //stress_.resize(numDof);
            linstress_.resize(numMech);
            std::fill(linstress_.begin(),linstress_.end(),0.0);
            outputstress_.resize(numDof);
            std::fill(outputstress_.begin(),outputstress_.end(),0.0);
            strain_.resize(numMech);
            std::fill(strain_.begin(),strain_.end(),0.0);
            const auto& gv = this->mechGrid().leafGridView();
            displacement_.resize(gv.indexSet().size(3));
        };

        void setMaterial(const std::vector<std::shared_ptr<Opm::Elasticity::Material>>& materials){
            elacticitysolver_->setMaterial(materials);
        }
        void setMaterial(const std::vector<double>& ymodule,const std::vector<double>& pratio){
            if(mechMap_.isIdentity()){
                elacticitysolver_->setMaterial(ymodule,pratio);
                return;
            }
            std::vector<double> ymodule_mech, pratio_mech;
            mechMap_.restrict(ymodule, ymodule_mech);
            mechMap_.restrict(pratio, pratio_mech);
            elacticitysolver_->setMaterial(ymodule_mech,pratio_mech);
        }
        const Dune::FieldVector<double,3>& displacement(size_t vertexIndex) const{
            return displacement_[vertexIndex];
        }
        const double& mechPotentialForce(unsigned globalDofIdx) const
        {
            return mechPotentialForce_[globalDofIdx];
        }
        const double& pressure(unsigned globalDofIdx) const
        {
            return pressure_[globalDofIdx];
        }
        const double& mechPotentialTempForce(unsigned globalDofIdx) const
        {
            return mechPotentialTempForce_[globalDofIdx];
        }
        const double& mechPotentialPressForce(unsigned globalDofIdx) const
        {
            return mechPotentialPressForce_[globalDofIdx];
        }

        // Output forms of the potential forces: the stored quantities carry
        // the poroelastic amplification fac = (1-nu)/(1-2nu) used by the
        // equations (pcoeff = poelCoef*fac, tcoeff = thelCoef*fac); for
        // output/comparison the displacement potential without that factor
        // is the more natural quantity.  Keeping this as separate accessors
        // leaves the equations and the primary accessors untouched.
        double mechPotentialPressForceOutput(unsigned globalDofIdx) const
        {
            const double pratio = simulator_.problem().pRatio(globalDofIdx);
            const double fac = (1 - pratio) / (1 - 2 * pratio);
            return mechPotentialPressForce_[globalDofIdx] / fac;
        }

        double mechPotentialTempForceOutput(unsigned globalDofIdx) const
        {
            const double pratio = simulator_.problem().pRatio(globalDofIdx);
            const double fac = (1 - pratio) / (1 - 2 * pratio);
            return mechPotentialTempForce_[globalDofIdx] / fac;
        }
        const Dune::FieldVector<double,3> disp(size_t globalIdx,bool with_fracture = false) const{
            auto disp =  celldisplacement_[this->mechIdx(globalIdx)];
            if(fracHost_.includeFractureContributions() && with_fracture){
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    //size_t locglobalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracHost_.fractureModelPtr()){
                        disp += fracHost_.fractureModelPtr()->disp(obs);
                    }
                }
            }
            return disp;
        }

        const SymTensor delstress(size_t globalIdx) const{
	  Dune::FieldVector<double,6> delStress = this->effstress(globalIdx);
	  double effPress = this->mechPotentialForce(globalIdx);
	  for(int i=0; i < 3; ++i){
	    delStress[i] += effPress;
	    
	  }
	  return delStress;  
        }

        const SymTensor linstress(size_t globalIdx) const{
	  return linstress_[this->mechIdx(globalIdx)];
        }

        const SymTensor outputstress(size_t globalIdx) const{
	        return outputstress_[globalIdx];
        }

        const SymTensor effstress(size_t globalIdx) const{
	  // make stress in with positive with compression
	         return -1.0*linstress_[this->mechIdx(globalIdx)];
        }

        const SymTensor strain(size_t globalIdx,bool with_fracture = false) const{
            auto strain = strain_[this->mechIdx(globalIdx)];
            if(fracHost_.includeFractureContributions() && with_fracture){
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    //size_t locglobalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracHost_.fractureModelPtr()){
                        strain += fracHost_.fractureModelPtr()->strain(obs);
                    }
                }
            }
            return strain_[this->mechIdx(globalIdx)];
        }   

        const SymTensor stress(size_t globalIdx,bool with_fracture = false) const{
            Dune::FieldVector<double,6> effStress = this->effstress(globalIdx);
            effStress += simulator_.problem().initStress(globalIdx);
            double effPress = this->mechPotentialForce(globalIdx);
            for(int i=0; i < 3; ++i){
                effStress[i] += effPress;

            }
            if(fracHost_.includeFractureContributions() && with_fracture){
                assert(false);
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    //size_t loglobalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracHost_.fractureModelPtr()){
                        effStress += fracHost_.fractureModelPtr()->stress(obs);
                    }
                }
            }

            return effStress;
         }

         const SymTensor fractureStress(size_t globalIdx) const{
            // need to know effect of absolute pressure
	   // const auto& pratio = simulator_.problem.pRatio(globalIdx);
	   // double fac = (1-pratio)/(1-2*pratio);
	   // const auto& poelCoef = simulator_.problem.poelCoef(globalIdx);
	   // double pcoeff = poelCoef*fac;
	   const auto& iq = simulator_.model().intensiveQuantities(globalIdx,0);
	   const auto& fs = iq.fluidState();
       const auto& press = Toolbox::value(fs.pressure(referencePhaseIdx(fs)));
	   //
            Dune::FieldVector<double,6> fracStress = this->stress(globalIdx);
            // effStress += simulator_.problem().initStress(globalIdx);
            // double effPress = this->mechPotentialTempForce(globalIdx);
            // effPress += mechPotentialPressForceFracture_[globalIdx];
	    // effPress -= press;
            for(int i=0; i < 3; ++i){
                fracStress[i] -= press;
            }
            return fracStress;
         }

        // NB used in output should be eliminated

        double pressureDiff(unsigned dofIx) const{
            return mechPotentialForce_[dofIx];
        }

        // const double& disp(unsigned globalDofIdx, unsigned dim) const
        // {
        //     return celldisplacement_[globalDofIdx][dim];
            
        // }
        // const double stress(unsigned globalDofIdx, unsigned dim) const
        // {
        //     // not efficient
        //     auto stress = this->stress(globalDofIdx);
        //     return stress[dim];
        //     //return stress_[globalDofIdx][dim];
        // }

        // const double& delstress(unsigned globalDofIdx, unsigned dim) const
        // {
        //     return linstress_[globalDofIdx][dim];
        // }
        // const double& strain(unsigned globalDofIdx, unsigned dim) const
        // {
        //     return strain_[globalDofIdx][dim];
        // }


        // void setStress(const Dune::BlockVector<SymTensor >& stress){
        //     stress_ = stress;
        // }
        void makeDisplacement(const Opm::Elasticity::Vector& field) {
            // make displacement on all nodes used for output to vtk
            const auto& grid = this->mechGrid();
            const auto& gv = grid.leafGridView();
            int dim = 3;
            for (const auto& vertex : Dune::vertices(gv)){
                auto index = gv.indexSet().index(vertex);
                for(int k=0; k < dim; ++k){
                    displacement_[index][k] = field[index*dim+k];
                }
            }
            for (const auto& cell: elements(gv)){
                auto cellindex = gv.indexSet().index(cell);
                celldisplacement_[cellindex] = 0.0;
                const auto& vertices = Dune::subEntities(cell, Dune::Codim<Grid::dimension>{});
                for (const auto& vertex : vertices){
                    auto nodeidex = gv.indexSet().index(vertex);
                    for(int k=0; k < dim; ++k){
                        celldisplacement_[cellindex][k] += displacement_[nodeidex][k] ;
                    }
                }
                celldisplacement_[cellindex] /= vertices.size();
            }
        }

        std::string finalTimingSummary() const
        {
            std::ostringstream os;
            auto report_time = [&os](const char* label, const double value)
            {
                os << std::left << std::setw(28) << label
                   << std::right << std::setw(7) << std::fixed << std::setprecision(2) << value
                   << " s\n";
            };
            auto report_count = [&os](const char* label, const int value)
            {
                os << std::left << std::setw(28) << label
                   << std::right << std::setw(7) << value << "\n";
            };

            report_time("  Mechanical solve time:", total_mechanical_solve_time_seconds_);
            report_count("Overall Mechanical Solves:", elacticitysolver_->numSolves());
            report_count("Overall Mech Lin Iters:", elacticitysolver_->totalLinearIterations());

            if (simulator_.problem().hasFractures()) {
                const auto total_stats = fracHost_.fractureModelPtr() ? fracHost_.fractureModelPtr()->totalSolveStats() : FractureSolveStats{};
                report_time("  Fracture solve time:", total_stats.solve_time_seconds);
                report_count("Overall Fracture Solves:", total_stats.fractures_solved);
                report_count("Overall Fracture Nl Iters:", total_stats.nonlinear_iterations);
                report_count("Overall Fracture Lin Solves:", total_stats.linear_solves);
                report_count("Overall Fracture Lin Iters:", total_stats.linear_iterations);
            }

            return os.str();
        }
        
        void setFirstSolveTrue(){
            elacticitysolver_->resetOperator();
            first_solve_ = true;
        }
    private:
        bool first_solve_{true};
        bool write_system_{false};
        bool reduce_boundary_{false};
        double last_mechanical_solve_time_seconds_{0.0};
        double total_mechanical_solve_time_seconds_{0.0};
        Simulator& simulator_;

        Dune::BlockVector<Dune::FieldVector<double,1>> pressure_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForce_;
        // Only filled when the mechanics grid is a coarsening of the flow grid.
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForceMech_;
        MechFlowMap mechMap_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForceFracture_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialTempForce_;
        //Dune::BlockVector<Dune::FieldVector<double,1> > solution_;
        Dune::BlockVector<Dune::FieldVector<double,3> > celldisplacement_;
        Dune::BlockVector<Dune::FieldVector<double,3> > displacement_;
        //Dune::BlockVector<Dune::FieldVector<double,6> > stress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > linstress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > outputstress_;// used in to avoid trouble with initialstress_
        Dune::BlockVector<Dune::FieldVector<double,6> > strain_;
        //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_;
        std::unique_ptr<Opm::Elasticity::VemElasticitySolver<Grid>> elacticitysolver_;
        const MechGridContext* mechCtx_{nullptr};
        //
        FractureMechHost<TypeTag> fracHost_{simulator_};
    };
}



#endif
