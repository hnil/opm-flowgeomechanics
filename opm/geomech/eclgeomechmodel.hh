#ifndef OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#define OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>

#include <opm/geomech/FlowGeomechLinearSolverParameters.hpp>
#include <opm/geomech/FractureModel.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

namespace Opm{
    template<typename TypeTag>
    class EclGeoMechModel : public BaseAuxiliaryModule<TypeTag>
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
        enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;
    public:
        EclGeoMechModel(Simulator& simulator):
            //           Parent(simulator)
            first_solve_(true),
            write_system_(false),
            reduce_boundary_(false),
            simulator_(simulator),
            elacticitysolver_(simulator.vanguard().grid())
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
            std::cout << "Geomech begin iteration" << std::endl;
        }
        void endIteration(){
            //Parent::endIteration();
            std::cout << "Geomech end iteration" << std::endl;
        }
        void beginTimeStep(){
            //Parent::beginIteration();
            std::cout << "Geomech begin iteration" << std::endl;
        }
        void endTimeStep(){
            // always do post solve
            std::cout << "Geomech model endstimeStep" << std::endl;
            this->solveGeomechAndFracture();
        }
        void solveGeomechAndFracture(){
            //Parent::endIteration();
            const auto& problem = simulator_.problem();
            this->solveGeomechanics();
            if(problem.hasFractures()){
                this->solveFractures();
            }   
        }

       void solveFractures(){
            OPM_TIMEBLOCK(solveFractures);
            int reportStepIdx = simulator_.episodeIndex();
            const auto& schedule =  this->simulator_.vanguard().schedule();
            int end_step = schedule.size() - 1;
            bool no_seeds = schedule[end_step].wseed().empty();
            if(!no_seeds){
                //std::cout << "No fracture seeds found, on this step " << reportStepIdx << std::endl;
                std::cout << "Fracture seeds found, on this step " << std::endl;
            }else{
                std::cout << "No fracture seeds found, on this step "  << std::endl;
            }
            if(fracturemodel_){
                std::cout << "Fracture model already initialized, solving fractures using previous fractures" << std::endl;
            }

            if(!no_seeds && !fracturemodel_){
                    std::cout << "Fracture model not initialized, initializing now. report step" <<  reportStepIdx << std::endl;
                    const auto& problem = simulator_.problem();
                    //NB could probably be moved to some initialization
                    // let fracture contain all wells
                    Opm::PropertyTree param = problem.getFractureParam();
                    include_fracture_contributions_ = param.get<bool>("include_fracture_contributions");
                    //param.read("fractureparam.json");
                    //const auto& schedule =  this->simulator_.vanguard().schedule();
                    //int reportStepIdx = simulator_.episodeIndex();
                    // take all wells and perforations
                    //int end_step = schedule.size() - 1;
                    const std::vector<Opm::Well>& wells = problem.wellModel().getLocalWells(end_step);
                    //const std::vector<Opm::Well>& wells = schedule.getWells(reportStepIdx);
                    //const Opm::EclipseGrid& eclgrid = simulator_.vanguard().eclState().getInputGrid();
                    const auto& grid = simulator_.vanguard().grid();
                    std::string outputDir = Parameters::Get<Parameters::OutputDir>();
                    std::string caseName  = simulator_.vanguard().caseName();
                    param.put("outputdir", outputDir);
                    param.put("casename", caseName);
                    //
                    


                    fracturemodel_ = std::make_unique<FractureModel>(grid,
                                                                     wells,
                                                                     param
                        );
                    // not to get the reservoir properties along the well before initialising the well
                    // most important stress
                    fracturemodel_->updateReservoirWellProperties<TypeTag,Simulator>(simulator_);
                    // add fractures along the wells
                    //fracturemodel_->addFractures(schedule[reportStepIdx]);
                    fracturemodel_->addFractures(schedule[end_step]);

                    fracturemodel_->updateFractureReservoirCells();
                    fracturemodel_->initReservoirProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->updateReservoirAndWellProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->initFractureStates();
                }
                // get reservoir properties on fractures
                // simulator need
                if(fracturemodel_){
                    std::cout << "Frac modelfound, updating reservoir properties and solving fractures" << std::endl;
                    fracturemodel_->updateReservoirAndWellProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->solve<TypeTag, Simulator>(simulator_);
                }else{
                    std::cout << "Fracture model not initialized, not solving fractures" << std::endl;
                }
                // copy from apply action
       }
       
        
        void writeFractureSolution(){
            const auto& problem = simulator_.problem();
            if(problem.hasFractures() && fracturemodel_){
                // write first solution in standard format
                // this may ad some extra output of static variables
                //int reportStepIdx = simulator_.episodeIndex();
                if(first_output_){
                    //fracturemodel_->write(reportStepIdx);
                    // hack to get correct number of fracture output
                    fracturemodel_->writemulti(0.0);
                    first_output_ = false;
                } 
                double time = simulator_.time();
                fracturemodel_->writemulti(time);
            }

        }


        std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname){
            if(fracturemodel_){
                return fracturemodel_->getExtraWellIndices(wellname);
            }else{
                return std::vector<RuntimePerforation>();
            }
        }
      
        void updatePotentialForces(){
            std::cout << "Update Forces" << std::endl;
            size_t numDof = simulator_.model().numGridDof();
            const auto& problem = simulator_.problem();
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){
                const auto& iq = simulator_.model().intensiveQuantities(dofIdx,0);
                // pressure part
                const auto& fs = iq.fluidState();
                const auto& press = fs.pressure(waterPhaseIdx);
                // const auto& biotcoef = problem.biotCoef(dofIdx); //NB not used
                //thermal part
                //Properties::EnableTemperature
                const auto& poelCoef = problem.poelCoef(dofIdx);
                double pressure = Toolbox::value(press);
                double diffpress = (pressure - problem.initPressure(dofIdx));
		        const auto& pratio = problem.pRatio(dofIdx);
                double fac = (1-pratio)/(1-2*pratio);
                double pcoeff = poelCoef*fac;
		//assert pcoeff == biot
                pressure_[dofIdx] = pressure;
                mechPotentialForce_[dofIdx] = diffpress*pcoeff;
                mechPotentialPressForce_[dofIdx] = diffpress*pcoeff;
                assert(pcoeff<=1.0);
                mechPotentialPressForceFracture_[dofIdx] = diffpress*(1.0-pcoeff);
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();

                if(thermal_expansion){
                    OPM_TIMEBLOCK(addTermalParametersToMech);
                    const auto& temp = fs.temperature(waterPhaseIdx);//NB all phases have equal temperature
                    const auto& thelcoef = problem.thelCoef(dofIdx);
                    // const auto& termExpr = problem.termExpr(dofIdx); //NB not used
		    // tcoeff = (youngs*tempExp/(1-pratio))*fac;
                    double tcoeff = thelcoef*fac;
                    // assume difftemp = 0 for non termal runs
                    double difftemp = (Toolbox::value(temp) - problem.initTemperature(dofIdx));
                    mechPotentialForce_[dofIdx] += difftemp*tcoeff;
                    mechPotentialTempForce_[dofIdx] = difftemp*tcoeff;
                }
                //NB check sign !!
                mechPotentialForce_[dofIdx] *= 1.0;
            }
        }
        void setupMechSolver(){
                const auto& problem = simulator_.problem();
                Opm::PropertyTree param = problem.getFractureParam();
                reduce_boundary_ = param.get<bool>("reduce_boundary");
                OPM_TIMEBLOCK(SetupMechSolver);
                bool do_matrix = true;//assemble matrix
                bool do_vector = true;//assemble matrix
                // set boundary
                elacticitysolver_.setBodyForce(0.0);
                elacticitysolver_.fixNodes(problem.bcNodes());
                //
                elacticitysolver_.initForAssembly();
                elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector,reduce_boundary_);
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
                elacticitysolver_.setupSolver(prm);
                elacticitysolver_.comm()->communicator().barrier();
                first_solve_ = false;
                write_system_ = prm.get<int>("verbosity") > 10;
        }
        void writeMechSystem(){
            OPM_TIMEBLOCK(WriteMechSystem);
                    const auto& problem = simulator_.problem();
                    Opm::Helper::writeMechSystem(simulator_,
                    elacticitysolver_.A.getOperator(),
                    elacticitysolver_.A.getLoadVector(),
                    elacticitysolver_.comm());
                    {
                        int num_points = simulator_.vanguard().grid().size(3);
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
                                                elacticitysolver_.comm());
                    }
        }

        void calculateOutputQuantitiesMech(){
            OPM_TIMEBLOCK(CalculateOutputQuantitesMech);
            Opm::Elasticity::Vector field;
            const auto& grid = simulator_.vanguard().grid();
            const auto& gv = grid.leafGridView();
            static constexpr int dim = Grid::dimension;
            field.resize(grid.size(dim)*dim);
            if(reduce_boundary_){
                elacticitysolver_.expandSolution(field,elacticitysolver_.u);
            }else{
              assert(field.size() == elacticitysolver_.u.size());
              field =  elacticitysolver_.u;   
            }

            this->makeDisplacement(field);
            // update variables used for output to resinsight
            // NB TO DO
            {
            OPM_TIMEBLOCK(calculateStress);
            elacticitysolver_.calculateStressPrecomputed(field);
            elacticitysolver_.calculateStrainPrecomputed(field);
            }
            const auto& linstress = elacticitysolver_.stress();
            const auto& linstrain = elacticitysolver_.strain();

            for (const auto& cell: elements(gv)){
                auto cellindex = simulator_.problem().elementMapper().index(cell);
                // add initial stress
                assert(cellindex == gv.indexSet().index(cell));
                //auto cellindex2 = gv.indexSet().index(cell);
                //stress_[cellindex] = linstress[cellindex];
                strain_[cellindex] = linstrain[cellindex];
                linstress_[cellindex] = linstress[cellindex];
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
                Dune::storeMatrixMarket(elacticitysolver_.A.getOperator(), "A.mtx");
                Dune::storeMatrixMarket(elacticitysolver_.A.getLoadVector(), "b.mtx");
                Dune::storeMatrixMarket(elacticitysolver_.u, "u.mtx");
                Dune::storeMatrixMarket(field, "field.mtx");
                Dune::storeMatrixMarket(mechPotentialForce_, "pressforce.mtx");
            }
        }


        void setupAndUpdateGemechanics(){
            OPM_TIMEBLOCK(endTimeStepMech);
            this->updatePotentialForces();
            // for now assemble and set up solver her
            const auto& problem = simulator_.problem();
            if(first_solve_){
                this->setupMechSolver();
            }
            //else
            {
                // reset the rhs even in the first iteration maybe bug in rhs for reduce_boundary=false;
                OPM_TIMEBLOCK(AssembleRhs);
                // need "static boundary conditions is changing"
                //bool do_matrix = false;//assemble matrix
                //bool do_vector = true;//assemble matrix
                //elacticitysolver_.A.initForAssembly();
                //elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector);
                // need precomputed divgrad operator
                elacticitysolver_.updateRhsWithGrad(mechPotentialForce_);
            }    
        }
        void solveGeomechanics(){
            setupAndUpdateGemechanics();
            {
                OPM_TIMEBLOCK(SolveMechanicalSystem);
                elacticitysolver_.solve();
                if(write_system_){
                    this->writeMechSystem();
                }
            }
            this->calculateOutputQuantitiesMech();// and properties used for fracturing
        }
        template<class Serializer>
        void serializeOp(Serializer&)
        {
            //serializer(tracerConcentration_);
            //serializer(wellTracerRate_);
         }

        // used in eclproblemgeomech
        void init(bool /*restart*/){
            std::cout << "Geomech init" << std::endl;
            size_t numDof = simulator_.model().numGridDof();
            pressure_.resize(numDof);
            mechPotentialForce_.resize(numDof);
            mechPotentialTempForce_.resize(numDof);
            mechPotentialPressForce_.resize(numDof);
            mechPotentialPressForceFracture_.resize(numDof);
            // hopefully temperature and pressure initilized
            celldisplacement_.resize(numDof);
            std::fill(celldisplacement_.begin(),celldisplacement_.end(),0.0);
            //stress_.resize(numDof);
            linstress_.resize(numDof);
            std::fill(linstress_.begin(),linstress_.end(),0.0);
            strain_.resize(numDof);
            std::fill(strain_.begin(),strain_.end(),0.0);
            const auto& gv = simulator_.vanguard().grid().leafGridView();
            displacement_.resize(gv.indexSet().size(3));
        };

        void setMaterial(const std::vector<std::shared_ptr<Opm::Elasticity::Material>>& materials){
            elacticitysolver_.setMaterial(materials);
        }
        void setMaterial(const std::vector<double>& ymodule,const std::vector<double>& pratio){
            elacticitysolver_.setMaterial(ymodule,pratio);
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
        const Dune::FieldVector<double,3> disp(size_t globalIdx,bool with_fracture = false) const{
            auto disp =  celldisplacement_[globalIdx];
            if(include_fracture_contributions_ && with_fracture){
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    size_t globalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracturemodel_){
                        disp += fracturemodel_->disp(obs);
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
	  return linstress_[globalIdx];
        }

        const SymTensor effstress(size_t globalIdx) const{
	  // make stress in with positive with compression
	  return -1.0*linstress_[globalIdx];
        }

        const SymTensor strain(size_t globalIdx,bool with_fracture = false) const{
            auto strain = strain_[globalIdx];
            if(include_fracture_contributions_ && with_fracture){
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    size_t globalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracturemodel_){
                        strain += fracturemodel_->strain(obs);
                    }
                }
            }
            return strain_[globalIdx];
        }   

        const SymTensor stress(size_t globalIdx,bool with_fracture = false) const{
            Dune::FieldVector<double,6> effStress = this->effstress(globalIdx);
            effStress += simulator_.problem().initStress(globalIdx);
            double effPress = this->mechPotentialForce(globalIdx);
            for(int i=0; i < 3; ++i){
                effStress[i] += effPress;

            }
            if(include_fracture_contributions_ && with_fracture){
                for(auto& elem: Dune::elements(simulator_.vanguard().grid().leafGridView())){
                    size_t globalIdx = simulator_.problem().elementMapper().index(elem);
                    auto geom = elem.geometry();
                    auto center = geom.center();
                    Dune::FieldVector<double,3> obs = {center[0],center[1],center[2]};
                    // check if this is correct stress
                    if(fracturemodel_){
                        effStress += fracturemodel_->stress(obs);
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
	   const auto& press = Toolbox::value(fs.pressure(waterPhaseIdx));
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
            const auto& grid = simulator_.vanguard().grid();
            const auto& gv = grid.leafGridView();
            int dim = 3;
            for (const auto& vertex : Dune::vertices(gv)){
                auto index = gv.indexSet().index(vertex);
                for(int k=0; k < dim; ++k){
                    displacement_[index][k] = field[index*dim+k];
                }
            }
            for (const auto& cell: elements(gv)){
                auto cellindex = simulator_.problem().elementMapper().index(cell);
                assert(cellindex== gv.indexSet().index(cell));
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

        bool fractureModelActive() const{
            if(!fracturemodel_){                
                return false;
            }else{
                return true;
            }           
        }
        const FractureModel& fractureModel() const{
            if(!fracturemodel_){
                std::cout << "Fracture model not initialized, returning nullptr" << std::endl;
                throw std::runtime_error("Fracture model not initialized");
            }
            return *fracturemodel_;
        }
    private:
        bool first_solve_{true};
        bool first_output_{true};
        bool write_system_{false};
        bool reduce_boundary_{false};
        bool include_fracture_contributions_{false};
        Simulator& simulator_;

        Dune::BlockVector<Dune::FieldVector<double,1>> pressure_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForceFracture_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialTempForce_;
        //Dune::BlockVector<Dune::FieldVector<double,1> > solution_;
        Dune::BlockVector<Dune::FieldVector<double,3> > celldisplacement_;
        Dune::BlockVector<Dune::FieldVector<double,3> > displacement_;
        //Dune::BlockVector<Dune::FieldVector<double,6> > stress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > linstress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > strain_;
        //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_;
        Opm::Elasticity::VemElasticitySolver<Grid> elacticitysolver_;
        //
        std::unique_ptr<FractureModel> fracturemodel_;
    };
}



#endif
