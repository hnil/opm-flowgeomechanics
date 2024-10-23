#ifndef OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#define OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
//#include <opm/elasticity/elasticity_preconditioners.hpp>
//#include <opm/elasticity/elasticity_upscale.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>

#include <opm/geomech/FlowGeomechLinearSolverParameters.hpp>
#include <opm/geomech/FractureModel.hpp>
//#include <opm/geomech/ElasticitySolverUpscale.hpp>
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
            //Parent::endIteration();
            const auto& problem = simulator_.problem();
            this->solveGeomechanics();
            if(problem.hasFractures()){
                if(!fracturemodel_){
                    //NB could probably be moved to some initialization
                    // let fracture contain all wells
                    Opm::PropertyTree param = problem.getFractureParam();
                    //param.read("fractureparam.json");
                    const auto& schedule =  this->simulator_.vanguard().schedule();
                    int reportStepIdx = simulator_.episodeIndex();
                    const std::vector<Opm::Well>& wells = schedule.getWells(reportStepIdx);
                    const Opm::EclipseGrid& eclgrid = simulator_.vanguard().eclState().getInputGrid();
                    const auto& grid = simulator_.vanguard().grid();
                    std::string outputDir = Parameters::Get<Parameters::OutputDir>();
                    std::string caseName  = simulator_.vanguard().caseName();
                    param.put("outputdir", outputDir);
                    param.put("casename", caseName);

                    fracturemodel_ = std::make_unique<FractureModel>(grid,
                                                                     wells,
                                                                     eclgrid,
                                                                     param,
                                                                     /*default fracture*/false
                        );
                    // not to get the reservoir properties along the well before initialising the well
                    // most important stress
                    fracturemodel_->updateReservoirWellProperties<TypeTag,Simulator>(simulator_);
                    // add fractures along the wells
                    fracturemodel_->addFractures();

                    fracturemodel_->updateFractureReservoirCells(grid,eclgrid);
                    fracturemodel_->initReservoirProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->updateReservoirProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->initFractureStates();
                }
                // get reservoir properties on fractures
                // simulator need
                fracturemodel_->updateReservoirProperties<TypeTag,Simulator>(simulator_);
                fracturemodel_->solve();
                // copy from apply action
            }
        }

        void writeFractureSolution(){
            const auto& problem = simulator_.problem();
            if(problem.hasFractures()){
                // write first solution in standard format
                int reportStepIdx = simulator_.episodeIndex();
                if(reportStepIdx==1){
                    fracturemodel_->write(reportStepIdx);
                }
                double time = simulator_.time();
                fracturemodel_->writemulti(time);
            }

        }

        std::vector<std::tuple<int,double,double>> getExtraWellIndices(std::string wellname){
            return fracturemodel_->getExtraWellIndices(wellname);
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
                double diffpress = (Toolbox::value(press) - problem.initPressure(dofIdx));
		const auto& pratio = problem.pRatio(dofIdx);
                double fac = (1-pratio)/(1-2*pratio);
                double pcoeff = poelCoef*fac;
		//assert pcoeff == biot
                mechPotentialForce_[dofIdx] = diffpress*pcoeff;
                mechPotentialPressForce_[dofIdx] = diffpress*pcoeff;
                assert(pcoeff<1.0);
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

        void solveGeomechanics(){
            OPM_TIMEBLOCK(endTimeStepMech);
            this->updatePotentialForces();
            // for now assemble and set up solver her
            const auto& problem = simulator_.problem();
            if(first_solve_){
                OPM_TIMEBLOCK(SetupMechSolver);
                bool do_matrix = true;//assemble matrix
                bool do_vector = true;//assemble matrix
                // set boundary
                elacticitysolver_.setBodyForce(0.0);
                elacticitysolver_.fixNodes(problem.bcNodes());
                //
                elacticitysolver_.initForAssembly();
                elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector);
                FlowLinearSolverParametersGeoMech p;
                p.init<TypeTag>();
                // Print parameters to PRT/DBG logs.

                PropertyTree prm = setupPropertyTree(p, true, true);
                if (p.linear_solver_print_json_definition_) {
                    std::ostringstream os;
                    os << "Property tree for mech linear solver:\n";
                    prm.write_json(os, true);
                    OpmLog::note(os.str());
                }
                elacticitysolver_.setupSolver(prm);
                first_solve_ = false;
            }else{
                OPM_TIMEBLOCK(AssembleRhs);

                // need "static boundary conditions is changing"
                //bool do_matrix = false;//assemble matrix
                //bool do_vector = true;//assemble matrix
                //elacticitysolver_.A.initForAssembly();
                //elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector);

                // need precomputed divgrad operator
                elacticitysolver_.updateRhsWithGrad(mechPotentialForce_);
            }

            {
                OPM_TIMEBLOCK(SolveMechanicalSystem);
                elacticitysolver_.solve();
            }
            {
            OPM_TIMEBLOCK(CalculateOutputQuantitesMech);
            Opm::Elasticity::Vector field;
            const auto& grid = simulator_.vanguard().grid();
            const auto& gv = grid.leafGridView();
            static constexpr int dim = Grid::dimension;
            field.resize(grid.size(dim)*dim);
            elacticitysolver_.expandSolution(field,elacticitysolver_.u);

            this->makeDisplacement(field);
            // update variables used for output to resinsight
            // NB TO DO
            {
            OPM_TIMEBLOCK(calculateStress);
            elacticitysolver_.calculateStress(true);
            elacticitysolver_.calculateStrain(true);
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
        const double& mechPotentialTempForce(unsigned globalDofIdx) const
        {
            return mechPotentialTempForce_[globalDofIdx];
        }
        const double& mechPotentialPressForce(unsigned globalDofIdx) const
        {
            return mechPotentialPressForce_[globalDofIdx];
        }
        const Dune::FieldVector<double,3>& disp(size_t globalIdx) const{
            return celldisplacement_[globalIdx];
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

        const SymTensor& strain(size_t globalIdx) const{
            return strain_[globalIdx];
        }

        const SymTensor stress(size_t globalIdx) const{
            Dune::FieldVector<double,6> effStress = this->effstress(globalIdx);
            effStress += simulator_.problem().initStress(globalIdx);
            double effPress = this->mechPotentialForce(globalIdx);
            for(int i=0; i < 3; ++i){
                effStress[i] += effPress;

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

        const double& disp(unsigned globalDofIdx, unsigned dim) const
        {
            return celldisplacement_[globalDofIdx][dim];
        }
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
        const FractureModel& fractureModel() const{return *fracturemodel_;}
    private:
        bool first_solve_;
        Simulator& simulator_;
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
