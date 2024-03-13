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

//#include <opm/geomech/ElasticitySolverUpscale.hpp>
namespace Opm{
    template<typename TypeTag>
    class EclGeoMechModel : public BaseAuxiliaryModule<TypeTag>
    {
        //using Parent = BaseAuxiliaryModule<TypeTag>;
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
            OPM_TIMEBLOCK(endTimeStepMech);
            std::cout << "Geomech end iteration" << std::endl;
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
                double fac = (1-2*pratio)/(1-pratio);
                double pcoeff = poelCoef*fac;
                mechPotentialForce_[dofIdx] = diffpress*pcoeff;
                bool thermal_expansion = getPropValue<TypeTag, Properties::EnableEnergy>();

                if(thermal_expansion){
                    OPM_TIMEBLOCK(addTermalParametersToMech);
                    const auto& temp = fs.temperature(waterPhaseIdx);//NB all phases have equal temperature
                    const auto& thelcoef = problem.thelCoef(dofIdx);
                    // const auto& termExpr = problem.termExpr(dofIdx); //NB not used

                    double tcoeff = thelcoef*fac;//+ youngs*tempExp;

                    // assume difftemp = 0 for non termal runs
                    double difftemp = (Toolbox::value(temp) - problem.initTemperature(dofIdx));
                    mechPotentialForce_[dofIdx] += difftemp*tcoeff;
                    mechPotentialTempForce_[dofIdx] = difftemp*tcoeff;
                }
                //NB check sign !!
                mechPotentialForce_[dofIdx] *= 1.0;
            }
            // for now assemble and set up solver her

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
                stress_[cellindex] = linstress[cellindex];
                strain_[cellindex] = linstrain[cellindex];
                delstress_[cellindex] = linstress[cellindex];
            }
            size_t lsdim = 6;
            for(size_t i = 0; i < stress_.size(); ++i){
                for(size_t j = 0; j < lsdim; ++j){
                    stress_[i][j] += problem.initStress(i,j);
                }
            }
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
            celldisplacement_.resize(numDof);
            stress_.resize(numDof);
            delstress_.resize(numDof);
            strain_.resize(numDof);
            const auto& gv = simulator_.vanguard().grid().leafGridView();
            displacement_.resize(gv.indexSet().size(3));
        };
        double pressureDiff(unsigned dofIx) const{
            return mechPotentialForce_[dofIx];
        }
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
        const double& disp(unsigned globalDofIdx, unsigned dim) const
        {
            return celldisplacement_[globalDofIdx][dim];
        }
        const double& stress(unsigned globalDofIdx, unsigned dim) const
        {
            return stress_[globalDofIdx][dim];
        }
        const double& delstress(unsigned globalDofIdx, unsigned dim) const
        {
            return delstress_[globalDofIdx][dim];
        }
        const double& strain(unsigned globalDofIdx, unsigned dim) const
        {
            return strain_[globalDofIdx][dim];
        }

        void setStress(const Dune::BlockVector<Dune::FieldVector<double,6> >& stress){
            stress_ = stress;
        }
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

    private:
        bool first_solve_;
        Simulator& simulator_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialPressForce_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialTempForce_;
        //Dune::BlockVector<Dune::FieldVector<double,1> > solution_;
        Dune::BlockVector<Dune::FieldVector<double,3> > celldisplacement_;
        Dune::BlockVector<Dune::FieldVector<double,3> > displacement_;
        Dune::BlockVector<Dune::FieldVector<double,6> > stress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > delstress_;//NB is also stored in esolver
        Dune::BlockVector<Dune::FieldVector<double,6> > strain_;
        //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_;
        Opm::Elasticity::VemElasticitySolver<Grid> elacticitysolver_;
    };
}



#endif
