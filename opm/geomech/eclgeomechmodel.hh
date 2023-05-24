#ifndef OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#define OPM_ECLPROBLEM_GEOMECH_MODEL_HH
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
//#include <opm/elasticity/elasticity_preconditioners.hpp>
//#include <opm/elasticity/elasticity_upscale.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>
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
        void postSolve(GlobalEqVector& deltaX){
            std::cout << "Geomech dummy PostSolve Aux" << std::endl;
        }
        void addNeighbors(std::vector<NeighborSet>& neighbors) const{
            std::cout << "Geomech add neigbors" << std::endl;
        };
        void applyInitial(){
            std::cout << "Geomech applyInitial" << std::endl;
        };
        unsigned numDofs() const{return 0;};
        void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual){
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
            std::cout << "Geomech end iteration" << std::endl;
            size_t numDof = simulator_.model().numGridDof();
            const auto& problem = simulator_.problem();
            for(size_t dofIdx=0; dofIdx < numDof; ++dofIdx){    
                const auto& iq = simulator_.model().intensiveQuantities(dofIdx,0);
                const auto& fs = iq.fluidState();
                const auto& press = fs.pressure(waterPhaseIdx);
                const auto& biotcoef = problem.biotCoef(dofIdx);
                mechPotentialForce_[dofIdx] = (Toolbox::value(press) - problem.initPressure(dofIdx))*biotcoef;
            }
            // for now assemble and set up solver her
            
            if(first_solve_){
                bool do_matrix = true;//assemble matrix
                bool do_vector = true;//assemble matrix
                // set boundary
                elacticitysolver_.setBodyForce(0.0);
                elacticitysolver_.fixNodes(problem.fixedNodes()); 
                //
                elacticitysolver_.initForAssembly();
                elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector);
                Opm::PropertyTree prm("mechsolver.json");
                elacticitysolver_.setupSolver(prm);
                first_solve_ = false;
            }else{
                bool do_matrix = false;//assemble matrix
                bool do_vector = true;//assemble matrix
                //elacticitysolver_.A.initForAssembly();
                elacticitysolver_.assemble(mechPotentialForce_, do_matrix, do_vector);
            }    
            
            elacticitysolver_.solve();
            Opm::Elasticity::Vector field;
            const auto& grid = simulator_.vanguard().grid();
            const auto& gv = grid.leafGridView();
            static constexpr int dim = Grid::dimension;
            field.resize(grid.size(dim)*dim);
            elacticitysolver_.expandSolution(field,elacticitysolver_.u);
            
            this->makeDisplacement(field);
            // update variables used for output to resinsight
            // NB TO DO
            elacticitysolver_.calculateStress();
            stress_ = elacticitysolver_.stress();
            
            
            bool verbose = true;
            if(verbose){
                // debug output to matrixmaket format
                Dune::storeMatrixMarket(elacticitysolver_.A.getOperator(), "A.mtx");
                Dune::storeMatrixMarket(elacticitysolver_.A.getLoadVector(), "b.mtx");
                Dune::storeMatrixMarket(elacticitysolver_.u, "u.mtx");
                Dune::storeMatrixMarket(field, "field.mtx");
                Dune::storeMatrixMarket(mechPotentialForce_, "pressforce.mtx");
            }
            
        }       
        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            //serializer(tracerConcentration_);
            //serializer(wellTracerRate_);           
         }

        // used in eclproblemgeomech
        void init(bool restart){
            std::cout << "Geomech init" << std::endl;
            size_t numDof = simulator_.model().numGridDof();
            mechPotentialForce_.resize(numDof);
            celldisplacement_.resize(numDof);
            stress_.resize(numDof);
            const auto& gv = simulator_.vanguard().grid().leafGridView();
            displacement_.resize(gv.indexSet().size(3));
            celldisplacement_.resize(gv.indexSet().size(3));
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
        const double mechPotentialForce(unsigned globalDofIdx) const
        {
            return mechPotentialForce_[globalDofIdx];
        }
        const double disp(unsigned globalDofIdx, unsigned dim) const
        {
            return celldisplacement_[globalDofIdx][dim];
        }
        const double stress(unsigned globalDofIdx, unsigned dim) const
        {
            return stress_[globalDofIdx][dim];
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
                auto cellindex = gv.indexSet().index(cell);
                const auto& vertices = Dune::subEntities(cell, Dune::Codim<Grid::dimension>{});
                for (const auto& vertex : vertices){
                    auto nodeidex = gv.indexSet().index(vertex);
                    for(int k=0; k < dim; ++k){
                        celldisplacement_[cellindex][nodeidex] += displacement_[nodeidex][k] ;
                    }
                }
                celldisplacement_[cellindex] /= vertices.size(); 
            }
        }
        
    private:
        bool first_solve_;
        Simulator& simulator_;
        Dune::BlockVector<Dune::FieldVector<double,1>> mechPotentialForce_;
        //Dune::BlockVector<Dune::FieldVector<double,1> > solution_;
        Dune::BlockVector<Dune::FieldVector<double,3> > celldisplacement_;
        Dune::BlockVector<Dune::FieldVector<double,3> > displacement_;
        Dune::BlockVector<Dune::FieldVector<double,6> > stress_;//NB is also stored in esolver
        //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_;
        Opm::Elasticity::VemElasticitySolver<Grid> elacticitysolver_;
    };
}



#endif
