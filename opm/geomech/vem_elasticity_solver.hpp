//==============================================================================
//!
//! \file elasticity_upscale.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class
//!
//==============================================================================
#ifndef VEM_ELASTICITY_SOLVER_HPP_
#define VEM_ELASTICITY_SOLVER_HPP_


#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/TimingMacros.hpp>
#include <dune/common/fmatrix.hh>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/elasticity/shapefunctions.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/elasticity/asmhandler.hpp>
#include <opm/elasticity/boundarygrid.hh>
#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/elasticity_preconditioners.hpp>
#include <opm/elasticity/logutils.hpp>
#include <opm/elasticity/materials.hh>
#include <opm/elasticity/matrixops.hpp>
#include <opm/elasticity/meshcolorizer.hpp>
#include <opm/elasticity/mpc.hh>
#include <opm/elasticity/mortar_schur.hpp>
#include <opm/elasticity/mortar_utils.hpp>
#include <opm/elasticity/mortar_evaluator.hpp>
#include <opm/elasticity/mortar_schur_precond.hpp>
#include <opm/elasticity/uzawa_solver.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Schedule/BCProp.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
//#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
namespace Opm {
namespace Elasticity {


//! \brief The main driver class
template<class GridType>
class VemElasticitySolver
{
  public:
    //! \brief Dimension of our grid
    static const int dim = GridType::dimension;

    //! \brief A basic number
    typedef typename GridType::LeafGridView::ctype ctype;

    //! \brief A vectorial node value
    typedef Dune::FieldVector<double,dim> NodeValue;

    //! \brief A global coordinate
    typedef typename GridType::LeafGridView::template Codim<1>::Geometry::GlobalCoordinate GlobalCoordinate;

    //! \brief A set of indices
    typedef typename GridType::LeafGridView::IndexSet LeafIndexSet;

    //! \brief An iterator over grid cells
    typedef typename GridType::LeafGridView::template Codim<0>::Iterator LeafIterator;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > Matrix;
    //! \brief The linear operator
    ASMHandler<GridType> A;//NB NB names have to change

    //! \brief The solution vectors
    Vector u;//NB NB names have to change
    //! \brief The load vectors
    //Vector rhs;//NB NB names have to change

    //! \brief Main constructor
    //! \param[in] gv_ The grid to operate on
    //! \param[in] tol_ The tolerance to use when deciding whether or not a coordinate falls on a plane/line/point. \sa tol
    //! \param[in] Escale_ A scale value for E-moduluses to avoid numerical issues
    //! \param[in] file The eclipse grid file
    //! \param[in] rocklist If not blank, file is a rocklist
    //! \param[in] verbose If true, give verbose output
    VemElasticitySolver(const GridType& grid)
        : A(grid),grid_(grid)
    {
    }

    void setMaterial(const std::vector<std::shared_ptr<Material>>& materials_){
        materials = materials_;
    }
    void setMaterial(const std::vector<double>& ymodule,const std::vector<double>& pratio){
        ymodule_ = ymodule;
        pratio_ = pratio;
    }
    //! \brief Find boundary coordinates
    //! \param[out] min The miminum coordinates of the grid
    //! \param[out] max The maximum coordinates of the grid
    void initForAssembly(){}
    void fixNodes(const std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes){
        //void fixNodesVem(const std::vector<size_t>& fixed_nodes){
        std::vector<int> fixed_dof_ixs;
        std::vector<double> fixed_dof_values;
        for (int i = 0; i < int(bc_nodes.size()); ++i){
            int node = std::get<0>(bc_nodes[i]);
            const auto& mechbcvalue = std::get<1>(bc_nodes[i]);
            const auto& mask = mechbcvalue.fixeddir;
            const auto& disp = mechbcvalue.disp;
            for(int m = 0; m < 3; ++m){
                if(mask[m]){
                    fixed_dof_ixs.insert(fixed_dof_ixs.end(), {3*node + m});
                    fixed_dof_values.insert(fixed_dof_values.end(),{disp[m]});
                }
                //fixed_dof_ixs.insert(fixed_dof_ixs.end(), {3*node, 3*node+1, 3*node+2});
            }
        }
        const int num_fixed_dofs = (int)fixed_dof_ixs.size();


        dirichlet_={num_fixed_dofs, fixed_dof_ixs, fixed_dof_values};
    }
    void expandSolution(Vector& result, const Vector& u){
        //vector<double> displacements_full(prob.points.size(), nan("1"));
         // cout << prob.fixed_dof_ixs.size() << endl;
        std::fill(result.begin(),result.end(),nan("1"));
        const auto& fixed_dof_values = std::get<2>(dirichlet_);
        const auto& fixed_dof_ixs = std::get<1>(dirichlet_);
        //const auto& num_dof_dofs = std::get<0>(dirichlet_);
        for (int i = 0; i != int(fixed_dof_ixs.size()); ++i)
             result[fixed_dof_ixs[i]] = fixed_dof_values[i];

         for (int i = 0, cur_ix = 0; i != int(result.size()); ++i){
             if (std::isnan(result[i][0])) {
                 result[i] = u[cur_ix++];
             }
         }
    }
    //! \brief Assemble (optionally) stiffness matrix A and load vector
    //! \param[in] loadcase The strain load case. Set to -1 to skip
    //! \param[in] matrix Whether or not to assemble the matrix
    void assemble(const Vector& pressure, bool matrix,bool vector);


    //! \brief Solve Au = b for u
    //! \param[in] loadcase The load case to solve
    void solve();

    // //! \param[in] params The linear solver parameters
    void setupSolver(const Opm::PropertyTree& prm){
        OPM_TIMEBLOCK(setupLinearSolver);
        // bool parallel=false;
        // if(parallel){
        //     OPM_THROW(std::runtime_error,"Parallel for mechanics not implemented");
        //     const std::function<Vector()> weightsCalculator;
        //     std::size_t pressureIndex;
        //     //NB NB ! need other operator in parallel
        //     using ParOperatorType = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, Comm>;
        //     auto pop = std::make_unique<ParOperatorType>(A.getOperator(), comm);
        //     using FlexibleSolverType = Dune::FlexibleSolver<ParOperatorType>;
        //     tsolver = std::make_shared<FlexibleSolverType>(*pop, comm, prm,
        //                                                    weightsCalculator,
        //                                                     pressureIndex);
        // }else{
            std::size_t pressureIndex;//Dummy
            const std::function<Vector()> weightsCalculator;//Dummy
            using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
            auto sop = std::make_unique<SeqOperatorType>(A.getOperator());
            using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
            auto tsolver = std::make_unique<FlexibleSolverType>(*sop, prm,
                                                                weightsCalculator,
                                                                pressureIndex);
            sop_ = std::move(sop);
            tsolver_ = std::move(tsolver);
            //NB names have to change.
            //b = A.getLoadVector();
            //}
    }
    template<int comp>
    void averageStress(Dune::BlockVector<Dune::FieldVector<ctype,comp>>& sigmacells,
                       const Vector& uarg);
    void setBodyForce(double gravity){
        OPM_TIMEBLOCK(setBodyFrorce);
        const int num_cells = grid_.leafGridView().size(0); // entities of codim 0
        // assemble the mechanical system
        body_force_.resize(3*num_cells, 0.0);
        for(int i=0; i<num_cells; ++i){
            body_force_[3*i+2] = 2000*gravity;
        }
    }

    static void makeDuneMatrixCompressed(const std::vector<std::tuple<int, int, double>>& A_entries, Matrix& mat){
        OPM_TIMEBLOCK(makeDuneMatrixCompressed);
        // // hopefully consisten with matrix
        // // should be moved to initialization
        // std::vector< std::set<int> > rows;
        // int nrows=0;
        // int ncols=0;
        // {
        // OPM_TIMEBLOCK(findRowStructure);
        // for(auto matel: A_entries){
        //     int i = std::get<0>(matel);
        //     int j = std::get<1>(matel);
        //     nrows = std::max(nrows,i);
        //     ncols = std::max(ncols,j);
        // }
        // nrows = nrows+1;
        // ncols = ncols+1;
        // //assert(nrows==ncols);
        // rows.resize(nrows);
        // for(auto matel: A_entries){
        //     int i = std::get<0>(matel);
        //     int j = std::get<1>(matel);
        //     rows[i].insert(j);
        // }
        // }
        // // set up matrix structure
        // {
        // OPM_TIMEBLOCK(makeMatrix);
        // MatrixOps::fromAdjacency(mat, rows, nrows, ncols);
        // }
        {
            OPM_TIMEBLOCK(buildMatrixImplicite);
            //build mode need to be implicite and corredect matrix settings has to be done
            //mat = 0;
            for (auto matel: A_entries) {
                int i = std::get<0>(matel);
                int j = std::get<1>(matel);
                mat.entry(i,j)=0;
            }
            //Dune::CompressionStatistics<Matrix::size_type> stats = mat.compress();
            mat.compress();
            //std::cout << stats;
        }
        {
        OPM_TIMEBLOCK(setMatrixValues);
        mat = 0;
        for(auto matel:A_entries){
            int i = std::get<0>(matel);
            int j = std::get<1>(matel);
            double val = std::get<2>(matel);
            mat[i][j] += val;
        }
        }
    }
    void updateRhsWithGrad(const Vector& mechpot){
        OPM_TIMEBLOCK(updateRhsWithGrad);
         Vector& b = A.getLoadVector();
         b = 0;
         b.resize(rhs_force_.size());
         divmat_.mv(mechpot,b);
        // end initialization
        for(int i=0; i < rhs_force_.size(); ++i){
                b[i] += rhs_force_[i];
        }
    }

    void makeDuneSystemMatrix(const std::vector<std::tuple<int, int, double>>& A_entries){
        OPM_TIMEBLOCK(makeDuneSystemMatrix);
        // hopefully consisten with matrix
        // should be moved to initialization
        // std::vector< std::set<int> > rows;
        // int nrows=0;
        // int ncols=0;

        // {
        // OPM_TIMEBLOCK(makeStructure);
        // for (auto matel: A_entries) {
        //     int i = std::get<0>(matel);
        //     int j = std::get<1>(matel);
        //     nrows = std::max(nrows,i);
        //     ncols = std::max(ncols,j);
        // }
        // nrows = nrows+1;
        // ncols = ncols+1;
        // //assert(nrows==ncols);
        // rows.resize(nrows);
        // for(const auto& matel: A_entries){
        //     int i = std::get<0>(matel);
        //     int j = std::get<1>(matel);
        //     rows[i].insert(j);
        // }
        // }
        // {
        // OPM_TIMEBLOCK(matrixFromAdjacency);
        // MatrixOps::fromAdjacency(MAT, rows, nrows, ncols);
        // }
        {
            OPM_TIMEBLOCK(buildMatrixImplicite);
            int ncols = 0;
            int nrows = 0;
            for (auto matel: A_entries) {
                int i = std::get<0>(matel);
                int j = std::get<1>(matel);
                nrows = std::max(nrows,i);
                ncols = std::max(ncols,j);
            }
            nrows = nrows+1;
            ncols = ncols+1;
            Matrix& MAT =this->A.getOperator();
            MAT = 0;
            MAT.setBuildMode(Matrix::implicit);
            MAT.setImplicitBuildModeParameters (81, 0.4);
            MAT.setSize(nrows, ncols);
            makeDuneMatrixCompressed(A_entries,MAT);
            //    for (auto matel: A_entries) {
            //         int i = std::get<0>(matel);
            //         int j = std::get<1>(matel);
            //         MAT.entry(i,j)=0;
            //     }
            //     Dune::CompressionStatistics<Matrix::size_type> stats = MAT.compress();
            //     //std::cout << stats;
            // }
            // {
            // OPM_TIMEBLOCK(insertValues);
            // MAT = 0;
            // for(const auto& matel:A_entries){
            //     int i = std::get<0>(matel);
            //     int j = std::get<1>(matel);
            //     double val = std::get<2>(matel);
            //     A.addMatElement(i,j, val);// += val;
            // }
        }
    }
    void calculateStress(bool precalculated);
    void calculateStrain(bool precalculated);
    //const std::vector<std::array<double,6>>& stress(){return stress_;}
    const Dune::BlockVector<Dune::FieldVector<ctype,6>>& stress() const{return stress_;}
    const Dune::BlockVector<Dune::FieldVector<ctype,6>>& strain() const{return strain_;}

  private:
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;

    using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    std::unique_ptr<SeqOperatorType> sop_;

        //! \brief Linear solver
    typedef std::unique_ptr<Dune::InverseOperator<Vector, Vector> > SolverPtr;
    SolverPtr tsolver_;

    //! \brief An iterator over grid vertices
    typedef typename GridType::LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;

    //! \brief A reference to our grid
    const GridType& grid_;
    int num_cells_;
    std::vector<double> coords_;
    std::vector<int> num_cell_faces_, num_face_corners_,face_corners_;
    std::vector<int> idx_free_;

    std::vector<double> rhs_force_;
    //! \brief Vector holding material parameters for each active grid cell
    std::vector< std::shared_ptr<Material> > materials;
    std::vector<double> ymodule_;
    std::vector<double> pratio_;
    std::vector<double> body_force_;//= set_body_force(num_cells, bfcase);
    std::tuple<int, std::vector<int>, std::vector<double>> dirichlet_; //set_dirichlet(coords, dircase);

    Dune::BlockVector<Dune::FieldVector<ctype,6>> stress_;
    Dune::BlockVector<Dune::FieldVector<ctype,6>> strain_;
    Matrix stressmat_; // from all dofs (not eliminating bc) to cell
    Matrix strainmat_;
    Matrix divmat_;   // from cell pressure to active dofs
    //std::vector<std::array<double,6>> stress_;

};

}} // namespace Opm, Elasticity

#include "vem_elasticity_solver_impl.hpp"

#endif
