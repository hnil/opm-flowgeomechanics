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
#ifndef ELASTICITY_SOLVER_HPP_
#define ELASTICITY_SOLVER_HPP_


#include <opm/common/utility/platform_dependent/disable_warnings.h>

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
#include <opm/simulators/linalg/FlexibleSolver.hpp>
//#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
namespace Opm {
namespace Elasticity {


//! \brief The main driver class
template<class GridType>
class ElasticitySolver
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
    ElasticitySolver(const GridType& gv_)
        :  A(gv_), gv(gv_), tol(1e-10), Escale(1.0), E(gv_), verbose(true),
           color(gv_)
    {
    }
    ElasticitySolver(const GridType& gv_, ctype tol_, ctype Escale_, 
                     bool verbose_)
        :  A(gv_), gv(gv_), tol(tol_), Escale(Escale_), E(gv_), verbose(verbose_),
           color(gv_)
    {
    }
    
    void setMaterial(const std::vector<std::shared_ptr<Material>>& materials_){
        materials = materials_;
    }
    void setMaterial(const std::vector<double>& ymodule,const std::vector<double>& pratio){
        materials.resize(0);
        for(size_t i=0; i < ymodule.size(); ++i){
            using IsoMat = Opm::Elasticity::Isotropic;
            if(pratio[i]>0.5 || pratio[i] < 0){
                OPM_THROW(std::runtime_error,"Pratio not valid");
            }
            materials.push_back(std::make_shared<IsoMat>(i,ymodule[i],pratio[i]));
        }
    }
    //! \brief Find boundary coordinates
    //! \param[out] min The miminum coordinates of the grid
    //! \param[out] max The maximum coordinates of the grid
    void findBoundaries(double* min, double* max);
    void initForAssembly(){A.initForAssembly();};
    void fixNodes(const std::vector<size_t>& fixed_nodes);
    void expandSolution(Vector& result, const Vector& u){A.expandSolution(result,u);};
    //! \brief Assemble (optionally) stiffness matrix A and load vector
    //! \param[in] loadcase The strain load case. Set to -1 to skip
    //! \param[in] matrix Whether or not to assemble the matrix
    void assemble(const Vector& pressure, bool matrix,bool vector);
    
    //! \brief Solve Au = b for u
    //! \param[in] loadcase The load case to solve
    void solve();

    // //! \param[in] params The linear solver parameters
    void setupSolver(const Opm::PropertyTree& prm){
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
  private:
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;

    using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    //using AbstractPrecondType = Dune::PreconditionerWithUpdate<Vector, Vector>;
    //std::unique_ptr<AbstractPrecondType> sop_;
    std::unique_ptr<SeqOperatorType> sop_;

        //! \brief Linear solver
    typedef std::unique_ptr<Dune::InverseOperator<Vector, Vector> > SolverPtr;
    SolverPtr tsolver_;

    //! \brief Matrix adaptor for the elasticity block
    //std::shared_ptr<Operator> op;

    
    //! \brief An iterator over grid vertices
    typedef typename GridType::LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;

    //! \brief A reference to our grid
    const GridType& gv;

    //! \brief Tolerance used to decide whether or not a coordinate falls on a plane/line/point.
    ctype tol;

    //! \brief Minimum E-modulus (scaling factor)
    ctype Escale;

    //! \brief Minimum real E for materials
    ctype Emin;

    //! \brief Vector holding material parameters for each active grid cell
    std::vector< std::shared_ptr<Material> > materials;

    //! \brief Extract the vertices on a given face
    //! \param[in] dir The direction of the face normal
    //! \param[in] coord The coordinate of the face plane
    //! \returns A vector holding the matching vertices
    std::vector<BoundaryGrid::Vertex> extractFace(Direction dir, ctype coord);

    //! \brief An enumeration used to indicate which side to extract from a cube
    enum SIDE {
      LEFT,
      RIGHT
    };

    //! \brief Extract a quad grid over a given face
    //! \param[in] dir The direction of the face normal
    //! \param[in] coord the coordinate of the face plance
    //! \param[in] side Extract left or right side
    //! \param[in] dc If true, order vertices in dune convention
    //! \returns A quad grid spanning the face
    BoundaryGrid extractMasterFace(Direction dir, ctype coord,
                                   SIDE side=LEFT, bool dc=false);

    //! \brief Find and establish master/slave grids (MPC)
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    void determineSideFaces(const double* min, const double* max);


    //! \brief Fix the DOFs in a given point on the grid
    //! \param[in] dir The coordinate direction to fix in
    //! \param[in] coord The coordinates of the node to fix
    //! \param[in] value The values to fix the given DOFs to
    void fixPoint(Direction dir, GlobalCoordinate coord,
                  const NodeValue& value = NodeValue(0));

    //! \brief Fix the DOFs in a given line on the grid
    //! \param[in] dir The coordinate direction to fix in
    //! \param[in] x The first coordinate of the line
    //! \param[in] y The second coordinate of the line
    //! \param[in] value The values to fix the given DOFs to
    void fixLine(Direction dir, ctype x, ctype y,
                 const NodeValue& value = NodeValue(0));
    

    //! \brief Elasticity helper class
    Elasticity<GridType> E;

    //! \brief Verbose output
    bool verbose;

    //! \brief Mesh colorizer used with multithreaded assembly
    MeshColorizer<GridType> color;
};

}} // namespace Opm, Elasticity

#include "elasticity_solver_impl.hpp"

#endif
