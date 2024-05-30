#include <dune/common/exceptions.hh>
#include <dune/istl/solvercategory.hh>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

#include <dune/istl/matrixmarket.hh>

//@@ there must be a more correct way to ensure UMFpack is included here
#define HAVE_SUITESPARSE_UMFPACK 1
#include <dune/istl/umfpack.hh>

//#include "opm/geomech/Fracture.hpp"
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/istl/io.hh> // needed for printSparseMatrix??


#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include "opm/geomech/DiscreteDisplacement.hpp"


namespace {
template<class M, class X, class Y> class ReducedMatrixAdapter; // forward declaration
  
using Grid = Dune::FoamGrid<2, 3>;
using Vector = Dune::BlockVector<double>;
using HTrans = std::tuple<size_t,size_t, double, double>;
  
using FullMatrix = Dune::DynamicMatrix<double>; //Dune::Matrix<double, std::allocator<double> >;
using SparseMatrix = Dune::BCRSMatrix<double>;
using EquationSystem = std::tuple<std::shared_ptr<FullMatrix>,
                                  std::shared_ptr<SparseMatrix>,
                                  std::vector<HTrans>>;

using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
using RMAdapter = ReducedMatrixAdapter<SparseMatrix, Vector, Vector>;
// ----------------------------------------------------------------------------
std::unique_ptr<const Grid> readgrid(const char* const name)
// ----------------------------------------------------------------------------  
{
  std::unique_ptr<const Grid> result;
  try {
    result = Dune::GmshReader<Grid>::read(name); // unique_ptr
  } catch (const Dune::Exception& e) {
    std::cout << "Unable to read file: " << name << std::endl;
  }
  return result;
}

// ----------------------------------------------------------------------------
template<int N>
std::array<size_t, N> n_closest(const Grid& grid)
// ----------------------------------------------------------------------------
{
  using Elem = std::tuple<size_t, double>;
  std::vector<Elem> distances;
  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto center = element.geometry().center();
    distances.push_back({eIdx, center.two_norm()});
  }

  std::sort(distances.begin(), distances.end(),
            [](Elem& a, Elem& b) {return std::get<1>(a) < std::get<1>(b);});
  
  std::array<size_t, N> result;
  for (int i = 0; i != N; ++i)
    result[i] = std::get<0>(distances[i]);
  
  return result;
}

// ----------------------------------------------------------------------------  
std::shared_ptr<FullMatrix> computeApertureMatrix(const Grid& G,
                                                  const double young,
                                                  const double poisson)
// ----------------------------------------------------------------------------  
{
  const int nc = G.leafGridView().size(0);
  std::shared_ptr<FullMatrix> result(new FullMatrix(nc, nc));
  ddm::assembleMatrix(*result, young, poisson, G);

  *result *= -1;
  return result;
}

// ----------------------------------------------------------------------------
std::vector<HTrans> computeHTrans(const Grid& grid)
// ----------------------------------------------------------------------------
{
  std::vector<HTrans> result;
  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto geom = element.geometry();
    for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
      if (is.boundary())
        continue;
      const int nIdx = mapper.index(is.outside());
      const auto igeom = is.geometry();
      if( eIdx < nIdx) {
        const auto iscenter = igeom.center();
        const auto ncenter = is.outside().geometry().center() - iscenter;
        const auto ecenter = geom.center() - iscenter;
        const double area = igeom.volume();
        const double h1 = ecenter.two_norm() / area;
        const double h2 = ncenter.two_norm() / area;

        result.push_back({nIdx, eIdx, h1, h2});
      }
    }
  }
  return result;
}

  
// ----------------------------------------------------------------------------
std::tuple<std::shared_ptr<SparseMatrix>, std::vector<HTrans>>
computePressureMatrix(const Grid& G)
// ----------------------------------------------------------------------------
{
  // initialize the sparsity pattern of a pressure matrix, with all entries set to zero.

  std::vector<HTrans> hvec = computeHTrans(G);

  const size_t nc = G.leafGridView().size(0);
  std::shared_ptr<SparseMatrix> mat(new SparseMatrix());
  mat->setBuildMode(SparseMatrix::implicit);
  mat->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  mat->setSize(nc, nc);

  for (const auto& elem : hvec) {
    const size_t ix(std::get<0>(elem));
    const size_t jx(std::get<1>(elem));
    mat->entry(ix, jx) = 0.0;
    mat->entry(jx, ix) = 0.0;
    mat->entry(jx, jx) = 0.0;
    mat->entry(ix, ix) = 0.0;
  }
  mat->compress();
  
  return {mat, hvec}; 
}
  
  
// ----------------------------------------------------------------------------
EquationSystem computeEquationSystem(const Grid& G,
                                     const double young, const double poisson)
// ----------------------------------------------------------------------------
{
  return std::tuple_cat(std::make_tuple(computeApertureMatrix(G, young, poisson)),
                        computePressureMatrix(G));
}

// ----------------------------------------------------------------------------
template<class M, class X, class Y>
class ReducedMatrixAdapter : public Dune::LinearOperator<X, Y> 
// ----------------------------------------------------------------------------
{
public:
  ReducedMatrixAdapter(const M& mat, const std::vector<size_t> elim) :
    mat_(mat), elim_(sorted_(elim)), keep_(nonelim_(elim_, mat.M())),
    tmpX_(mat.M()), tmpY_(mat.M()) {}
  virtual ~ReducedMatrixAdapter() {}

  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;

  virtual void apply(const X& x, Y& y) const
  {
    expand(x, tmpX_); // expand x into tmpX_
    mat_.mv(tmpX_, tmpY_);
    contract(tmpY_, y); //contract tmpY_ into y
  }

  virtual void applyscaleadd(field_type alpha, const X& x, Y& y) const
  {
    expand(x, tmpX_); // expand x into tmpX_
    mat_.usmv(alpha, tmpX_, tmpY_);
    contract(tmpY_, y); //contract tmpY_ into y    
  }
  
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

  void expand(const X& source, X& target) const {
    target.resize(mat_.M());
    target = 0;
    for (size_t i = 0; i != source.N(); ++i)
      target[keep_[i]] = source[i];
  }

  void contract(const Y& source, Y& target) const {
    target.resize(keep_.size());
    for (size_t i = 0; i != keep_.size(); ++i)
      target[i] = source[keep_[i]];
  }
  
private:

  std::vector<size_t> sorted_(const std::vector<size_t>& vec)
  {
    std::vector<size_t> result(vec);
    std::sort(result.begin(), result.end());
    return result;
  }

  std::vector<size_t> nonelim_(const std::vector<size_t>& elim, size_t nc)
  {
    // return vector with all incides that are _not_ found in elim
    std::vector<size_t> all_ixs(nc);
    std::iota(all_ixs.begin(), all_ixs.end(), 0);
    std::vector<size_t> nonelim;
    std::set_difference(all_ixs.begin(), all_ixs.end(), elim.begin(), elim.end(),
                        std::back_inserter(nonelim));
    return nonelim;
  }
      
  const M& mat_;
  const std::vector<size_t> elim_;
  const std::vector<size_t> keep_;
  mutable Y tmpY_;
  mutable X tmpX_;
};

// ----------------------------------------------------------------------------  
template<typename T> inline T hmean(const T a, const T b) {return T(1) / ( T(1)/a + T(1)/b);};
// ----------------------------------------------------------------------------
  
// ----------------------------------------------------------------------------
  void updateTrans(SparseMatrix& mat, const std::vector<HTrans>& htransvec, const Vector& h)
// ----------------------------------------------------------------------------
{
  for (const auto& e : htransvec) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e);
    assert(i != j);
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = h[i];
    const double h2 = h[j];

    const double trans1 = h1 * h1 * t1 / 12.0;
    const double trans2 = h2 * h2 * t2 / 12.0;
    
    const double T = hmean(trans1, trans2);

    // update diagonal
    mat[i][i] += T;
    mat[j][j] += T;

    // update off-diagonal terms
    mat[i][j] -= T;
    mat[j][i] -= T;
    
  }
}

// ----------------------------------------------------------------------------
  void solveFixedBHP(Vector& p, Vector& h, const EquationSystem& eqsys,
                     const std::vector<size_t> wellcells,
                     const double bhp)
// ----------------------------------------------------------------------------
{
  // reduce pressure system
  //Dune::MatrixAdapter<FullMatrix, Vector, Vector> MA_h(*std::get<0>(eqsys));
  const FullMatrix&          M_h       = *std::get<0>(eqsys);
  SparseMatrix&              M_p       = *std::get<1>(eqsys);
  const std::vector<HTrans>& htransvec =  std::get<2>(eqsys);

  const RMAdapter MA_p(M_p, wellcells);
  
  // initialize pressure
  p = 0;
  for (auto w : wellcells) p[w] = bhp;

  // solve for aperture
  M_h.solve(h, p); // solve for aperture (h) given pressure (p)

  // update pressure matrix entries
  updateTrans(M_p, htransvec, h);
  Vector rhs_full(p.N());
  M_p.mv(p, rhs_full); // only imposed values of p should be nonzero here.
  Vector rhs;
  MA_p.contract(rhs_full, rhs); // remove entries corresponding to eliminated equations
  
  // solve for pressure

  using PressureOperatorType = Dune::MatrixAdapter<SparseMatrix, Vector, Vector>;
  using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;


  
  
  const auto prm = Opm::setupPropertyTree(Opm::FlowLinearSolverParameters(), true, true);

  auto op = PressureOperatorType(M_p);
  auto psolver = FlexibleSolverType(op, prm, std::function<Vector()>(), size_t(0));



  // auto psolver = Dune::FlexibleSolver<const RMAdapter>(MA_p, prm,
  //                                                      std::function<Vector()>(),
  //                                                      size_t(0));
  Dune::InverseOperatorResult iores;
  psolver.apply(p, rhs, iores);
                                                         
}
  
}; // end anonymous namespace

// ============================================================================
int main(int varnum, char** vararg)
// ============================================================================
{
  // check argument list
  if (varnum < 2) {
    std::cout << "Grid file name must be given as argument." << std::endl;
    return -1;
  }
  
  // read grid
  const auto grid = readgrid(vararg[1]);
  if (!grid)
    return -1;

  // identify well cells
  const auto wellcells = n_closest<2>(*grid);

  std::cout << "Identified wellcells: ";
  std::copy(wellcells.begin(), wellcells.end(), std::ostream_iterator<size_t>(std::cout, " "));
  std::cout << std::endl;
  
  // compute equation system
  const double young = 1e9; // Young's modulus
  const double poisson = 0.25; // Poisson's ratio
  const auto eqsys = computeEquationSystem(*grid, young, poisson);
  
  // solve fixed pressure system
  const size_t nc = grid->leafGridView().size(0);
  Vector pressure(nc), aperture(nc);
  pressure = 0;
  const double bhp = 1e7; // in Pascal
  solveFixedBHP(pressure, aperture, eqsys,
                std::vector<size_t>(wellcells.begin(), wellcells.end()), bhp);

  // solve one injection timestep


  // write output
  Dune::VTKWriter<Grid::LeafGridView> vtkwriter(grid->leafGridView(), Dune::VTK::nonconforming);
  vtkwriter.addCellData(aperture, "aperture");
  vtkwriter.addCellData(pressure, "pressure");
  vtkwriter.write("output");
  
}; // end main
