#include <dune/common/exceptions.hh>
#include <dune/istl/solvercategory.hh>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
//#include <dune/istl/paamg/amg.hh>

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

using FullMatrix = Dune::DynamicMatrix<double>;
using SparseMatrix = Dune::BCRSMatrix<double>;
using EquationSystem = std::tuple<std::shared_ptr<FullMatrix>,    // aperture matrix
                                  std::shared_ptr<SparseMatrix>,  // pressure matrix
                                  std::vector<HTrans>,            // inter-cell transmissibility factors
                                  std::vector<double>>;           // leakoff factors

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
  std::tuple<std::vector<HTrans>, std::vector<double>>
  computeHTrans(const Grid& grid, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  std::vector<HTrans> htransvec;
  std::vector<double> leakvec(grid.leafGridView().size(0), 0.0);

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

        htransvec.push_back({nIdx, eIdx, h1, h2});
      }
    }
    leakvec[eIdx] = geom.volume() * leakoff_fac;
  }
  return {htransvec, leakvec};
}


// ----------------------------------------------------------------------------
std::tuple<std::shared_ptr<SparseMatrix>, std::vector<HTrans>, std::vector<double>>
computePressureMatrix(const Grid& G, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  // initialize the sparsity pattern of a pressure matrix, with all entries set to zero.

  std::tuple<std::vector<HTrans>,
             std::vector<double>> hvec = computeHTrans(G, leakoff_fac);

  const size_t nc = G.leafGridView().size(0);
  std::shared_ptr<SparseMatrix> mat(new SparseMatrix());
  mat->setBuildMode(SparseMatrix::implicit);
  mat->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  mat->setSize(nc, nc);

  for (const auto& elem : std::get<0>(hvec)) {
    const size_t ix(std::get<0>(elem));
    const size_t jx(std::get<1>(elem));
    mat->entry(ix, jx) = 0.0;
    mat->entry(jx, ix) = 0.0;
    mat->entry(jx, jx) = 0.0;
    mat->entry(ix, ix) = 0.0;
  }
  mat->compress();

  return std::tuple_cat(std::make_tuple(mat), hvec);
}


// ----------------------------------------------------------------------------
EquationSystem computeEquationSystem(const Grid& G,
                                     const double young, const double poisson,
                                     const double leakoff_fac)
// ----------------------------------------------------------------------------
{
  return std::tuple_cat(std::make_tuple(computeApertureMatrix(G, young, poisson)),
                        computePressureMatrix(G, leakoff_fac));
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
    expand(y, tmpY_); // expand y into tmpY_
    mat_.usmv(alpha, tmpX_, tmpY_);
    contract(tmpY_, y); //contract tmpY_ into y
  }

  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::sequential;
  }

  void expand(const X& source, X& target) const {
    target.resize(fullSize());
    target = 0;
    for (size_t i = 0; i != reducedSize(); ++i)
      target[keep_[i]] = source[i];
  }

  void contract(const Y& source, Y& target) const {
    target.resize(reducedSize());
    for (size_t i = 0; i != reducedSize(); ++i)
      target[i] = source[keep_[i]];
  }

  size_t fullSize() const {return mat_.M();}
  size_t reducedSize() const {return keep_.size();}

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
  mutable X tmpX_;
  mutable Y tmpY_;
};

// ----------------------------------------------------------------------------
template<typename T> inline T hmean(const T a, const T b) {return T(1) / ( T(1)/a + T(1)/b);};
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void updateTrans(SparseMatrix& mat, const std::vector<HTrans>& htransvec,
                 const Vector& h, const std::vector<double>& leakvec)
// ----------------------------------------------------------------------------
{
  mat = 0; // reset matrix before starting to fill in values

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

  // add-in leakage term on diagonal
  for (size_t ix = 0; ix != mat.N(); ++ix)
    mat[ix][ix] += leakvec[ix];

}

// ----------------------------------------------------------------------------
void solveCoupled(Vector& p, Vector& h, const EquationSystem& eqsys,
                  const std::vector<size_t> bhpcells, const double bhp,
                  const std::vector<size_t> ratecells, const double rate,
                  const double convergence_tol=1e-4, const int max_nonlin_iter=400)
// ----------------------------------------------------------------------------
{
  // NB: there should be no overlap between the cells in 'bhpcells' and 'ratecells'
  // reduce pressure system

  //Dune::MatrixAdapter<FullMatrix, Vector, Vector> MA_h(*std::get<0>(eqsys));
  const FullMatrix&          M_h       = *std::get<0>(eqsys);
  SparseMatrix&              M_p       = *std::get<1>(eqsys);
  const std::vector<HTrans>& htransvec =  std::get<2>(eqsys);
  const std::vector<double>& leakvec   =  std::get<3>(eqsys);

  RMAdapter MA_p(M_p, bhpcells); // allows reducing the system without creating new matrices

  // define convergence criterion
  auto converged = [&](const Vector& v1, const Vector& v2) {
    auto diff = v1; diff -= v2; // diff = v1 - v2
    const double delta = diff.infinity_norm();
    const double denom = std::max(v1.infinity_norm(), v2.infinity_norm());
    std::cout << "Delta: " << delta << " Denom: " << denom;
    std::cout << " Delta/denom: " << delta/denom << std::endl;
    return delta/denom < convergence_tol;
  };

  // helper function to set vector values
  auto setvals = [](Vector& v, const std::vector<size_t>& ixs, const double val) {
    for (auto i : ixs) v[i] = val;
  };

  // initialize pressure
  p = 0;
  if (bhpcells.size() > 0)
    setvals(p, bhpcells, bhp);
  else
    p = 1e6; // we need something to get started without h being 0.

  // solve for aperture
  M_h.solve(h, p); // solve for aperture (h) given pressure (p)

  Vector htmp(h.N()); htmp = 0;
  Vector rhs_full_rate(p.N()), rhs_full(p.N());
  Vector rhs;
  Dune::InverseOperatorResult iores;
  rhs_full_rate = 0; setvals(rhs_full_rate, ratecells, rate);

  int i;
  for (i = 0; i != max_nonlin_iter; ++i) {

    // update pressure matrix entries
    updateTrans(M_p, htransvec, h, leakvec);

    // determine right-hand side modifications from eliminated degrees of freedom

    M_p.mv(p, rhs_full); // only imposed values of p should be nonzero here.
    rhs_full *= -1;
    rhs_full += rhs_full_rate;

    MA_p.contract(rhs_full, rhs); // remove entries corresponding to eliminated equations

    // solve for pressure
    Dune::Richardson<Vector, Vector> precond(1); // "no" preconditioner
    //Dune::SeqJac<SparseMatrix, Vector, Vector> precond(M_p, 3, 0.3);
    // using Smoother = Dune::SeqSSOR<SparseMatrix, Vector, Vector>;
    // Dune::Amg::SmootherTraits<Smoother>::Arguments smootherArgs;
    // smootherArgs.iterations = 3;
    // smootherArgs.relaxationFactor = 1;

    // auto criterion =
    //   Dune::Amg::CoarsenCriterion<Dune::Amg::UnSymmetricCriterion<SparseMatrix,
    //                                                               Dune::Amg::FirstDiagonal>>(15, 50);
    // criterion.setDefaultValuesIsotropic(2);
    // criterion.setAlpha(.67);
    // criterion.setBeta(1.0e-4);
    // criterion.setGamma(1);
    // criterion.setDebugLevel(2);
    // using AMG = Dune::Amg::AMG<RMAdapter, Vector, Smoother>;
    // AMG precond(MA_p, criterion, smootherArgs);

    auto psolver = Dune::CGSolver<Vector>(MA_p,
                                          precond,
                                          1e-7, // desired residual reduction factor
                                          100, // max number of iterations
                                          1); // verbose

    Vector p_reduced(MA_p.reducedSize());
    psolver.apply(p_reduced, rhs, iores);
    MA_p.expand(p_reduced, p);
    setvals(p, bhpcells, bhp);

  // solve for aperture again
    M_h.solve(htmp, p); // solve for aperture (h) given pressure (p)

    if (converged(h, htmp))
      break;

    h = htmp;

  }

  if (i == max_nonlin_iter)
    std::cout << "Warning, did not converge in max number of nonlinear iterations." << std::endl;
  else
    std::cout << "Converged in: " << i << " iterations." << std::endl;
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
  const double leakoff_fac = 1e-6; //1e-13; // a bit heuristic; conceptually rock perm divided by distance
  const auto eqsys = computeEquationSystem(*grid, young, poisson, leakoff_fac);

  // solve fixed pressure system
  const size_t nc = grid->leafGridView().size(0);
  Vector pressure(nc), aperture(nc);
  pressure = 0;
  const double bhp = 1e7; // in Pascal
  const double rate = -0.1;

  // solve fixed pressure system
  // solveCoupled(pressure, aperture, eqsys,
  //              std::vector<size_t>(wellcells.begin(), wellcells.end()), bhp,
  //              std::vector<size_t>(), rate);

  // solve fixed rate system
  solveCoupled(pressure, aperture, eqsys,
               std::vector<size_t>(), bhp,
               std::vector<size_t>(wellcells.begin(), wellcells.end()), rate);


  // write output
  Dune::VTKWriter<Grid::LeafGridView> vtkwriter(grid->leafGridView(), Dune::VTK::nonconforming);
  vtkwriter.addCellData(aperture, "aperture");
  vtkwriter.addCellData(pressure, "pressure");
  vtkwriter.write("output");

}; // end main
