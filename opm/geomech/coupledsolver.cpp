#include "coupledsolver.hpp"
#include <opm/common/ErrorMacros.hpp>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <functional>

namespace { // anonymous namespace

using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;
using SystemMatrix =
  Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                             Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
using Dune::Indices::_0;
using Dune::Indices::_1;

// ----------------------------------------------------------------------------
void update_flow_system(FMatrix& M, FMatrix& C,
                        const VectorHP& hp,
                        const double rate,
                        const std::vector<size_t> ratecells,
                        const double bhp,
                        const std::vector<size_t> bhpcells)
// ----------------------------------------------------------------------------  
{
  OPM_THROW(std::runtime_error,"update_flow_system unimplemented.");
  // @ unimplemented for now
}

// ----------------------------------------------------------------------------  
SMatrix makeIdentity(size_t num, double fac = 1)
// ----------------------------------------------------------------------------
{
  //build a sparse matrix to represent the identity matrix of a given size
  SMatrix M;
  M.setBuildMode(SMatrix::implicit);
  M.setImplicitBuildModeParameters(1, 0.1);
  M.setSize(num, num);
  for (size_t i = 0; i != num; ++i)
    M.entry(i,i) = fac;
  M.compress();
  return M;
}
  
// ----------------------------------------------------------------------------
bool nonlinear_iteration(const SMatrix& A, // mech. matrix
                         FMatrix& M, // flow matrix, will be updated (depends on current aperture)
                         FMatrix& C, // coupling matrix, will be updated
                         const VectorHP& hp,
                         VectorHP& dhp, // increment (to be computed)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<size_t> ratecells,
                         const double bhp,
                         const std::vector<size_t> bhpcells,
                         const std::function<bool(const VectorHP&)>& convergence_test)
// ----------------------------------------------------------------------------
{
  update_flow_system(M, C, hp, rate, ratecells, bhp, bhpcells);

  const SystemMatrix
  
  const auto negI = makeIdentity(A.N(), -1);
  
  
  
  return true; // @@ still unimplemented
}
   
}; // end anonymous namespace

namespace Opm
{

// ----------------------------------------------------------------------------
void solve_fully_coupled(ResVector& pvec, // output: fracture pressure
                         ResVector& hvec, // output: aperture
                         const SMatrix& pmat, // pressure matrix (sparse)
                         const FMatrix& amat, // aperture matrix (full)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<size_t> ratecells,
                         const double bhp,
                         const std::vector<size_t> bhpcells,
                         const int max_nonlin_iter,
                         const double conv_tol)
// ----------------------------------------------------------------------------  
{
  // We consider the nonlinear system:
  // | A          -I |  | h |   | 0 |   (mechanics equation)
  // |               |  |   | = |   |
  // | C(h, p)  M(h) |  | p |   | q |   ( flow equation)

  // note below that C and M have same structure as 'amat'
  SystemMatrix S { { pmat, makeIdentity(pmat.N(), -1) }, // | A  -I|
                   { amat, amat } };                     // | C   M|  
  
  VectorHP dhp; // update
  VectorHP hp {hvec, pvec}; // aperture h and pressure p

  const double h_tol = (bhpcells.size() > 0) ? conv_tol * bhp : conv_tol;
  const double p_tol = conv_tol;

  std::function<bool(const VectorHP&)> conv_test = [&](const VectorHP& res) {
    return res[_0].infinity_norm() < h_tol && res[_1].infinity_norm() < p_tol;
  };
  
  int iter = 0;
  while ( ! nonlinear_iteration(S, hp, dhp, htrans, rate,
                                ratecells, bhp, bhpcells, conv_test) &&
          iter++ < max_nonlin_iter)
    hp += dhp;

  if (iter == max_nonlin_iter)
    std::cout << "Warning: failed to converge." << std::endl; // @@
  
  hvec = hp[_0];
  pvec = hp[_1];
}
  
}; // end namespace Opm
