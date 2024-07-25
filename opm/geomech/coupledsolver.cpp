#include <functional>
#include <iostream>
#include <fstream>

#include "coupledsolver.hpp"
#include <opm/common/ErrorMacros.hpp>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

namespace { // anonymous namespace

using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;
using SystemMatrix =
  Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                             Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
using Dune::Indices::_0;
using Dune::Indices::_1;

// ----------------------------------------------------------------------------
//template<class MatrixType> void dump_matrix(const MatrixType& m, const char* const name)
void dump_matrix(const FMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
  const size_t rows = m.N();
  const size_t cols = m.M();
  std::ofstream os(name);
  for (size_t i = 0; i != rows; ++i)
    for (size_t j = 0; j != cols; ++j)
      os << m[i][j] << ((j == cols-1) ? "\n" : " ");
  os.close();
}

// ----------------------------------------------------------------------------
//template<>
void dump_matrix(const SMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
  std::ofstream os(name);
  for (auto rowIt = m.begin(); rowIt != m.end(); ++ rowIt) {
    const size_t i = rowIt.index();
    for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
      const size_t j = colIt.index();
      os << i+1 << " " << j+1 << " " << m[i][j] << "\n";
    }
  }
  os.close();
}

// ----------------------------------------------------------------------------
ResVector v0(const VectorHP& v)
// ----------------------------------------------------------------------------
{
  return v[_0];
}

// ----------------------------------------------------------------------------
ResVector v1(const VectorHP& v)
// ----------------------------------------------------------------------------
{
  return v[_1];
}
  
// ----------------------------------------------------------------------------
void dump_vector(const ResVector& v, const char* const name)
// ----------------------------------------------------------------------------
{
  std::ofstream os(name);
  for (size_t i = 0; i != v.size(); ++i)
    os << v[i] << "\n";
  os.close();
}
  
// ============================================================================
  
// ----------------------------------------------------------------------------  
template<typename T> inline T hmean(const T a, const T b) 
// ----------------------------------------------------------------------------
{
  T tmp = T(1) / ( T(1)/a + T(1)/b);
  return std::isnan(tmp) ? 0.0 : tmp;
}
  
// ----------------------------------------------------------------------------
void update_M_and_rhs(SMatrix& M, ResVector& rhs,
                      const std::vector<Htrans>& htrans,
                      const VectorHP& hp,
                      const double rate,
                      const std::vector<size_t>& ratecells,
                      const double bhp,
                      const std::vector<double>& leakoff_fac,
                      const std::vector<int> remapping)
// ----------------------------------------------------------------------------  
{
  // update M and rhs_reduct_
  M = 0; rhs = 0;
  for (const auto& e : htrans) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e);
    assert(i != j);
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = hp[_0][i];
    const double h2 = hp[_0][j];
    
    const double trans1 = h1 * h1 * t1 / 12.0;
    const double trans2 = h2 * h2 * t2 / 12.0;
    
    const double T = hmean(trans1, trans2);
    
    // remapped indices
    const int ir = remapping[i];
    const int jr = remapping[j];
    
    // update diagonal
    if (ir >= 0) M[ir][ir] += T; 
    if (jr >= 0) M[jr][jr] += T; 
    
    // update off-diagonal terms, and the reduction terms of rhs_reduct_
    if (ir >= 0) {if (jr >= 0) M[ir][jr] -= T; else rhs[ir] += T * bhp; }
    if (jr >= 0) {if (ir >= 0) M[jr][ir] -= T; else rhs[jr] += T * bhp; }
  }
  
  // add leakage terms
  for (size_t ix = 0; ix != leakoff_fac.size(); ++ix)
    if (remapping[ix] >= 0)
      M[remapping[ix]][remapping[ix]] += leakoff_fac[ix];
    
  // add source terms
  for (const auto& ix : ratecells)
    if (remapping[ix] >= 0)
      rhs[remapping[ix]] = rate;
}

// ----------------------------------------------------------------------------
void update_C(SMatrix& C, const std::vector<Htrans>& htrans, const VectorHP& hp,
              const std::vector<int> remapping)
// ----------------------------------------------------------------------------  
{
  C = 0;
  // intermediary matrices for computing the partial derivatives
  for (const auto& e : htrans) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e); assert(i != j);
    
    const int ir = remapping[i];
    const int jr = remapping[j];
    
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = hp[_0][i];
    const double h2 = hp[_0][j];
    const double p1 = hp[_1][i];
    const double p2 = hp[_1][j];
    
    const double q   = h1 * h1 * h2 * h2 * t1 * t2; // numerator
    const double d1q = 2 * h1 * h2 * h2 * t1 * t2;
    const double d2q = h1 * h1 * 2 * h2 * t1 * t2;
    
    const double r   = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
    const double d1r = 24 * h1 * t1;
    const double d2r = 24 * h2 * t2;
    
    const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
    const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);
    
    // diagonal elements
    if (ir >= 0) C[ir][i] += dTdh1 * (p1-p2);
    if (jr >= 0) C[jr][j] += dTdh2 * (p2-p1); 
    
    // off-diagonal elements
    if (ir >= 0) C[ir][j] += dTdh2 * (p1-p2);
    if (jr >= 0) C[jr][i] += dTdh1 * (p2-p1); 
    
  }
}

// ----------------------------------------------------------------------------
const std::vector<int> compute_remapping(const std::vector<size_t>& eliminated,
                                         const std::vector<size_t>& kept,
                                         const size_t num_cells)
// ----------------------------------------------------------------------------  
{
  std::vector<int> result(num_cells);
  for (size_t i = 0; i != eliminated.size(); ++i) result[eliminated[i]] = -i-1;
  for (size_t i = 0; i != kept.size(); ++i) result[kept[i]] = i;
  return result;
}
  
// ----------------------------------------------------------------------------
void update_flow_system(SMatrix& M, SMatrix& C, // flow and coupling matrices
                        ResVector& rhs, // right-hand side (sources/sinks)
                        const std::vector<Htrans>& htrans, 
                        const VectorHP& hp,
                        const double rate,
                        const std::vector<size_t>& ratecells,
                        const double bhp,
                        const std::vector<size_t>& eliminated_cells, 
                        const std::vector<size_t>& kept_cells,
                        const std::vector<double>& leakfac)
// ----------------------------------------------------------------------------  
{
  // compute remapping, in case there have been eliminated pressure values
  const std::vector<int> remap = compute_remapping(eliminated_cells, kept_cells, M.M());

  // update the output arguments: the matrices 'M' and 'C', and the vector 'rhs'
  update_M_and_rhs(M, rhs, htrans, hp, rate, ratecells, bhp, leakfac, remap);
  update_C(C, htrans, hp, remap);
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
template<typename Mat> ResVector diagvec(const Mat& M)
// ----------------------------------------------------------------------------
{
  ResVector res(M.M());
  for (size_t i = 0; i != res.size(); ++i)
    res[i] = M[i][i];
  return res;
}
  
// ----------------------------------------------------------------------------
class TailoredPrecondDiag : public Dune::Preconditioner<VectorHP, VectorHP>
// ----------------------------------------------------------------------------
{
public:
  TailoredPrecondDiag(const SystemMatrix& S) :
    A_diag_(diagvec(S[_0][_0])), M_diag_(diagvec(S[_1][_1])) {}
  virtual void apply(VectorHP& v, const VectorHP& d) {
    for (size_t i = 0; i != A_diag_.size(); ++i) v[_0][i] = d[_0][i]/A_diag_[i];
    for (size_t i = 0; i != M_diag_.size(); ++i) v[_1][i] = d[_1][i]/M_diag_[i];
  };
  virtual void post(VectorHP& v) {};
  virtual void pre(VectorHP& x, VectorHP& b) {};
  virtual Dune::SolverCategory::Category category() const {return Dune::SolverCategory::sequential;}
private:
  const ResVector A_diag_;
  const ResVector M_diag_;
};
  
// ----------------------------------------------------------------------------
double estimate_step_fac(const VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------  
{
  // estimate what might be a safe step size to avoid exiting the convergence radius
  const double f1 = dx[_0].infinity_norm() / x[_0].infinity_norm();
  const double f2 = dx[_1].infinity_norm() / x[_1].infinity_norm();
  const double fmax = std::max(f1, f2);
  const double threshold = 0.95;
  const double fac_min = 1e-2; 
  return (fmax < threshold) ? 1.0 : std::max(threshold / fmax, fac_min);
}
  
// ----------------------------------------------------------------------------
bool nonlinear_iteration(SystemMatrix& S, // will be updated (depends on current aperture)
                         const SMatrix& negI_rhs, // needed to setup rhs, in case of eliminated vars
                         const VectorHP& hp,
                         VectorHP& dhp, // increment (to be computed)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<size_t>& ratecells,
                         const double bhp,
                         const std::vector<size_t>& eliminated_cells, // sorted bhp_cells
                         const std::vector<size_t>& kept_cells, 
                         const std::vector<double>& leakoff_fac,
                         const std::function<bool(const VectorHP&)>& convergence_test)
// ----------------------------------------------------------------------------
{
  VectorHP rhs(hp); // same size as hp, but only zeros
  rhs = 0;
  
  auto& M = S[_1][_1]; // flow matrix (reference; will be modified)
  auto& C = S[_1][_0]; // coupling matrix (reference; will be modified)

  // --- update matrix entries based on the latest iteration of the solution vector ---
  update_flow_system(M, C, rhs[_1], htrans, hp, rate, ratecells,
                     bhp, eliminated_cells, kept_cells, leakoff_fac); 

  // --- update right-hand side ---
  if (negI_rhs.M() > 0) {
    // move impact of eliminated pressure values over to rhs
    ResVector tmp(negI_rhs.M()); tmp = -bhp;
    negI_rhs.mv(tmp, rhs[_0]);
  }
  auto S0 = S;
  S0[_1][_0] = 0; // the C matrix should be zero here
  S0.mmv(hp, rhs); // rhs = rhs - S0 * hp_reduced

  if (convergence_test(rhs))
    return true; // system converged
  
  // --- solve system equations --- 
  const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S); 
  TailoredPrecondDiag precond(S); // cannot be 'const' due to BiCGstabsolver interface
  Dune::InverseOperatorResult iores; // cannot be 'const' due to BiCGstabsolver interface
  auto psolver = Dune::BiCGSTABSolver<VectorHP>(S_linop,
                                                precond,
                                                1e-20, // desired rhs reduction factor
                                                200, // max number of iterations
                                                1); // verbose
  psolver.apply(dhp, rhs, iores); // NB: will modify 'rhs'

  std::cout << "hp:  " << hp[_0].infinity_norm() << " " << hp[_1].infinity_norm() << std::endl;
  std::cout << "dhp: " << dhp[_0].infinity_norm() << " " << dhp[_1].infinity_norm() << std::endl;
  
  // the following is a heuristic way to limit stepsize to stay within convergence radius
  const double step_fac = estimate_step_fac(hp, dhp);
  std::cout << "fac: " << step_fac << std::endl;
  dhp *= step_fac;

  return false;
}

// ----------------------------------------------------------------------------
template <typename T>
std::vector<T> range_vec(T startval, T num_elem)
// ----------------------------------------------------------------------------  
{
  std::vector<T> result(num_elem);
  std::iota(result.begin(), result.end(), startval);
  return result;
}

// ----------------------------------------------------------------------------
  std::vector<size_t> remaining(const size_t N, const std::vector<size_t>& elim)
// ----------------------------------------------------------------------------  
{
  std::vector<size_t> all = range_vec(size_t(0), N);
  std::vector<size_t> res;
  std::set_difference(all.begin(), all.end(), elim.begin(), elim.end(),
                      std::back_inserter(res));
  return res;
}
  
// ----------------------------------------------------------------------------
// Note: 'keep' should be sorted. The matrix 'M' will be reduced, and any
// columns extracted will be collected and returned. (Rows will just be
// removed).
SMatrix reduceMatrix(SMatrix& M, bool rows, bool cols, const std::vector<size_t> elim)
// ----------------------------------------------------------------------------
{
  assert(cols || rows);
  const size_t total_size = cols ? M.M() : M.N();
  const std::vector<size_t> keep = remaining(total_size, elim);

  SMatrix Mreduced, reduced_cols;

  Mreduced.setBuildMode(SMatrix::implicit);
  Mreduced.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  Mreduced.setSize(rows ? keep.size() : M.N(),
                   cols ? keep.size() : M.M());

  reduced_cols.setBuildMode(SMatrix::implicit);
  reduced_cols.setImplicitBuildModeParameters(3*6-3-2, 0.4);
  reduced_cols.setSize(Mreduced.N(), 
                       cols ? elim.size() : 0);

  std::vector<int> mapping(keep.size() + elim.size(), 0);
  for (size_t i = 0; i != keep.size(); ++i)  mapping[keep[i]] = i;
  for (size_t i = 0; i != elim.size(); ++i)  mapping[elim[i]] = -(i+1);
  
  for (auto rowIt = M.begin(); rowIt != M.end(); ++ rowIt) {
    const size_t i = rowIt.index();
    if (rows && mapping[i] < 0)
      continue; // this row should be eliminated, nothing further to do
    
    for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
      const size_t j = colIt.index();
      
      if (cols && mapping[j] < 0)
        // this column should be eliminated
        reduced_cols.entry(rows ? mapping[i] : i, -mapping[i] - 1) = M[i][j];
      else 
        // this column should be kept
        Mreduced.entry(rows ? mapping[i] : i,
                       cols ? mapping[j] : j) = M[i][j];
    }
  }
  
  Mreduced.compress();
  reduced_cols.compress();
  std::swap(M, Mreduced);
  return reduced_cols;
}

// ----------------------------------------------------------------------------  
ResVector reduceVector(const ResVector& v, const std::vector<size_t>& keep)
// ----------------------------------------------------------------------------  
{
  ResVector res(keep.size());
  size_t pos = 0;
  for (auto it : keep)
    res[pos++] =v[it];
  return res;
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
                         const std::vector<size_t>& ratecells,
                         const double bhp,
                         const std::vector<size_t>& bhpcells,
                         const std::vector<double>& leakfac,  
                         const int max_nonlin_iter,
                         const double conv_tol)
// ----------------------------------------------------------------------------  
{
  // We consider the nonlinear system:
  // | A          -I |  | h |   | 0 |   (mechanics equation)
  // |               |  |   | = |   |
  // | C(h, p)  M(h) |  | p |   | q |   ( flow equation)

  // note below that C and M have same structure as 'amat'
  auto eliminated(bhpcells); std::sort(eliminated.begin(), eliminated.end());
  auto kept = remaining(pmat.N(), eliminated);
  auto C(pmat), M(pmat);
  auto negI = makeIdentity(amat.N(), -1);
  reduceMatrix(C, true, false, eliminated);
  reduceMatrix(M, true, true, eliminated);
  const auto negI_rhs = reduceMatrix(negI, false, true, eliminated);
  
  SystemMatrix S { { amat, negI }, // | A  -I|
                   { C,    M } };  // | C   M|  
  
  VectorHP hp {hvec, reduceVector(pvec, kept)}; // aperture h and pressure p
  VectorHP dhp = hp; dhp = 0; // update vector has same shape as hp

  const double h_tol = (bhpcells.size() > 0) ? conv_tol * bhp : conv_tol;
  const double p_tol = conv_tol;

  std::function<bool(const VectorHP&)> conv_test = [&](const VectorHP& res) {
    return res[_0].infinity_norm() < h_tol && res[_1].infinity_norm() < p_tol;
  };
  
  int iter = 0;
  while ( ! nonlinear_iteration(S, negI_rhs, hp, dhp, htrans, rate,
                                ratecells, bhp, eliminated, kept, leakfac, conv_test) &&
          iter++ < max_nonlin_iter)
    hp += dhp;

  if (iter == max_nonlin_iter)
    std::cout << "Warning: failed to converge." << std::endl; // @@
  
  hvec = hp[_0];
  for (size_t i = 0; i != kept.size(); ++i)
    pvec[kept[i]] = hp[_1][i];
}
  
}; // end namespace Opm
