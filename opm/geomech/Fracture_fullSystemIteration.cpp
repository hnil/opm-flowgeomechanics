#include <tuple>
#include <vector>
#include <functional>
#include <iostream>
#include <fstream>

#include <dune/common/indices.hh> // needed for _0, _1, etc.
          // 
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/multitypeblockvector.hh>

#include "config.h"
#include <opm/geomech/Fracture.hpp>
#include <opm/common/ErrorMacros.hpp>



namespace {
// ========================== Convenience definitions ==========================
using Dune::Indices::_0;
using Dune::Indices::_1;

  //using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using ResVector = Opm::Fracture::Vector; 
using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

  // using Krull = Dune::MultiTypeBlockVector<Dune::BlockVector<double>>;
  // Dune::MultiTypeBlockVector<Dune::BlockVector<double, std::allocator<double>>> dill;
  // Krull tull(2);  
  
using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
using FMatrix = Dune::DynamicMatrix<double>;                       // full matrix

using SystemMatrix =
  Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                             Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;

using Htrans = std::tuple<size_t, size_t, double, double>;

// ============================ Debugging functions ============================
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
void dump_vector(const ResVector& v, const char* const name)
// ----------------------------------------------------------------------------
{
  if (!name)
    for (size_t i = 0; i != v.size(); ++i)
      std::cout << v[i] << "\n";
  if (!name) return;
    
  std::ofstream os(name);
  for (size_t i = 0; i != v.size(); ++i)
    os << v[i] << "\n";
  os.close();
}

// ----------------------------------------------------------------------------
void dump_vector(const VectorHP& v, const char* const name1, const char* const name2)
// ----------------------------------------------------------------------------
{
  dump_vector(v[_0], name1);
  dump_vector(v[_1], name2);
}
  
// ============================= Helper functions =============================

// ----------------------------------------------------------------------------  
SMatrix makeIdentity(size_t num, double right_padding = 0, double fac = 1)
// ----------------------------------------------------------------------------
{
  //build a sparse matrix to represent the identity matrix of a given size
  SMatrix M;
  M.setBuildMode(SMatrix::implicit);
  M.setImplicitBuildModeParameters(1, 0.1);
  M.setSize(num, num + right_padding);
  for (size_t i = 0; i != num; ++i)
    M.entry(i,i) = fac;
  M.compress();
  return M;
}

// ----------------------------------------------------------------------------
std::unique_ptr<Opm::Fracture::Matrix>
initCouplingMatrixSparsity(const std::vector<Htrans>& htrans, size_t num_cells, size_t num_wells)
// ----------------------------------------------------------------------------
{
  auto C = std::make_unique<Opm::Fracture::Matrix>(num_cells + num_wells,
                                                   num_cells,
                                                   4, 0.4,
                                                   Opm::Fracture::Matrix::implicit);
  for (const auto& e : htrans) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e); assert(i != j);
    C->entry(i, j) = 0;
    C->entry(j, i) = 0;
    C->entry(i,i) = 0;
    C->entry(j,j) = 0;
  }
  C->compress();
  return C;
}
  
// ----------------------------------------------------------------------------
void updateCouplingMatrix(std::unique_ptr<Opm::Fracture::Matrix>& Cptr,
                          const size_t num_cells,
                          const size_t num_wells,
                          const std::vector<Htrans>& htrans,
                          const ResVector& pressure,
                          const ResVector& aperture)
{
  // create C if not done already
  if (!Cptr) Cptr = initCouplingMatrixSparsity(htrans, num_cells, num_wells);
  // if (!Cptr) Cptr = std::make_unique<Opm::Fracture::Matrix>(*M); // copy sparsity of M

  // @@ NB: If the implementation of the pressure matrix changes (i.e. `assemblePressure()`),
  //        the code below might need to be updated accordingly as well.
  auto& C = *Cptr;
  C = 0;
  for (const auto& e : htrans) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e); assert(i != j);
    
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = aperture[i];
    const double h2 = aperture[j];
    const double p1 = pressure[i];
    const double p2 = pressure[j];
    
    const double q   = (h1 * h1 * h1) * (h2 * h2 * h2) * (t1 * t2); // numerator
    const double d1q = 3 * (h1 * h1) * (h2 * h2 * h2) * (t1 * t2);
    const double d2q = 3 * (h1 * h1 * h1) * (h2 * h2) * (t1 * t2);
    
    const double r   = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
    const double d1r = 36 * (h1 * h1) * t1;
    const double d2r = 36 * (h2 * h2) * t2;
    
    const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
    const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);

    const double krull = 1; // 1e4; // @@ Not sure if this should be removed?
    // diagonal elements
    C[i][i] += dTdh1 * (p1-p2) *krull;
    C[j][j] += dTdh2 * (p2-p1) *krull;
    
    // off-diagonal elements
    C[i][j] += dTdh2 * (p1-p2) *krull;
    C[j][i] += dTdh1 * (p2-p1) *krull; 
    
  }
}

// ----------------------------------------------------------------------------  
inline bool convergence_test(const VectorHP& res, const double tol_flow, double tol_mech)
// ----------------------------------------------------------------------------
{
    std::cout << "Residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
    std::cout << "Residual norm[1] is " << res[_1].infinity_norm()<<std::endl;

    std::cout << "tol mech is: " << tol_mech << std::endl;
    std::cout << "tol flow is: " << tol_flow << std::endl;
  return res[_0].infinity_norm() < tol_mech && res[_1].infinity_norm() < tol_flow;
  //return res.infinity_norm() < tol;
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


  
}; // end anonymous namespace

namespace Opm
{

// ----------------------------------------------------------------------------
bool Fracture::fullSystemIteration(const double tol)
// ----------------------------------------------------------------------------
{
  // update pressure matrix with the current values of `fracture_width_` and
  // `fracture_pressure_` 
  assemblePressure(); // update pressure matrix
  setSource();         // update right-hand side of pressure system;

  // update the coupling matrix (possibly create it if not already initialized)
  updateCouplingMatrix(coupling_matrix_,
                       pressure_matrix_->N() - numWellEquations(), // num cells
                       numWellEquations(),                         // num wells
                       htrans_, fracture_pressure_, fracture_width_);

  const auto& A = fractureMatrix(); // will be created if not already existing
  const auto& M = *pressure_matrix_;
  const auto& C = *coupling_matrix_;
  const auto I = makeIdentity(A.N(), numWellEquations());

  // system Jacobian (with cross term)  @@ should S be included as a member variable of Fracture?
  SystemMatrix S { {A, I},   // mechanics system (since A is negative, we leave I positive here)
                   {C, M} }; // flow system

  // system equations
  SystemMatrix S0 = S; S0[_1][_0] = 0; // the equations themselves have no cross term

  VectorHP x {fracture_width_, fracture_pressure_};
  VectorHP dx = x; dx = 0; // gradient of 'x' (which we aim to compute below)

  // initialize right-hand aide
  VectorHP rhs {x}; rhs = 0; // same size as x, but initially just with zeros

  //rhs[_0] = 0; // @@ we currently do not use rhs_width_ here, but include p through the I block
  normalFractureTraction(rhs[_0], false); // right-hand side equals the normal fracture traction
  rhs[_1] = rhs_pressure_; // should have been updated in call to `assemblePressure` above

  std::cout << "rhs[0]: " << rhs[_0].infinity_norm() << std::endl;
  std::cout << "pressure: ";
  std::copy(fracture_pressure_.begin(), fracture_pressure_.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  
  S0.mmv(x, rhs); // rhs = rhs - S0 * x;   (we are working in the tanget plane)

  // check if system is already at a converged state (in which case we return immediately)
  //
  // for convergence, we use 'tol' directly for the mechanics system (where
  // residuals are expected to scale with pressure), and 'tol * ||M||' for the
  // flow system (where residuals scale with M*p)
  if (convergence_test(rhs, tol * M.infinity_norm(), tol)) 
    return true;
                         
  // solve system equations
  const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S); 
  TailoredPrecondDiag precond(S); // cannot be 'const' due to BiCGstabsolver interface
  Dune::InverseOperatorResult iores; // cannot be 'const' due to BiCGstabsolver interface
  auto psolver = Dune::BiCGSTABSolver<VectorHP>(S_linop,
                                                precond,
                                                1e-10, //1e-20, // desired rhs reduction factor
                                                200, // max number of iterations
                                                1); // verbose
  psolver.apply(dx, rhs, iores); // NB: will modify 'rhs'

  std::cout << "x:  " << x[_0].infinity_norm() << " " << x[_1].infinity_norm() << std::endl;
  std::cout << "dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm() << std::endl;
  
  // the following is a heuristic way to limit stepsize to stay within convergence radius
  const double step_fac = estimate_step_fac(x, dx);
  std::cout << "fac: " << step_fac << std::endl;

  dx *= step_fac;
  x += dx;

  // copying modified variables back to member variables
  fracture_width_ = x[_0];
  fracture_pressure_ = x[_1];
  
  return false;
}
    
}; // end namespace Opm
