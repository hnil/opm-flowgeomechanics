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
//#include <opm/models/blackoil/blackoilmodel.hh>
//#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/common/ErrorMacros.hpp>



// namespace {
// // ========================== Convenience definitions ==========================
// using Dune::Indices::_0;
// using Dune::Indices::_1;

// using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
//   //using ResVector = Opm::Fracture::Vector; @@
// using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

//   using Krull = Dune::MultiTypeBlockVector<Dune::BlockVector<double>>;
//   Dune::MultiTypeBlockVector<Dune::BlockVector<double, std::allocator<double>>> dill;
//   Krull tull(2);  
  
// using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
// using FMatrix = Dune::DynamicMatrix<double>;                       // full matrix

// using SystemMatrix =
//   Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
//                              Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;

//   SystemMatrix krull;
// using Htrans = std::tuple<size_t, size_t, double, double>;
  
// // ============================= Helper functions =============================

// // ----------------------------------------------------------------------------  
// SMatrix makeIdentity(size_t num, double fac = 1)
// // ----------------------------------------------------------------------------
// {
//   //build a sparse matrix to represent the identity matrix of a given size
//   SMatrix M;
//   M.setBuildMode(SMatrix::implicit);
//   M.setImplicitBuildModeParameters(1, 0.1);
//   M.setSize(num, num);
//   for (size_t i = 0; i != num; ++i)
//     M.entry(i,i) = fac;
//   M.compress();
//   return M;
// }

// // ----------------------------------------------------------------------------
// void updateCouplingMatrix(std::unique_ptr<Opm::Fracture::Matrix>& Cptr,
//                           const std::unique_ptr<Opm::Fracture::Matrix>& M, 
//                           const std::vector<Htrans>& htrans,
//                           const ResVector& pressure,
//                           const ResVector& aperture)
// {
//   // create C if not done already
//   if (!Cptr) Cptr = std::make_unique<Opm::Fracture::Matrix>(*M); // copy sparsity of M

//   // @@ NB: If the implementation of the pressure matrix changes (i.e. `assemblePressure()`),
//   //        the code below might need to be updated accordingly as well.
//   auto& C = *Cptr;
//   C = 0;
//   for (const auto& e : htrans) {
//     const size_t i = std::get<0>(e);
//     const size_t j = std::get<1>(e); assert(i != j);
    
//     const double t1 = std::get<2>(e);
//     const double t2 = std::get<3>(e);
//     const double h1 = aperture[i];
//     const double h2 = aperture[j];
//     const double p1 = pressure[i];
//     const double p2 = pressure[j];
    
//     const double q   = h1 * h1 * h2 * h2 * t1 * t2; // numerator
//     const double d1q = 2 * h1 * h2 * h2 * t1 * t2;
//     const double d2q = h1 * h1 * 2 * h2 * t1 * t2;
    
//     const double r   = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
//     const double d1r = 24 * h1 * t1;
//     const double d2r = 24 * h2 * t2;
    
//     const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
//     const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);
    
//     // diagonal elements
//     C[i][i] += dTdh1 * (p1-p2);
//     C[j][j] += dTdh2 * (p2-p1); 
    
//     // off-diagonal elements
//     C[i][j] += dTdh2 * (p1-p2);
//     C[j][i] += dTdh1 * (p2-p1); 
    
//   }
  
// }

// // ----------------------------------------------------------------------------  
// inline bool convergence_test(const VectorHP& res, const double tol)
// // ----------------------------------------------------------------------------
// {
//   return true; //@@@@ res.infinity_norm() < tol;
// }

// }; // end anonymous namespace

namespace Opm
{

// ----------------------------------------------------------------------------
bool Fracture::fullSystemIteration(const double tol)
// ----------------------------------------------------------------------------
{
  // // update pressure matrix (fracture matrix is constant and does not need updating)
  // assemblePressure(); 
  // updateCouplingMatrix(coupling_matrix_, pressure_matrix_,
  //                      htrans_, fracture_pressure_, fracture_width_);

  // const auto& A = *fracture_matrix_;
  // const auto& M = *pressure_matrix_;
  // const auto& C = *coupling_matrix_;
  // const auto I = makeIdentity(A.N());

  // // system Jacobian (with cross term)  @@ should S be included as a member variable of Fracture?
  // SystemMatrix S { {A, I},   // mechanics system
  //                  {C, M} }; // flow system

  // // system equations
  // SystemMatrix S0 = S; S0[_1][_0] = 0; // the equations themselves have no cross term

  // VectorHP x {fracture_width_, fracture_pressure_};
  // vectorHP dx = x; dx = 0; // gradient of 'x' (which we aim to compute

  // // initialize right-hand aide
  // VectorHP rhs {x}; rhs = 0; // same size as x, but initially just with zeros

  // rhs[_0] = 0; // @@ we currently do not use rhs_width_ here, but include p through the I block
  // rhs[_1] = rhs_pressure_; // should have been updated in call to `assemblePressure` above

  // S0.mmv(x, rhs); // rhs = rhs - S0 * x;   (we are working in the tanget plane)

  // // check if system is already at a converged state (in which case we return immediately)
  // if (convergence_test(rhs, tol))
  //   return true;
                         
  // // solve system equations
  // const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S); 
  // TailoredPrecondDiag precond(S); // cannot be 'const' due to BiCGstabsolver interface
  // Dune::InverseOperatorResult iores; // cannot be 'const' due to BiCGstabsolver interface
  // auto psolver = Dune::BiCGSTABSolver<VectorHP>(S_linop,
  //                                               precond,
  //                                               1e-20, // desired rhs reduction factor
  //                                               200, // max number of iterations
  //                                               1); // verbose
  // psolver.apply(dx, rhs, iores); // NB: will modify 'rhs'

  // std::cout << "x:  " << x[_0].infinity_norm() << " " << x[_1].infinity_norm() << std::endl;
  // std::cout << "dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm() << std::endl;
  
  // // the following is a heuristic way to limit stepsize to stay within convergence radius
  // const double step_fac = estimate_step_fac(x, dx);
  // std::cout << "fac: " << step_fac << std::endl;

  // dx *= step_fac;
  // x += dx;

  // // copying modified variables back to member variables
  // fracture_width_ = x[_0];
  // fracture_pressure_ = x[_1];
  
  return false;
}
    
}; // end namespace Opm
