#pragma once

#include <tuple>
#include <vector>

#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix

namespace // anonymous namespace
{ 
  using Htrans = std::tuple<size_t, size_t, double, double>;
  using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
  using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;   // sparse matrix
  using FMatrix = Dune::DynamicMatrix<double>; // full matrix

};  // end anonoymous namespace



namespace Opm
{
  
void solve_fully_coupled(ResVector& pvec, // output: fracture pressure
                         ResVector& hvec, // output: aperture
                         const SMatrix& pmat, // pressure matrix (sparse)
                         const FMatrix& amat, // aperture matrix (full)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<size_t> ratecells,
                         const double bhp,
                         const std::vector<size_t> bhpcells,
                         const int max_nonlin_iter = 100,
                         const double conv_tol=1e-5);

}; // end namespace Opm
