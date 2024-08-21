#include <array>
#include <vector>
#include <algorithm>

//#include <dune/common/densematrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include "param_interior.hpp"

namespace {
  
// ----------------------------------------------------------------------------
std::array<double, 3> barycentric(const double* const p1,
                                  const double* const p2,
                                  const double* const p3,
                                  const double* q)
// ----------------------------------------------------------------------------
{
  // Compute barycentric coordinates of q in terms of triangle corners p1, p2,
  // and p3 NB: If q is outside the triangle, at least one of the barycentric
  // coordinates will be negative.
  //
  // solve | a  b | |t1   | r
  //       |      | |   = |
  //       | c  d | |t2   | s

  const double a = p1[0] - p3[0];
  const double b = p2[0] - p3[0];
  const double c = p1[1] - p3[1];
  const double d = p2[1] - p3[1];
  const double r = q[0] - p3[0];
  const double s = q[1] - p3[1];

  const double denom = a *d - b * c;

  const double t1 = (r * d - s * b) / denom;
  const double t2 = (s * a - r * c) / denom;

  return std::array<double, 3> {t1, t2, 1-t1-t2};  // barycentric coordinates of q
  
}

// ----------------------------------------------------------------------------
std::vector<double> partition_of_unity(const std::vector<double>& bpoints,
                                       const std::vector<double>& ipoints)
// ----------------------------------------------------------------------------
{
  const size_t N = bpoints.size() / 2;
  const size_t M = ipoints.size() / 2;

  std::vector<double> result(N*M);

  for (size_t i = 0; i != N; ++i) {// loop over boundary points
    for (size_t i2 = 1; i2 != N-1; ++i2) {// loop over neighbor boundary points
      for (size_t j = 0; j != M; ++j) { // loop over interior points

        const size_t p2ix = (i + i2) % N;
        const size_t p3ix = (p2ix + 1) % N;
        const auto alpha = barycentric(&bpoints[2*i],
                                       &bpoints[2 * p2ix],
                                       &bpoints[2*p3ix],
                                       &ipoints[2*j]);
        if (alpha[0] >= 0 &&
            alpha[1] >= 0 &&
            alpha[2] >= 0) {
          result[j * N + i] = alpha[0];
        }
      }
    }
  }

  // normalize
  for (size_t j = 0; j != M; ++j) { // loop over ipoints
    double sum = 0;
    for (size_t i = 0; i != N; ++i) { // loop over bpoints
      sum += result[j * N + i];
    }
    for (size_t i = 0; i != N; ++i) {
      result[j * N + i] /= sum;
    }
  }
  
  return result; 
}

// ----------------------------------------------------------------------------
void compute_parameters(const std::vector<double>& bpoints,
                        const double* const point,
                        const double* const punity,
                        double* target)
// ----------------------------------------------------------------------------  
{
  const size_t N = bpoints.size() / 2;
  
  Dune::DynamicMatrix<double> A(N+3, N+3);
  Dune::DynamicVector<double> b(N+3); b = 0;
  A = 0; b = 0; // initialize elements to zero  

  // matrix derives from lagrangian function, and should be on form:
  // | 2*I                , (bpoints * punity), punity |
  // | (bpoints * punity)', 0                 , 0      |
  // | punity'            , 0                 , 0      |

  // right hand side should be:
  //     | 2
  // b = | point_x
  //     | point_y
  //     | 1
  
  for (size_t row = 0; row != N; ++row) {
    A[row][row] = 2;
    A[row][N]   = bpoints[2*row]   * punity[row]; // x coordinate
    A[row][N+1] = bpoints[2*row+1] * punity[row]; // y coordinate
    A[row][N+2] = punity[row];

    A[N][row]   = A[row][N];
    A[N+1][row] = A[row][N+1];
    A[N+2][row] = A[row][N+2];

    b[row] = 2;
  }

  b[N]   = point[0];
  b[N+1] = point[1];
  b[N+2] = 1;

  Dune::DynamicVector<double> res(N+3); res = 0;

  A.solve(res, b);

  for (size_t i = 0; i != N; ++i)
    target[i] = res[i] * punity[i];
    
}
  
}; // end anonymous namespace

namespace Opm{

// ----------------------------------------------------------------------------
void parametrize_interior_2D(const std::vector<double>& bpoints,
                             const std::vector<double>& ipoints,
                             std::vector<double>& result)
// ----------------------------------------------------------------------------
{
  
  const int num_bpoints = bpoints.size() / 2; // number of boundary points
  const  int num_ipoints = ipoints.size() / 2; // number of interior points

  const std::vector<double> punity = partition_of_unity(bpoints, ipoints);

  result.resize(num_bpoints * num_ipoints);
  fill(result.begin(), result.end(), 0);
  
  for (int i = 0; i != num_ipoints; ++i) 
    // compute parameterization for interior point 'i'
    compute_parameters(bpoints, &ipoints[2*i],
                       &punity[num_bpoints * i], &result[num_bpoints * i]);
}


  
}; // end namespace Opm
