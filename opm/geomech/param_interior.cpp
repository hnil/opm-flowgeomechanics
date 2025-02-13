#include <array>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <assert.h>

//#include <dune/common/densematrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include "param_interior.hpp"

using namespace Opm;

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

// ----------------------------------------------------------------------------  
Axis3D determine_suitable_axis_2D(const std::vector<double>& p3d)
// ----------------------------------------------------------------------------  
{
  const size_t N = p3d.size() / 3;
  assert(N >= 3);
  using Entry = std::tuple<double, size_t>;
  
  // find boundingbox (as well as indices to the points that determined its bounds
  std::array<Entry, 6> bbox {{ {p3d[0], 0}, {p3d[1], 0}, {p3d[2], 0},
                               {p3d[0], 0}, {p3d[1], 0}, {p3d[2], 0}}};
  
  for (size_t i = 1; i != N; ++i) {
    for (int dim = 0; dim != 3; ++dim) {
      bbox[dim]     = std::get<0>(bbox[dim])     < p3d[3*i+dim] ? bbox[dim]     : Entry {p3d[3*i+dim], i}; // min
      bbox[dim + 3] = std::get<0>(bbox[dim + 3]) > p3d[3*i+dim] ? bbox[dim + 3] : Entry {p3d[3*i+dim], i}; // max
    }
  }
  
  std::array<double, 3> origin {0, 0, 0};

  // find center point (which should lie on the sought-after 2d plane)
  for (int dim = 0; dim != 3; ++dim) 
    origin[dim] = (std::get<0>(bbox[dim]) + std::get<0>(bbox[dim+3])) / 2;

  // make function for normalizing an 3D vector represented as an array<double, 3>:
  auto normalize = [](const std::array<double, 3>& arg) {
    double norm = std::sqrt(arg[0] * arg[0] + arg[1] * arg[1] + arg[2] * arg[2]);
    return std::array<double, 3>{arg[0] / norm, arg[1] / norm, arg[2] / norm};
  };

  // find two good axes, start out by identify the dimension along which the
  // bounding box has largest extent
  int ix_longest = 0;
  double l = 0;
  for (int dim = 0; dim != 3; ++dim) {
    double tmp = std::get<0>(bbox[dim+3]) - std::get<0>(bbox[dim]);
    ix_longest = tmp > l ? dim : ix_longest;
    l  = tmp > l ? tmp : l;
  }

  // decide on first axis (which we store in the variable 'ax1')
  const size_t pix = std::get<1>(bbox[ix_longest+3]); // point defining the upper bound of the bounding box in the chosen direction
  std::array<double, 3> ax1(normalize( {p3d[3 * pix]     - origin[0],
                                        p3d[3 * pix + 1] - origin[1],
                                        p3d[3 * pix + 2] - origin[2]}));

  // identify second axis by among the five remaining points defining the bounding box
  std::array<double, 3> ax2(ax1); // start out with overlapping axes, which
                                  // obviously has to change as we search for a
                                  // better pick
  double scal_prod = 1;  // ax1 . ax2
  for (int i = 0; i != 6; ++i) {
    if (i != ix_longest + 3) { // skip the direction that was used to construct 'ax1'
      const size_t point_ix = std::get<1>(bbox[i]);
      auto cand(normalize( {p3d[3 * point_ix]     - origin[0],
                            p3d[3 * point_ix + 1] - origin[1],
                            p3d[3 * point_ix + 2] - origin[2]}));
      const double cand_scal_prod = cand[0] * ax1[0] + cand[1] * ax1[1] + cand[2] * ax1[2];
      if (std::abs(cand_scal_prod) < std::abs(scal_prod)) {
        scal_prod = cand_scal_prod;
        ax2 = cand;
      }
    }
  }
  // axes ax1 and ax2 are now chosen to avoid degeneracy. Now, we modify ax2 to assure perpendicularity
  for (int dim = 0; dim != 3; ++dim) 
    ax2[dim] -= ax1[dim] * scal_prod;
  ax2 = normalize(ax2);

  return {origin, ax1, ax2};
}

  
}; // end anonymous namespace


  
namespace Opm{

// ----------------------------------------------------------------------------
  Axis3D project_to_2D(const std::vector<double>& p3d,
                       std::vector<double>& p2d)
// ----------------------------------------------------------------------------
{
  // determine suitable axis
  const Axis3D ax = determine_suitable_axis_2D(p3d);

  const size_t N = p3d.size() / 3;
  p2d = std::vector<double>(2*N, 0);

  // computing 2D points by projecting on (normalized) axes
  for (size_t i = 0; i != N; ++i) // loop over points
    for (size_t d2 = 0; d2 != 2; ++d2) // loop over 2D
      for  (size_t d3 = 0; d3 != 3; ++d3) // loop over 3D
        p2d[2 * i + d2] += (p3d[3*i + d3] - ax[0][d3]) * ax[1 + d2][d3];
  
  return ax;
}

// ----------------------------------------------------------------------------
void lift_to_3D(const std::vector<double>& p2d, const Axis3D& ax, std::vector<double>& p3d)
// ----------------------------------------------------------------------------
{
  const size_t N = p2d.size() / 2;
  p3d = std::vector<double>(3*N, 0);
        
  for (size_t i = 0; i != N; ++i) {// loop over points
    for (size_t d3 = 0; d3 != 3; ++d3) {// loop over 3D
      p3d[3 * i + d3] = ax[0][d3];
      for  (size_t d2 = 0; d2 != 2; ++d2) {// loop over 2D
        p3d[3 * i + d3] += p2d[2*i + d2] * ax[1 + d2][d3];
      }
    }
  }
}
  
// ----------------------------------------------------------------------------
// redistribute a set of 2D points on a loop so that they are equally
// distributed around the loop.  The first point is the "anchor point", which will
// remain in place.
void redistribute_2D(const std::vector<double>& points, std::vector<double>& result)
// ----------------------------------------------------------------------------
{
  // calculate total length around the loop
  const size_t N = points.size() / 2;
  double total_length = 0;
  for (size_t i = 0; i < N; ++i) {
    const size_t j = (i + 1) % N;
    total_length += sqrt((points[2*i] - points[2*j]) * (points[2*i] - points[2*j]) +
                         (points[2*i+1] - points[2*j+1]) * (points[2*i+1] - points[2*j+1]));
  }

  // calculate the length of each segment
  const double segment_length = total_length / N;

  // redistribute the points to lie on the loop, but equally spaced
  result.resize(N * 2);
  size_t ix_next = 1;
  std::array<double, 2> last_corner { points[0], points[1] };
  std::copy(last_corner.begin(), last_corner.end(), result.begin()); // first point is fixed

  // define distance function to use in loop below
  auto dist = [](const double* p1, const double* p2) {
    return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                (p1[1] - p2[1]) * (p1[1] - p2[1]));
  };
  
  // placing the other N-1 points equally spaced
  double remaining_seglength = segment_length;
  size_t count = 1;
  while (count < N) {
    const double len = dist(&last_corner[0], &points[2 * ix_next]);
    if (len >= remaining_seglength) {
      const double alpha = remaining_seglength / len;
      for (int i = 0; i != 2; ++i) {
        last_corner[i] = (1-alpha) * last_corner[i] + alpha * points[2*ix_next + i];
      }
      std::copy(last_corner.begin(), last_corner.end(), result.begin() + 2*count);
      remaining_seglength = segment_length;
      count++;
    } else {
      remaining_seglength -= len;
      last_corner[0] = points[2*ix_next];
      last_corner[1] = points[2*ix_next + 1];
      ix_next = (ix_next + 1) % N;
    }
  }
}

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
