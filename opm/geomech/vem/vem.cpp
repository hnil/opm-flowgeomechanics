/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include "vem.hpp"
#include <opm/common/OpmLog/OpmLog.hpp>
#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace vem;

namespace vem
{

namespace {
FaceQuadrature face_quadrature_ = FaceQuadrature::KITE;
}
void set_face_quadrature(FaceQuadrature q) { face_quadrature_ = q; }
FaceQuadrature face_quadrature() { return face_quadrature_; }

// ============================================================================
// Very basic helpers (since we here avoid having a defined point class)
// ============================================================================

// ----------------------------------------------------------------------------
// Compute an arbitrary linear combination of two ND-points, and write result
// to target
// ('plc' = point linear combination)
template <int Dim>
void
plc(const double* p1, const double* const p2, const double fac1, const double fac2, double* target)
// ----------------------------------------------------------------------------
{
    transform(
        p1, p1 + Dim, p2, target, [fac1, fac2](double d1, double d2) { return d1 * fac1 + d2 * fac2; });
}

// ----------------------------------------------------------------------------
// Compute an arbitrary linear combination of two ND-points, and return
// result as array
// ('plc' = point linear combination)
template <int Dim>
array<double, Dim>
plc(const double* const p1, const double* const p2, const double fac1, const double fac2)
// ----------------------------------------------------------------------------
{
    array<double, Dim> result;
    plc<Dim>(p1, p2, fac1, fac2, &result[0]);
    return result;
}

// ----------------------------------------------------------------------------
// Compute sum of two points, return result as array
template <int Dim>
array<double, Dim>
pointsum(const double* const p1, const double* const p2)
{
    return plc<Dim>(p1, p2, 1.0, 1.0);
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Compute sum of two points, write result to target
template <int Dim>
void
pointsum(const double* const p1, const double* const p2, double* target)
{
    plc<Dim>(p1, p2, 1.0, 1.0, target);
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Compute difference of two points, return result as array
template <int Dim>
array<double, Dim>
pointdiff(const double* const p1, const double* const p2)
{
    return plc<Dim>(p1, p2, 1.0, -1.0);
}

// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Compute difference of two points, write result to target
template <int Dim>
void
pointdiff(const double* const p1, const double* const p2, double* target)
{
    plc<Dim>(p1, p2, 1.0, -1.0, target);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Compute mean of two points, return result as array
template <int Dim>
array<double, Dim>
pointmean(const double* const p1, const double* const p2)
{
    return plc<Dim>(p1, p2, 0.5, 0.5);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Compute mean of two points, write result to target
template <int Dim>
void
pointmean(const double* const p1, const double* const p2, double* target)
{
    plc<Dim>(p1, p2, 0.5, 0.5, target);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Compute L2-norm of a N-dimensional vector
template <int Dim>
double
norm(const double* const v)
// ----------------------------------------------------------------------------
{
    const double sumsquared
        = accumulate(v, v + Dim, 0.0, [](double acc, double el) { return acc + el * el; });
    return sqrt(sumsquared);
}

// ----------------------------------------------------------------------------
// Compute L2-norm of a N-dimensional vector
template <int Dim>
double
norm(const array<double, Dim>& arr)
{
    return norm<Dim>(&arr[0]);
};
// ----------------------------------------------------------------------------


// ============================================================================
// Other helper functions
// ============================================================================

// ----------------------------------------------------------------------------
// Compute centroid (i.e. geometric center) of (approximately) planar face
// embedded in 2D or 2D space, return result as array
template <int Dim>
array<double, Dim> face_centroid(const double* const points, const int num_points);
template <>
array<double, 2>
face_centroid<2>(const double* const points, const int num_points)
{
    return centroid_2D(points, num_points);
}
template <>
array<double, 3>
face_centroid<3>(const double* const points, const int num_points)
{
    return centroid_2D_3D(points, num_points);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Compute determinant of 3D matrix, given pointer to its three rows (or columns)
double
determinant3D(const double* const c1, const double* const c2, const double* const c3)
// ----------------------------------------------------------------------------
{
    return c1[0] * (c2[1] * c3[2] - c2[2] * c3[1]) - c1[1] * (c2[0] * c3[2] - c2[2] * c3[0])
        + c1[2] * (c2[0] * c3[1] - c2[1] * c3[0]);
};

// ----------------------------------------------------------------------------
// Compute the "average" point of a number of points in N-dimensional space,
// defined as the coordinate mean of the points.  Note that this point is not
// necessarily the same as the geometric centroid.
template <int Dim>
array<double, Dim>
point_average(const double* const corners, const int num_corners)
// ----------------------------------------------------------------------------
{
    array<double, Dim> result;
    fill(result.begin(), result.end(), 0.0);

    // accumulate coordinate values
    for (int i = 0; i != num_corners; ++i)
        pointsum<Dim>(&result[0], &corners[Dim * i], &result[0]);

    // divide by number of corners
    for (auto& c : result)
        c /= num_corners;

    return result;
}

// ----------------------------------------------------------------------------
// Pick a selection of points from those listed in 'pts', and return them
// consecutively in a vector
template <int Dim>
vector<double>
pick_points(const double* const pts, const int* const p_ixs, int num_points)
// ----------------------------------------------------------------------------
{
    vector<double> result;
    result.reserve(Dim * num_points);

    for (int i = 0; i != num_points; ++i)
        for (int d = 0; d != Dim; ++d)
            result.push_back(pts[p_ixs[i] * Dim + d]);

    return result;
}

// ----------------------------------------------------------------------------
// Compute the area of a triangle embedded in N-dimensional space, given
// pointers to the coordinates of its three corners.  Heron's formula is used,
// so there is no notion of signed area.
template <int Dim>
double
triarea(const double* const c1, const double* const c2, const double* const c3)
// ----------------------------------------------------------------------------
{
    // computing area of triangle using Heron's formula
    const double l1 = norm<Dim>(pointdiff<Dim>(c2, c1));
    const double l2 = norm<Dim>(pointdiff<Dim>(c3, c2));
    const double l3 = norm<Dim>(pointdiff<Dim>(c1, c3));

    const double s = 0.5 * (l1 + l2 + l3);

    return sqrt(s * (s - l1) * (s - l2) * (s - l3));
}

// ----------------------------------------------------------------------------
// Computes the normal of a triangle embedded in 3D space.  The normal
// returned is scaled by the area of the triangle.
// (This function only makes sense in 3D, so no template on dimension used.)
array<double, 3>
trinormal(const double* const c1, const double* const c2, const double* const c3)
// ----------------------------------------------------------------------------
{
    const array<double, 3> v1 {c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]};
    const array<double, 3> v2 {c3[0] - c1[0], c3[1] - c1[1], c3[2] - c1[2]};

    // return the cross product
    return array<double, 3> {
        v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]};
}

// corner_mask flags the face nodes that are grid corners; null, or none flagged, means all are
bool
has_hanging_nodes(const char* const corner_mask, const int n)
{
    return corner_mask && std::find(corner_mask, corner_mask + n, 0) != corner_mask + n
        && std::find(corner_mask, corner_mask + n, 1) != corner_mask + n;
}

// Mean of a face's corners, the fan centre; fv indexes corners.
array<double, 3>
fan_centre(const double* const corners, const int* const fv, const int m, const char* const corner_mask)
{
    const char* const mask = has_hanging_nodes(corner_mask, m) ? corner_mask : nullptr;
    int nc = 0;
    for (int k = 0; k != m; ++k)
        nc += !mask || mask[k];
    array<double, 3> c {0, 0, 0};
    for (int k = 0; k != m; ++k)
        if (!mask || mask[k])
            for (int d = 0; d != 3; ++d)
                c[d] += corners[3 * fv[k] + d] / nc;
    return c;
}

// ----------------------------------------------------------------------------
// Tessellates a face defined by a set of corner points into triangles.
//
// @tparam Dim The dimension of the space in which the points lie.
// @param corners Pointer to the array of corner points, where each point is a
//        sequence of Dim values. The points must be ordered consecutively.
// @param num_corners The number of corner points that define the face.
// @param skip_if_tri If true and the face is already a triangle, returns the face
//        unchanged. If false, still performs the computation to split the face
//        into triangles.
// @param centroid If a pointer that is not a nullpointer is given here, the
//        3-component face centroid will be returned here
// @return A vector of vectors, where each inner vector represents the coordinates
//         of a triangle that together tessellate the face. Each inner vector has
//         3*Dim values: the first Dim values are the coordinates of the centroid of
//         the face, and the remaining 2*Dim values are the coordinates of the
//         corner points of the triangle.
template <int Dim>
vector<vector<double>>
tessellate_face(const double* const corners,
                const int num_corners,
                bool skip_if_tri = true,
                double* centroid = nullptr,
                const char* const corner_mask = nullptr)
// ----------------------------------------------------------------------------
{
    if (skip_if_tri && num_corners == 3)
        // the face is already a triangle - no need to split it up
        return {vector<double>(corners, corners + Dim * num_corners)};

    // the centre comes from the face's grid corners, not its hanging nodes
    vector<double> cpts;
    if (has_hanging_nodes(corner_mask, num_corners))
        for (int c = 0; c != num_corners; ++c)
            if (corner_mask[c])
                cpts.insert(cpts.end(), &corners[Dim * c], &corners[Dim * (c + 1)]);
    const array<double, Dim> center = cpts.empty()
        ? face_centroid<Dim>(corners, num_corners)
        : face_centroid<Dim>(cpts.data(), int(cpts.size()) / Dim);

    if (centroid != nullptr)
        copy(center.begin(), center.end(), centroid);

    vector<vector<double>> result;
    for (int c = 0; c != num_corners; ++c) {
        const int cnext = (c + 1) % num_corners;

        // midpoint between current and next corner point (along edge)
        const array<double, Dim> midpt = pointmean<Dim>(&corners[Dim * c], &corners[Dim * cnext]);
        // insert first triangle
        vector<double> tri(center.begin(), center.end());
        tri.insert(tri.end(), &corners[Dim * c], &corners[Dim * (c + 1)]);
        tri.insert(tri.end(), &midpt[0], &midpt[Dim]);
        result.push_back(tri);

        // insert second triangle
        tri.clear();
        tri.insert(tri.end(), center.begin(), center.end());
        tri.insert(tri.end(), &midpt[0], &midpt[Dim]);
        tri.insert(tri.end(), &corners[Dim * cnext], &corners[Dim * (cnext + 1)]);
        result.push_back(tri);
    }
    return result;
}

// ----------------------------------------------------------------------------
// Tessellates a face defined by its corner points into triangles.
// If the face is already a triangle and `skip_if_tri` is set to `true`,
// then the original triangle is returned as-is.
//
// @tparam Dim Dimensionality of the space (2 or 3)
// @param corners Vector of corner points of the face.
// @param skip_if_tri If set to `true` and the face is already a triangle,
//                    the original triangle is returned as-is.
// @return A vector of vectors, where each inner vector contains the coordinates
//         of one triangle in the tessellation. The vertices of each triangle are
//         listed in counterclockwise order.
template <int Dim>
vector<vector<double>>
tessellate_face(const vector<double>& corners,
                bool skip_if_tri = true,
                double* centroid = nullptr,
                const char* const corner_mask = nullptr)
// ----------------------------------------------------------------------------
{
    return tessellate_face<Dim>(&corners[0], corners.size() / Dim, skip_if_tri, centroid, corner_mask);
}

// ----------------------------------------------------------------------------
template <int Dim>
double
face_integral_impl(const double* const corners,
                   const int num_corners,
                   const double* corner_values,
                   const char* const corner_mask = nullptr)
// ----------------------------------------------------------------------------
{
    const auto tris = tessellate_face<Dim>(&corners[0], num_corners, false, nullptr, corner_mask);

    assert(int(tris.size()) == 2 * num_corners);

    double result = 0;

    for (int i = 0; i != num_corners; ++i) {
        const auto& t1 = tris[2 * i];
        const auto& t2 = tris[2 * i + 1];

        const int inext = (i + 1) % num_corners;
        const double fval1 = corner_values ? corner_values[i] : 1.0;
        const double fval2 = corner_values ? corner_values[inext] : 1.0;

        result += triarea<Dim>(&t1[0], &t1[Dim], &t1[2 * Dim]) * fval1;
        result += triarea<Dim>(&t2[0], &t2[Dim], &t2[2 * Dim]) * fval2;
    }
    return result;
}

// ----------------------------------------------------------------------------
// Computes the matrix multiplication between two matrices and stores the result
// in a given array.
//
// @param data1 Pointer to the first matrix data.
// @param r1 Number of rows in the first matrix.
// @param c1 Number of columns in the first matrix.
// @param transposed1 Flag indicating whether to transpose the first matrix.
// @param data2 Pointer to the second matrix data.
// @param r2 Number of rows in the second matrix.
// @param c2 Number of columns in the second matrix.
// @param transposed2 Flag indicating whether to transpose the second matrix.
// @param result Pointer to the memory location where the result will be stored.
// @param fac Optional scaling factor to apply to the result. Default value is 1.
void
matmul(const double* data1,
       const int r1,
       const int c1,
       bool transposed1,
       const double* data2,
       const int r2,
       const int c2,
       bool transposed2,
       double* result,
       const double fac = 1)
// ----------------------------------------------------------------------------
{
    const pair<int, int> dim1 = transposed1 ? make_pair(c1, r1) : make_pair(r1, c1);
    const pair<int, int> dim2 = transposed2 ? make_pair(c2, r2) : make_pair(r2, c2);
    const pair<int, int> stride1 = transposed1 ? make_pair(1, c1) : make_pair(c1, 1);
    const pair<int, int> stride2 = transposed2 ? make_pair(1, c2) : make_pair(c2, 1);

    if (dim1.second != dim2.first)
        throw invalid_argument("Matrices are not compatible for multiplication.");

    const int num_elements = dim1.first * dim2.second;
    fill(result, result + num_elements, 0);

    for (int r = 0; r != dim1.first; ++r)
        for (int c = 0; c != dim2.second; ++c)
            for (int k = 0; k != dim1.second; ++k)
                result[r * dim2.second + c] += data1[r * stride1.first + k * stride1.second]
                    * data2[k * stride2.first + c * stride2.second];

    // Multiply matrix by a scalar factor, if requested
    if (fac != 1.0)
        for_each(result, result + num_elements, [fac](double& e) { e *= fac; });
}

// ----------------------------------------------------------------------------
// Computes the matrix multiplication between two matrices and returns the result as a vector.
//
// @param data1 pointer to the data of the first matrix
// @param r1 number of rows of the first matrix
// @param c1 number of columns of the first matrix
// @param transposed1 flag indicating whether the first matrix should be transposed before multiplication
// @param data2 pointer to the data of the second matrix
// @param r2 number of rows of the second matrix
// @param c2 number of columns of the second matrix
// @param transposed2 flag indicating whether the second matrix should be transposed before
// multiplication
// @param fac scaling factor applied to the multiplication (default is 1)
//
// @return a vector containing the result of the multiplication
vector<double>
matmul(const double* data1,
       const int r1,
       const int c1,
       bool transposed1,
       const double* data2,
       const int r2,
       const int c2,
       bool transposed2,
       const double fac = 1)
// ----------------------------------------------------------------------------
{
    vector<double> result((transposed1 ? c1 : r1) * (transposed2 ? r2 : c2));
    matmul(data1, r1, c1, transposed1, data2, r2, c2, transposed2, &result[0], fac);
    return result;
}

// ----------------------------------------------------------------------------
// compute the trace of a (full) matrix A, whose elements are stored in a vector
// (row-major or column-major order does not really matter here).
double
trace(const vector<double>& A)
// ----------------------------------------------------------------------------
{
    const auto N = std::size_t(std::sqrt(A.size()));
    assert(N * N == A.size()); // should be a square matrix

    double result = 0;
    for (std::size_t i = 0; i < N; ++i)
        result += A[i + N * i];

    return result;
}

// ----------------------------------------------------------------------------
// extract the diagonal elements of a matrix
vector<double>
diag_elems(const vector<double>& A)
// ----------------------------------------------------------------------------
{
    const int N = int(sqrt(A.size()));
    vector<double> result(N, 0);
    for (int i = 0; i != N; ++i)
        result[i] = A[i * N + i];
    return result;
}



// ----------------------------------------------------------------------------
// Compute the trace of the inverse matrix.  The function only works correctly
// for positive symmetric matrices
double
invtrace(const vector<double>& A)
// ----------------------------------------------------------------------------
{
    const auto diag = diag_elems(A);
    double result = 0;
    for (int i = 0; i != (int)diag.size(); ++i)
        result += 1.0 / diag[i];
    return result;
}

// ----------------------------------------------------------------------------
// Identity matrix of size N, returned as a vector of lenght N^2 containing.
// its elements (row-major or column-major order is equivalent here).
// The matrix is multiplied by the number 'fac' before being returned.
vector<double>
identity_matrix(const double fac, const int N)
// ----------------------------------------------------------------------------
{
    // create zero matrix
    vector<double> result(N * N, 0);

    // fill in diagonal elements with the value 'fac'
    for (int i = 0; i != N; ++i)
        result[i + N * i] = fac;

    // return the result
    return result;
}

// ----------------------------------------------------------------------------
// Compute the VEM stability term for a given 2D or 3D element.
// This is a very simple term that is described in (Gain, 2014)
// DOI:10.1016/j.cma.2014.05.005
vector<double>
compute_S(const std::vector<double>& Nc,
          const std::vector<double>& D,
          const int num_corners,
          const double volume,
          const int dim,
          const StabilityChoice stability_choice)
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);
    assert(stability_choice == SIMPLE || stability_choice == HARMONIC);
    const int r = dim * num_corners;
    const int c = dim == 2 ? 3 : 6;
    const auto NtN = matmul(&Nc[0], r, c, true, &Nc[0], r, c, false);

    const double gamma = 1.0 / 9.0; // since trace is about 9 times E
    const double alpha = stability_choice == SIMPLE ? volume * trace(D) / trace(NtN)
                                                    : // from Gain et al. 2014
        gamma*volume * trace(D) * invtrace(NtN); // Andersen et al. 2017

    return identity_matrix(alpha, dim * num_corners);
}

// ----------------------------------------------------------------------------
vector<double>
compute_S_D_recipe(const std::vector<double>& EWcDWct, const int dofs, const double volume)
// ----------------------------------------------------------------------------
{
    vector<double> result(dofs * dofs, 0);

    const double h = cbrt(volume); // diameter of cell scales with cube root of volume
    for (int i = 0; i != dofs; ++i)
        result[i + dofs * i] = EWcDWct[i + dofs * i]; // extract diagonal from matrix

    return result;
}

// ----------------------------------------------------------------------------
// Eigendecomposition of a small (dim <= 3) symmetric matrix by cyclic Jacobi
// rotations.  On return, `A` holds the eigenvalues on its diagonal and the
// columns of `V` are the corresponding eigenvectors.
void
small_sym_eig(const int dim, double A[3][3], double V[3][3])
// ----------------------------------------------------------------------------
{
    for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 3; ++j)
            V[i][j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep != 50; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (int p = 0; p != dim; ++p) {
            diag += A[p][p] * A[p][p];
            for (int q = p + 1; q != dim; ++q)
                off += A[p][q] * A[p][q];
        }
        if (off <= 1e-30 * (diag + 1e-300))
            break;
        for (int p = 0; p != dim; ++p)
            for (int q = p + 1; q != dim; ++q) {
                if (A[p][q] == 0.0)
                    continue;
                const double theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
                const double t = (theta >= 0 ? 1.0 : -1.0)
                    / (fabs(theta) + sqrt(theta * theta + 1));
                const double c = 1.0 / sqrt(t * t + 1);
                const double s = t * c;
                for (int k = 0; k != dim; ++k) {
                    const double akp = A[k][p], akq = A[k][q];
                    A[k][p] = c * akp - s * akq;
                    A[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k != dim; ++k) {
                    const double apk = A[p][k], aqk = A[q][k];
                    A[p][k] = c * apk - s * aqk;
                    A[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k != dim; ++k) {
                    const double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
            }
    }
}

// ----------------------------------------------------------------------------
// Anisotropic stabilization floor *tensor* of an element, computed in the
// principal frame of its corner point cloud (so it does not assume the cell to
// be aligned with the coordinate axes):
//
//   F = Escale * volume * sum_k  v_k v_k^T / h_k^2,
//
// where v_k are the principal axes (eigenvectors of the corner second-moment
// tensor) and h_k = 2*sqrt(lambda_k) the corresponding principal lengths (equal
// to the edge lengths for a box, in any orientation).  For an axis-aligned cell
// this reduces to diag(Escale * volume / h_d^2), i.e. the natural per-direction
// scale of the consistency-matrix diagonal (tr(D)/9 * cbrt(volume) for a cube).
array<double, 9>
compute_aniso_floor_tensor(const vector<double>& corners,
                           const int num_corners,
                           const int dim,
                           const double Escale,
                           const double volume)
// ----------------------------------------------------------------------------
{
    array<double, 3> centroid {0.0, 0.0, 0.0};
    for (int i = 0; i != num_corners; ++i)
        for (int d = 0; d != dim; ++d)
            centroid[d] += corners[dim * i + d];
    for (int d = 0; d != dim; ++d)
        centroid[d] /= num_corners;

    double M[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i = 0; i != num_corners; ++i)
        for (int d1 = 0; d1 != dim; ++d1)
            for (int d2 = 0; d2 != dim; ++d2)
                M[d1][d2] += (corners[dim * i + d1] - centroid[d1])
                    * (corners[dim * i + d2] - centroid[d2]) / num_corners;

    double V[3][3];
    small_sym_eig(dim, M, V);

    // principal lengths, guarded against (near-)degenerate point clouds
    double h[3], hmax = 0.0;
    for (int k = 0; k != dim; ++k) {
        h[k] = 2.0 * sqrt(max(M[k][k], 0.0)); // box with edge L: 2*sqrt(L^2/4) = L
        hmax = max(hmax, h[k]);
    }
    for (int k = 0; k != dim; ++k)
        h[k] = max(h[k], 1e-8 * hmax);

    array<double, 9> F {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int k = 0; k != dim; ++k) {
        const double w = Escale * volume / (h[k] * h[k]);
        for (int d1 = 0; d1 != dim; ++d1)
            for (int d2 = 0; d2 != dim; ++d2)
                F[3 * d1 + d2] += w * V[d1][k] * V[d2][k];
    }
    return F;
}

// ----------------------------------------------------------------------------
// Anisotropic stabilization: the floor tensor F above replaces the isotropic
// cbrt(volume)-based floor of D_RECIPE, which over-stabilizes the soft direction
// of high-aspect-ratio cells.  Two variants:
//  - use_consistency_diag (ANISO_DIAG): diagonal S, each dof floored by the
//    corresponding diagonal entry of F (frame-aware floor; the consistency
//    diagonal itself is naturally expressed in the coordinate frame);
//  - otherwise (ANISO_HARMONIC): block-diagonal S with the full tensor F as the
//    per-node block, which makes the stabilization frame *invariant* (the
//    element stiffness transforms exactly under rigid rotations of the cell).
vector<double>
compute_S_aniso(const std::vector<double>& EWcDWct,
                const std::vector<double>& D,
                const vector<double>& corners,
                const int num_nodes,
                const double volume,
                const int dim,
                const bool use_consistency_diag)
// ----------------------------------------------------------------------------
{
    const int dofs = dim * num_nodes;
    // tr(D) ~ 9E in 3D and ~4E in 2D (with the doubled shear entries of this D)
    const double Escale = trace(D) / ((dim == 3) ? 9.0 : 4.0);
    const auto F = compute_aniso_floor_tensor(corners, num_nodes, dim, Escale, volume);

    vector<double> result(dofs * dofs, 0.0);
    if (use_consistency_diag) {
        for (int i = 0; i != dofs; ++i) {
            const int d = i % dim;
            result[i + dofs * i] = max(EWcDWct[i + dofs * i], F[3 * d + d]);
        }
    } else {
        for (int node = 0; node != num_nodes; ++node)
            for (int d1 = 0; d1 != dim; ++d1)
                for (int d2 = 0; d2 != dim; ++d2)
                    result[(node * dim + d1) * dofs + node * dim + d2] = F[3 * d1 + d2];
    }
    return result;
}

std::vector<double> getStabilityMatrix(vem::StabilityChoice stability_choice,const std::vector<double>&  EWcDWct,const std::vector<double>& Nc,const std::vector<double>& D,
                                       const std::vector<double>& corners,
                                       int num_nodes,double volume,int  dim,double totdim,
                                       const std::vector<double>& node_weights = {}){
    std::vector<double> S;
    // FEM choices reaching this point means the cell could not be treated with Q1
    // elements; use the anisotropic diagonal stabilization for the VEM fallback.
    if(is_fem_choice(stability_choice)){
        stability_choice = ANISO_DIAG;
    }
    if (is_face_weighted_choice(stability_choice)) {
        const auto base = stability_choice == FACE_D_RECIPE ? D_RECIPE
            : stability_choice == FACE_ANISO_HARMONIC       ? ANISO_HARMONIC
                                                            : SIMPLE;
        const bool weighted = int(node_weights.size()) == num_nodes;
        if (base == D_RECIPE) {
            // only the isotropic floor is rescaled; the consistency diagonal is already area-weighted
            S = compute_S_D_recipe(EWcDWct, totdim, volume);
            const double val = trace(D) / 9.0 * cbrt(volume);
            for (int i = 0; i != totdim; ++i) {
                const double w = weighted ? node_weights[i / dim] : 1.0;
                S[i + totdim * i] = max(w * val, S[i + totdim * i]);
            }
            return S;
        }
        S = getStabilityMatrix(base, EWcDWct, Nc, D, corners, num_nodes, volume, dim, totdim);
        if (weighted) {
            for (int a = 0; a != totdim; ++a)
                for (int b = 0; b != totdim; ++b)
                    S[a + totdim * b] *= sqrt(node_weights[a / dim] * node_weights[b / dim]);
        }
        return S;
    }
    if (stability_choice == TRACE || stability_choice == SMOOTH_HOURGLASS)
        throw logic_error("TRACE/SMOOTH_HOURGLASS stabilization needs the face data");
    if(stability_choice == ANISO_DIAG){
        S = compute_S_aniso(EWcDWct, D, corners, num_nodes, volume, dim, true);
    } else if(stability_choice == ANISO_HARMONIC){
        S = compute_S_aniso(EWcDWct, D, corners, num_nodes, volume, dim, false);
    } else if(stability_choice == D_RECIPE){
        const auto SD= compute_S_D_recipe(EWcDWct, totdim, volume);
        S=SD;
        double val = trace(D)/9.0;
        val *= cbrt(volume);
        for (int i = 0; i != totdim; ++i){
            int ind = i + totdim * i;
            S[ind]= max(val,S[ind]);
        }
    } else if((stability_choice == SIMPLE) || (stability_choice == HARMONIC)){
        S  = compute_S(Nc, D, num_nodes, volume, dim, stability_choice);
    } else if( stability_choice == ALLMAX){
        const auto SD= compute_S_D_recipe(EWcDWct, totdim, volume);
        const auto SS = compute_S(Nc, D, num_nodes, volume, dim, SIMPLE);
        const auto SH = compute_S(Nc, D, num_nodes, volume, dim, HARMONIC);
        S = SD;
        for(int i=0;i< int(S.size());++i){
            S[i]= max(S[i],SS[i]);
            S[i]= max(S[i],SH[i]);
        }
    }else if( stability_choice == ALLMIN){
        const auto SD= compute_S_D_recipe(EWcDWct, totdim, volume);
        const auto SS = compute_S(Nc, D, num_nodes, volume, dim, SIMPLE);
        const auto SH = compute_S(Nc, D, num_nodes, volume, dim, HARMONIC);
        S = SD;
        for(int i=0;i< int(S.size());++i){
            S[i]= min(S[i],SS[i]);
            S[i]= min(S[i],SH[i]);
        }
    }else if( stability_choice == EXPERIMENTAL){
        double val = trace(D)/9.0;
        val *= 4.0;
        S = identity_matrix(val, totdim);
    }else{
        throw invalid_argument("Unknown stability choice in final_assembly");
    }
    return S;
}

// ----------------------------------------------------------------------------
// Trace stabilization  kappa * mu * sum_F |E|/|F| int_F |grad_F w|^2, w piecewise linear on the
// face fan around the vertex mean.  Scales with the face, not the cell, so short edges and small
// faces keep their stiffness; mu only, so no volumetric locking from S.
constexpr double trace_stab_factor = 0.25;

vector<double>
compute_S_trace_3D(const double* const corners,
                   const int num_corners,
                   const int* const faces,
                   const int* const num_face_edges,
                   const int num_faces,
                   const double volume,
                   const double young,
                   const double poisson,
                   const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    const int n = 3 * num_corners;
    const double mu = young / (2 * (1 + poisson));
    vector<double> S(n * n, 0.0);
    for (int f = 0, pos = 0; f != num_faces; pos += num_face_edges[f++]) {
        const int m = num_face_edges[f];
        const int* const fv = &faces[pos];
        const char* const fmask = corner_mask && has_hanging_nodes(corner_mask + pos, m) ? corner_mask + pos : nullptr;
        const array<double, 3> c = fan_centre(corners, fv, m, fmask);
        int nc = 0;
        for (int k = 0; k != m; ++k)
            nc += !fmask || fmask[k];
        double area = 0;
        vector<array<double, 9>> kloc(m); // per fan triangle (centre, a, b)
        for (int i = 0; i != m; ++i) {
            const double* pa = &corners[3 * fv[i]];
            const double* pb = &corners[3 * fv[(i + 1) % m]];
            const double* P[3] = {&c[0], pa, pb};
            const array<double, 3> u {pa[0] - c[0], pa[1] - c[1], pa[2] - c[2]};
            const array<double, 3> v {pb[0] - c[0], pb[1] - c[1], pb[2] - c[2]};
            const double A = 0.5 * norm<3>(array<double, 3> {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]});
            area += A;
            for (int a = 0; a != 3; ++a)
                for (int b = 0; b != 3; ++b) {
                    // P1 stiffness: (e_a . e_b) / (4A), e_a the edge opposite vertex a
                    double ea[3], eb[3];
                    for (int d = 0; d != 3; ++d) {
                        ea[d] = P[(a + 2) % 3][d] - P[(a + 1) % 3][d];
                        eb[d] = P[(b + 2) % 3][d] - P[(b + 1) % 3][d];
                    }
                    kloc[i][3 * a + b] = A > 0 ? (ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]) / (4 * A) : 0.0;
                }
        }
        if (area <= 0)
            continue;
        const double coef = trace_stab_factor * mu * volume / area;
        // scalar face matrix on the face vertices, the centre eliminated as the vertex mean
        vector<double> Kf(m * m, 0.0);
        for (int i = 0; i != m; ++i) {
            const int loc[3] = {-1, i, (i + 1) % m};
            for (int a = 0; a != 3; ++a)
                for (int b = 0; b != 3; ++b) {
                    const double k = coef * kloc[i][3 * a + b];
                    for (int r = 0; r != m; ++r) {
                        const double wa = loc[a] < 0 ? (!fmask || fmask[r] ? 1.0 / nc : 0.0)
                                                     : (r == loc[a] ? 1.0 : 0.0);
                        if (wa == 0.0)
                            continue;
                        for (int s2 = 0; s2 != m; ++s2) {
                            const double wb = loc[b] < 0 ? (!fmask || fmask[s2] ? 1.0 / nc : 0.0)
                                                         : (s2 == loc[b] ? 1.0 : 0.0);
                            Kf[r * m + s2] += k * wa * wb;
                        }
                    }
                }
        }
        for (int r = 0; r != m; ++r)
            for (int s2 = 0; s2 != m; ++s2)
                for (int d = 0; d != 3; ++d)
                    S[(3 * fv[r] + d) * n + 3 * fv[s2] + d] += Kf[r * m + s2];
    }
    return S;
}

// ----------------------------------------------------------------------------
// 2D version: kappa * mu * sum_e |E|/|e| int_e |d_t w|^2
vector<double>
compute_S_trace_2D(const double* const corners, const int num_corners, const double area,
                   const double young, const double poisson)
// ----------------------------------------------------------------------------
{
    const int n = 2 * num_corners;
    const double mu = young / (2 * (1 + poisson));
    vector<double> S(n * n, 0.0);
    for (int i = 0; i != num_corners; ++i) {
        const int j = (i + 1) % num_corners;
        const double le = hypot(corners[2 * j] - corners[2 * i], corners[2 * j + 1] - corners[2 * i + 1]);
        if (le <= 0)
            continue;
        const double k = trace_stab_factor * mu * area / (le * le);
        for (int d = 0; d != 2; ++d) {
            S[(2 * i + d) * n + 2 * i + d] += k;
            S[(2 * j + d) * n + 2 * j + d] += k;
            S[(2 * i + d) * n + 2 * j + d] -= k;
            S[(2 * j + d) * n + 2 * i + d] -= k;
        }
    }
    return S;
}

// ----------------------------------------------------------------------------
// Assemble the VEM stiffness matrix for a single element, based on a number
// of intermediary matrices and values, as defined in (Gain, 2014)
// DOI:10.1016/j.cma.2014.05.005
void
final_assembly(const vector<double>& Wc,
               const vector<double>& D,
               const vector<double>& Nc,
               const vector<double>& ImP,
               const vector<double>& corners,
               const StabilityChoice stability_choice,
               const double volume,
               const double num_nodes,
               const double dim,
               double* target,
               const vector<double>* S_given = nullptr,
               const vector<double>& node_weights = {})
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);
    const int lsdim = (dim == 2) ? 3 : 6; // dimension of "linear strain space"
    const int totdim = dim * num_nodes; // total number of unknowns

    // compute stiffness matrix components
    const auto DWct = matmul(&D[0], lsdim, lsdim, false, &Wc[0], totdim, lsdim, true);
    const auto EWcDWct = matmul(&Wc[0], totdim, lsdim, false, &DWct[0], lsdim, totdim, false, volume);

    std::vector<double> S = S_given ? *S_given
        : getStabilityMatrix(stability_choice, EWcDWct, Nc, D, corners, num_nodes, volume, dim, totdim, node_weights);


    const auto SImP = matmul(&S[0], totdim, totdim, false, &ImP[0], totdim, totdim, false);
    const auto ImpSImp = matmul(&ImP[0], totdim, totdim, true, &SImP[0], totdim, totdim, false);
    assert(EWcDWct.size() == ImpSImp.size());

    // add conformance and stability term, and write result to target
    transform(EWcDWct.begin(), EWcDWct.end(), ImpSImp.begin(), target, [](double a, double b) {
        return a + b;
    });
}

// ----------------------------------------------------------------------------
template <class T>
void
append(vector<T>& target, const vector<T>& elems)
// ----------------------------------------------------------------------------
{
    target.insert(target.end(), elems.begin(), elems.end());
}

// ----------------------------------------------------------------------------
// Create the sub-matrix associated with a particular 2D node; used to assemble
// Nr, Nc, Wr and Wc, which all have the same basic structure
vector<double>
matentry_2D(const double e1, const double e2, const double e3, const double e4)
// ----------------------------------------------------------------------------
{
    // return the matrix [e1,  0, e2;
    //                     0, e1, e4]
    // in row-major order
    return vector<double> {e1, 0, e2, 0, e3, e4};
}

// ----------------------------------------------------------------------------
// Create the sub-matrix associated with a particular 3D node; used to assemble
// Nr, Nc, Wr and Wc, which all have the same basic structure
vector<double>
matentry_3D(const double e1,
            const double e2,
            const double e3,
            const double e4,
            const double e5,
            const double e6,
            const double e7,
            const double e8,
            const double e9)
// ----------------------------------------------------------------------------
{

    return vector<double> {e1, 0, 0, e2, 0, e3, 0, e4, 0, e5, e6, 0, 0, 0, e7, 0, e8, e9};
}


// ----------------------------------------------------------------------------
// Compute and distribute the impact of a 2D force applied to a given 2D element
// area ("volume") among its nodes, such that the full impact of the force is
// equal to the sum of the forces applied at the nodes.
void
compute_bodyforce_2D(const double* const points,
                     const int* cell_corners, // corners of the 2D element
                     const int number_cell_faces, // number of corners
                     const double* bforce, // 2D body force
                     double* target)
// ----------------------------------------------------------------------------
{
    auto coords = pick_points<2>(points, cell_corners, number_cell_faces);

    const auto tris = tessellate_face<2>(&coords[0], number_cell_faces, false);
    assert(int(tris.size()) == 2 * number_cell_faces);

    // note: in the list of tesselated faces, the triangles associated with node 'i'
    // are tris[i-1] and tris[i].  Hence the somewhat complicated indexing below.
    for (int c = 0; c != 2 * number_cell_faces; ++c)
        for (int d = 0; d != 2; ++d)
            target[2 * ((c / 2 + c % 2) % number_cell_faces) + d]
                += bforce[d] * triarea<2>(&tris[c][0], &tris[c][2], &tris[c][4]);
}

// ----------------------------------------------------------------------------
// Compute and distribute the impact of a 2D force applied to a given edge
// such that the full impact of the force is equal to the sum of the forces
// appliedat its two end nodes.
array<double, 4>
compute_applied_forces_2D(
    const double* const points, const int n1, const int n2, const double fx, const double fy)
// ----------------------------------------------------------------------------
{
    const double hL = norm<2>(pointdiff<2>(&points[2 * n1], &points[2 * n2])) / 2; // 1/2 L

    // return {f_x on n1, f_y on n1, f_x on n2, f_y on n2}
    return {hL * fx, hL * fy, hL * fx, hL * fy};
}

// ----------------------------------------------------------------------------
// Compute and distribute the impact of a 3D force applied to a given cell face
// such that the full impact of the force is equal to the sum of the forces
// applied to its corner points.  (See the documentation of `assemble_mech_system_3D`
// in `vem.hpp` for an explanation of the function arguments).
// The result is  written directly into the coresponding global positions in the
// right-hand-side vector (`b_global`).
void
compute_applied_forces_3D(const double* const points,
                          const int* const num_face_corners,
                          const int* const face_corners,
                          const int num_neumann_faces,
                          const int* neumann_faces,
                          const double* const neumann_forces,
                          double* b_global,
                          const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    for (int nf = 0; nf != num_neumann_faces; ++nf) {
        const int nfc = num_face_corners[neumann_faces[nf]];
        const int ix_start = accumulate(num_face_corners, num_face_corners + neumann_faces[nf], 0);
        const auto face_corners_loc = pick_points<3>(points, &face_corners[ix_start], nfc);
        const auto tris = tessellate_face<3>(&face_corners_loc[0], nfc, false, nullptr,
                                             corner_mask ? corner_mask + ix_start : nullptr);
        for (int c = 0; c != 2 * nfc; ++c) { // two relevant triangles per corner
            const double area = triarea<3>(&tris[c][0], &tris[c][3], &tris[c][6]);
            const int corner = (c / 2 + c % 2) % nfc;
            for (int d = 0; d != 3; ++d)
                b_global[3 * face_corners[ix_start + corner] + d] += neumann_forces[3 * nf + d] * area;
        }
    }
}

// ----------------------------------------------------------------------------
// Compute and distribute the impact of a 3D force applied to a given 3D element
// volume among its corner nodes, such that the full impact of the force
// is equal to the sum of the forces applied at the nodes.
// To avoid potential multiple costly evaluations, the the computation requires
// that the centroid of the element is provided by the user.
// (See the documentation of `assemble_mech_system_3D` in `vem.hpp` for an
// explanation of the function arguments).
// The result is written directly into the corresponding global positions in the
// right-had side vector (`b_global`).
void
compute_bodyforce_3D(const double* const points,
                     const int* const face_corners,
                     const int* const num_face_corners,
                     const int num_faces,
                     const double* const centroid,
                     const double* const bforce, // 3D body force
                     double* b_global, // global right-hand side
                     const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    int fcorner_start = 0;
    for (int f = 0; f != num_faces; ++f) {
        const int nfc = num_face_corners[f];
        const auto tris = tessellate_face<3>(pick_points<3>(points, face_corners + fcorner_start, nfc), false,
                                             nullptr, corner_mask ? corner_mask + fcorner_start : nullptr);
        assert(int(tris.size()) == nfc * 2);
        for (int c = 0; c != nfc; ++c) {
            // the following two triangles in the tesselation pertain to corner 'c'.
            const auto t1 = tris[2 * c];
            const auto t2 = (c == 0) ? tris.back() : tris[2 * c - 1];
            const double vol1 = tetrahedron_volume(&t1[0], &t1[3], &t1[6], &centroid[0]);
            const double vol2 = tetrahedron_volume(&t2[0], &t2[3], &t2[6], &centroid[0]);
            for (int d = 0; d != 3; ++d) {
                b_global[3 * face_corners[fcorner_start + c] + d] += (vol1 + vol2) * bforce[d];
            }
        }
        fcorner_start += nfc;
    }
}

// ----------------------------------------------------------------------------
// Eliminate degrees of freedom associated with Dirichlet boundary conditions.
// The original system (A and b) will be modified/reduced accordingly.
void
reduce_system(vector<tuple<int, int, double>>& A,
              vector<double>& b,
              const int num_fixed_dofs,
              const int* const fixed_dof_ixs,
              const double* const fixed_dof_values)
// ----------------------------------------------------------------------------
{
    // check input
    if (!is_sorted(fixed_dof_ixs, fixed_dof_ixs + num_fixed_dofs))
        throw invalid_argument("The indices of fixed degrees of freedom must be "
                               "provided in ascending order.");

    // sort entries of A into columns

    sort(A.begin(), A.end(), [](const auto& aa, const auto& bb) { return get<1>(aa) < get<1>(bb); });

    // eliminate columns associated with fixed dofs (move value to b)
    auto cur_end = A.begin();
    for (int i = 0; i != num_fixed_dofs; ++i) {
        const int dof = fixed_dof_ixs[i];
        // find the range of the elements whose column index equal 'dof'
        auto start = lower_bound(
            cur_end, A.end(), dof, [](const auto& a, int value) { return get<1>(a) < value; });
        if (start != A.end() && get<1>(*start) == dof) {
            cur_end = find_if(start, A.end(), [dof](const auto& a) { return get<1>(a) != dof; });
            for (; start != cur_end; ++start)
                b[get<0>(*start)] -= get<2>(*start) * fixed_dof_values[i];
        }
    }

    // determine renumbering
    vector<int> keep(b.size());
    iota(keep.begin(), keep.end(), 0);
    vector<int> renum;
    renum.reserve(b.size());

    set_difference(
        keep.begin(), keep.end(), fixed_dof_ixs, fixed_dof_ixs + num_fixed_dofs, back_inserter(renum));

    const int discard_flag = b.size() + 1;
    vector<int> renum_inv(b.size(), discard_flag);
    for (int i = 0; i != int(renum.size()); ++i)
        renum_inv[renum[i]] = i;

    // eliminate rows associated with fixed dofs and renumber
    A.erase(remove_if(A.begin(),
                      A.end(),
                      [&discard_flag, &renum_inv](const tuple<int, int, double>& el) {
                          return renum_inv[get<0>(el)] == discard_flag
                              || renum_inv[get<1>(el)] == discard_flag;
                      }),
            A.end());
    transform(A.begin(), A.end(), A.begin(), [&renum_inv](const tuple<int, int, double>& el) {
        return tuple<int, int, double> {renum_inv[get<0>(el)], renum_inv[get<1>(el)], get<2>(el)};
    });

    for (int i = 0; i != int(renum.size()); ++i)
        b[i] = b[renum[i]];

    b.resize(renum.size());
}

void
set_boundary_conditions(vector<tuple<int, int, double>>& A,
                        vector<double>& b,
                        const int num_fixed_dofs,
                        const int* const fixed_dof_ixs,
                        const double* const fixed_dof_values)
// ----------------------------------------------------------------------------
{
    // check input
    if (!is_sorted(fixed_dof_ixs, fixed_dof_ixs + num_fixed_dofs))
        throw invalid_argument("The indices of fixed degrees of freedom must be "
                               "provided in ascending order.");

    // sort entries of A into columns

    sort(A.begin(), A.end(), [](const auto& aa, const auto& bb) { return get<1>(aa) < get<1>(bb); });

    // eliminate columns associated with fixed dofs (move value to b)
    // not changing the matrix, just the right hand side??
    auto cur_end = A.begin();
    for (int i = 0; i != num_fixed_dofs; ++i) {
        const int dof = fixed_dof_ixs[i];
        // find the range of the elements whose column index equal 'dof'
        auto start = lower_bound(
            cur_end, A.end(), dof, [](const auto& a, int value) { return get<1>(a) < value; });
        if (start != A.end() && get<1>(*start) == dof) {
            cur_end = find_if(start, A.end(), [dof](const auto& a) { return get<1>(a) != dof; });
            for (; start != cur_end; ++start) {
                b[get<0>(*start)] -= get<2>(*start) * fixed_dof_values[i];
                get<2>(*start) = 0; // set matrix element to zero
            }
        }
    }


    // determine renumbering
    vector<int> keep(b.size());
    iota(keep.begin(), keep.end(), 0);
    vector<int> renum;
    renum.reserve(b.size());

    set_difference(
        keep.begin(), keep.end(), fixed_dof_ixs, fixed_dof_ixs + num_fixed_dofs, back_inserter(renum));

    const int discard_flag = b.size() + 1;
    vector<int> renum_inv(b.size(), discard_flag);
    for (int i = 0; i != int(renum.size()); ++i)
        renum_inv[renum[i]] = i;

    // eliminate rows associated with fixed dofs and renumber
    for (auto& el : A) {
        bool is_diag = get<0>(el) == get<1>(el);
        bool is_dof = (renum_inv[get<0>(el)] == discard_flag) || (renum_inv[get<1>(el)] == discard_flag);
        if (is_dof) {
            if (is_diag) {
                get<2>(el) = 1;
            } else {
                get<2>(el) = 0;
            }
        }
    }
    for (int i = 0; i != num_fixed_dofs; ++i) {
        const int dof = fixed_dof_ixs[i];
        b[dof] = fixed_dof_values[i];
    }

}

// ----------------------------------------------------------------------------
// Create a local indexing of the subset of nodes that are referenced by the
// given faces (which refer to the global indexing of the nodes).
// Write this local indexing into the `indexing` output argument, and return a
// map from global to local indexing.
map<int, int>
global_to_local_indexing(const int* const faces,
                         const int* const num_face_edges,
                         const int num_faces,
                         vector<int>& indexing)
// ----------------------------------------------------------------------------
{
    // determine number of references to vertices
    const int faces_len = accumulate(num_face_edges, num_face_edges + num_faces, 0);

    // by creating a set, we keep only unique (global) node indices, sorted.
    const set<int> unique_elements(faces, faces + faces_len);

    // copying the sorted list of unique indices to the 'indexing' variable,
    // which is an output variable of this function
    indexing = vector<int>(unique_elements.begin(), unique_elements.end());

    // return a map to facilitate global-to-local indexing
    map<int, int> reindex;
    for (int i = 0; i != int(indexing.size()); ++i)
        reindex[indexing[i]] = i;

    return reindex;
}


// ----------------------------------------------------------------------------
// Compute the q-values that are related to virtual basis function integrals
// over element faces in 2D.  These are used as intermediary values when
// assembling the 2D stiffness matrix for a given element.  See (Gain, 2014)
// DOI:10.1016/j.cma.2014.05.005
// It is expected that the corners are listed in anti-clockwise order around the
// element.
vector<double>
compute_q_2D(const double* const corners, const int num_corners)
// ----------------------------------------------------------------------------
{
    vector<double> result(num_corners * 2, 0); // matrix dimension: num_corners x 2

    // (half of) the common factor outside the integral: 1 / (2 |E|)
    const double fac = 1.0 / (4 * element_volume_2D(corners, num_corners));

    // loop over faces (edges) and compute \phi integrated over the edge, and add
    // it to the q-functions associated with the indicent nodes
    for (int i = 0, inext = 1; i != num_corners; ++i, inext = (i + 1) % num_corners) {
        const array<double, 2> scaled_normal {(corners[2 * inext + 1] - corners[2 * i + 1]),
                                              -(corners[2 * inext + 0] - corners[2 * i + 0])};

        for (int d = 0; d != 2; ++d) {
            result[i * 2 + d] += fac * scaled_normal[d];
            result[inext * 2 + d] += fac * scaled_normal[d];
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
// Vector areas of the fan triangles (vertex mean, p_i, p_i+1) of a face, oriented
// along the given outward normal; 3*m values.
vector<double>
fan_vector_areas(const double* const corners,
                 const int* const fcorners,
                 const int m,
                 const double* const normal,
                 const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    const array<double, 3> c = fan_centre(corners, fcorners, m, corner_mask);
    vector<double> res(3 * m);
    double dot = 0;
    for (int i = 0; i != m; ++i) {
        const double* p = &corners[3 * fcorners[i]];
        const double* q = &corners[3 * fcorners[(i + 1) % m]];
        const array<double, 3> u {p[0] - c[0], p[1] - c[1], p[2] - c[2]};
        const array<double, 3> v {q[0] - c[0], q[1] - c[1], q[2] - c[2]};
        res[3 * i + 0] = 0.5 * (u[1] * v[2] - u[2] * v[1]);
        res[3 * i + 1] = 0.5 * (u[2] * v[0] - u[0] * v[2]);
        res[3 * i + 2] = 0.5 * (u[0] * v[1] - u[1] * v[0]);
        for (int d = 0; d != 3; ++d)
            dot += res[3 * i + d] * normal[d];
    }
    if (dot < 0)
        for (auto& x : res)
            x = -x;
    return res;
}

// ----------------------------------------------------------------------------
// Volume enclosed by the fan triangulation of the faces (divergence theorem).
double
fan_volume(const double* const corners,
           const int num_corners,
           const int* const faces,
           const int* const num_face_edges,
           const int num_faces,
           const vector<double>& outward_normals,
           const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    const array<double, 3> o = point_average<3>(corners, num_corners);
    double vol = 0;
    for (int f = 0, pos = 0; f != num_faces; pos += num_face_edges[f++]) {
        const int m = num_face_edges[f];
        const char* const fmask = corner_mask && has_hanging_nodes(corner_mask + pos, m) ? corner_mask + pos : nullptr;
        const auto tris = fan_vector_areas(corners, &faces[pos], m, &outward_normals[3 * f], fmask);
        const array<double, 3> c = fan_centre(corners, &faces[pos], m, fmask);
        for (int i = 0; i != m; ++i) {
            const double* p = &corners[3 * faces[pos + i]];
            const double* q = &corners[3 * faces[pos + (i + 1) % m]];
            for (int d = 0; d != 3; ++d)
                vol += tris[3 * i + d] * (c[d] + p[d] + q[d] - 3 * o[d]) / 9.0;
        }
    }
    return vol;
}

// ----------------------------------------------------------------------------
// Each node's share of the cell boundary, int_dE phi_i, normalised to mean one.
vector<double>
compute_boundary_weights_3D(const double* const corners,
                            const int num_corners,
                            const int* const faces,
                            const int* const num_face_edges,
                            const int num_faces)
// ----------------------------------------------------------------------------
{
    vector<double> result(num_corners, 0);
    int cur_vertex_ix = 0;
    for (int f = 0; f != num_faces; ++f) {
        const auto facecorners = pick_points<3>(corners, &faces[cur_vertex_ix], num_face_edges[f]);
        vector<double> cvals(num_face_edges[f], 0);
        for (int e = 0; e != num_face_edges[f]; ++e, ++cur_vertex_ix) {
            fill(cvals.begin(), cvals.end(), 0);
            cvals[e] = 1;
            result[faces[cur_vertex_ix]] += face_integral(&facecorners[0], num_face_edges[f], 3, &cvals[0]);
        }
    }
    const double mean = accumulate(result.begin(), result.end(), 0.0) / num_corners;
    for (auto& w : result)
        w /= mean;
    return result;
}

// ----------------------------------------------------------------------------
// Compute the q-values that are related to virtual basis function integrals
// over element faces in 3D.  These are used as intermediary values when
// assembling the 2D stiffness matrix for a given element.  See (Gain, 2014)
// DOI:10.1016/j.cma.2014.05.005
vector<double>
compute_q_3D(const double* const corners,
             const int num_corners,
             const int* const faces,
             const int* const num_face_edges,
             const int num_faces,
             const double volume,
             const vector<double>& outward_normals,
             const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    vector<double> result(num_corners * 3, 0); // the vector of q

    // common factor outside the integral: 1 / (2 |E|)
    const double fac = 1.0 / (2 * volume);

    if (face_quadrature_ == FaceQuadrature::TRIANGULATED) {
        // P1 on the fan (corner mean, p_i, p_i+1): int_tri phi n = A n / 3 per triangle vertex,
        // the centre carrying 1/nc of each face corner
        int pos = 0;
        for (int f = 0; f != num_faces; pos += num_face_edges[f++]) {
            const int m = num_face_edges[f];
            const char* const fmask = corner_mask && has_hanging_nodes(corner_mask + pos, m) ? corner_mask + pos : nullptr;
            const auto tris = fan_vector_areas(corners, &faces[pos], m, &outward_normals[3 * f], fmask);
            int nc = 0;
            for (int k = 0; k != m; ++k)
                nc += !fmask || fmask[k];
            for (int i = 0; i != m; ++i) {
                const int a = faces[pos + i], b = faces[pos + (i + 1) % m];
                for (int d = 0; d != 3; ++d) {
                    const double w = fac * tris[3 * i + d] / 3.0;
                    result[3 * a + d] += w;
                    result[3 * b + d] += w;
                    for (int k = 0; k != m; ++k)
                        if (!fmask || fmask[k])
                            result[3 * faces[pos + k] + d] += w / nc;
                }
            }
        }
        return result;
    }

    // loop over faces, an accumulate q-values for the associated vertices
    int cur_vertex_ix = 0;
    for (int f = 0; f != num_faces; ++f) {
        const auto facecorners = pick_points<3>(corners, &faces[cur_vertex_ix], num_face_edges[f]);
        const double* const normal = &outward_normals[3 * f];
        vector<double> cvals(num_face_edges[f], 0);
        for (int e = 0; e != num_face_edges[f]; ++e, ++cur_vertex_ix) {
            fill(cvals.begin(), cvals.end(), 0);
            cvals[e] = 1; // indicate which basis function we are integrating

            const double phi = face_integral(&facecorners[0], num_face_edges[f], 3, &cvals[0],
                                             corner_mask ? corner_mask + cur_vertex_ix - e : nullptr);
            for (int d = 0; d != 3; ++d) {
                result[3 * faces[cur_vertex_ix] + d] += fac * phi * normal[d];
            }
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Nr matrix for a given 2D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Nr_2D(const double* const corners, const int num_corners)

// ----------------------------------------------------------------------------
{
    const array<double, 2> midpoint = point_average<2>(corners, num_corners);

    vector<double> result;
    result.reserve(6 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_2D(1, corners[2 * i + 1] - midpoint[1], 1, -(corners[2 * i] - midpoint[0])));
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Nr matrix for a given 3D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Nr_3D(const double* const corners, const int num_corners)
// ----------------------------------------------------------------------------
{
    const array<double, 3> midpoint = point_average<3>(corners, num_corners);
    vector<double> result;
    result.reserve(18 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_3D(1,
                           corners[3 * i + 1] - midpoint[1],
                           -(corners[3 * i + 2] - midpoint[2]),
                           1,
                           -(corners[3 * i] - midpoint[0]),
                           corners[3 * i + 2] - midpoint[2],
                           1,
                           -(corners[3 * i + 1] - midpoint[1]),
                           corners[3 * i] - midpoint[0]));

    return result;
}

// ----------------------------------------------------------------------------
// Computes the Nc matrix for a given 3D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Nc_3D(const double* const corners, const int num_corners)
// ----------------------------------------------------------------------------
{
    const array<double, 3> mid = point_average<3>(corners, num_corners);
    vector<double> result;
    result.reserve(18 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_3D(corners[3 * i] - mid[0],
                           corners[3 * i + 1] - mid[1],
                           corners[3 * i + 2] - mid[2],
                           corners[3 * i + 1] - mid[1],
                           corners[3 * i] - mid[0],
                           corners[3 * i + 2] - mid[2],
                           corners[3 * i + 2] - mid[2],
                           corners[3 * i + 1] - mid[1],
                           corners[3 * i] - mid[0]));
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Nc matrix for a given 2D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Nc_2D(const double* const corners, const int num_corners)
// ----------------------------------------------------------------------------
{
    const array<double, 2> midpoint = point_average<2>(corners, num_corners);
    vector<double> result;
    result.reserve(6 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_2D(corners[2 * i] - midpoint[0],
                           corners[2 * i + 1] - midpoint[1],
                           corners[2 * i + 1] - midpoint[1],
                           corners[2 * i] - midpoint[0]));
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Wc matrix for a given 3D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Wr_3D(const vector<double>& q)
// ----------------------------------------------------------------------------
{
    const int num_corners = q.size() / 3;
    const double ncinv = 1.0 / num_corners;
    assert(num_corners * 3 == int(q.size()));
    vector<double> result;
    result.reserve(18 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_3D(ncinv,
                           q[3 * i + 1],
                           -q[3 * i + 2],
                           ncinv,
                           -q[3 * i],
                           q[3 * i + 2],
                           ncinv,
                           -q[3 * i + 1],
                           q[3 * i]));
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Wr matrix for a given 2D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Wr_2D(const vector<double>& q)
// ----------------------------------------------------------------------------
{
    const int num_corners = q.size() / 2;
    assert(num_corners * 2 == int(q.size())); // should have even number of elements

    vector<double> result;
    result.reserve(6 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result, matentry_2D(1.0 / num_corners, q[2 * i + 1], 1.0 / num_corners, -q[2 * i]));

    return result;
}

// ----------------------------------------------------------------------------
// Computes the Wc matrix for a given 3D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Wc_3D(const vector<double>& q)
// ----------------------------------------------------------------------------
{
    const int num_corners = q.size() / 3;
    assert(num_corners * 3 == int(q.size()));
    vector<double> result;
    result.reserve(18 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result,
               matentry_3D(2 * q[3 * i],
                           q[3 * i + 1],
                           q[3 * i + 2],
                           2 * q[3 * i + 1],
                           q[3 * i],
                           q[3 * i + 2],
                           2 * q[3 * i + 2],
                           q[3 * i + 1],
                           q[3 * i]));
    return result;
}

// ----------------------------------------------------------------------------
// Computes the Wc matrix for a given 2D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
vector<double>
compute_Wc_2D(const vector<double>& q)
// ----------------------------------------------------------------------------
{
    const int num_corners = q.size() / 2;
    assert(num_corners * 2 == int(q.size())); // should have even number of elements

    vector<double> result;
    result.reserve(6 * num_corners);

    for (int i = 0; i != num_corners; ++i)
        append(result, matentry_2D(2 * q[2 * i], q[2 * i + 1], 2 * q[2 * i + 1], q[2 * i]));

    return result;
}

// ----------------------------------------------------------------------------
// Computes the D matrix for a given 2D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
// The D matrix is obtained from the 3x3 elasticity tensor matrix (in Voigt)
// notation by multiplying the last row and last column by 2 (so the lower
// right element ends up being multiplied by 4 in total).
vector<double>
compute_D_2D(const double young, const double poisson)
// ----------------------------------------------------------------------------
{
    const double fac = young / (1 + poisson) / (1 - 2 * poisson);

    // Note that element (3,3) is four times higher than what may be sometimes
    // be seen in literature.   This is because of using half the value for shear
    // elements when using voigt notation.
    vector<double> result {
        1 - poisson, poisson, 0, poisson, 1 - poisson, 0, 0, 0, 2 * (1 - 2 * poisson)};

    // multiply by the constant factor
    for_each(result.begin(), result.end(), [fac](double& d) { d *= fac; });

    return result;
}

// ----------------------------------------------------------------------------
// Computes the D matrix for a given 3D element
// (See Gain 2014, DOI:10.1016/j.cma.2014.05.005)
// The D matrix is obtained from the 6x6 elasticity tensor matrix (in Voigt)
// notation by multiplying the last three rows and last three columns by 2
// (so the lower right 3x3 element block ends up being multiplied by 4 in total).
vector<double>
compute_D_3D(const double young, const double poisson)
// ----------------------------------------------------------------------------
{
    const double fac = young / (1 + poisson) / (1 - 2 * poisson);

    // Note that element (3,3) is four times higher than what may be sometimes
    // be seen in literature.   This is because of using half the value for shear
    // elements when using voigt notation.
    vector<double> result {1 - poisson,           poisson, poisson, 0, 0, 0,       poisson,
                           1 - poisson,           poisson, 0,       0, 0, poisson, poisson,
                           1 - poisson,           0,       0,       0, 0, 0,       0,
                           2 * (1 - 2 * poisson), 0,       0,       0, 0, 0,       0,
                           2 * (1 - 2 * poisson), 0,       0,       0, 0, 0,       0,
                           2 * (1 - 2 * poisson)};

    // multiply by the constant factor
    for_each(result.begin(), result.end(), [fac](double& d) { d *= fac; });

    return result;
}

// ----------------------------------------------------------------------------
// Compute the matrix corresponding to (I - P), as described in
// (Gain, 2014) DOI:10.1016/j.cma.2014.05.005.
vector<double>
compute_ImP(const vector<double>& Nr,
            const vector<double>& Nc,
            const vector<double>& Wr,
            const vector<double>& Wc,
            const int dim)
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);

    const int N = (dim == 2) ? Nc.size() / 6 : Nc.size() / 18;
    const int r = dim * N;
    const int c = (dim == 2) ? 3 : 6;

    const auto Pr = matmul(&Nr[0], r, c, false, &Wr[0], r, c, true);
    const auto Pc = matmul(&Nc[0], r, c, false, &Wc[0], r, c, true);

    // We want to compute I-P, which equals I-(Pr+Pc)

    vector<double> result = identity_matrix(1, N * dim);

    assert(result.size() == Pr.size());

    for (int i = 0; i != int(result.size()); ++i)
        result[i] -= (Pr[i] + Pc[i]);

    return result;
}

// ----------------------------------------------------------------------------
bool
is_behind_face(const array<double, 3>& point,
               const double* const face_normal,
               const double* const face_centroid)
// ----------------------------------------------------------------------------
{
    const double tol = 1e-13;
    const array<double, 3> v {
        face_centroid[0] - point[0], face_centroid[1] - point[1], face_centroid[2] - point[2]};
    return (v[0] * face_normal[0] + v[1] * face_normal[1] + v[2] * face_normal[2]) + tol > 0;
}


// ----------------------------------------------------------------------------
bool
is_star_point(const array<double, 3>& point,
              const vector<double>& face_normals,
              const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
    const int N = (int)face_normals.size() / 3;
    for (int i = 0; i != N; ++i)
        if (!is_behind_face(point, &face_normals[3 * i], &face_centroids[3 * i]))
            return false;
    return true;
}

// ----------------------------------------------------------------------------
// NB, the face normals must be _unit_ normals (not scaled by face area)
// in order for this algorithm to work as intended.
//
// Assumes planar faces; warped faces may need tessellating into triangles first.
array<double, 3>
identify_star_point(const array<double, 3>& point,
                    const vector<double>& face_normals,
                    const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
    const int N = (int)face_normals.size() / 3;
    const int max_iter = 20 * N; // presumably enough for most realistic cases?

    array<double, 3> result(point);
    int count = 0;
    for (int i = 0; i != max_iter; ++i) {
        const int f_ix = i % N;
        if (!is_behind_face(result, &face_normals[3 * f_ix], &face_centroids[3 * f_ix])) {
            count = 0;
            //  projecting the point onto the plane defined by the normal and the centroid
            const double proj = (result[0] - face_centroids[3 * f_ix]) * face_normals[3 * f_ix]
                + (result[1] - face_centroids[3 * f_ix + 1]) * face_normals[3 * f_ix + 1]
                + (result[2] - face_centroids[3 * f_ix + 2]) * face_normals[3 * f_ix + 2];
            for (int ii = 0; ii != 3; ++ii)
                result[ii] -= 1.1 * proj * face_normals[3 * f_ix + ii]; // move slightly past plane
        }
        if (++count == N)
            break;
    }
    if (count != N) {
        std::ostringstream os;
        os << "Unable to find a star point for cell. faces=" << N
           << " start=(" << point[0] << ',' << point[1] << ',' << point[2]
           << ") last=(" << result[0] << ',' << result[1] << ',' << result[2] << ')';
        throw runtime_error(os.str());
    }

    return result;
}

// ----------------------------------------------------------------------------
// Point deepest behind all face planes (Chebyshev centre): max t subject to
// n_f.(x - c_f) + t <= 0, solved exactly by enumerating the LP's vertices.
array<double, 3>
chebyshev_star_point(const vector<double>& face_normals,
                     const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
    const int N = (int)face_normals.size() / 3;
    array<double, 3> lo {face_centroids[0], face_centroids[1], face_centroids[2]};
    array<double, 3> hi = lo;
    for (int f = 1; f != N; ++f)
        for (int d = 0; d != 3; ++d) {
            lo[d] = min(lo[d], face_centroids[3 * f + d]);
            hi[d] = max(hi[d], face_centroids[3 * f + d]);
        }
    const array<double, 3> mid {(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2};
    const double scale = max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-300});

    // rows a.(x - mid, t) <= b in units of scale; faces first, then a bounding box
    vector<array<double, 5>> rows;
    for (int f = 0; f != N; ++f) {
        const double* n = &face_normals[3 * f];
        double b = 0;
        for (int d = 0; d != 3; ++d)
            b += n[d] * (face_centroids[3 * f + d] - mid[d]) / scale;
        rows.push_back({n[0], n[1], n[2], 1.0, b});
    }
    for (int d = 0; d != 3; ++d) {
        array<double, 5> r {0, 0, 0, 0, 1.0};
        r[d] = 1;
        rows.push_back(r);
        r[d] = -1;
        rows.push_back(r);
    }

    const int M = (int)rows.size();
    double best_t = -numeric_limits<double>::infinity();
    array<double, 4> best {0, 0, 0, 0};
    for (int i = 0; i < M; ++i)
        for (int j = i + 1; j < M; ++j)
            for (int k = j + 1; k < M; ++k)
                for (int l = k + 1; l < M; ++l) {
                    // solve the 4x4 system of the active rows by partial pivoting
                    double A[4][5];
                    const int act[4] = {i, j, k, l};
                    for (int r = 0; r != 4; ++r)
                        for (int c = 0; c != 5; ++c)
                            A[r][c] = rows[act[r]][c];
                    bool singular = false;
                    for (int c = 0; c != 4 && !singular; ++c) {
                        int p = c;
                        for (int r = c + 1; r != 4; ++r)
                            if (abs(A[r][c]) > abs(A[p][c]))
                                p = r;
                        if (abs(A[p][c]) < 1e-12) {
                            singular = true;
                            break;
                        }
                        swap(A[c], A[p]);
                        for (int r = 0; r != 4; ++r)
                            if (r != c) {
                                const double fac = A[r][c] / A[c][c];
                                for (int cc = c; cc != 5; ++cc)
                                    A[r][cc] -= fac * A[c][cc];
                            }
                    }
                    if (singular)
                        continue;
                    array<double, 4> x;
                    for (int r = 0; r != 4; ++r)
                        x[r] = A[r][4] / A[r][r];
                    if (x[3] <= best_t)
                        continue;
                    bool feasible = true;
                    for (const auto& row : rows)
                        if (row[0] * x[0] + row[1] * x[1] + row[2] * x[2] + row[3] * x[3] > row[4] + 1e-10) {
                            feasible = false;
                            break;
                        }
                    if (feasible) {
                        best_t = x[3];
                        best = x;
                    }
                }

    const array<double, 3> result {mid[0] + scale * best[0], mid[1] + scale * best[1], mid[2] + scale * best[2]};
    if (!(best_t > 1e-9) || !is_star_point(result, face_normals, face_centroids)) {
        std::ostringstream os;
        os << "Cell is not star-shaped: faces=" << N << " depth=" << best_t * scale;
        throw runtime_error(os.str());
    }
    return result;
}

// ----------------------------------------------------------------------------
void
compute_face_geometry(const vector<double>& points, double* normal, double* centroid,
                      const char* const corner_mask = nullptr)
// ----------------------------------------------------------------------------
{
    // compute tessellation and face centroid
    const auto face_tessellation
        = tessellate_face<3>(&points[0], (int)points.size() / 3, false, centroid, corner_mask);

    // compute (normalized) face normal
    fill(normal, normal + 3, 0);

    for (auto tri : face_tessellation) {
        const auto tn = trinormal(&tri[0], &tri[3], &tri[6]);
        for (int d = 0; d != 3; ++d)
            normal[d] += tn[d];
    }

    // normalize
    const double n = norm<3>(&normal[0]);
    for (int d = 0; d != 3; ++d)
        normal[d] /= n;
}

// ----------------------------------------------------------------------------
vector<int>
extract_local_faces(const int* const faces, const int* const num_face_edges, const int num_faces)
// ----------------------------------------------------------------------------
{
    // Picks the locally relevant face corner indices
    vector<int> result(faces, faces + accumulate(num_face_edges, num_face_edges + num_faces, 0));
    return result;
}

// ----------------------------------------------------------------------------
bool
inward_pointing_normals(const vector<double>& normals, const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
    // Allows for non-planar faces; may still fail in pathological cases.

    // compute the "centroid" of the polyhedron
    const int N = (int)normals.size() / 3;
    const array<double, 3> mean_point = point_average<3>(&face_centroids[0], N);

    // compute the distance vectors from the polyhedron centroid to the face centroids
    vector<array<double, 3>> dists(N);
    for (int i = 0; i != N; ++i)
        for (int d = 0; d != 3; ++d)
            dists[i][d] = face_centroids[3 * i + d] - mean_point[d];

    // add up the scalar products between distance vectors and face normals
    double sum = 0;
    for (int i = 0; i != N; ++i)
        for (int d = 0; d != 3; ++d)
            sum += dists[i][d] * normals[3 * i + d];

    return sum < 0;
}



// ----------------------------------------------------------------------------
// Compute key parts of cell geometry, including a consistent set of outward
// normals, the cell volume, centroid and a point for which the cell is
// star-shaped (which may or may not be the centroid)
void
compute_cell_geometry(const double* points,
                      const int num_points,
                      const int* const faces,
                      const int* const num_face_edges,
                      const int num_faces,
                      vector<double>& outward_normals,
                      vector<double>& face_centroids,
                      array<double, 3>& cell_centroid,
                      array<double, 3>& star_point,
                      double& volume,
                      const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    // ensure faces are consistently oriented so that normals all point outwards
    // or inwards
    vector<int> faces2 = extract_local_faces(faces, num_face_edges, num_faces);

    outward_normals.resize(num_faces * 3);
    face_centroids.resize(num_faces * 3);
    for (int f = 0, faces_offset = 0; f != num_faces; faces_offset += num_face_edges[f++])
        compute_face_geometry(pick_points<3>(points, &faces2[faces_offset], num_face_edges[f]),
                              &outward_normals[f * 3],
                              &face_centroids[f * 3],
                              corner_mask ? corner_mask + faces_offset : nullptr);
    assert(face_centroids.size() == static_cast<std::size_t>(num_faces * 3));


    // ensure normals are pointing out of, rather than into, polyhedron
    if (inward_pointing_normals(outward_normals, face_centroids)) {
        Opm::OpmLog::error("Inward pointing normals try to fix them");
        assert(false);
    }

    // identify a star point (usually, mean_point qualifies, but not necessarily)
    bool signed_tets = false;
    try {
        star_point = identify_star_point(point_average<3>(points, num_points),
                                         outward_normals, face_centroids);
    } catch (const runtime_error&) {
        // A cell merged from several others has many faces of very different
        // size, and the average of its corners is pulled towards the finely
        // divided ones, sometimes out of the cell. The face centroids are
        // spread evenly over its surface.
        const int nf = (int)face_centroids.size()/3;
        array<double, 3> mid {0, 0, 0};
        for (int f = 0; f != nf; ++f)
            for (int d = 0; d != 3; ++d)
                mid[d] += face_centroids[3*f + d]/nf;
        try {
            star_point = identify_star_point(mid, outward_normals, face_centroids);
        } catch (const runtime_error&) {
            try {
                // the projection iteration can cycle when the kernel is a thin sliver
                star_point = chebyshev_star_point(outward_normals, face_centroids);
            } catch (const runtime_error&) {
                // Warped faces can leave no point behind every face plane. Any apex
                // gives the exact volume and centroid with signed tetrahedra.
                star_point = mid;
                signed_tets = true;
                static bool warned = false;
                if (!warned) {
                    warned = true;
                    Opm::OpmLog::warning("VEM: cell without a star point (warped faces); "
                                         "using signed tetrahedra from the face-centroid mean");
                }
            }
        }
    }

    // compute cell centroid and volume
    volume = 0;
    cell_centroid = {0, 0, 0};

    // loop over faces, split each face into triangles, and accumulate the volumes
    // defined by these triangles and the inside point
    for (int f = 0, fpos = 0; f != num_faces; fpos += num_face_edges[f++]) {
        for (auto tri : tessellate_face<3>(pick_points<3>(points, &faces[fpos], num_face_edges[f]), true, nullptr,
                                           corner_mask ? corner_mask + fpos : nullptr)) {
            tri.insert(tri.end(), &star_point[0], &star_point[0] + 3);
            double tvol = tetrahedron_volume(&tri[0], &tri[3], &tri[6], &tri[9]);
            if (signed_tets) {
                // Positive when the apex is behind the face, by the outward normal.
                const double* n = &outward_normals[3 * f];
                double side = 0.0;
                for (int d = 0; d != 3; ++d)
                    side += n[d] * ((tri[d] + tri[3 + d] + tri[6 + d]) / 3.0 - star_point[d]);
                if (side < 0)
                    tvol = -tvol;
            }
            const auto tet_centroid = point_average<3>(&tri[0], 4);
            volume += tvol;
            for (int d = 0; d != 3; ++d)
                cell_centroid[d] += tvol * tet_centroid[d];
        }
    }

    for (int d = 0; d != 3; ++d)
        cell_centroid[d] /= volume;

    // the strain projection divides by the volume its own face triangulation encloses
    if (face_quadrature_ == FaceQuadrature::TRIANGULATED)
        volume = fan_volume(points, num_points, faces, num_face_edges, num_faces, outward_normals, corner_mask);

    // if cell_centroid is also a star_point, use that as star_point instead
    if (is_star_point(cell_centroid, outward_normals, face_centroids))
        star_point = cell_centroid;
}




}; // end anonymous namespace

// ============================================================================
namespace vem
{
// ============================================================================

// ----------------------------------------------------------------------------
void
potential_gradient_force_3D(const double* const points,
                            const int num_cells,
                            const int* const num_cell_faces, // cell faces per cell
                            const int* const num_face_corners, // corners per cellface
                            const int* const face_corners,
                            const double* const field,
                            vector<double>& fgrad,
                            vector<tuple<int, int, double>>& div,
                            bool get_matrix,
                            const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    // initializing result vector
    const int tot_num_faces = accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    const int tot_num_fcorners = accumulate(num_face_corners, num_face_corners + tot_num_faces, 0);
    const int tot_num_nodes = *max_element(face_corners, face_corners + tot_num_fcorners) + 1;
    fgrad = vector<double>(3 * tot_num_nodes, 0.0);

    // The contribution from a given cell to the gradient of pressure evaluated at
    // one of its corners can be obtained from vem by considering the first three componnets
    // of `volume * P * Wc`, where `P = [p p p 0 0 0]`.  Therefore, only the three first columns
    // of Wc matter, which is `diag([2*q1, 2*q2, 2*q3])`.  Moreover, `q_i = 1/(2*volume) z_i`,
    // where z_i is the face integral of the phi_i basis function of the node in question.  As such,
    // we can compute the contribution of this cell to the gradient of pressure in the given node
    // directly as `p * z_i`, and bypass the volume computation entirely.

    // Below, we compute the contribution from all corners of all faces to the global 'fgrad' vector.

    vector<int> indexing;
    int cur_fcor_start = 0;
    int cur_cface_start = 0;
    for (int cell = 0; cell != num_cells; ++cell) {
        const auto reindex = global_to_local_indexing(&face_corners[cur_fcor_start],
                                                      &num_face_corners[cur_cface_start],
                                                      num_cell_faces[cell],
                                                      indexing);
        const int num_corners = int(indexing.size());
        const auto corners_loc = pick_points<3>(points, &indexing[0], num_corners);

        const int tot_num_cellface_corners
            = accumulate(&num_face_corners[cur_cface_start],
                         &num_face_corners[cur_cface_start] + num_cell_faces[cell],
                         0);

        vector<int> faces_loc(tot_num_cellface_corners);
        transform(&face_corners[cur_fcor_start],
                  &face_corners[cur_fcor_start] + tot_num_cellface_corners,
                  faces_loc.begin(),
                  [&reindex](const int I) { return reindex.find(I)->second; });

        double volume; // to be computed below
        vector<double> outward_normals, face_centroids; // to be computed below
        array<double, 3> star_point, cell_centroid; // to be computed below
        compute_cell_geometry(&corners_loc[0],
                              num_corners,
                              &faces_loc[0],
                              &num_face_corners[cur_cface_start],
                              num_cell_faces[cell],
                              outward_normals,
                              face_centroids,
                              cell_centroid,
                              star_point,
                              volume,
                              corner_mask ? corner_mask + cur_fcor_start : nullptr);

        // Since computing q involves dividing by volume, we are really computing q * volume
        // in the below call, by passing the value '1' for volume.  In other words, we are computing
        // (1/2) z_i.
        const auto qv = compute_q_3D(&corners_loc[0],
                                     num_corners,
                                     &faces_loc[0],
                                     &num_face_corners[cur_cface_start],
                                     num_cell_faces[cell],
                                     1,
                                     outward_normals,
                                     corner_mask ? corner_mask + cur_fcor_start : nullptr);

        // fill in entries in global fgrad vector
        for (int c = 0; c != num_corners; ++c) {
            for (int d = 0; d != 3; ++d) {
                fgrad[3 * indexing[c] + d] += 2 * field[cell] * qv[3 * c + d];
                if (get_matrix) {
                    int I = 3 * indexing[c] + d;
                    int J = cell;
                    double val = 2 * qv[3 * c + d];
                    div.push_back(tuple<int, int, double> {I, J, val});
                }
            }
        }

        // keep track of our current position in the `face_corners` array
        cur_fcor_start += tot_num_cellface_corners;
        cur_cface_start += num_cell_faces[cell];
    }
}


// ============================================================================
// ================= Q1 FEM element routines and hex detection ================
// ============================================================================

namespace {

// evaluate the gradients (w.r.t. reference coordinates) of the 8 trilinear shape
// functions of the Dune reference hexahedron at reference point xi
void
hex_shape_gradients(const double* const xi, double grad[8][3])
{
    for (int k = 0; k != 8; ++k) {
        const int c[3] = {k & 1, (k >> 1) & 1, (k >> 2) & 1};
        for (int d = 0; d != 3; ++d) {
            double g = (c[d] == 1) ? 1.0 : -1.0;
            for (int e = 0; e != 3; ++e) {
                if (e == d)
                    continue;
                g *= (c[e] == 1) ? xi[e] : 1.0 - xi[e];
            }
            grad[k][d] = g;
        }
    }
}

// compute physical shape function gradients g[8][3] and the Jacobian determinant of
// the trilinear map at reference point xi; returns detJ (may be <= 0 for invalid cells)
double
hex_physical_gradients(const double* const corners, const double* const xi, double g[8][3])
{
    double gref[8][3];
    hex_shape_gradients(xi, gref);

    // J[i][j] = d x_i / d xi_j
    double J[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int k = 0; k != 8; ++k)
        for (int i = 0; i != 3; ++i)
            for (int j = 0; j != 3; ++j)
                J[i][j] += corners[3 * k + i] * gref[k][j];

    const double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
        - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
        + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    if (detJ == 0.0)
        return 0.0;

    // inverse of J
    double Jinv[3][3];
    Jinv[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / detJ;
    Jinv[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / detJ;
    Jinv[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ;
    Jinv[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / detJ;
    Jinv[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ;
    Jinv[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / detJ;
    Jinv[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ;
    Jinv[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / detJ;
    Jinv[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / detJ;

    // physical gradient: g = J^{-T} gref
    for (int k = 0; k != 8; ++k)
        for (int i = 0; i != 3; ++i)
            g[k][i] = Jinv[0][i] * gref[k][0] + Jinv[1][i] * gref[k][1] + Jinv[2][i] * gref[k][2];

    return detJ;
}

} // anonymous namespace

// ----------------------------------------------------------------------------
void
stiffness_matrix_fem_hex_3D(const double* const corners,
                            const double young,
                            const double poisson,
                            const bool bbar,
                            double* target)
// ----------------------------------------------------------------------------
{
    const double lambda = young * poisson / ((1 + poisson) * (1 - 2 * poisson));
    const double mu = young / (2 * (1 + poisson));

    // 2x2x2 Gauss rule on the unit cube
    const double gp[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
    const double w = 1.0 / 8.0;

    // first pass: physical gradients and integration weights at each Gauss point,
    // and (for B-bar) the volume-averaged dilatational row
    double g[8][8][3]; // [gauss point][node][direction]
    double dV[8];
    double bbar_row[24] = {0};
    double volume = 0.0;
    int q = 0;
    for (int kz = 0; kz != 2; ++kz)
        for (int ky = 0; ky != 2; ++ky)
            for (int kx = 0; kx != 2; ++kx, ++q) {
                const double xi[3] = {gp[kx], gp[ky], gp[kz]};
                const double detJ = hex_physical_gradients(corners, xi, g[q]);
                dV[q] = fabs(detJ) * w;
                volume += dV[q];
                if (bbar)
                    for (int k = 0; k != 8; ++k)
                        for (int d = 0; d != 3; ++d)
                            bbar_row[3 * k + d] += g[q][k][d] * dV[q];
            }
    if (bbar)
        for (int i = 0; i != 24; ++i)
            bbar_row[i] /= volume;

    // second pass: K = sum_q dV_q * B^T C B, with the dilatational part of B replaced
    // by its volume average if B-bar is requested (mean dilatation, Hughes 1980)
    fill(target, target + 24 * 24, 0.0);
    for (q = 0; q != 8; ++q) {
        // build B (6 x 24), Voigt order (xx, yy, zz, xy, yz, zx), engineering shear
        double B[6][24] = {{0}};
        for (int k = 0; k != 8; ++k) {
            const double gx = g[q][k][0], gy = g[q][k][1], gz = g[q][k][2];
            B[0][3 * k + 0] = gx;
            B[1][3 * k + 1] = gy;
            B[2][3 * k + 2] = gz;
            B[3][3 * k + 0] = gy;
            B[3][3 * k + 1] = gx;
            B[4][3 * k + 1] = gz;
            B[4][3 * k + 2] = gy;
            B[5][3 * k + 0] = gz;
            B[5][3 * k + 2] = gx;
        }
        if (bbar) {
            // replace the dilatational part (1/3 of the column gradient in each normal
            // strain row) by its volume average
            for (int k = 0; k != 8; ++k)
                for (int d = 0; d != 3; ++d) {
                    const int col = 3 * k + d;
                    const double corr = (bbar_row[col] - g[q][k][d]) / 3.0;
                    for (int r = 0; r != 3; ++r)
                        B[r][col] += corr;
                }
        }
        // CB = C * B
        double CB[6][24];
        for (int j = 0; j != 24; ++j) {
            const double tr = B[0][j] + B[1][j] + B[2][j];
            CB[0][j] = lambda * tr + 2 * mu * B[0][j];
            CB[1][j] = lambda * tr + 2 * mu * B[1][j];
            CB[2][j] = lambda * tr + 2 * mu * B[2][j];
            CB[3][j] = mu * B[3][j];
            CB[4][j] = mu * B[4][j];
            CB[5][j] = mu * B[5][j];
        }
        for (int i = 0; i != 24; ++i)
            for (int j = 0; j != 24; ++j) {
                double s = 0.0;
                for (int r = 0; r != 6; ++r)
                    s += B[r][i] * CB[r][j];
                target[i * 24 + j] += s * dV[q];
            }
    }
}

// ----------------------------------------------------------------------------
void
stiffness_matrix_fem_quad_2D(const double* const corners,
                             const double young,
                             const double poisson,
                             const bool bbar,
                             double* target)
// ----------------------------------------------------------------------------
{
    // plane strain
    const double lambda = young * poisson / ((1 + poisson) * (1 - 2 * poisson));
    const double mu = young / (2 * (1 + poisson));

    const double gp[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
    const double w = 1.0 / 4.0;

    double g[4][4][2]; // [gauss point][node][direction]
    double dA[4];
    double bbar_row[8] = {0};
    double area = 0.0;
    int q = 0;
    for (int ky = 0; ky != 2; ++ky)
        for (int kx = 0; kx != 2; ++kx, ++q) {
            const double xi[2] = {gp[kx], gp[ky]};
            // bilinear shape gradients on the Dune reference square
            double gref[4][2];
            for (int k = 0; k != 4; ++k) {
                const int c[2] = {k & 1, (k >> 1) & 1};
                gref[k][0] = (c[0] == 1 ? 1.0 : -1.0) * (c[1] == 1 ? xi[1] : 1.0 - xi[1]);
                gref[k][1] = (c[1] == 1 ? 1.0 : -1.0) * (c[0] == 1 ? xi[0] : 1.0 - xi[0]);
            }
            double J[2][2] = {{0, 0}, {0, 0}};
            for (int k = 0; k != 4; ++k)
                for (int i = 0; i != 2; ++i)
                    for (int j = 0; j != 2; ++j)
                        J[i][j] += corners[2 * k + i] * gref[k][j];
            const double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            const double Jinv[2][2]
                = {{J[1][1] / detJ, -J[0][1] / detJ}, {-J[1][0] / detJ, J[0][0] / detJ}};
            for (int k = 0; k != 4; ++k)
                for (int i = 0; i != 2; ++i)
                    g[q][k][i] = Jinv[0][i] * gref[k][0] + Jinv[1][i] * gref[k][1];
            dA[q] = fabs(detJ) * w;
            area += dA[q];
            if (bbar)
                for (int k = 0; k != 4; ++k)
                    for (int d = 0; d != 2; ++d)
                        bbar_row[2 * k + d] += g[q][k][d] * dA[q];
        }
    if (bbar)
        for (int i = 0; i != 8; ++i)
            bbar_row[i] /= area;

    fill(target, target + 64, 0.0);
    for (q = 0; q != 4; ++q) {
        // B (3 x 8), Voigt order (xx, yy, xy), engineering shear
        double B[3][8] = {{0}};
        for (int k = 0; k != 4; ++k) {
            const double gx = g[q][k][0], gy = g[q][k][1];
            B[0][2 * k + 0] = gx;
            B[1][2 * k + 1] = gy;
            B[2][2 * k + 0] = gy;
            B[2][2 * k + 1] = gx;
        }
        if (bbar) {
            // in-plane mean dilatation: replace the dilatational part (1/2 of the
            // column gradient in each normal strain row) by its area average
            for (int k = 0; k != 4; ++k)
                for (int d = 0; d != 2; ++d) {
                    const int col = 2 * k + d;
                    const double corr = (bbar_row[col] - g[q][k][d]) / 2.0;
                    for (int r = 0; r != 2; ++r)
                        B[r][col] += corr;
                }
        }
        double CB[3][8];
        for (int j = 0; j != 8; ++j) {
            const double tr = B[0][j] + B[1][j];
            CB[0][j] = lambda * tr + 2 * mu * B[0][j];
            CB[1][j] = lambda * tr + 2 * mu * B[1][j];
            CB[2][j] = mu * B[2][j];
        }
        for (int i = 0; i != 8; ++i)
            for (int j = 0; j != 8; ++j) {
                double s = 0.0;
                for (int r = 0; r != 3; ++r)
                    s += B[r][i] * CB[r][j];
                target[i * 8 + j] += s * dA[q];
            }
    }
}

// ----------------------------------------------------------------------------
bool
hex_corner_ordering(const double* const corners,
                    const int num_corners,
                    const int* const faces,
                    const int* const num_face_edges,
                    const int num_faces,
                    array<int, 8>& order)
// ----------------------------------------------------------------------------
{
    if (num_corners != 8 || num_faces != 6)
        return false;
    for (int f = 0; f != 6; ++f)
        if (num_face_edges[f] != 4)
            return false;

    // A combinatorial hex needs a closed quad surface: 12 distinct edges, each in exactly 2
    // faces, every vertex of degree 3. Counts of corners/quads/edges alone do not exclude
    // degenerate non-manifold face complexes from collapsed corner-point cells.
    map<pair<int, int>, int> edge_count;
    for (int f = 0, pos = 0; f != 6; pos += num_face_edges[f++])
        for (int i = 0; i != 4; ++i) {
            const int u = faces[pos + i];
            const int v = faces[pos + (i + 1) % 4];
            if (u == v)
                return false; // collapsed face edge
            edge_count[{min(u, v), max(u, v)}] += 1;
        }
    if (int(edge_count.size()) != 12)
        return false;
    vector<int> degree(num_corners, 0);
    set<pair<int, int>> edges;
    for (const auto& [e, cnt] : edge_count) {
        if (cnt != 2)
            return false; // edge not shared by exactly two faces: not a closed surface
        degree[e.first] += 1;
        degree[e.second] += 1;
        edges.insert(e);
    }
    for (int v = 0; v != num_corners; ++v)
        if (degree[v] != 3)
            return false;

    // use the first face as 'bottom'; find for each of its corners the unique
    // neighbor outside the bottom face (the corresponding 'top' corner)
    const int b[4] = {faces[0], faces[1], faces[2], faces[3]};
    const set<int> bottom(b, b + 4);
    int top[4];
    for (int i = 0; i != 4; ++i) {
        int found = -1;
        for (const auto& e : edges) {
            const int other = (e.first == b[i]) ? e.second : (e.second == b[i]) ? e.first : -1;
            if (other >= 0 && bottom.count(other) == 0) {
                if (found >= 0)
                    return false; // more than one neighbor above: not a hex
                found = other;
            }
        }
        if (found < 0)
            return false;
        top[i] = found;
    }
    // the four top corners must be distinct and complete the set of 8
    if (set<int>(top, top + 4).size() != 4)
        return false;

    // Dune reference hex: bottom quad visited cyclically as vertices 0,1,3,2
    order[0] = b[0];
    order[1] = b[1];
    order[3] = b[2];
    order[2] = b[3];
    order[4] = top[0];
    order[5] = top[1];
    order[7] = top[2];
    order[6] = top[3];

    // check the Jacobian sign at the element center; if negative the bottom face was
    // oriented the other way -- mirror the ordering in the x reference direction
    array<double, 24> hc;
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            hc[3 * k + d] = corners[3 * order[k] + d];
    double g[8][3];
    const double center[3] = {0.5, 0.5, 0.5};
    double detJ = hex_physical_gradients(&hc[0], center, g);
    if (detJ < 0) {
        swap(order[0], order[1]);
        swap(order[2], order[3]);
        swap(order[4], order[5]);
        swap(order[6], order[7]);
        for (int k = 0; k != 8; ++k)
            for (int d = 0; d != 3; ++d)
                hc[3 * k + d] = corners[3 * order[k] + d];
    }

    // finally, require a positive Jacobian at all quadrature points *and* at the
    // corners -- badly warped or degenerate cells should be treated with VEM
    // instead.  (Positivity at finitely many points does not strictly guarantee
    // global injectivity of the trilinear map, but Gauss points + corners is the
    // standard practical criterion and covers what the quadrature actually sees.)
    const double gpt[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
    for (int kz = 0; kz != 2; ++kz)
        for (int ky = 0; ky != 2; ++ky)
            for (int kx = 0; kx != 2; ++kx) {
                const double xi[3] = {gpt[kx], gpt[ky], gpt[kz]};
                if (hex_physical_gradients(&hc[0], xi, g) <= 0)
                    return false;
                const double xc[3] = {double(kx), double(ky), double(kz)};
                if (hex_physical_gradients(&hc[0], xc, g) <= 0)
                    return false;
            }
    return true;
}

// ----------------------------------------------------------------------------
void
diagonal_scale_system(vector<tuple<int, int, double>>& A_entries,
                      vector<double>& b,
                      vector<double>& scale)
// ----------------------------------------------------------------------------
{
    const int n = int(b.size());
    vector<double> diag(n, 0.0);
    for (const auto& e : A_entries)
        if (get<0>(e) == get<1>(e))
            diag[get<0>(e)] += get<2>(e);

    scale.assign(n, 1.0);
    for (int i = 0; i != n; ++i)
        if (diag[i] > 0)
            scale[i] = 1.0 / sqrt(diag[i]);

    for (auto& e : A_entries)
        get<2>(e) *= scale[get<0>(e)] * scale[get<1>(e)];
    for (int i = 0; i != n; ++i)
        b[i] *= scale[i];
}

// ----------------------------------------------------------------------------
int
assemble_mech_system_2D(const double* const points,
                        const int num_cells,
                        const int* const num_cell_faces,
                        const int* const cell_corners,
                        const double* const young,
                        const double* const poisson,
                        const double* const body_force, // 2 x number of cells
                        const int num_fixed_dofs, // dirichlet
                        const int* const fixed_dof_ixs,
                        const double* const fixed_dof_values,
                        const int num_neumann_faces, // neumann
                        const int* const neumann_faces,
                        const double* const neumann_forces, // 2 * number of neumann faces
                        vector<tuple<int, int, double>>& A_entries,
                        vector<double>& b,
                        const StabilityChoice stability_choice,
                        bool reduce_boundary)
// ----------------------------------------------------------------------------
{
    // preliminary computations
    const int tot_num_cell_faces = accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    const int num_points = *max_element(cell_corners, cell_corners + tot_num_cell_faces) + 1;


    // assemble full system matrix
    A_entries.clear();
    b.resize(num_points * 2);
    fill(b.begin(), b.end(), 0);

    // loop over cells and assemble system matrix
    vector<double> loc;
    for (int c = 0, corner_ixs = 0; c != num_cells; corner_ixs += num_cell_faces[c++]) {

        // computing local stiffness matrix, and writing its entries into the global
        // system matrix
        const int ncf = num_cell_faces[c];
        assemble_stiffness_matrix_2D(
            points, &cell_corners[corner_ixs], ncf, young[c], poisson[c], stability_choice, loc);

        for (int i = 0; i != 2 * ncf; ++i) {
            for (int j = 0; j != 2 * ncf; ++j) {
                const int I = 2 * cell_corners[corner_ixs + (i / 2)] + i % 2;
                const int J = 2 * cell_corners[corner_ixs + (j / 2)] + j % 2;
                const double val = loc[i * (2 * ncf) + j];
                A_entries.push_back(tuple<int, int, double> {I, J, val});
            }
        }

        // add contribution to right hand side from body forces
        fill(loc.begin(), loc.begin() + 2 * ncf, 0);
        compute_bodyforce_2D(points, &cell_corners[corner_ixs], ncf, &body_force[2 * c], &loc[0]);
        for (int i = 0; i != ncf * 2; ++i)
            b[2 * cell_corners[corner_ixs + (i / 2)] + i % 2] += loc[i];
    }

    // add contribution to right hand side from applied forces
    for (int f = 0; f != num_neumann_faces; ++f) {
        const auto fvals = compute_applied_forces_2D(points,
                                                     neumann_faces[2 * f],
                                                     neumann_faces[2 * f + 1],
                                                     neumann_forces[2 * f],
                                                     neumann_forces[2 * f + 1]);
        b[2 * neumann_faces[2 * f]] += fvals[0];
        b[2 * neumann_faces[2 * f] + 1] += fvals[1];
        b[2 * neumann_faces[2 * f + 1]] += fvals[2];
        b[2 * neumann_faces[2 * f + 1] + 1] += fvals[3];
    }

    // reduce system by eliminating Dirichlet degrees of freedom, and returning
    if (reduce_boundary) {
        reduce_system(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    } else {
        set_boundary_conditions(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    }
    return b.size();
};

// ----------------------------------------------------------------------------
int
assemble_mech_system_3D(const double* const points,
                        const int num_cells,
                        const int* const num_cell_faces, // cell faces per cell
                        const int* const num_face_corners, // corners per face
                        const int* const face_corners,
                        const double* const young,
                        const double* const poisson,
                        const double* const body_force, // 3 * number of cells
                        const int num_fixed_dofs, // dirichlet
                        const int* const fixed_dof_ixs, // indices must be sorted
                        const double* const fixed_dof_values,
                        const int num_neumann_faces,
                        const int* const neumann_faces,
                        const double* const neumann_forces, // 3 * number of neumann faces
                        std::vector<std::tuple<int, int, double>>& A_entries,
                        std::vector<double>& b,
                        const StabilityChoice stability_choice,
                        bool reduce_boundary,
                        const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    // preliminary computations
    const int tot_num_cell_faces = accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    const int tot_num_face_corners
        = accumulate(num_face_corners, num_face_corners + tot_num_cell_faces, 0);
    const int num_points = *max_element(face_corners, face_corners + tot_num_face_corners) + 1;

    // assemble full system matrix
    A_entries.clear();
    b.resize(num_points * 3);
    fill(b.begin(), b.end(), 0);

    // loop over cells and assemble system matrix
    vector<int> loc_indexing;
    vector<double> loc; // use as local 'scratch' vector
    array<double, 3> centroid; // will contain the centroid for the currently treated cell
    int cf_ix = 0; // index to first cell face for the current cell
    int fcorners_start = 0; // index to first cell corner for current cell

    for (int c = 0; c != num_cells; ++c) {
        // computing local stiffness matrix, writing its entries into the global matrix
        assemble_stiffness_matrix_3D(points,
                                     &face_corners[fcorners_start],
                                     &num_face_corners[cf_ix],
                                     num_cell_faces[c],
                                     young[c],
                                     poisson[c],
                                     stability_choice,
                                     centroid,
                                     loc_indexing,
                                     loc,
                                     nullptr,
                                     corner_mask ? corner_mask + fcorners_start : nullptr);
        const int ncv = int(loc_indexing.size()); // number of cell vertices
        for (int i = 0; i != 3 * ncv; ++i) {
            for (int j = 0; j != 3 * ncv; ++j) {
                // using integer division to go from 'i' and 'j' to corresponding
                // local point indices, which are subsequently converted to global indexing
                const int I = 3 * loc_indexing[i / 3] + i % 3;
                const int J = 3 * loc_indexing[j / 3] + j % 3;
                const double val = loc[i * 3 * ncv + j]; // local matrix entry
                A_entries.push_back(tuple<int, int, double> {I, J, val});
            }
        }
        // add contribution to right hand side from body forces.  These are
        // written directly into the global right-hand side vector
        compute_bodyforce_3D(points,
                             &face_corners[fcorners_start],
                             &num_face_corners[cf_ix],
                             num_cell_faces[c],
                             &centroid[0],
                             &body_force[3 * c],
                             &b[0],
                             corner_mask ? corner_mask + fcorners_start : nullptr);
        // increment local indexing
        fcorners_start
            += accumulate(&num_face_corners[cf_ix], &num_face_corners[cf_ix + num_cell_faces[c]], 0);
        cf_ix += num_cell_faces[c];
    }

    // add contribution to right hand side from applied forces.  Values are written
    // directly into the global right-hand-side vector
    compute_applied_forces_3D(
        points, num_face_corners, face_corners, num_neumann_faces, neumann_faces, neumann_forces, &b[0],
        corner_mask);

    if (reduce_boundary) {
        reduce_system(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    } else {
        set_boundary_conditions(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    }


    return b.size();
}

void
compute_stress_3D(const double* const points,
                  const int num_cells,
                  const int* const num_cell_faces, // cell faces per cell
                  const int* const num_face_corners, // corners per face
                  const int* const face_corners,
                  const double* const young,
                  const double* const poisson,
                  const std::vector<double>& disp,
                  std::vector<std::array<double, 6>>& stress,
                  vector<tuple<int, int, double>>& stressmat,
                  bool do_matrix,
                  bool do_stress,
                  vem::StabilityChoice stability_choice,
                  bool stab_on_stress,
                  const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    int cf_ix = 0; // index to first cell face for the current cell
    int fcorners_start = 0; // index to first cell corner for current cell

    for (int c = 0; c != num_cells; ++c) {
        // computing local stiffness matrix, writing its entries into the global matrix
        calculate_stress_3D_local(points,
                                  &face_corners[fcorners_start],
                                  &num_face_corners[cf_ix],
                                  num_cell_faces[c],
                                  young[c],
                                  poisson[c],
                                  disp, // global displacement
                                  c,
                                  stress[c],
                                  stressmat,
                                  do_matrix,
                                  do_stress,
                                  corner_mask ? corner_mask + fcorners_start : nullptr);
        fcorners_start
            += accumulate(&num_face_corners[cf_ix], &num_face_corners[cf_ix + num_cell_faces[c]], 0);
        cf_ix += num_cell_faces[c];
    }
}



// ============================================================================
// ===== SMOOTH_HOURGLASS: smooth-field (assumed-strain) hourglass control =====
// ============================================================================

namespace {

// S is computed with lambda capped at this Poisson ratio, so it stays mu-scaled as nu -> 1/2.
constexpr double smooth_hg_nu_cap = 0.3;
// Numerical-degeneracy guard only: the energy stays bounded on near-pyramids (ratio 1e-3).
constexpr double smooth_hg_min_jac_ratio = 1e-8;

double
capped_lambda(const double young, const double poisson)
{
    const double mu = young / (2 * (1 + poisson));
    const double nu = min(poisson, smooth_hg_nu_cap);
    return min(young * poisson / ((1 + poisson) * (1 - 2 * poisson)), 2 * mu * nu / (1 - 2 * nu));
}

// Cholesky in place (lower), n x n row-major; false if not positive definite.
bool
small_cholesky(vector<double>& A, const int n)
{
    for (int j = 0; j != n; ++j) {
        double d = A[j * n + j];
        for (int k = 0; k != j; ++k)
            d -= A[j * n + k] * A[j * n + k];
        if (!(d > 0))
            return false;
        A[j * n + j] = sqrt(d);
        for (int i = j + 1; i != n; ++i) {
            double s = A[i * n + j];
            for (int k = 0; k != j; ++k)
                s -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / A[j * n + j];
        }
    }
    return true;
}

// Energy of the minimum-energy field in Q1 + span{xi_d (1 - xi_d)} (Wilson bubbles, per
// component) through the 8 nodal values; hc in Dune corner order. Exact in bending and twist on boxes.
// remove_linear: drop the energy of linear fields, which the consistency term already carries.
bool
smooth_hex_energy(const double* const hc, const double lambda, const double mu, const bool remove_linear,
                  double* Ks)
{
    constexpr int nq = 24, nb = 9, nt = nq + nb;
    const double gp[3] = {0.5 - 0.5 * sqrt(0.6), 0.5, 0.5 + 0.5 * sqrt(0.6)};
    const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};
    vector<double> A(nt * nt, 0.0);
    double jmin = 1e300, jmax = 0;
    for (int a = 0; a != 3; ++a)
        for (int b = 0; b != 3; ++b)
            for (int c = 0; c != 3; ++c) {
                const double xi[3] = {gp[a], gp[b], gp[c]};
                double gref[8][3];
                hex_shape_gradients(xi, gref);
                double J[3][3] = {{0}};
                for (int k = 0; k != 8; ++k)
                    for (int i = 0; i != 3; ++i)
                        for (int j = 0; j != 3; ++j)
                            J[i][j] += hc[3 * k + i] * gref[k][j];
                const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                    - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                    + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
                jmin = min(jmin, det);
                jmax = max(jmax, det);
                if (!(det > 0))
                    return false;
                double Ji[3][3]; // J^{-1}
                Ji[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
                Ji[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det;
                Ji[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
                Ji[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det;
                Ji[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
                Ji[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
                Ji[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
                Ji[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det;
                Ji[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
                // physical gradients of the 8 shape functions and the 3 bubbles 4 xi_d (1 - xi_d)
                double g[11][3];
                for (int m = 0; m != 11; ++m) {
                    double r[3] = {0, 0, 0};
                    if (m < 8)
                        for (int e = 0; e != 3; ++e)
                            r[e] = gref[m][e];
                    else
                        r[m - 8] = 4 * (1 - 2 * xi[m - 8]);
                    for (int i = 0; i != 3; ++i)
                        g[m][i] = Ji[0][i] * r[0] + Ji[1][i] * r[1] + Ji[2][i] * r[2];
                }
                // dof I: function m = I / 3 (bubbles after the nodes), component I % 3
                double B[6][nt];
                for (int I = 0; I != nt; ++I) {
                    const int m = I < nq ? I / 3 : 8 + (I - nq) / 3, d = I % 3;
                    for (int r = 0; r != 6; ++r)
                        B[r][I] = 0;
                    B[d][I] = g[m][d];
                    // engineering shear rows xy, yz, zx
                    if (d == 0) { B[3][I] = g[m][1]; B[5][I] = g[m][2]; }
                    if (d == 1) { B[3][I] = g[m][0]; B[4][I] = g[m][2]; }
                    if (d == 2) { B[4][I] = g[m][1]; B[5][I] = g[m][0]; }
                }
                const double w = gw[a] * gw[b] * gw[c] * det;
                for (int I = 0; I != nt; ++I) {
                    const double trI = B[0][I] + B[1][I] + B[2][I];
                    double CB[6];
                    for (int r = 0; r != 3; ++r)
                        CB[r] = lambda * trI + 2 * mu * B[r][I];
                    for (int r = 3; r != 6; ++r)
                        CB[r] = mu * B[r][I];
                    for (int Jx = 0; Jx != nt; ++Jx) {
                        double s = 0;
                        for (int r = 0; r != 6; ++r)
                            s += CB[r] * B[r][Jx];
                        A[I * nt + Jx] += w * s;
                    }
                }
            }
    if (jmin < smooth_hg_min_jac_ratio * jmax)
        return false;
    // condense the bubbles: Ks = Aqq - Aqb Abb^-1 Abq
    vector<double> L(nb * nb);
    for (int i = 0; i != nb; ++i)
        for (int j = 0; j != nb; ++j)
            L[i * nb + j] = A[(nq + i) * nt + nq + j];
    if (!small_cholesky(L, nb))
        return false;
    vector<double> X(nb * nq); // Abb^-1 Abq
    for (int col = 0; col != nq; ++col) {
        double y[nb];
        for (int i = 0; i != nb; ++i) {
            double s = A[(nq + i) * nt + col];
            for (int k = 0; k != i; ++k)
                s -= L[i * nb + k] * y[k];
            y[i] = s / L[i * nb + i];
        }
        for (int i = nb - 1; i >= 0; --i) {
            double s = y[i];
            for (int k = i + 1; k != nb; ++k)
                s -= L[k * nb + i] * y[k];
            y[i] = s / L[i * nb + i];
        }
        for (int i = 0; i != nb; ++i)
            X[i * nq + col] = y[i];
    }
    for (int i = 0; i != nq; ++i)
        for (int j = 0; j != nq; ++j) {
            double s = A[i * nt + j];
            for (int k = 0; k != nb; ++k)
                s -= A[i * nt + nq + k] * X[k * nq + j];
            Ks[i * nq + j] = s;
        }
    if (!remove_linear)
        return true;
    // Ks <- Ks - Ks L (L^T Ks L)^-1 L^T Ks, L the 6 linear strain modes at the corners
    double xc[3] = {0, 0, 0}, h = 0;
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            xc[d] += hc[3 * k + d] / 8;
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            h = max(h, fabs(hc[3 * k + d] - xc[d]));
    const int sym[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}};
    vector<double> Lm(nq * 6, 0.0);
    for (int k = 0; k != 8; ++k)
        for (int m = 0; m != 6; ++m) {
            const int a = sym[m][0], b = sym[m][1];
            Lm[(3 * k + a) * 6 + m] += (hc[3 * k + b] - xc[b]) / h;
            if (a != b)
                Lm[(3 * k + b) * 6 + m] += (hc[3 * k + a] - xc[a]) / h;
        }
    vector<double> KL(nq * 6, 0.0), M(36, 0.0);
    for (int i = 0; i != nq; ++i)
        for (int j = 0; j != nq; ++j)
            for (int m = 0; m != 6; ++m)
                KL[i * 6 + m] += Ks[i * nq + j] * Lm[j * 6 + m];
    for (int i = 0; i != nq; ++i)
        for (int m = 0; m != 6; ++m)
            for (int l = 0; l != 6; ++l)
                M[m * 6 + l] += Lm[i * 6 + m] * KL[i * 6 + l];
    if (!small_cholesky(M, 6))
        return false;
    for (int i = 0; i != nq; ++i) {
        double y[6]; // M^-1 (KL)^T column i
        for (int r = 0; r != 6; ++r) {
            double s = KL[i * 6 + r];
            for (int k = 0; k != r; ++k)
                s -= M[r * 6 + k] * y[k];
            y[r] = s / M[r * 6 + r];
        }
        for (int r = 5; r >= 0; --r) {
            double s = y[r];
            for (int k = r + 1; k != 6; ++k)
                s -= M[k * 6 + r] * y[k];
            y[r] = s / M[r * 6 + r];
        }
        for (int j = 0; j != nq; ++j) {
            double s = 0;
            for (int r = 0; r != 6; ++r)
                s += KL[j * 6 + r] * y[r];
            Ks[j * nq + i] -= s;
        }
    }
    return true;
}

// 2D analogue on a quad (Dune corner order): Q4 + 2 bubbles per component, plane strain.
bool
smooth_quad_energy(const double* const qc, const double lambda, const double mu, double* Ks)
{
    constexpr int nq = 8, nb = 4, nt = nq + nb;
    const double gp[3] = {0.5 - 0.5 * sqrt(0.6), 0.5, 0.5 + 0.5 * sqrt(0.6)};
    const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};
    vector<double> A(nt * nt, 0.0);
    double jmin = 1e300, jmax = 0;
    for (int a = 0; a != 3; ++a)
        for (int b = 0; b != 3; ++b) {
            const double xi[2] = {gp[a], gp[b]};
            double gref[4][2];
            for (int k = 0; k != 4; ++k) {
                const int c[2] = {k & 1, (k >> 1) & 1};
                gref[k][0] = (c[0] ? 1.0 : -1.0) * (c[1] ? xi[1] : 1.0 - xi[1]);
                gref[k][1] = (c[1] ? 1.0 : -1.0) * (c[0] ? xi[0] : 1.0 - xi[0]);
            }
            double J[2][2] = {{0, 0}, {0, 0}};
            for (int k = 0; k != 4; ++k)
                for (int i = 0; i != 2; ++i)
                    for (int j = 0; j != 2; ++j)
                        J[i][j] += qc[2 * k + i] * gref[k][j];
            const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            jmin = min(jmin, det);
            jmax = max(jmax, det);
            if (!(det > 0))
                return false;
            const double Ji[2][2] = {{J[1][1] / det, -J[0][1] / det}, {-J[1][0] / det, J[0][0] / det}};
            double g[6][2];
            for (int m = 0; m != 6; ++m) {
                double r[2] = {0, 0};
                if (m < 4) { r[0] = gref[m][0]; r[1] = gref[m][1]; }
                else r[m - 4] = 4 * (1 - 2 * xi[m - 4]);
                for (int i = 0; i != 2; ++i)
                    g[m][i] = Ji[0][i] * r[0] + Ji[1][i] * r[1];
            }
            double B[3][nt];
            for (int I = 0; I != nt; ++I) {
                const int m = I < nq ? I / 2 : 4 + (I - nq) / 2, d = I % 2;
                B[0][I] = d == 0 ? g[m][0] : 0;
                B[1][I] = d == 1 ? g[m][1] : 0;
                B[2][I] = g[m][1 - d];
            }
            const double w = gw[a] * gw[b] * det;
            for (int I = 0; I != nt; ++I)
                for (int Jx = 0; Jx != nt; ++Jx) {
                    const double tI = B[0][I] + B[1][I], tJ = B[0][Jx] + B[1][Jx];
                    A[I * nt + Jx] += w * (lambda * tI * tJ + 2 * mu * (B[0][I] * B[0][Jx] + B[1][I] * B[1][Jx])
                                           + mu * B[2][I] * B[2][Jx]);
                }
        }
    if (jmin < smooth_hg_min_jac_ratio * jmax)
        return false;
    vector<double> L(nb * nb);
    for (int i = 0; i != nb; ++i)
        for (int j = 0; j != nb; ++j)
            L[i * nb + j] = A[(nq + i) * nt + nq + j];
    if (!small_cholesky(L, nb))
        return false;
    for (int i = 0; i != nq; ++i)
        for (int j = 0; j != nq; ++j) {
            // Abb^-1 Abq column j by forward/back substitution
            double y[nb];
            for (int r = 0; r != nb; ++r) {
                double s = A[(nq + r) * nt + j];
                for (int k = 0; k != r; ++k)
                    s -= L[r * nb + k] * y[k];
                y[r] = s / L[r * nb + r];
            }
            for (int r = nb - 1; r >= 0; --r) {
                double s = y[r];
                for (int k = r + 1; k != nb; ++k)
                    s -= L[k * nb + r] * y[k];
                y[r] = s / L[r * nb + r];
            }
            double s = A[i * nt + j];
            for (int k = 0; k != nb; ++k)
                s -= A[i * nt + nq + k] * y[k];
            Ks[i * nq + j] = s;
        }
    return true;
}

double
hex_shape_value(const int k, const double* const xi)
{
    double v = 1;
    for (int e = 0; e != 3; ++e)
        v *= ((k >> e) & 1) ? xi[e] : 1 - xi[e];
    return v;
}

// Reference coordinates of every node of a hex with extra nodes, by Newton on the trilinear map of
// its 8 corners (Dune order); false if a node falls outside the hex.
bool
extra_node_xi(const double* const corners,
              const int num_corners,
              const array<int, 8>& order,
              vector<array<double, 3>>& xi)
{
    array<double, 24> hc;
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            hc[3 * k + d] = corners[3 * order[k] + d];
    // reference coordinates of the extra nodes by Newton on the trilinear map
    double h = 0;
    for (int k = 1; k != 8; ++k)
        h = max(h, hypot(hc[3 * k] - hc[0], hc[3 * k + 1] - hc[1], hc[3 * k + 2] - hc[2]));
    xi.assign(num_corners, {0.5, 0.5, 0.5});
    for (int k = 0; k != 8; ++k)
        xi[order[k]] = {double(k & 1), double((k >> 1) & 1), double((k >> 2) & 1)};
    for (int v = 0; v != num_corners; ++v) {
        if (find(order.begin(), order.end(), v) != order.end())
            continue;
        auto& r = xi[v];
        double res = 1e300;
        for (int it = 0; it != 30 && res > 1e-10 * h; ++it) {
            double x[3] = {0, 0, 0}, gref[8][3], J[3][3] = {{0}};
            hex_shape_gradients(&r[0], gref);
            for (int k = 0; k != 8; ++k) {
                const double N = hex_shape_value(k, &r[0]);
                for (int i = 0; i != 3; ++i) {
                    x[i] += hc[3 * k + i] * N;
                    for (int j = 0; j != 3; ++j)
                        J[i][j] += hc[3 * k + i] * gref[k][j];
                }
            }
            double f[3];
            for (int i = 0; i != 3; ++i)
                f[i] = corners[3 * v + i] - x[i];
            res = norm<3>(array<double, 3> {f[0], f[1], f[2]});
            const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
            if (!(fabs(det) > 0))
                return false;
            // Cramer's rule for J dxi = f
            for (int c = 0; c != 3; ++c) {
                double M[3][3];
                for (int i = 0; i != 3; ++i)
                    for (int j = 0; j != 3; ++j)
                        M[i][j] = j == c ? f[i] : J[i][j];
                r[c] += (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                         - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                         + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
            }
        }
        if (res > 1e-6 * h)
            return false;
        for (auto& c : r) {
            if (c < -0.05 || c > 1.05)
                return false;
            c = min(1.0, max(0.0, c));
        }
    }
    return true;
}

// Recognise a polyhedron as a hexahedron with extra nodes on its edges/faces (faults, hanging
// nodes): faces are merged across edges with dihedral angle < 30 degrees. On success, order holds
// the 8 corners in Dune order and xi the reference coordinates of every node.
bool
logical_hex(const double* const corners,
            const int num_corners,
            const int* const faces,
            const int* const num_face_edges,
            const int num_faces,
            array<int, 8>& order,
            vector<array<double, 3>>& xi)
{
    vector<array<double, 3>> nrm(num_faces);
    vector<int> start(num_faces + 1, 0);
    for (int f = 0; f != num_faces; ++f)
        start[f + 1] = start[f] + num_face_edges[f];
    map<pair<int, int>, vector<int>> edge_faces;
    for (int f = 0; f != num_faces; ++f) {
        array<double, 3> nn {0, 0, 0}; // Newell normal
        const int m = num_face_edges[f];
        for (int i = 0; i != m; ++i) {
            const double* p = &corners[3 * faces[start[f] + i]];
            const double* q = &corners[3 * faces[start[f] + (i + 1) % m]];
            nn[0] += (p[1] - q[1]) * (p[2] + q[2]);
            nn[1] += (p[2] - q[2]) * (p[0] + q[0]);
            nn[2] += (p[0] - q[0]) * (p[1] + q[1]);
            const int u = faces[start[f] + i], v = faces[start[f] + (i + 1) % m];
            if (u != v)
                edge_faces[{min(u, v), max(u, v)}].push_back(f);
        }
        const double l = norm<3>(nn);
        if (!(l > 0))
            return false;
        for (auto& x : nn)
            x /= l;
        nrm[f] = nn;
    }
    vector<int> root(num_faces);
    iota(root.begin(), root.end(), 0);
    const auto find = [&](int f) { while (root[f] != f) f = root[f] = root[root[f]]; return f; };
    const double cos_tol = cos(M_PI / 6);
    for (const auto& [e, fs] : edge_faces) {
        if (fs.size() != 2)
            return false;
        const auto& a = nrm[fs[0]];
        const auto& b = nrm[fs[1]];
        if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] > cos_tol)
            root[find(fs[0])] = find(fs[1]);
    }
    map<int, int> lface; // root -> logical face id
    for (int f = 0; f != num_faces; ++f)
        lface.emplace(find(f), int(lface.size()));
    if (lface.size() != 6)
        return false;
    vector<set<int>> node_lf(num_corners);
    for (int f = 0; f != num_faces; ++f)
        for (int i = 0; i != num_face_edges[f]; ++i)
            node_lf[faces[start[f] + i]].insert(lface[find(f)]);
    vector<int> cid(num_corners, -1), cnode;
    for (int v = 0; v != num_corners; ++v) {
        if (node_lf[v].size() > 3)
            return false;
        if (node_lf[v].size() == 3) {
            cid[v] = int(cnode.size());
            cnode.push_back(v);
        }
    }
    if (cnode.size() != 8)
        return false;
    // corners of each logical face, in boundary order
    vector<int> quads;
    for (const auto& [r, id] : lface) {
        map<int, int> next;
        int nbnd = 0;
        for (int f = 0; f != num_faces; ++f) {
            if (find(f) != r)
                continue;
            const int m = num_face_edges[f];
            for (int i = 0; i != m; ++i) {
                const int u = faces[start[f] + i], v = faces[start[f] + (i + 1) % m];
                if (u == v)
                    continue;
                const auto& fs = edge_faces[{min(u, v), max(u, v)}];
                const int other = fs[0] == f ? fs[1] : fs[0];
                if (find(other) == r)
                    continue;
                if (!next.emplace(u, v).second)
                    return false;
                ++nbnd;
            }
        }
        int v0 = -1;
        for (const auto& [u, v] : next)
            if (cid[u] >= 0) { v0 = u; break; }
        if (v0 < 0)
            return false;
        vector<int> qc;
        int v = v0, steps = 0;
        do {
            if (cid[v] >= 0)
                qc.push_back(cid[v]);
            const auto it = next.find(v);
            if (it == next.end())
                return false;
            v = it->second;
        } while (v != v0 && ++steps <= nbnd);
        if (v != v0 || steps + 1 != nbnd || qc.size() != 4)
            return false;
        quads.insert(quads.end(), qc.begin(), qc.end());
    }
    vector<double> cc(24);
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            cc[3 * k + d] = corners[3 * cnode[k] + d];
    const int nfe4[6] = {4, 4, 4, 4, 4, 4};
    array<int, 8> ord;
    if (!hex_corner_ordering(&cc[0], 8, &quads[0], nfe4, 6, ord))
        return false;
    for (int k = 0; k != 8; ++k)
        order[k] = cnode[ord[k]];
    return extra_node_xi(corners, num_corners, order, xi);
}

// S = C^T K_smooth C + R^T S_trace R, C the corner selection and R = I - (trilinear interpolation
// of the corners); R = 0 on plain hexes. Pure TRACE when no well-conditioned corner hex exists.
// hex_corners: the cell's 8 corners (local node ids, Dune order) if the grid knows them.
vector<double>
compute_S_smooth_hg_3D(const double* const corners,
                       const int num_corners,
                       const int* const faces,
                       const int* const num_face_edges,
                       const int num_faces,
                       const double volume,
                       const double young,
                       const double poisson,
                       const int* const hex_corners,
                       const char* const corner_mask)
{
    const auto trace = [&] {
        return compute_S_trace_3D(corners, num_corners, faces, num_face_edges, num_faces, volume, young, poisson,
                                  corner_mask);
    };
    array<int, 8> order;
    vector<array<double, 3>> xi;
    bool plain = false;
    if (hex_corners) {
        // collapsed corners share a node: a degenerate hex whose tied corners add up below
        copy(hex_corners, hex_corners + 8, order.begin());
        const int distinct = int(set<int>(order.begin(), order.end()).size());
        if (distinct < 4)
            return trace();
        array<double, 24> hc;
        for (int k = 0; k != 8; ++k)
            for (int d = 0; d != 3; ++d)
                hc[3 * k + d] = corners[3 * order[k] + d];
        double g[8][3];
        const double centre[3] = {0.5, 0.5, 0.5};
        if (hex_physical_gradients(&hc[0], centre, g) < 0) // left-handed: mirror in x
            for (int k = 0; k != 8; k += 2)
                swap(order[k], order[k + 1]);
        plain = num_corners == distinct;
        if (!plain && !extra_node_xi(corners, num_corners, order, xi))
            return trace();
    } else {
        plain = hex_corner_ordering(corners, num_corners, faces, num_face_edges, num_faces, order);
        if (!plain && !logical_hex(corners, num_corners, faces, num_face_edges, num_faces, order, xi))
            return trace();
    }
    array<double, 24> hc;
    for (int k = 0; k != 8; ++k)
        for (int d = 0; d != 3; ++d)
            hc[3 * k + d] = corners[3 * order[k] + d];
    const double mu = young / (2 * (1 + poisson));
    vector<double> Ks(24 * 24);
    // with extra nodes, P mixes them into the corner values: keep S off the linear part
    if (!smooth_hex_energy(&hc[0], capped_lambda(young, poisson), mu, !plain, &Ks[0]))
        return trace();
    const int n = 3 * num_corners;
    vector<double> S(n * n, 0.0);
    for (int a = 0; a != 24; ++a)
        for (int b = 0; b != 24; ++b)
            S[(3 * order[a / 3] + a % 3) * n + 3 * order[b / 3] + b % 3] += Ks[a * 24 + b];
    if (plain)
        return S;
    // scalar R (num_corners^2); S_trace acts identically on each component
    const int m = num_corners;
    vector<double> R(m * m, 0.0);
    vector<bool> is_corner(m, false);
    for (int k = 0; k != 8; ++k)
        is_corner[order[k]] = true;
    for (int v = 0; v != m; ++v) {
        if (is_corner[v])
            continue;
        R[v * m + v] = 1;
        for (int k = 0; k != 8; ++k)
            R[v * m + order[k]] -= hex_shape_value(k, &xi[v][0]);
    }
    const auto St = trace();
    vector<double> Ts(m * m, 0.0), TR(m * m, 0.0);
    for (int i = 0; i != m; ++i)
        for (int j = 0; j != m; ++j)
            Ts[i * m + j] = St[(3 * i) * n + 3 * j];
    for (int i = 0; i != m; ++i)
        for (int k = 0; k != m; ++k)
            if (Ts[i * m + k] != 0)
                for (int j = 0; j != m; ++j)
                    TR[i * m + j] += Ts[i * m + k] * R[k * m + j];
    for (int i = 0; i != m; ++i)
        for (int k = 0; k != m; ++k)
            if (R[k * m + i] != 0)
                for (int j = 0; j != m; ++j) {
                    const double v = R[k * m + i] * TR[k * m + j];
                    for (int d = 0; d != 3; ++d)
                        S[(3 * i + d) * n + 3 * j + d] += v;
                }
    return S;
}

// 2D: smooth quad energy on quadrilaterals (perimeter order), TRACE on other polygons.
vector<double>
compute_S_smooth_hg_2D(const double* const corners, const int num_corners, const double area,
                       const double young, const double poisson)
{
    if (num_corners == 4) {
        // Dune corner order, mirrored for clockwise input so the Jacobian is positive
        const double* c = corners;
        const bool ccw = (c[2] - c[0]) * (c[7] - c[1]) - (c[3] - c[1]) * (c[6] - c[0]) > 0;
        const int p[4] = {ccw ? 0 : 1, ccw ? 1 : 0, ccw ? 3 : 2, ccw ? 2 : 3};
        array<double, 8> qc;
        for (int k = 0; k != 4; ++k)
            for (int d = 0; d != 2; ++d)
                qc[2 * k + d] = corners[2 * p[k] + d];
        const double mu = young / (2 * (1 + poisson));
        vector<double> Ks(64);
        if (smooth_quad_energy(&qc[0], capped_lambda(young, poisson), mu, &Ks[0])) {
            vector<double> S(64);
            for (int a = 0; a != 8; ++a)
                for (int b = 0; b != 8; ++b)
                    S[(2 * p[a / 2] + a % 2) * 8 + 2 * p[b / 2] + b % 2] = Ks[a * 8 + b];
            return S;
        }
    }
    return compute_S_trace_2D(corners, num_corners, area, young, poisson);
}

} // anonymous namespace

// ----------------------------------------------------------------------------
void
assemble_stiffness_matrix_2D(const double* const points,
                             const int* const corner_ixs,
                             const int num_corners,
                             const double young,
                             const double poisson,
                             const StabilityChoice stability_choice,
                             vector<double>& target)
// ----------------------------------------------------------------------------
{
    // collect all affected points consecutively in one vector
    const auto corners = pick_points<2>(points, corner_ixs, num_corners);

    const double area = element_volume_2D(&corners[0], num_corners);

    // Q1 FEM requested: applicable whenever the cell is a quadrilateral.  The corners
    // arrive in (counter)clockwise perimeter order; the Dune reference square visits
    // them in the order 0,1,3,2, so corner 2 and 3 must be swapped.
    if (is_fem_choice(stability_choice) && num_corners == 4) {
        array<double, 8> qc;
        const int perm[4] = {0, 1, 3, 2};
        for (int k = 0; k != 4; ++k) {
            qc[2 * k] = corners[2 * perm[k]];
            qc[2 * k + 1] = corners[2 * perm[k] + 1];
        }
        target.resize(64);
        vector<double> Kfem(64);
        stiffness_matrix_fem_quad_2D(&qc[0], young, poisson, stability_choice == FEM_BBAR, &Kfem[0]);
        // permute back to perimeter corner ordering
        for (int a = 0; a != 8; ++a)
            for (int b = 0; b != 8; ++b) {
                const int i = 2 * perm[a / 2] + a % 2;
                const int j = 2 * perm[b / 2] + b % 2;
                target[i * 8 + j] = Kfem[a * 8 + b];
            }
        return;
    }

    // compute all intermediary matrices
    const auto q = compute_q_2D(&corners[0], num_corners);
    const auto Nr = compute_Nr_2D(&corners[0], num_corners);
    const auto Nc = compute_Nc_2D(&corners[0], num_corners);
    const auto Wr = compute_Wr_2D(q);
    const auto Wc = compute_Wc_2D(q);
    const auto D = compute_D_2D(young, poisson);
    const auto ImP = compute_ImP(Nr, Nc, Wr, Wc, 2);

    // do the final assembly of matrices, and write result to target
    target.resize(pow(2 * num_corners, 2));
    if (stability_choice == TRACE || stability_choice == SMOOTH_HOURGLASS) {
        const auto S = stability_choice == TRACE
            ? compute_S_trace_2D(&corners[0], num_corners, area, young, poisson)
            : compute_S_smooth_hg_2D(&corners[0], num_corners, area, young, poisson);
        final_assembly(Wc, D, Nc, ImP, corners, stability_choice, area, num_corners, 2, &target[0], &S);
        return;
    }
    final_assembly(Wc, D, Nc, ImP, corners, stability_choice, area, num_corners, 2, &target[0]);
}

// ----------------------------------------------------------------------------
void
assemble_stiffness_matrix_3D(const double* const points,
                             const int* const faces,
                             const int* const num_face_edges,
                             const int num_faces,
                             const double young,
                             const double poisson,
                             const StabilityChoice stability_choice,
                             array<double, 3>& centroid,
                             vector<int>& indexing,
                             vector<double>& target,
                             const int* const hex_corners,
                             const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    // compute mapping between local vertex indices (0, 1, ... Ne) to
    // indexing in the global 'points' vector (0, 1, .....N)
    const auto reindex = global_to_local_indexing(faces, num_face_edges, num_faces, indexing);
    const int num_corners = int(indexing.size());
    const int num_face_entries = accumulate(num_face_edges, num_face_edges + num_faces, 0);

    // the given corners are used only if all of them are nodes of the cell
    array<int, 8> hex_loc;
    bool have_hex = hex_corners != nullptr;
    for (int k = 0; have_hex && k != 8; ++k) {
        const auto it = reindex.find(hex_corners[k]);
        if (it != reindex.end()) {
            hex_loc[k] = it->second;
            continue;
        }
        // a collapsed corner may be a separate point at the same place as a cell node
        const double* p = &points[3 * hex_corners[k]];
        have_hex = false;
        for (int v = 0; v != num_corners && !have_hex; ++v) {
            const double* q = &points[3 * indexing[v]];
            if (fabs(p[0] - q[0]) + fabs(p[1] - q[1]) + fabs(p[2] - q[2]) <= 1e-9 * (1 + fabs(q[2]))) {
                hex_loc[k] = v;
                have_hex = true;
            }
        }
    }

    // make local list of corner coordinates, and a locally indexed version of
    // 'faces'
    const auto corners_loc = pick_points<3>(points, &indexing[0], num_corners);
    vector<int> faces_loc(num_face_entries);
    transform(faces, faces + num_face_entries, faces_loc.begin(), [&reindex](const int I) {
        return reindex.find(I)->second;
    });

    double volume; // computed in 'compute_cell_geometry' below
    vector<double> outward_normals, face_centroids; // computed in 'compute_cell_geometry' below
    array<double, 3> star_point; // computed in 'compute_cell_geometry' below
    compute_cell_geometry(&corners_loc[0],
                          num_corners,
                          &faces_loc[0],
                          num_face_edges,
                          num_faces,
                          outward_normals,
                          face_centroids,
                          centroid,
                          star_point,
                          volume,
                          corner_mask);

    // Q1 FEM requested: applicable whenever the cell is a topological hexahedron.
    // Non-hex cells fall through to VEM (with ANISO_DIAG stabilization, see
    // getStabilityMatrix).
    if (is_fem_choice(stability_choice)) {
        array<int, 8> order;
        if (hex_corner_ordering(
                &corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces, order)) {
            array<double, 24> hc;
            for (int k = 0; k != 8; ++k)
                for (int d = 0; d != 3; ++d)
                    hc[3 * k + d] = corners_loc[3 * order[k] + d];
            vector<double> Kfem(24 * 24);
            stiffness_matrix_fem_hex_3D(
                &hc[0], young, poisson, stability_choice == FEM_BBAR, &Kfem[0]);
            // scatter from reference-hex dof ordering back to local corner ordering
            target.assign(pow(3 * num_corners, 2), 0.0);
            const int n = 3 * num_corners;
            for (int a = 0; a != 24; ++a)
                for (int b = 0; b != 24; ++b) {
                    const int i = 3 * order[a / 3] + a % 3;
                    const int j = 3 * order[b / 3] + b % 3;
                    target[i * n + j] = Kfem[a * 24 + b];
                }
            return;
        }
    }

    // compute all intermediary matrices
    const auto q = compute_q_3D(&corners_loc[0],
                                int(corners_loc.size() / 3),
                                &faces_loc[0],
                                num_face_edges,
                                num_faces,
                                volume,
                                outward_normals,
                                corner_mask);
    const auto Nr = compute_Nr_3D(&corners_loc[0], num_corners);
    const auto Nc = compute_Nc_3D(&corners_loc[0], num_corners);
    const auto Wr = compute_Wr_3D(q);
    const auto Wc = compute_Wc_3D(q);
    const auto D = compute_D_3D(young, poisson);
    const auto ImP = compute_ImP(Nr, Nc, Wr, Wc, 3);

    // do the final assembly of matrices, and write result to target
    target.resize(pow(3 * num_corners, 2));
    if (stability_choice == TRACE || stability_choice == SMOOTH_HOURGLASS) {
        const auto S = stability_choice == TRACE
            ? compute_S_trace_3D(&corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces,
                                 volume, young, poisson, corner_mask)
            : compute_S_smooth_hg_3D(&corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces,
                                     volume, young, poisson, have_hex ? &hex_loc[0] : nullptr, corner_mask);
        final_assembly(Wc, D, Nc, ImP, corners_loc, stability_choice, volume, num_corners, 3, &target[0], &S);
        return;
    }
    const auto node_weights = is_face_weighted_choice(stability_choice)
        ? compute_boundary_weights_3D(&corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces)
        : vector<double>{};
    final_assembly(Wc, D, Nc, ImP, corners_loc, stability_choice, volume, num_corners, 3, &target[0], nullptr,
                   node_weights);
}

void
calculate_stress_3D_local(const double* const points,
                          const int* const faces,
                          const int* const num_face_edges,
                          const int num_faces,
                          const double young,
                          const double poisson,
                          const std::vector<double>& disp, /// global displacement
                          const int cell,
                          std::array<double, 6>& stress,
                          vector<tuple<int, int, double>>& stressmat,
                          bool do_matrix,
                          bool do_stress,
                          const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    array<double, 3> centroid;
    vector<int> indexing;
    vector<double> target;
    // compute mapping between local vertex indices (0, 1, ... Ne) to
    // indexing in the global 'points' vector (0, 1, .....N)
    const auto reindex = global_to_local_indexing(faces, num_face_edges, num_faces, indexing);
    const int num_corners = int(indexing.size());
    const int num_face_entries = accumulate(num_face_edges, num_face_edges + num_faces, 0);

    // make local list of corner coordinates, and a locally indexed version of
    // 'faces'
    const auto corners_loc = pick_points<3>(points, &indexing[0], num_corners);
    vector<int> faces_loc(num_face_entries);
    transform(faces, faces + num_face_entries, faces_loc.begin(), [&reindex](const int I) {
        return reindex.find(I)->second;
    });

    double volume; // computed in 'compute_cell_geometry' below
    vector<double> outward_normals, face_centroids; // computed in 'compute_cell_geometry' below
    array<double, 3> star_point; // computed in 'compute_cell_geometry' below
    compute_cell_geometry(&corners_loc[0],
                          num_corners,
                          &faces_loc[0],
                          num_face_edges,
                          num_faces,
                          outward_normals,
                          face_centroids,
                          centroid,
                          star_point,
                          volume,
                          corner_mask);

    // compute all intermediary matrices
    const auto q = compute_q_3D(&corners_loc[0],
                                int(corners_loc.size() / 3),
                                &faces_loc[0],
                                num_face_edges,
                                num_faces,
                                volume,
                                outward_normals,
                                corner_mask);

    const int dim = 3;
    const int lsdim = (dim == 2) ? 3 : 6; // dimension of "linear strain space"
    const int totdim = dim * num_corners; // total number of unknowns
    const auto Wc = compute_Wc_3D(q);
    std::vector<double> matrix(lsdim * totdim, 0.0);
    if (do_stress) {
        const auto D = compute_D_3D(young, poisson);
        // compute stiffness matrix components
        const auto DWct = matmul(&D[0], lsdim, lsdim, false, &Wc[0], totdim, lsdim, true);
        matrix = DWct;
    } else {
        for (int i = 0; i < lsdim; ++i) {
            for (int j = 0; j < totdim; ++j) {
                matrix[i * totdim + j] = Wc[j * lsdim + i];
            }
        }
    }

    std::vector<double> local_disp;
    local_disp.resize(indexing.size() * 3);
    std::vector<int> global_index;
    global_index.resize(indexing.size() * 3);
    for (size_t i = 0; i < indexing.size(); ++i) {
        for (size_t d = 0; d != 3; ++d) {
            local_disp[3 * i + d] = disp[3 * indexing[i] + d];
            global_index[3 * i + d] = 3 * indexing[i] + d;
        }
    }
    for (int i = 0; i < lsdim; ++i) {
        for (int j = 0; j < totdim; ++j) {
            stress[i] += matrix[i * totdim + j] * local_disp[j];
            if (do_matrix) {
                int I = lsdim * cell + i;
                int J = global_index[j];
                double val = matrix[i * totdim + j];
                if (i > 2 && do_stress) {
                    val /= 2.0;
                }
                stressmat.push_back(tuple<int, int, double> {I, J, val});
            }
        }
    }
    if (do_stress) {
        for (int i = 3; i < lsdim; ++i) {
            stress[i] /= 2;
        }
    }
}
// ----------------------------------------------------------------------------
void
matprint(const double* data, const int r, const int c, bool transposed, const double zthreshold)
// ----------------------------------------------------------------------------
{
    const pair<int, int> dim = transposed ? make_pair(c, r) : make_pair(r, c);
    const pair<int, int> stride = transposed ? make_pair(1, c) : make_pair(c, 1);

    for (int i = 0; i != dim.first; ++i)
        for (int j = 0; j != dim.second; ++j)
            cout << setw(12)
                 << ((abs(data[i * stride.first + j * stride.second]) <= zthreshold) ? fixed
                                                                                     : scientific)
                 << ((abs(data[i * stride.first + j * stride.second]) <= zthreshold) ? setprecision(0)
                                                                                     : setprecision(2))
                 << (abs(data[i * stride.first + j * stride.second]) <= zthreshold
                         ? 0
                         : data[i * stride.first + j * stride.second])
                 << (j == dim.second - 1 ? "\n" : "");
}

// ----------------------------------------------------------------------------
vector<double>
sparse2full(const vector<tuple<int, int, double>>& nz, const int r, const int c)
// ----------------------------------------------------------------------------
{
    vector<double> result(r * c, 0);
    for (auto it = nz.begin(); it != nz.end(); ++it) {
        const int ix = get<0>(*it) * r + get<1>(*it);
        assert(ix < int(result.size()));
        result[ix] += get<2>(*it);
    }
    return result;
}



// ----------------------------------------------------------------------------
// Return unsigned tetrahedron volume
double
tetrahedron_volume(const double* const p1,
                   const double* const p2,
                   const double* const p3,
                   const double* const p4)
// ----------------------------------------------------------------------------
{
    const auto v1 = pointdiff<3>(p1, p4);
    const auto v2 = pointdiff<3>(p2, p4);
    const auto v3 = pointdiff<3>(p3, p4);

    return abs(determinant3D(&v1[0], &v2[0], &v3[0])) / 6.0;
}

// ----------------------------------------------------------------------------
array<double, 2>
centroid_2D(const double* const points, const int num_points)
// ----------------------------------------------------------------------------
{
    // To find the centroid (C_x, C_y) we use the formula
    // C_d = 1/6A * \sum_{i=0}^{n-1} (d_i + d_{i+1}) (x_i y_{i+1} - x_{i+1} y_i)
    // where 'A' denotes area

    // We need the signed area, so cannot use the element_volume_2D function.
    // Instead, we compute the area alongside with the centroid coordinates below

    array<double, 2> result {0, 0};
    double area = 0;
    for (int i = 0, inext = 1; i != num_points; ++i, inext = (i + 1) % num_points) {
        // x_i * y_{i+1} - x_{i+1} * y_i
        double fac = (points[2 * i] * points[2 * inext + 1]) - (points[2 * inext] * points[2 * i + 1]);
        area += 0.5 * fac;
        for (int d = 0; d != 2; ++d)
            result[d] += (points[2 * i + d] + points[2 * inext + d]) * fac;
    }

    result[0] /= 6 * area;
    result[1] /= 6 * area;

    return result;
}

// ---------------------------------------------------------------------------
// May fail for strongly non-convex faces
array<double, 3>
centroid_2D_3D(const double* const points, const int num_points)
// ----------------------------------------------------------------------------
{
    // Compute centroid of planar polygon embedded in 3D space.

    const array<double, 3> inside_point = point_average<3>(points, num_points);

    array<double, 3> result {0, 0, 0};
    double area = 0;
    for (int i = 0, inext = 1; i != num_points; ++i, inext = (i + 1) % num_points) {
        const double tri_area = triarea<3>(&points[3 * i], &points[3 * inext], &inside_point[0]);
        area += tri_area;
        for (int d = 0; d != 3; ++d)
            result[d] += tri_area * (points[3 * i + d] + points[3 * inext + d] + inside_point[d]) / 3;
    }

    for (int d = 0; d != 3; ++d)
        result[d] /= area;

    return result;
}

// ----------------------------------------------------------------------------
// requires that corner points are already ordered consecutively
double
face_integral(const double* const corners,
              const int num_corners,
              const int dim,
              const double* const corner_values,
              const char* const corner_mask)
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);

    return dim == 2 ? face_integral_impl<2>(corners, num_corners, corner_values, corner_mask)
                    : face_integral_impl<3>(corners, num_corners, corner_values, corner_mask);
}

// ----------------------------------------------------------------------------
double
element_volume_2D(const double* const points, const int* const faces, const int num_faces)
// ----------------------------------------------------------------------------
{
    const auto corners = pick_points<2>(points, faces, num_faces);
    return face_integral(&corners[0], num_faces, 2);
}

// ----------------------------------------------------------------------------
double
element_volume_2D(const double* const corners, const int num_faces)
// ----------------------------------------------------------------------------
{
    return face_integral(corners, num_faces, 2);
}

// ----------------------------------------------------------------------------
// throw away these functions... only included to allow direct testing of intermediary results
vector<double>
pick_points_2D(const double* pts, const int* const p_ixs, int num_points)
// ----------------------------------------------------------------------------
{
    return pick_points<2>(pts, p_ixs, num_points);
}

// ----------------------------------------------------------------------------
vector<double>
pick_points_3D(const double* pts, const int* const p_ixs, int num_points)
// ----------------------------------------------------------------------------
{
    return pick_points<3>(pts, p_ixs, num_points);
}


}; // end namespace vem
