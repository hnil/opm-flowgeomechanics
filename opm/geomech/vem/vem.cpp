#include "vem.hpp"

#include <opm/porsol/common/blas_lapack.hpp>

#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator> //@@ only for debugging
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace std;
using namespace vem;

namespace
{

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
    transform(p1, p1 + Dim, p2, target, [fac1, fac2](double d1, double d2) { return d1 * fac1 + d2 * fac2; });
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
    const double sumsquared = accumulate(v, v + Dim, 0.0, [](double acc, double el) { return acc + el * el; });
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
    return array<double, 3> {v1[1] * v2[2] - v1[2] * v2[1],
                             v1[2] * v2[0] - v1[0] * v2[2],
                             v1[0] * v2[1] - v1[1] * v2[0]};
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
tessellate_face(const double* const corners, const int num_corners,
                bool skip_if_tri = true, double* centroid = nullptr)
// ----------------------------------------------------------------------------
{
    if (skip_if_tri && num_corners == 3)
        // the face is already a triangle - no need to split it up
        return {vector<double>(corners, corners + Dim * num_corners)};

    // compute face centroid
    const array<double, Dim> center = face_centroid<Dim>(corners, num_corners);

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
tessellate_face(const vector<double>& corners, bool skip_if_tri = true)
// ----------------------------------------------------------------------------
{
  return tessellate_face<Dim>(&corners[0], corners.size() / Dim, skip_if_tri);
}

// ----------------------------------------------------------------------------
template <int Dim>
double
face_integral_impl(const double* const corners, const int num_corners, const double* corner_values)
// ----------------------------------------------------------------------------
{
    const auto tris = tessellate_face<Dim>(&corners[0], num_corners, false);

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
    if ((transposed1 ? r1 : c1) != (transposed2 ? c2 : r2))
        throw invalid_argument("Matrices are not compatible for multiplication.");

    double null = 0.0;
    char T = 'T';
    char N = 'N';
    DGEMM(transposed1 ? &T : &N,
          transposed2 ? &T : &N,
          &r1, &c1, transposed1 ? &r1 : &c1,
          &fac, data1, &r1,
          data2, &r2, &null, result, &r1);
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
// @param transposed2 flag indicating whether the second matrix should be transposed before multiplication
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
double trace(const vector<double>& A)
// ----------------------------------------------------------------------------
{
    const int N = int(sqrt(A.size()));
    assert(N * N == A.size()); // should be a square matrix

    double result = 0;
    for (int i = 0; i != N; ++i)
        result += A[i + N * i];

    return result;
}

// ----------------------------------------------------------------------------
// extract the diagonal elements of a matrix
vector<double> diag_elems(const vector<double>& A)
// ----------------------------------------------------------------------------
{
  const int N = int(sqrt(A.size()));
  vector<double> result(N, 0);
  for (int i = 0; i != N; ++i)
    result[i] = A[i*N+i];
  return result;
}



// ----------------------------------------------------------------------------
// Compute the trace of the inverse matrix.  The function only works correctly
// for positive symmetric matrices
double invtrace(const vector<double>& A)
// ----------------------------------------------------------------------------
{
  const auto diag = diag_elems(A);
  double result = 0;
  for (int i = 0; i != (int)diag.size(); ++i)
    result += 1.0/diag[i];
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
// @@TODO: find a better estimate for S
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

    // const double alpha = 1;  // @@ use this line for comparison with vemmech
    const double alpha = stability_choice == SIMPLE ?
      volume * trace(D) / trace(NtN) :             // from Gain et al. 2014
      (1/9.0) * volume * trace(D) * invtrace(NtN); // Andersen et al. 2017

    return identity_matrix(alpha, dim * num_corners);
}

// ----------------------------------------------------------------------------
vector<double> compute_S_D_recipe(const std::vector<double>& EWcDWct,
                                  const int dofs,
                                  const double volume)
// ----------------------------------------------------------------------------
{
  vector<double> result(dofs * dofs, 0);

  const double h = cbrt(volume); // diameter of cell scales with cube root of volume
  for (int i = 0; i != dofs; ++i)
    result[i + dofs*i] = max(h, EWcDWct[i + dofs * i]); // extract diagonal from matrix

  return result;
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
               const StabilityChoice stability_choice,
               const double volume,
               const double num_nodes,
               const double dim,
               double* target)
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);
    const int lsdim = (dim == 2) ? 3 : 6; // dimension of "linear strain space"
    const int totdim = dim * num_nodes; // total number of unknowns

    // compute stiffness matrix components
    const auto DWct = matmul(&D[0], lsdim, lsdim, false, &Wc[0], totdim, lsdim, true);
    const auto EWcDWct = matmul(&Wc[0], totdim, lsdim, false, &DWct[0], lsdim, totdim, false, volume);

    const auto S = stability_choice == D_RECIPE ?
      compute_S_D_recipe(EWcDWct, totdim, volume) :
      compute_S(Nc, D, num_nodes, volume, dim, stability_choice);

    const auto SImP = matmul(&S[0], totdim, totdim, false, &ImP[0], totdim, totdim, false);
    const auto ImpSImp = matmul(&ImP[0], totdim, totdim, true, &SImP[0], totdim, totdim, false);
    assert(EWcDWct.size() == ImpSImp.size());

    // add conformance and stability term, and write result to target
    transform(EWcDWct.begin(), EWcDWct.end(), ImpSImp.begin(), target, [](double a, double b) { return a + b; });
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
    return vector<double> {e1, 0, e2,
                           0, e3, e4};
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

    return vector<double> {e1, 0, 0, e2, 0, e3,
                           0, e4, 0, e5, e6, 0,
                           0, 0, e7, 0, e8, e9};
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
    // @@ for the time being, do nothing
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
compute_applied_forces_2D(const double* const points, const int n1, const int n2, const double fx, const double fy)
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
                          double* b_global)
// ----------------------------------------------------------------------------
{
    for (int nf = 0; nf != num_neumann_faces; ++nf) {
        const int nfc = num_face_corners[neumann_faces[nf]];
        const int ix_start = accumulate(num_face_corners, num_face_corners + neumann_faces[nf], 0);
        const auto face_corners_loc = pick_points<3>(points, &face_corners[ix_start], nfc);
        const auto tris = tessellate_face<3>(&face_corners_loc[0], nfc, false);
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
                     double* b_global) // global right-hand side
// ----------------------------------------------------------------------------
{
    int fcorner_start = 0;
    for (int f = 0; f != num_faces; ++f) {
        const int nfc = num_face_corners[f];
        const auto tris = tessellate_face<3>(pick_points<3>(points, face_corners + fcorner_start, nfc), false);
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
    cout << "Reducing system: moving elements to right hand side" << endl;

    sort(A.begin(), A.end(), [](const auto& aa, const auto& bb) { return get<1>(aa) < get<1>(bb); });

    // eliminate columns associated with fixed dofs (move value to b)
    auto cur_end = A.begin();
    for (int i = 0; i != num_fixed_dofs; ++i) {
        const int dof = fixed_dof_ixs[i];
        // find the range of the elements whose column index equal 'dof'
        auto start = lower_bound(cur_end, A.end(), dof, [](const auto& a, int value) { return get<1>(a) < value; });
        if (start != A.end() && get<1>(*start) == dof) {
            cur_end = find_if(start, A.end(), [dof](const auto& a) { return get<1>(a) != dof; });
            for (; start != cur_end; ++start)
                b[get<0>(*start)] -= get<2>(*start) * fixed_dof_values[i];
        }
    }

    // determine renumbering
    cout << "Reducing system: determining renumbering" << endl;
    vector<int> keep(b.size());
    iota(keep.begin(), keep.end(), 0);
    vector<int> renum;
    renum.reserve(b.size());

    set_difference(keep.begin(), keep.end(), fixed_dof_ixs, fixed_dof_ixs + num_fixed_dofs, back_inserter(renum));

    const int discard_flag = b.size() + 1;
    vector<int> renum_inv(b.size(), discard_flag);
    for (int i = 0; i != int(renum.size()); ++i)
        renum_inv[renum[i]] = i;

    // eliminate rows associated with fixed dofs and renumber
    cout << "Reducing system: eliminating entries" << endl;
    A.erase(remove_if(A.begin(),
                      A.end(),
                      [&discard_flag, &renum_inv](const tuple<int, int, double>& el) {
                          return renum_inv[get<0>(el)] == discard_flag || renum_inv[get<1>(el)] == discard_flag;
                      }),
            A.end());
    transform(A.begin(), A.end(), A.begin(), [&renum_inv](const tuple<int, int, double>& el) {
        return tuple<int, int, double> {renum_inv[get<0>(el)], renum_inv[get<1>(el)], get<2>(el)};
    });

    for (int i = 0; i != int(renum.size()); ++i)
        b[i] = b[renum[i]];

    b.resize(renum.size());
    cout << "Reducing system: FINISHED" << endl;
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
             const vector<double>& outward_normals)
// ----------------------------------------------------------------------------
{
    vector<double> result(num_corners * 3, 0); // the vector of q

    // common factor outside the integral: 1 / (2 |E|)
    const double fac = 1.0 / (2 * volume);

    // loop over faces, an accumulate q-values for the associated vertices
    int cur_vertex_ix = 0;
    for (int f = 0; f != num_faces; ++f) {
        const auto facecorners = pick_points<3>(corners, &faces[cur_vertex_ix], num_face_edges[f]);
        const double* const normal = &outward_normals[3*f];
        vector<double> cvals(num_face_edges[f], 0);
        for (int e = 0; e != num_face_edges[f]; ++e, ++cur_vertex_ix) {
            fill(cvals.begin(), cvals.end(), 0);
            cvals[e] = 1; // indicate which basis function we are integrating

            const double phi = face_integral(&facecorners[0], num_face_edges[f], 3, &cvals[0]);
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
        append(result, matentry_2D(1, corners[2 * i + 1] - midpoint[1], 1, -(corners[2 * i] - midpoint[0])));
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
               matentry_3D(ncinv, q[3 * i + 1], -q[3 * i + 2],
                           ncinv, -q[3 * i], q[3 * i + 2],
                           ncinv, -q[3 * i + 1], q[3 * i]));
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
               matentry_3D(2 * q[3 * i],     q[3 * i + 1], q[3 * i + 2],
                           2 * q[3 * i + 1], q[3 * i],     q[3 * i + 2],
                           2 * q[3 * i + 2], q[3 * i + 1], q[3 * i]));
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
    vector<double> result {1 - poisson, poisson,     0,
                           poisson,     1 - poisson, 0,
                           0,           0,           2 * (1 - 2 * poisson)};

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
    vector<double> result {
      1 - poisson,     poisson,     poisson,                     0,                     0,                     0,
      poisson,     1 - poisson,     poisson,                     0,                     0,                     0,
      poisson,         poisson, 1 - poisson,                     0,                     0,                     0,
      0,                     0,           0, 2 * (1 - 2 * poisson),                     0,                     0,
      0,                     0,           0,                     0, 2 * (1 - 2 * poisson),                     0,
      0,                     0,           0,                     0,                     0, 2 * (1 - 2 * poisson)};

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
bool is_behind_face(const array<double, 3>& point,
                    const double* const face_normal,
                    const double* const face_centroid)
// ----------------------------------------------------------------------------
{
  const double tol = 1e-13;
  const array<double, 3> v {face_centroid[0] - point[0],
                            face_centroid[1] - point[1],
                            face_centroid[2] - point[2]};
  return (v[0] * face_normal[0] +
          v[1] * face_normal[1] +
          v[2] * face_normal[2]) + tol > 0;
}


// ----------------------------------------------------------------------------
bool is_star_point(const array<double, 3>& point,
                   const vector<double>& face_normals,
                   const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
  const int N = (int)face_normals.size() / 3;
  for (int i = 0; i != N; ++i)
    if (!is_behind_face(point, &face_normals[3*i], &face_centroids[3*i]))
      return false;
  return true;
}

// ----------------------------------------------------------------------------
// NB, the face normals must be _unit_ normals (not scaled by face area)
// in order for this algorithm to work as intended.
//
// @@ NB: The algorithm assumes planar faces, which is not always guaranteed in
// reality.  If it proves that the obtained star point creates inconsistencies
// in some settings, one might have to tessellate the faces into planar triangles first
array<double, 3> identify_star_point(const array<double, 3>& point,
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
      std::cout << "Averga point not insize cell " << std::endl;
      // projecting the point onto the plane defined by the normal and the centroid
      const double proj = (result[0] - face_centroids[3*f_ix])   * face_normals[3*f_ix] +
                          (result[1] - face_centroids[3*f_ix+1]) * face_normals[3*f_ix+1] +
                          (result[2] - face_centroids[3*f_ix+2]) * face_normals[3*f_ix+2];
      for (int ii = 0; ii != 3; ++ii)
        result[ii] -= proj * face_normals[3*f_ix+ii];
    }
    if (++count == N)
      break;
  }
  if (count != N)
      throw runtime_error("Unable to find a star point for cell.");
  return result;
}

// ----------------------------------------------------------------------------
void compute_face_geometry(const vector<double>& points, double* normal, double* centroid)
// ----------------------------------------------------------------------------
{
  // compute tessellation and face centroid
  const auto face_tessellation = tessellate_face<3>(&points[0],
                                                    (int)points.size()/3, false, centroid);

  // compute (normalized) face normal
  fill(normal, normal+3, 0);

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
vector<int> consistent_face_ordering(const int* const faces,
                                     const int* const num_face_edges,
                                     const int num_faces)
// ----------------------------------------------------------------------------
{
  vector<int> result(faces, faces + accumulate(num_face_edges, num_face_edges + num_faces, 0));
  //return result;
  vector<int> sorted_facenodes(result); // we need to sort the facenodes for quicker set operations
  for (int f = 0, fpos=0; f != num_faces; fpos += num_face_edges[f++])
    sort(&sorted_facenodes[fpos], &sorted_facenodes[fpos] + num_face_edges[f]);
  vector<int*> face_ptrs(1, &result[0]);
  vector<int*> sorted_fnode_ptrs(1, &sorted_facenodes[0]);
  for (int i = 1; i != num_faces; ++i) {
    face_ptrs.push_back(face_ptrs[i-1] + num_face_edges[i-1]);
    sorted_fnode_ptrs.push_back(sorted_fnode_ptrs[i-1] + num_face_edges[i-1]);
  }
  vector<int> ordered(num_faces, false), completed(num_faces, false);
  ordered[0] = true;

  while (!all_of(ordered.begin(), ordered.end(), [](auto x){return x;})) {
    // identfying the reference face to use this time
    int ref_face = -1;
    for (ref_face = 0; ref_face != num_faces; ++ref_face)
      if (ordered[ref_face] && !completed[ref_face])
        break;

    // check which faces that can be oriented based on the current reference face
    vector<int> isect;
    for (int i = 0; i != num_faces; ++i, isect.clear()) {
      if (!ordered[i]) {
        //std::cout << "Fix ordering of face" << std::endl;
        set_intersection(sorted_fnode_ptrs[ref_face], sorted_fnode_ptrs[ref_face] + num_face_edges[ref_face],
                         sorted_fnode_ptrs[i], sorted_fnode_ptrs[i] + num_face_edges[i], back_inserter(isect));
        if (isect.size() == 2 ) {


          assert(isect.size() > 1); // should be at least two corners to a shared edge
          const auto ref1 = find(face_ptrs[ref_face], face_ptrs[ref_face] + num_face_edges[ref_face], isect[0]);
          const auto ref2 = find(face_ptrs[ref_face], face_ptrs[ref_face] + num_face_edges[ref_face], isect[1]);
          const auto i1 = find(face_ptrs[i], face_ptrs[i] + num_face_edges[i], isect[0]);
          const auto i2 = find(face_ptrs[i], face_ptrs[i] + num_face_edges[i], isect[1]);
          if(!(i1-i2) == 1){
              bool is_ok = false;
              if( i1 == (face_ptrs[i]) ){
                  if( i2 == face_ptrs[i]+num_face_edges[i] ){
                      is_ok = true;
                  }
              }
              if( i2 == face_ptrs[i] ){
                  if( i1 == face_ptrs[i]+num_face_edges[i] ){
                      is_ok = true;
                  }
              }
              if(!is_ok){
                  std::cout << "strange intersection of faces" << std::endl;
                  assert(false);
              }
          }

          const bool orient1 = (ref2 == ref1+1 || ref2 == face_ptrs[ref_face]);  // order is ref1, ref2 if true
          const bool orient2 = (i2 == i1 + 1 || i2 == face_ptrs[i]);            // order is i1, i2 if true
          if (orient1 == orient2)
            // faces have the same orientation.  Flip the orientation of the non-reference face
            reverse(face_ptrs[i], face_ptrs[i] + num_face_edges[i]);

          // at this point, we should be assured that face 'i' is oriented consistently with face 'ref'
          ordered[i] = true;
        }
        else if(!isect.empty()){
            std::cout << "Face-face intersection with not 2 nodes" << std::endl;
        }
      }
    }
    completed[ref_face] = true;
  }
  return result;
}

// ----------------------------------------------------------------------------
bool inward_pointing_normals(const vector<double>& normals,
                             const vector<double>& face_centroids)
// ----------------------------------------------------------------------------
{
  const int N = (int)normals.size()/3;
  const double tol = 1e-9;

  for (int d = 0; d != 3; ++d) {
    // usually, it should be enough to check for d==0, but there may be pathological
    // cases where we would have to resort to other coordinate directions
    int max_ix = 0, min_ix = 0;
    double max_val = face_centroids[d], min_val = max_val;
    for (int i = 0; i != N; ++i) {
      if (face_centroids[3*i+d] < min_val) {
        min_val = face_centroids[3*i+d];
        min_ix = i;
      } else if (face_centroids[3*i+d] > max_val) {
        max_val = face_centroids[3*i+d];
        max_ix = i;
      }
    }
    if (fabs(normals[3 * max_ix + d]) > tol)
      return normals[3 * max_ix + d] < 0;
    else if (fabs(normals[3 * min_ix + d]) > tol)
      return normals[3 * min_ix + d] > 0;
  }

  // if we got here, there is something peculiar about the input data
  throw runtime_error("Unable to determine orientation of normals in element.");
}

// ----------------------------------------------------------------------------
// Compute key parts of cell geometry, including a consistent set of outward
// normals, the cell volume, centroid and a point for which the cell is
// star-shaped (which may or may not be the centroid)
void compute_cell_geometry(const double* points,
                           const int num_points,
                           const int* const faces,
                           const int* const num_face_edges,
                           const int num_faces,
                           vector<double>& outward_normals,
                           vector<double>& face_centroids,
                           array<double, 3>& cell_centroid,
                           array<double, 3>& star_point,
                           double& volume)
// ----------------------------------------------------------------------------
{
  // ensure faces are consistently oriented so that normals all point outwards
  // or inwards
  vector<int> faces2 = consistent_face_ordering(faces, num_face_edges, num_faces);

  outward_normals.resize(num_faces * 3);
  face_centroids.resize(num_faces * 3);
  for (int f = 0, faces_offset = 0; f != num_faces; faces_offset += num_face_edges[f++])
    compute_face_geometry(pick_points<3>(points, &faces2[faces_offset], num_face_edges[f]),
                          &outward_normals[f*3], &face_centroids[f*3]);
  assert(face_centroids.size() == num_faces*3);


  // ensure normals are pointing out of, rather than into, polyhedron
  if (inward_pointing_normals(outward_normals, face_centroids))
    for_each(outward_normals.begin(), outward_normals.end(), [](double& d) {d *= -1;});
  assert(face_centroids.size() == num_faces*3);
  //std::cout << "Size # "<< face_centroids.size() << std::endl;
  // identify a star point (usually, mean_point qualifies, but not necessarily)
  star_point = identify_star_point(point_average<3>(points, num_points),
                                   outward_normals, face_centroids);

  // compute cell centroid and volume
  volume = 0;
  cell_centroid = {0, 0, 0};

  // loop over faces, split each face into triangles, and accumulate the volumes
  // defined by these triangles and the inside point
  for (int f = 0, fpos = 0; f != num_faces; fpos += num_face_edges[f++]) {
    for (auto tri : tessellate_face<3>(pick_points<3>(points, &faces[fpos], num_face_edges[f]))) {
      tri.insert(tri.end(), &star_point[0], &star_point[0] + 3);
      const double tvol = tetrahedron_volume(&tri[0], &tri[3], &tri[6], &tri[9]);
      const auto tet_centroid = point_average<3>(&tri[0], 4);
      volume += tvol;
      for (int d = 0; d != 3; ++d)
        cell_centroid[d] += tvol * tet_centroid[d];
    }
  }

  for (int d = 0; d != 3; ++d)
    cell_centroid[d] /= volume;

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
void potential_gradient_force_3D(const double* const points,
                                 const int num_cells,
                                 const int* const num_cell_faces, // cell faces per cell
                                 const int* const num_face_corners, // corners per cellface
                                 const int* const face_corners,
                                 const double* const field,
                                 vector<double>& fgrad,
                                 vector<tuple<int, int, double>>& div,
                                 bool get_matrix
    )
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
                                                  num_cell_faces[cell], indexing);
    const int num_corners = int(indexing.size());
    const auto corners_loc = pick_points<3>(points, &indexing[0], num_corners);

    const int tot_num_cellface_corners =
      accumulate(&num_face_corners[cur_cface_start],
                 &num_face_corners[cur_cface_start] + num_cell_faces[cell], 0);

    vector<int> faces_loc(tot_num_cellface_corners);
    transform(&face_corners[cur_fcor_start],
              &face_corners[cur_fcor_start] + tot_num_cellface_corners,
              faces_loc.begin(), [&reindex] (const int I) { return reindex.find(I)->second;});

    double volume;                                  // to be computed below
    vector<double> outward_normals, face_centroids; // to be computed below
    array<double, 3> star_point, cell_centroid;     // to be computed below
    compute_cell_geometry(&corners_loc[0], num_corners, &faces_loc[0],
                          &num_face_corners[cur_cface_start], num_cell_faces[cell],
                          outward_normals, face_centroids, cell_centroid, star_point, volume);

    // Since computing q involves dividing by volume, we are really computing q * volume
    // in the below call, by passing the value '1' for volume.  In other words, we are computing
    // (1/2) z_i.
    const auto qv = compute_q_3D(&corners_loc[0], num_corners, &faces_loc[0],
                                 &num_face_corners[cur_cface_start], num_cell_faces[cell],
                                 1, outward_normals);

    // fill in entries in global fgrad vector
    for (int c = 0; c != num_corners; ++c){
        for (int d = 0; d != 3; ++d){
            fgrad[3 * indexing[c] + d] += 2 * field[cell] * qv[3 * c + d];
            if(get_matrix){
                int I = 3 * indexing[c] + d;
                int J = cell;
                double val = 2*qv[3 * c + d];
                div.push_back(tuple<int, int, double> {I, J, val});
            }
        }
    }

    // keep track of our current position in the `face_corners` array
    cur_fcor_start += tot_num_cellface_corners;
    cur_cface_start += num_cell_faces[cell];
  }
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
                        const StabilityChoice stability_choice)
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
        assemble_stiffness_matrix_2D(points, &cell_corners[corner_ixs], ncf, young[c],
                                     poisson[c], stability_choice, loc);

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
        const auto fvals = compute_applied_forces_2D(
            points, neumann_faces[2 * f], neumann_faces[2 * f + 1], neumann_forces[2 * f], neumann_forces[2 * f + 1]);
        b[2 * neumann_faces[2 * f]] += fvals[0];
        b[2 * neumann_faces[2 * f] + 1] += fvals[1];
        b[2 * neumann_faces[2 * f + 1]] += fvals[2];
        b[2 * neumann_faces[2 * f + 1] + 1] += fvals[3];
    }

    // reduce system by eliminating Dirichlet degrees of freedom, and returning
    reduce_system(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);

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
                        const StabilityChoice stability_choice)
// ----------------------------------------------------------------------------
{
    // preliminary computations
    const int tot_num_cell_faces = accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    const int tot_num_face_corners = accumulate(num_face_corners, num_face_corners + tot_num_cell_faces, 0);
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

    cout << "Starting assembly" << endl;
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
                                     loc);
        const int ncv = int(loc_indexing.size()); // number of cell vertices
        for (int i = 0; i != 3 * ncv; ++i) {
            for (int j = 0; j != 3 * ncv; ++j) {
                // using integer division to go from 'i' and 'j' to corresponding
                // local point indices, which are subsequently converted to global indexing
                const int I = 3 * loc_indexing[i / 3] + i % 3;
                const int J = 3 * loc_indexing[j / 3] + j % 3;
                const double val = loc[i * 3 * ncv + j]; // local matrix entry
                A_entries.push_back(tuple<int, int, double> {I, J, val});
                // const auto it = A_entrymap.find(make_pair(I, J));
                // (it != A_entrymap.end()) ? it->second += val :
                //                            A_entrymap[make_pair(I, J)] = val;
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
                             &b[0]);
        // increment local indexing
        fcorners_start += accumulate(&num_face_corners[cf_ix], &num_face_corners[cf_ix + num_cell_faces[c]], 0);
        cf_ix += num_cell_faces[c];
    }
    cout << "Applying forces" << endl;

    // add contribution to right hand side from applied forces.  Values are written
    // directly into the global right-hand-side vector
    compute_applied_forces_3D(
        points, num_face_corners, face_corners, num_neumann_faces, neumann_faces, neumann_forces, &b[0]);

    cout << "Reducing system" << endl;
    reduce_system(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);

    cout << "Finished assembly.  Returning." << endl;

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
                  //const double* const body_force, // 3 * number of cells
                  //const int num_fixed_dofs, // dirichlet
                  //const int* const fixed_dof_ixs, // indices must be sorted
                  //const double* const fixed_dof_values,
                  //const int num_neumann_faces,
                  //const int* const neumann_faces,
                  //const double* const neumann_forces, // 3 * number of neumann faces
                  const std::vector<double>& disp,
                  std::vector<std::array<double,6>>& stress,
                  const StabilityChoice stability_choice,
                  vector<tuple<int, int, double>>& stressmat,
                  bool do_matrix,
                  bool do_stress
    )
// ----------------------------------------------------------------------------
{
    // preliminary computations
    //const int tot_num_cell_faces = accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    //const int tot_num_face_corners = accumulate(num_face_corners, num_face_corners + tot_num_cell_faces, 0);
    //const int num_points = *max_element(face_corners, face_corners + tot_num_face_corners) + 1;

    // assemble full system matrix
    // loop over cells and assemble system matrix
    // vector<int> loc_indexing;
    // vector<double> loc; // use as local 'scratch' vector
    // array<double, 3> centroid; // will contain the centroid for the currently treated cell
    int cf_ix = 0; // index to first cell face for the current cell
    int fcorners_start = 0; // index to first cell corner for current cell

    cout << "Starting assembly" << endl;
    for (int c = 0; c != num_cells; ++c) {
        // computing local stiffness matrix, writing its entries into the global matrix
        calculate_stress_3D_local(points,
                                  &face_corners[fcorners_start],
                                  &num_face_corners[cf_ix],
                                  num_cell_faces[c],
                                  young[c],
                                  poisson[c],
                                  stability_choice,
                                  disp,// global displacement
                                  c,
                                  stress[c],
                                  stressmat,
                                  do_matrix,
                                  do_stress
            );
        fcorners_start += accumulate(&num_face_corners[cf_ix], &num_face_corners[cf_ix + num_cell_faces[c]], 0);
        cf_ix += num_cell_faces[c];
    }
}



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

    // compute all intermediary matrices
    const auto q = compute_q_2D(&corners[0], num_corners);
    const auto Nr = compute_Nr_2D(&corners[0], num_corners);
    const auto Nc = compute_Nc_2D(&corners[0], num_corners);
    const auto Wr = compute_Wr_2D(q);
    const auto Wc = compute_Wc_2D(q);
    const auto D = compute_D_2D(young, poisson);
    //const auto S = compute_S(Nc, D, num_corners, area, 2);
    const auto ImP = compute_ImP(Nr, Nc, Wr, Wc, 2);

    // do the final assembly of matrices, and write result to target
    target.resize(pow(2 * num_corners, 2));
    final_assembly(Wc, D, Nc, ImP, stability_choice, area, num_corners, 2, &target[0]);
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
                             vector<double>& target)
// ----------------------------------------------------------------------------
{
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

    double volume;                                    // computed in 'compute_cell_geometry' below
    vector<double> outward_normals, face_centroids;   // computed in 'compute_cell_geometry' below
    array<double, 3> star_point;                      // computed in 'compute_cell_geometry' below
    compute_cell_geometry(&corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces,
                          outward_normals, face_centroids, centroid, star_point, volume);

    // // compute all intermediary matrices
    const auto q = compute_q_3D(&corners_loc[0], int(corners_loc.size() / 3), &faces_loc[0],
                                num_face_edges, num_faces, volume, outward_normals);
    const auto Nr = compute_Nr_3D(&corners_loc[0], num_corners);
    const auto Nc = compute_Nc_3D(&corners_loc[0], num_corners);
    const auto Wr = compute_Wr_3D(q);
    const auto Wc = compute_Wc_3D(q);
    const auto D = compute_D_3D(young, poisson);
    //const auto S = compute_S(Nc, D, num_corners, volume, 3);
    const auto ImP = compute_ImP(Nr, Nc, Wr, Wc, 3);

    // // do the final assembly of matrices, and write result to target
    target.resize(pow(3 * num_corners, 2));
    final_assembly(Wc, D, Nc, ImP, stability_choice, volume, num_corners, 3, &target[0]);
}

void
calculate_stress_3D_local(const double* const points,
                          const int* const faces,
                          const int* const num_face_edges,
                          const int num_faces,
                          const double young,
                          const double poisson,
                          const StabilityChoice stability_choice,
                          const std::vector<double>& disp,/// global displacement
                          const int cell,
                          std::array<double,6>& stress,
                          vector<tuple<int, int, double>>& stressmat,
                          bool do_matrix,
                          bool do_stress
    )
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

    double volume;                                    // computed in 'compute_cell_geometry' below
    vector<double> outward_normals, face_centroids;   // computed in 'compute_cell_geometry' below
    array<double, 3> star_point;                      // computed in 'compute_cell_geometry' below
    compute_cell_geometry(&corners_loc[0], num_corners, &faces_loc[0], num_face_edges, num_faces,
                          outward_normals, face_centroids, centroid, star_point, volume);

    // // compute all intermediary matrices
    const auto q = compute_q_3D(&corners_loc[0], int(corners_loc.size() / 3), &faces_loc[0],
                                num_face_edges, num_faces, volume, outward_normals);
    //const auto Nr = compute_Nr_3D(&corners_loc[0], num_corners);
    //const auto Nc = compute_Nc_3D(&corners_loc[0], num_corners);
    //const auto Wr = compute_Wr_3D(q);

    const int dim = 3;
    //assert(dim == 2 || dim == 3);
    const int lsdim = (dim == 2) ? 3 : 6; // dimension of "linear strain space"
    const int totdim = dim * num_corners; // total number of unknowns
    const auto Wc = compute_Wc_3D(q);
    std::vector<double> matrix(lsdim*totdim,0.0);
    if(do_stress){
        const auto D = compute_D_3D(young, poisson);
        // compute stiffness matrix components
        const auto DWct = matmul(&D[0], lsdim, lsdim, false, &Wc[0], totdim, lsdim, true);
        matrix = DWct;
    }else{
        for(int i = 0; i < lsdim; ++i){
            for(int j = 0; j < totdim; ++ j){
                matrix[i*totdim+j] =Wc[j*lsdim+i];
            }
        }
    }

    std::vector<double> local_disp;
    local_disp.resize(indexing.size()*3);
    std::vector<int> global_index;
    global_index.resize(indexing.size()*3);
    for(size_t i = 0; i < indexing.size(); ++i){
        for (size_t d = 0; d != 3; ++d){
            local_disp[3*i+ d] = disp[3*indexing[i] + d];
            global_index[3*i+ d] = 3*indexing[i] + d;
        }
    }
    for(int i=0; i < lsdim; ++i){
        for(int j=0; j < totdim; ++j){
            stress[i] += matrix[i*totdim + j]*local_disp[j];
            if(do_matrix){
                int I = lsdim*cell+i;
                int J = global_index[j];
                double val = matrix[i*totdim + j];
                if(i>2 && do_stress){
                    val /=2.0;
                }
                stressmat.push_back(tuple<int, int, double> {I, J, val});
            }
        }
    }
    if(do_stress){
        for(int i = 3; i<lsdim; ++i){
            stress[i] /=2;
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
        cout << setw(12) << ((abs(data[i * stride.first + j * stride.second]) <= zthreshold) ? fixed : scientific)
             << ((abs(data[i * stride.first + j * stride.second]) <= zthreshold) ? setprecision(0) : setprecision(2))
             << (abs(data[i * stride.first + j * stride.second]) <= zthreshold ? 0 : data[i * stride.first + j * stride.second])
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
tetrahedron_volume(const double* const p1, const double* const p2, const double* const p3, const double* const p4)
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
// @@ Might fail in case of strongly non-convex faces
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
face_integral(const double* const corners, const int num_corners, const int dim, const double* const corner_values)
// ----------------------------------------------------------------------------
{
    assert(dim == 2 || dim == 3);

    return dim == 2 ? face_integral_impl<2>(corners, num_corners, corner_values)
                    : face_integral_impl<3>(corners, num_corners, corner_values);
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
