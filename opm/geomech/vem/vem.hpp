#ifndef _VEM_HPP
#define _VEM_HPP

#include <array>
#include <tuple>
#include <vector>
namespace vem
{
  enum StabilityChoice {
    SIMPLE = 1,     // the basic stability term described in (Gain et al., 2014)
    HARMONIC = 2,   // modified stability term more robust for high aspect ratios (Andersen et al., 2017)
    D_RECIPE = 3    // another modified stability based on the diagonal of the consistency matrix
  };
  
// ============================================================================
// == Main functions for computing element stiffness matrices and assembling ==
// ==                      linear systems in 2D and 3D                       ==
// ============================================================================

// -----------------------------------------------------------------------------------------------
//  Assembles the global linear system for a 2D grid of virutal elements, with
//  Dirichlet and Neumann boundary conditions. The resulting system is represented
//  as the entries of a sparse matrix and a right-hand-side vector.
//
//  @param points             Pointer to an array of doubles, containing the x,y-coordinates
//                            of all points in the mesh.
//  @param num_cells          The number of cells in the mesh.
//  @param num_cell_faces     Pointer to an array of ints, containing the number of faces for
//                            each cell in the mesh.
//  @param cell_corners       Pointer to an array of ints, containing the indices of the
//                            corners of each cell in the `points` array. The corners for
//                            cell i are located at indices cell_corners[K],...,cell_corners[L]-1.
//                            where K is the sum of `num_cell_faces[0]` through `num_cell_faces[i-1]`
//                            and L equals `num_cell_faces[i]`.
//  @param young              Pointer to an array of doubles, containing the Young's modulus
//                            for each cell.
//  @param poisson            Pointer to an array of doubles, containing the Poisson ratio
//                            for each cell.
//  @param body_force         Pointer to an array of doubles, containing the body force
//                            (i.e., force per unit volume) for each cell in the mesh. The
//                            force for cell i is located at indices 2*i and 2*i+1.
//  @param num_fixed_dofs     The number of degrees of freedom restricted by Dirichlet conditions.
//                            Note that each node in the mesh has 2 degrees of freedom (x and y).
//  @param fixed_dof_ixs      Pointer to an array of ints, containing the indices of the
//                            fixed degrees of freedom (i.e., nodes coordinate components with
//                            Dirichlet boundary conditions) in the global system. The indices
//                            must be sorted in ascending order.
//  @param fixed_dof_values   Pointer to an array of doubles, containing the prescribed
//                            values for the fixed degrees of freedom.
//  @param num_neumann_faces  The number of Neumann boundary conditions.
//  @param neumann_faces      Pointer to an array of ints, containing the indices of the
//                            faces with Neumann boundary conditions.
//  @param neumann_forces     Pointer to an array of doubles, containing the force values
//                            for the Neumann boundary conditions. The force for face i is
//                            located at indices 2*i and 2*i+1.
//  @param A_entries          A vector of tuples, each containing the (i,j) indices and value
//                            for a non-zero entry in the sparse matrix.  There may be multiple
//                            entries for a given (i,j) coordinate, in which case these needs to
//                            be added up when constructing the actual, sparse matrix.
//                            (E.g. using the Eigen package, this can be achieved using the
//                            `setFromTriplet` member function of the sparse matrix class).
//  @param b                  The right-hand side vector of the linear system.
//
//  @return                   The total number of unknowns, i.e. the number of rows and columns
//                            in the final matrix (and the length of the right-hand side vector `b`).
//                            This number is equal to 3 times the number of points, minus the
//                            number of fixed degreees of freedom, `num_fixed_dofs`.
int assemble_mech_system_2D(const double* const points,
                            const int num_cells,
                            const int* const num_cell_faces,
                            const int* const cell_corners,
                            const double* const young,
                            const double* const poisson,
                            const double* const body_force, // 2 x number of cells
                            const int num_fixed_dofs, // dirichlet
                            const int* const fixed_dof_ixs, // indices must be sorted
                            const double* const fixed_dof_values,
                            const int num_neumann_faces, // neumann
                            const int* const neumann_faces,
                            const double* const neumann_forces, // 2 * number of neumann faces
                            std::vector<std::tuple<int, int, double>>& A_entries,
                            std::vector<double>& b,
                            const StabilityChoice stability_choice,
                            bool reduce_system);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Assembles the global linear system for a 3D grid of virtual elements, with
// Dirichlet and Neumann boundary conditions. The resulting system is represented
// as the entries of a sparse matrix and a right-hand-side vector.
//
// @param points             Pointer to an array of doubles, containing the x,y,z-coordinates
//                           of all points in the mesh.
// @param num_cells          The number of cells in the mesh.
// @param num_cell_faces     Pointer to an array of ints, of length `num_cells`, containing the
//                           number of faces for each cell in the mesh.
// @param num_face_corners   Pointer to an array of ints, containing the number of corners
//                           for each face in the mesh.  The length of this array equals
//                           the total number of cell-faces, i.e. the sum of
//                           `num_cell_faces[0]...[num_cells-1]`.  The number of face corners
//                           for the nth face of the kth cell is found at entry l, where
//                           l = sum(num_cell_faces[0]...[k-1]) + n.
// @param face_corners       Pointer to an array of ints, containing the indices of the
//                           corners of each face for each cell in the grid in the `points` array.
//                           The indices of corners for the nth face of the kth cell is found
//                           in `face_corners[m+n... m+n+p-1]`, where
//                           m = sum(num_face_corners[0 .. l-1]), and p = num_face_corners[l].
//                           Here, `l` is defined as above in the documentation of `num_face_corners`.
// @param young              Pointer to an array of doubles, containing the Young's modulus
//                           for each cell.
// @param poisson            Pointer to an array of doubles, containing the Poisson ratio
//                           for each cell.
// @param body_force         Pointer to an array of doubles, containing the body force
//                           (i.e., force per unit volume) for each cell in the mesh. The
//                           force for cell i is located at indices 3*i, 3*i+1 and 3*i+2.
// @param num_fixed_dofs     The number of degrees of freedom restricted by Dirichlet conditions.
//                           Note that each node in the mesh has 3 degrees of freedom (x, y and z).
// @param fixed_dof_ixs      Pointer to an array of ints, containing the indices of the
//                           fixed degrees of freedom (i.e., node coordinate components with
//                           Dirichlet boundary conditions) in the global system. The indices
//                           must be sorted in ascending order.
// @param fixed_dof_values   Pointer to an array of doubles, containing the prescribed
//                           values for the fixed degrees of freedom.
// @param num_neumann_faces  The number of Neumann boundary conditions.
// @param neumann_faces      Pointer to an array of ints, containing the indices of the
//                           faces with Neumann boundary conditions.
// @param neumann_forces     Pointer to an array of doubles, containing the force values
//                           for the Neumann boundary conditions. The force for face i is
//                           located at indices 3*i, 3*i+1 and 3*i+2.
// @param A_entries          A vector of tuples, each containing the (i,j) indices and value
//                           for a non-zero entry in the sparse matrix.  There may be multiple
//                           entries for a given (i,j) coordinate, in which case these need to
//                           be added up when constructing the actual sparse matrix.
//                           (E.g. using the Eigen package, this can be achieved using the
//                           `setFromTriplets` member function of the sparse matrix class).
// @param b                  The right-hand side vector of the linear system.
//
// @return                   The total number of unknowns, i.e., the number of rows and columns
//                           in the final matrix (and the length of the right-hand side vector `b`).
//                           This number is equal to 3 times the number of points, minus the
//                           number of fixed degrees of freedom, `num_fixed_dofs`.
int assemble_mech_system_3D(const double* const points,
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
                            bool reduce_system);
// -----------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------
// Assembles the stiffness matrix for a single 2D grid element, given the coordinates of
// its corners, its Young's modulus, and Poisson ratio. The assembled matrix is
// written to the target vector in row-major order.
//
// @param points        Pointer to an array of doubles, containing the x,y-coordinates
//                      of all points in the mesh.
// @param corner_ixs    Pointer to an array of ints, containing the indices of the
//                      corners of the element in the `points` array.
// @param num_corners   The number of corners in the element.
// @param young         The Young's modulus of the material.
// @param poisson       The Poisson ratio of the material.
// @param target        The vector that the assembled matrix will be written to, in
//                      row-major order. Its size will become (2*num_corners)^2.
void assemble_stiffness_matrix_2D(const double* const points,
                                  const int* const corner_ixs,
                                  const int num_corners,
                                  const double young,
                                  const double poisson,
                                  const StabilityChoice stability_choice,
                                  std::vector<double>& target);
// -----------------------------------------------------------------------------------------------




// -----------------------------------------------------------------------------------------------
// Assembles the stiffness matrix for a single 3D grid polyhedral element.
// The resulting stiffness matrix is stored in the output parameter `target`.
//
// @param points Pointer to an array of size 3*N, where N is the number of
//               points in the global mesh. The array stores the x, y, and z
//               coordinates of each point.
// @param faces Pointer to an array of size K, where K is the total number of
//              entries in the array of faces (i.e., the sum of all entries in
//              `num_face_edges`). The array stores the global indices of the
//              points that define each face.
// @param num_face_edges Pointer to an array of size `num_faces`, where `num_faces`
//                       is the number of polygonal faces in the element. The
//                       array stores the number of edges (or vertices) for
//                       each face.
// @param num_faces Number of polygonal faces in the element (and the number
//                  of entries in the array pointed to by `num_face_edges`).
// @param young Young's modulus of elasticity for the material.
// @param poisson Poisson's ratio for the material.
// @param centroid Output parameter that stores the geometric centroid of the element.
// @param indexing Output parameter that stores the local vertex indices in
//                 the global 'points' vector.
// @param target Output parameter that stores the resulting stiffness matrix
//               as a vector of size (3*num_corners)^2, where num_corners is the
//               number of corners in the element.  Matrix elements are stored
//               in row-major order (although it does not really matter, since
//               the matrix is symmetric).
void assemble_stiffness_matrix_3D(const double* const points,
                                  const int* const faces,
                                  const int* const num_face_edges,
                                  const int num_faces,
                                  const double young,
                                  const double poisson,
                                  const StabilityChoice stability_choice,
                                  std::array<double, 3>& centroid,
                                  std::vector<int>& indexing,
                                  std::vector<double>& target);
// -----------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------
// Computes the node-wise body force resulting from the gradient of a potential (e.g.
// for poromechanics, this would be the gradient of pressure.  The computed value represents
// the gradient of the potential in the node, scaled by the volumes of the neighboring cells.
//
// @param points Pointer to an array of size 3*N, where N is the number of
//               points in the global mesh. The array stores the x, y, and z
//               coordinates of each point.
// @param num_cells          The number of cells in the mesh.
// @param num_cell_faces     Pointer to an array of ints, of length `num_cells`, containing the
//                           number of faces for each cell in the mesh.
// @param num_face_corners   Pointer to an array of ints, containing the number of corners
//                           for each face in the mesh.  The length of this array equals
//                           the total number of cell-faces, i.e. the sum of
//                           `num_cell_faces[0]...[num_cells-1]`.  The number of face corners
//                           for the nth face of the kth cell is found at entry l, where
//                           l = sum(num_cell_faces[0]...[k-1]) + n.
// @param face_corners       Pointer to an array of ints, containing the indices of the
//                           corners of each face for each cell in the grid in the `points` array.
//                           The indices of corners for the nth face of the kth cell is found
//                           in `face_corners[m+n... m+n+p-1]`, where
//                           m = sum(num_face_corners[0 .. l-1]), and p = num_face_corners[l].
//                           Here, `l` is defined as above in the documentation of `num_face_corners`.
// @param field              The field representing the potential from which we want to compute the
//                           force.
// @param fgrad              Output parameter that stores the resulting, node-wise force field
//                           defined in each node by the local gradient of the field, scaled by
//                           the cell volumes of the cells neighboring the node.
void potential_gradient_force_3D(const double* const points,
                                 const int num_cells,
                                 const int* const num_cell_faces, // cell faces per cell
                                 const int* const num_face_corners, // corners per cellface
                                 const int* const face_corners,
                                 const double* const field,
                                 std::vector<double>& fgrad,
                                 std::vector<std::tuple<int, int, double>>& div,
                                 bool get_matrix);    

void
compute_stress_3D(const double* const points,
                  const int num_cells,
                  const int* const num_cell_faces, // cell faces per cell
                  const int* const num_face_corners, // corners per face
                  const int* const face_corners,
                  const double* const young,
                  const double* const poisson,
                  // const double* const body_force, // 3 * number of cells
                  // const int num_fixed_dofs, // dirichlet
                  // const int* const fixed_dof_ixs, // indices must be sorted
                  // const double* const fixed_dof_values,
                  // const int num_neumann_faces,
                  // const int* const neumann_faces,
                  // const double* const neumann_forces, // 3 * number of neumann faces
                  const std::vector<double>& disp,
                  std::vector<std::array<double,6>>& stress,
                  std::vector<std::tuple<int, int, double>>& stressmat,
                  bool do_matrix,
                  bool do_stress
    );    


void
calculate_stress_3D_local(const double* const points,
                          const int* const faces,
                          const int* const num_face_edges,
                          const int num_faces,
                          const double young,
                          const double poisson,
                          const std::vector<double>& disp,
                          const int cell,
                          std::array<double,6>& stress,
                          std::vector<std::tuple<int, int, double>>& stressmat,
                          bool do_matrix,
                          bool do_stress
    );    
// ============================================================================
// ============ Various utility functions for geometry computations ===========
// ============================================================================



// -----------------------------------------------------------------------------------------------
// Calculates the integral of a function over a face in 2D or 3D space, with an arbitrary
// number of corners, using a quadrature rule with a fixed number of points.
//
// @param corners        Pointer to an array of doubles, containing the (x,y) or (x,y,z)-coordinates
//                       of the corners of the face. The corners must be ordered consecutively
//                       around the perimeter of the face.
// @param num_corners    The number of corners of the face.
// @param dim            The dimension of the space in which the face is embedded (2 or 3).
// @param corner_values  (optional) Pointer to an array of doubles, containing the value of the
//                       integrand at each corner point. If this argument is nullptr, the integrand
//                       is assumed to be 1.
//
// @return               The integral of the integrand over the face, calculated using the simple,
//                       first order quadrature rule given by (Gain, 2014)
//                       DOI:10.1016/j.cma.2014.05.005
double face_integral(const double* const corners,
                     const int num_corners,
                     const int dim,
                     const double* const corner_values = nullptr);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Computes the volume of a 3D element defined by its vertices and faces.
//
// @param points           A pointer to an array of doubles containing the x,y,z-coordinates
//                         of all points in the element.
// @param num_points       The number of points in the element.
// @param faces            A pointer to an array of ints containing the indices of the
//                         corners that make up each face of the element. The corners of
//                         face i are located at indices faces[K],...,faces[L]-1, where
//                         K is the sum of `num_face_edges[0]` through `num_face_edges[i-1]`
//                         and L equals `num_face_edges[i]`.
// @param num_face_edges   A pointer to an array of ints containing the number of edges
//                         for each face in the element.
// @param num_faces        The number of faces in the element.
//
// @return                 The volume of the element.
double element_volume_3D(const double* const points,
                         const int num_points,
                         const int* const faces,
                         const int* const num_face_edges,
                         const int num_faces);
// -----------------------------------------------------------------------------------------------




// -----------------------------------------------------------------------------------------------
// Computes the volume of a 2D polygon, given the coordinates of its corner points. The points
// should be ordered in counterclockwise order. The function assumes that the polygon is
// non-self-intersecting and does not contain holes.
//
// @param corners Pointer to an array of doubles, containing the x,y-coordinates
//                of the corner points of the polygon. The corners should be ordered
//                in counterclockwise direction.
// @param num_corners The number of corner points in the polygon.
//
// @return The area of the polygon
double element_volume_2D(const double* const corners, const int num_faces);
// -----------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------
// Computes the geometric centroid of a polygon in 2D.
//
// @param corners Pointer to an array of doubles, containing the x,y-coordinates
//                of the corner points of the polygon. The corners should be ordered
//                in clockwise or counterclockwise direction.
// @param num_points: The number of points in the `points` array.
//
// @return: A std::array<double, 2> containing the coordinates of the centroid.
std::array<double, 2> centroid_2D(const double* const points, const int num_points);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Computes the geometric centroid of a 2D planar polygon embedded in 3D space
//
// @param corners Pointer to an array of doubles, containing the x,y,z-coordinates
//                of the corner points of the polygon. The corners should be ordered
//                in clockwise or counterclockwise direction.
// @param num_points: The number of points in the `points` array.
//
// @return: A std::array<double, 3> containing the 3D coordinates of the centroid.
std::array<double, 3> centroid_2D_3D(const double* const points, const int num_points);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Computes the 3D geometic centroid of a polyhedral cell defined by the given points and faces.
//
// @param points The array of 3D coordinates of the cell vertices.
// @param num_points The number of points in the array.
// @param faces The array of indices that define the faces of the cell.
// @param num_face_edges The array of numbers of edges for each face.
// @param num_faces The number of faces in the array.
// @return An array containing the 3D coordinates of the centroid of the cell.
std::array<double, 3> centroid_3D(const double* const points,
                                  const int num_points,
                                  const int* const faces,
                                  const int* const num_face_edges,
                                  const int num_faces);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Calculates the volume of a tetrahedron with the given corner points.
//
// @param p1 Pointer to the first corner point in 3D Cartesian space.
// @param p2 Pointer to the second corner point in 3D Cartesian space.
// @param p3 Pointer to the third corner point in 3D Cartesian space.
// @param p4 Pointer to the fourth corner point in 3D Cartesian space.
// @return The volume of the tetrahedron defined by the four corner points.
double
tetrahedron_volume(const double* const p1, const double* const p2, const double* const p3, const double* const p4);
// -----------------------------------------------------------------------------------------------


// ============================================================================
// ================== Other utility and diagnostic functions ==================
// ============================================================================



// -----------------------------------------------------------------------------------------------
// Returns a vector of selected points in 2D space, given an array of all points
// and an array of point indices.
//
// @param pts        The array of all points. Each point consists of two doubles,
//                   representing its x and y coordinates.
// @param p_ixs      The array of point indices, indicating which points from `pts`
//                   to include in the output vector.
// @param num_points The number of points to pick from `pts` based on `p_ixs`.
// @return           A vector of 2D points represented as doubles, where each
//                   consecutive pair of values correspond to the x and y
//                   coordinates of a point.
std::vector<double> pick_points_2D(const double* pts, const int* const p_ixs, int num_points);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Returns a vector of selected points in 3D space, given an array of all points
// and an array of point indices.
//
// @param pts        The array of all points. Each point consists of two doubles,
//                   representing its x, y and z coordinates.
// @param p_ixs      The array of point indices, indicating which points from
//                   `pts` to include in the output vector.
// @param num_points The number of points to pick from `pts` based on `p_ixs`.
// @return           A vector of 2D points represented as doubles, where each
//                   consecutive pair of values correspond to the x and y coordinates
//                   of a point.
std::vector<double> pick_points_3D(const double* pts, const int* const p_ixs, int num_points);
// -----------------------------------------------------------------------------------------------



// -----------------------------------------------------------------------------------------------
// Print the contents of a matrix to the console.  Elements in matrix are
// assumed to be stored in row-major order.
//
// @param data Pointer to the first element of the matrix.
// @param r The number of rows in the matrix.
// @param c The number of columns in the matrix.
// @param transposed If true, the matrix is printed transposed.
// @param ztreshold Values with absolute magnitude less than this value, will be rounded to zero
void matprint(const double* data, const int r, const int c, bool transposed, const double zthreshold=0);
// -----------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------
// Converts a sparse matrix represented as a vector of non-zero elements to a full matrix.
//
// @param nz A vector of tuples, where each tuple represents a non-zero element in the sparse matrix.
//           The first element of the tuple is the row index, the second element is the column index,
//           and the third element is the value.
// @param r The number of rows in the full matrix.
// @param c The number of columns in the full matrix.
//
// @return A vector representing the full matrix. The vector contains the elements of the matrix
//         in row-major order.
std::vector<double> sparse2full(const std::vector<std::tuple<int, int, double>>& nz, const int r, const int c);
// -----------------------------------------------------------------------------------------------
}; // end namespace vem

#endif
