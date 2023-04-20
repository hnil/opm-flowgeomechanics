#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/grid/yaspgrid.hh>
#include <stdexcept>
#include <iterator> // for ostream_iterator  (debugging purposes)
#include <numeric>
#include "opm/geomech/vem/vem.hpp"

#include <iostream>

using namespace std;

// ----------------------------------------------------------------------------
vector<double> set_body_force(const int num_cells, const string bfcase)
// ----------------------------------------------------------------------------
{
  const double gforce = 0.02;
  vector<double> body_force(3*num_cells, 0);
  if (bfcase == "none") 
    {} // do nothing 
  else if (bfcase == "down") 
    for (int i = 0; i != num_cells; ++i)
      body_force[3*i+(3-1)] = -gforce;
  else if (bfcase == "up") 
    for (int i = 0; i != num_cells; ++i)
      body_force[3*i+(3-1)] = gforce;
  else if (bfcase == "left") 
    for (int i = 0; i != num_cells; ++i)
      body_force[3*i] = -gforce;
  else if (bfcase == "right") 
    for (int i = 0; i != num_cells; ++i)
      body_force[3*i] = gforce;
  else 
    throw invalid_argument("Body force case must be 'none' or 'down'.");
  
  return body_force;
}

// ----------------------------------------------------------------------------
tuple<int, vector<int>, vector<double>> set_dirichlet(const vector<double>& coords,
                                                      const string& dircase)
// ----------------------------------------------------------------------------
{
  if (dircase != "fixed_bottom")
    throw runtime_error("Currently, only 'fixed_bottom' is supported as Dirichlet test case.");

  vector<int> fixed_dof_ixs;
  for (int i = 0; i != int(coords.size()); i += 3)
    if (coords[i+2] == 0) // z-coordinate is at bottom
      fixed_dof_ixs.insert(fixed_dof_ixs.end(), {i, i+1, i+2});

  const int num_fixed_dofs = (int)fixed_dof_ixs.size();
  const vector<double> fixed_dof_values(num_fixed_dofs, 0.0);

  return {num_fixed_dofs, fixed_dof_ixs, fixed_dof_values};
}

// ============================================================================
int main(int argc, char** argv)
// ============================================================================
{
  Dune::MPIHelper::instance(argc,argv); // to prevent compiler from complaining about MPI ??

  const int xres           = (argc > 1)  ? atoi(argv[1])    : 3;
  const int yres           = (argc > 2)  ? atoi(argv[2])    : 3;
  const int zres           = (argc > 3)  ? atoi(argv[3])    : 3;
  const double Lx          = (argc > 4)  ? atof(argv[4])    : 4.0;
  const double Ly          = (argc > 5)  ? atof(argv[5])    : 4.0;
  const double Lz          = (argc > 6)  ? atof(argv[6])    : 4.0;
  const double poisson_val = (argc > 7)  ? atof(argv[7])    : 0.25;
  const double young_val   = (argc > 8)  ? atof(argv[8])    : 1.0;
  const string bfcase      = (argc > 9)  ? string(argv[9])  : string("down");
  const string dircase     = (argc > 10) ? string(argv[10]) : string("fixed_bottom");
  
  const Dune::YaspGrid<3> grid({Lx, Ly, Lz}, {xres, yres, zres});

  const auto gv = grid.leafGridView();
  const auto& ixset = gv.indexSet();

  // make global point coordinate vector
  vector<double> coords;
  for (const auto& v : vertices(gv)) {
    const auto c = v.geometry().corner(0);
    coords.insert(coords.end(), c.begin(), c.end());
  }

  const int num_cells = gv.size(0); // entities of codim 0

  // count cell faces
  vector<int> num_cell_faces;
  for (const auto& c : elements(gv)) 
    num_cell_faces.push_back(Dune::subEntities(c, Dune::Codim<1>{}).size());
  const int tot_num_cfaces = accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

  // count face corners
  vector<int> num_face_corners;
  for (const auto& c : elements(gv))
    for (const auto& f : Dune::subEntities(c, Dune::Codim<1>{}))
      num_face_corners.push_back(f.geometry().corners());

  // establish all face corners
  vector<int> face_corners;
  for (const auto& c : elements(gv))  
    for (const auto& f: Dune::subEntities(c, Dune::Codim<1>{})) 
      for (int i = 0, f_ix = 0;
           f_ix != tot_num_cfaces && i != num_face_corners[f_ix];
           ++i, f_ix += (i == num_face_corners[f_ix]))
        face_corners.push_back(ixset.subIndex(f, i, 3));

  // correct order of nodes in each face, to ensure they are mentioned in
  // clockwise or counterclockwise order. @@ This is a hack that might only work
  // for 4-faces!
  for (int i = 0, j = 0; i != (int)num_face_corners.size(); j += num_face_corners[i++])
    swap(face_corners[j], face_corners[j+1]);

  // linear elastic parameters
  const vector<double> young(num_cells, young_val);
  const vector<double> poisson (num_cells, poisson_val);
  
  // body force
  const vector<double> body_force = set_body_force(num_cells, bfcase);

  // dirichlet boundary conditions
  const auto dirichlet = set_dirichlet(coords, dircase);
  const int num_fixed_dofs = get<0>(dirichlet);
  const vector<int>& fixed_dof_ixs = get<1>(dirichlet);
  const vector<double>& fixed_dof_values = get<2>(dirichlet);

  // neumann boundary conditions 
  const int num_neumann_faces = 0;
  
  // assemble the mechanical system
  vector<tuple<int, int, double>> A_entries;
  vector<double> b;
  const int numdof =
    vem::assemble_mech_system_3D(&coords[0], num_cells, &num_cell_faces[0], &num_face_corners[0],
                                 &face_corners[0], &young[0], &poisson[0], &body_force[0],
                                 num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                                 num_neumann_faces, nullptr, nullptr,
                                 A_entries, b);
  
  // view matrix

  cout << "System matrix is: " << endl;
  const auto M = vem::sparse2full(A_entries, numdof, numdof);
  vem::matprint(&M[0], numdof, numdof, false, 1e-10);
  
  cout << endl << endl;

  cout << "Right-hand side is: " << endl;
  copy(b.begin(), b.end(), ostream_iterator<double>(cout, " "));
  cout << endl << endl;

  cout << "Number of cells: " << num_cells << endl;
  
  return 0;
};
