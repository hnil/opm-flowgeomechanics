#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/grid/yaspgrid.hh>
#include <iterator> // for ostream_iterator
#include "opm/geomech/vem/vem.hpp"

#include <iostream>

using namespace std;

typedef Dune::YaspGrid<3>::LeafGridView GridView;
typedef Dune::Entity<0, 3, const Dune::YaspGrid<3>, Dune::YaspEntity> Cell;
typedef Dune::Entity<1, 3, const Dune::YaspGrid<3>, Dune::YaspEntity> Face;
typedef Dune::IndexSet<const Dune::YaspGrid<3>,
                       Dune::YaspIndexSet<const Dune::YaspGrid<3>, true>,
                       unsigned int, std::vector<Dune::GeometryType> > IndexSet;


// ----------------------------------------------------------------------------
vector<int> make_nodemap(const IndexSet& ixset, const Cell& cell)
// ----------------------------------------------------------------------------
{
  const int N = Dune::subEntities(cell, Dune::Codim<3>{}).size();
  vector<int> result;

  for (int i = 0; i != N; ++i)
    result.push_back(ixset.subIndex(cell, i, 3));

  return result;
}

// ----------------------------------------------------------------------------
vector<double> corner_coords(const Cell& cell)
// ----------------------------------------------------------------------------  
{
  const int N = Dune::subEntities(cell, Dune::Codim<3>{}).size();
  const auto cellgeo = cell.geometry();
  vector<double> result;
  
  for (int i = 0; i != N; ++i) {
    const auto cor = cellgeo.corner(i);
    result.insert(result.end(), cor.begin(), cor.end());
  }
  return result;
}

// ----------------------------------------------------------------------------
int main(int argc, char** argv)
// ----------------------------------------------------------------------------  
{
  Dune::MPIHelper::instance(argc,argv); // to prevent compiler from complaining about MPI ??
  
  using Grid = Dune::YaspGrid<3>;
  const Grid grid({4.0, 4.0, 4.0} , // physical dimensions
                  {2, 2, 2});       // resolution

  const auto gv = grid.leafGridView();
  const auto& ixset = gv.indexSet();

  const double young = 1;
  const double poisson = 0.25;
  
  for (const auto& c : elements(gv)) {
    
    // make local-to-global map of cell corners
    const vector<int> nodemap = make_nodemap(ixset, c);

    // define lookup function from local face corner indexing to local cell corner indexing
    const auto nodemap_lookup = [&nodemap, &ixset] (const Face& f, int i) {
      return int(find(nodemap.begin(), nodemap.end(), ixset.subIndex(f, i, 3))
                 - nodemap.begin());
    };
    
    // create 'faces' vector
    vector<int> faces;
    for(const auto& f : Dune::subEntities(c, Dune::Codim<1>{})) 
      for (int i = 0; i != f.geometry().corners(); ++i)
        faces.push_back(nodemap_lookup(f, i));

    // create 'num_face_corners' vector
    vector<int> num_face_corners;
    for(const auto& f : Dune::subEntities(c, Dune::Codim<1>{}))
      num_face_corners.push_back(f.geometry().corners());


    // correct order of nodes in each face, to ensure they are mentioned in
    // clockwise or counterclockwise order. @@ This is a hack that might only work
    // for 4-faces!
    for (int i = 0, j = 0; i != (int)num_face_corners.size(); j += num_face_corners[i++])
      swap(faces[j], faces[j+1]);

    
    // create local vector of point coordinates
    const vector<double> coords = corner_coords(c);

    // assemble element stiffness matrix
    array<double, 3> cell_centroid;
    vector<int> indexing;
    vector<double> target(24*24);
    const vem::StabilityChoice stability_choice = vem::D_RECIPE;
    vem::assemble_stiffness_matrix_3D(&coords[0], &faces[0], &num_face_corners[0],
                                      (int)num_face_corners.size(),
                                      young, poisson,
                                      stability_choice,
                                      cell_centroid, indexing, target);

    cout << "Centroid: " << cell_centroid[0] << " " << cell_centroid[1] << " " << cell_centroid[2] << endl;
    cout << "Indexing: " << endl;
    copy(indexing.begin(), indexing.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
    cout << "Matrix: " << endl;
    vem::matprint(&target[0], indexing.size()*3, indexing.size()*3, false, 1e-9);
    cout << endl;
    
    cout << "Faces are: " << endl;
    copy(faces.begin(), faces.end(), ostream_iterator<int>(cout, " "));
    cout << endl;

    cout << "Num face corners are: " << endl;
    copy(num_face_corners.begin(), num_face_corners.end(), ostream_iterator<int>(cout, " " ));
    cout << endl;

    cout << "Coordinates are: " << endl;
    for (std::size_t i = 0; i < coords.size(); i += 3)
      cout << coords[i] << " " << coords[i+1] << " " << coords[i+2] << endl;
    
  }
  
  return 0;
};
