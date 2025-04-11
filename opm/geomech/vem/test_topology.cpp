#include "topology.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

using namespace std;

int main(int varnum, char** vararg)
{
  ifstream is(vararg[1]);
  if (!is) {
    cerr << "Error opening file: " << vararg[1] << endl;
    return 1;
  }
  // read number of cells and number of cell faces per cell
  int num_cells; is >> num_cells;
  vector<int> num_cell_faces(num_cells);
  for (int i = 0; i < num_cells; ++i) { is >> num_cell_faces[i]; }
  const int total_num_cellfaces = 
    accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

  // read nodes
  int num_nodes; is >> num_nodes;
  vector<double> nodes(num_nodes * 3);
  for (int i = 0; i < num_nodes * 3; ++i) { is >> nodes[i]; }

  // read number of corners per cell face
  vector<int> num_face_corners(total_num_cellfaces);
  for (int i = 0; i < total_num_cellfaces; ++i) { is >> num_face_corners[i]; }

  // read face corners
  const int total_num_face_corners =
    accumulate(num_face_corners.begin(), num_face_corners.end(), 0);
  vector<int> face_corners(total_num_face_corners);
  for (int i = 0; i < total_num_face_corners; ++i) { is >> face_corners[i]; }

  // mention all cell face indices per default
  vector<size_t> cellface_ixs(total_num_cellfaces);
  for (size_t i = 0; i < total_num_cellfaces; ++i) {
    cellface_ixs[i] = i;
  }
  
  //const vector<size_t> cellface_ixs = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  
    // test compute centroids
  // const auto centroids =
  //   vem::cellface_centroids(&nodes[0],
  //                           cellface_ixs,
  //                           &num_face_corners[0],
  //                           &face_corners[0]);

  // for (const auto c : centroids)
  //   cout << c[0] << " " << c[1] << " " << c[2] << endl;

  // // test compute distances between centroids
  // auto dist = vem::mutual_distances(centroids);
  // sort(dist.begin(), dist.end(), [](tuple<double, vem::IndexPair>& a,
  //                                   tuple<double, vem::IndexPair>& b)
  // { return get<0>(a) < get<0>(b);});

  // for (const auto e : dist)
  //   cout << get<0>(e) << ": " << get<1>(e)[0] << " " << get<1>(e)[1] << endl;

  // test matching of cell faces

  // ----------------------------------------------------------------------------
  
  // const auto matching =
  //   vem::cellfaces_matching_faces(num_cells,
  //                                 &num_cell_faces[0],
  //                                 &num_face_corners[0],
  //                                 &face_corners[0]);

  // cout << "----- Matching cellfaces: ---" << endl;
  // for (const auto e : get<0>(matching))
  //   cout << e[0] << " " << e[1] << endl;

  // cout << "----- Outer cellfaces: ---" << endl;
  // for (const auto e : get<1>(matching))
  //   cout << e << " ";
  // cout << endl;
                                                 
  // cout << "----- Outer cellfaces coords: ---" << endl;
  // vector<size_t> face_corner_start(total_num_cellfaces);
  // face_corner_start[0] = 0;
  // for (size_t i = 1; i < total_num_cellfaces; ++i) {
  //   face_corner_start[i] = face_corner_start[i-1] + num_face_corners[i-1];
  // }
    
  // for (const auto e : get<1>(matching)) {
  //   const int nc = num_face_corners[e];
  //   const int start = face_corner_start[e];
  //   for (int i = 0; i < nc; ++i) {
  //     // get corner coordinates and output them
  //     cout << nodes[3 * face_corners[start + i]] << " "
  //          << nodes[3 * face_corners[start + i] + 1] << " "
  //          << nodes[3 * face_corners[start + i] + 2] << endl;
  //   }
  // }
  // cout << endl;
  // cout << "---- face corner numbers" << endl;
  // for (const auto e : get<1>(matching)) 
  //   cout << num_face_corners[e] << " ";
  // cout << endl;

// ----------------------------------------------------------------------------

  const auto tbfaces = 
    vem::identify_top_bottom_faces(&nodes[0],
                                   num_cells,
                                   &num_cell_faces[0],
                                   &num_face_corners[0],
                                   &face_corners[0]);

  //cout << tbfaces.size() << endl;

  // print all tbfaces to screen
  cout << "----- Top bottom faces: ---" << endl;
  for (const auto e : tbfaces) {
    cout << e[0] << " " << e[1] << endl;
  }
  
  return 0;
}
