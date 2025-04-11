#ifndef _VEM_TOPOLOGY_HPP
#define _VEM_TOPOLOGY_HPP

#include <array>
#include <vector>
#include <tuple>

namespace vem
{
  using IndexPair = std::array<int, 2>;
  
  // Identify cell faces with individual cells and their faces.  Each entry in
  // the return vector represents a cell face.  Its first element is the index
  // of the cell, and the second element is the index of the face.
  std::vector<IndexPair>
  cellfaces_cells_faces(const int num_cells,
                        const int* const num_cell_faces);

  // identify matching cell faces
  std::tuple<std::vector<IndexPair>, std::vector<int>>
  cellfaces_matching_faces(const int num_cells,
                           const int* const num_cell_faces,
                           const int* const num_face_corners,
                           const int* const face_corners);
  
  std::vector<std::array<double, 3>>
  cellface_centroids(const double* const coords, 
                     const std::vector<int>& cellface_ixs,
                     const int* const num_face_corners,
                     const int* const face_corners);
  
  std::vector<std::tuple<double, IndexPair>>
  mutual_distances(const std::vector<std::array<double, 3>>& points);

  // return vector with indices to the cellfaces considered to be 'top' and 'bottom'
  // for each cell.
  std::vector<std::array<int, 2>>
  identify_top_bottom_faces(const double* const coords,
                            const int num_cells,
                            const int* const num_cell_faces,
                            const int* const num_face_corners,
                            const int* const face_corners);

  
};

#endif
