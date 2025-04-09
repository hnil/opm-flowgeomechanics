#include "topology.hpp"
#include <map>
#include <algorithm>
#include <cmath>

using namespace std;

namespace {
  // ----------------------------------------------------------------------------
  struct FaceCorners
  // ----------------------------------------------------------------------------  
  {
    // Data
    vector<int> face_corners;
    
    // Methods
    FaceCorners(const int* start, const int num) : 
      face_corners(start, start + num)
    {
      sort(face_corners.begin(), face_corners.end());
    }

    size_t numCorners() const { return face_corners.size(); }

    // 'smaller than'-operator between two face corners
    bool operator<(const FaceCorners& other) const
    {
      if (numCorners() != other.numCorners())
        return numCorners() < other.numCorners();
      for (size_t i = 0; i != numCorners(); ++i)
        if (face_corners[i] != other.face_corners[i])
          return face_corners[i] < other.face_corners[i];
      return false;
    }
  };
};

namespace vem {
  // ----------------------------------------------------------------------------
  vector<IndexPair>
  cellfaces_cells_faces(const int num_cells,
                        const int* const num_cell_faces)
  // ----------------------------------------------------------------------------
  {
    vector<IndexPair> result;
      
    for (size_t cell = 0; cell != num_cells; ++cell)
      for (size_t face = 0; face != num_cell_faces[cell]; ++face)
        result.push_back({cell, face});
    return result;
  }


  // ----------------------------------------------------------------------------
  tuple<vector<IndexPair>, vector<size_t>>
  cellfaces_matching_faces(const int num_cells,
                           const int* const num_cell_faces,
                           const int* const num_face_corners,
                           const int* const face_corners)
  // ----------------------------------------------------------------------------    
  {
    tuple<vector<IndexPair>, vector<size_t>> result;
    map<FaceCorners, size_t> fc2cellface;
    
    const int* nfc_ptr = num_face_corners;
    const int* fc_ptr = face_corners;
    size_t cellface_ix = 0;

    // identify pairs of cellfaces with identical corners
    for (size_t cell = 0; cell != num_cells; ++cell) {
      for (size_t face = 0; face != num_cell_faces[cell]; ++face, ++cellface_ix) {
        const FaceCorners fc(fc_ptr, *nfc_ptr);
        fc_ptr += *nfc_ptr++; // increment fc_pointer and nfc_ptr
        auto it = fc2cellface.find(fc);
        if (it == fc2cellface.end()) {
          fc2cellface[fc] = cellface_ix;
        } else {
          get<0>(result).push_back({cellface_ix, it->second});
          fc2cellface.erase(it);
        }
      }
    }
    // create a vector of the remaining cellfaces
    get<1>(result).reserve(fc2cellface.size());
    for (auto it = fc2cellface.begin(); it != fc2cellface.end(); ++it) {
      get<1>(result).push_back(it->second);
    }
    return result;
  }

  // ----------------------------------------------------------------------------
  vector<array<double, 3>> cellface_centroids(const double* const coords, 
                                              const vector<size_t>& cellface_ixs,
                                              const int* const num_face_corners,
                                              const int* const face_corners)
  // ----------------------------------------------------------------------------    
  {
    vector<array<double, 3>> result(cellface_ixs.size(), {0.0, 0.0, 0.0});

    const auto max_cellface_ix =
      *max_element(cellface_ixs.begin(), cellface_ixs.end());
    
    vector<size_t> face_corner_starts(max_cellface_ix + 1, 0);
    for (size_t i = 1; i != face_corner_starts.size(); ++i)
      face_corner_starts[i] = face_corner_starts[i-1] + num_face_corners[i-1];

    for (size_t i = 0; i != cellface_ixs.size(); ++i) {
      const size_t cellface_ix = cellface_ixs[i];
      
      // compute "centroid" as means of face corners
      for (size_t c = 0; c != num_face_corners[cellface_ix]; ++c) {
        const size_t fc = face_corners[face_corner_starts[cellface_ix] + c];
        for (int d = 0; d !=3; ++d)
          result[i][d] += coords[3 * fc + d];
        
      }
      for (int d = 0; d != 3; ++d)
        result[i][d] /= num_face_corners[cellface_ix];
    }
    return result;
  }
  
  // ----------------------------------------------------------------------------
  vector<tuple<double, IndexPair>>
  mutual_distances(const vector<array<double, 3>>& points)
  // ----------------------------------------------------------------------------
  {
    vector<std::tuple<double, IndexPair>> result;
    const size_t num_points = points.size();
    for (size_t i = 0; i != num_points; ++i) {
      for (size_t j = i + 1; j != num_points; ++j) {
        const double dist = sqrt(pow(points[i][0] - points[j][0], 2) +
                                 pow(points[i][1] - points[j][1], 2) +
                                 pow(points[i][2] - points[j][2], 2));
        result.push_back({dist, {i, j}});
      }
    }
    return result;
  }

    
};
