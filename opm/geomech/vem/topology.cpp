#include "topology.hpp"
#include "vem.hpp"
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream> // @@ for debug/warning

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

  // ----------------------------------------------------------------------------
  bool same_face(tuple<double, double, int, int> f1,
                 tuple<double, double, int, int> f2,
                 const int* const num_face_corners,
                 const int* const face_corners)
  // ----------------------------------------------------------------------------    
  {
    const int f1_ix = get<2>(f1);
    const int f2_ix = get<2>(f2);
    const int fc1_start = get<3>(f1);
    const int fc2_start = get<3>(f2);
    
    
    if (num_face_corners[f1_ix] != num_face_corners[f2_ix])
      return false;

    for (int i = 0; i != num_face_corners[f1_ix]; ++i) 
      if (face_corners[fc1_start + i] != face_corners[fc2_start + i])
        return false;

      return true;
  }
  
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
        result.push_back({(int)cell, (int)face});
    return result;
  }


  // ----------------------------------------------------------------------------
  tuple<vector<IndexPair>, vector<int>>
  cellfaces_matching_faces(const int num_cells,
                           const int* const num_cell_faces,
                           const int* const num_face_corners,
                           const int* const face_corners)
  // ----------------------------------------------------------------------------    
  {
    tuple<vector<IndexPair>, vector<int>> result;
    map<FaceCorners, int> fc2cellface;
    
    const int* nfc_ptr = num_face_corners;
    const int* fc_ptr = face_corners;
    int cellface_ix = 0;

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
                                              const vector<int>& cellface_ixs,
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
    vector<tuple<double, IndexPair>> result;
    const size_t num_points = points.size();
    for (size_t i = 0; i != num_points; ++i) {
      for (size_t j = i + 1; j != num_points; ++j) {
        const double dist = sqrt(pow(points[i][0] - points[j][0], 2) +
                                 pow(points[i][1] - points[j][1], 2) +
                                 pow(points[i][2] - points[j][2], 2));
        result.push_back({dist, {int(i), int(j)}});
      }
    }
    return result;
  }


  // ----------------------------------------------------------------------------
  vector<array<int, 2>> identify_top_bottom_faces(const double* const coords,
                                                  const int num_cells,
                                                  const int* const num_cell_faces,
                                                  const int* const num_face_corners,
                                                  const int* const face_corners)
  // ----------------------------------------------------------------------------
  {
    // the top and bottom face for each cell are identified as the two faces with
    // largest area.  
    vector<array<int, 2>> result;
    result.reserve(num_cells);
    int cf_index = 0; // cell face index; will be counted steadily upward as we
                      // address each cell and its faces.
    int fcor_index = 0; // face corner index; will be counted steadily upward as we
                        // address each face and its corners.

    for (int cell = 0; cell != num_cells; ++cell) {

      // area, lowest z-value and cellface index (last int is the 'fcor_index', may be
      // removed if face duplicity check no longer needed)
      vector<tuple<double, double, int, int>> areas; 
      
      for (int i = 0; i != num_cell_faces[cell]; ++i, ++cf_index) {

        const vector<double> poly_corners = 
          pick_points_3D(coords, &face_corners[fcor_index], num_face_corners[cf_index]);

        // compute area
        const double area = face_integral(&poly_corners[0], num_face_corners[cf_index], 3);
        //const double area = element_volume_2D(&poly_corners[0], num_face_corners[cf_index]);

        // summed z-value
        double sum_zval = 0;
        for (int j = 0; j != num_face_corners[cf_index]; ++j) 
          sum_zval += poly_corners[3 * j + 2];

        areas.push_back({area, sum_zval, cf_index, fcor_index});
        fcor_index += num_face_corners[cf_index];        
      }

      // sort element in 'areas' according to their area (largest first)
      sort(areas.begin(), areas.end(), [](const auto& a, const auto& b) {
          return get<0>(a) > get<0>(b);
      });
      // get the two largest areas, the one with the lowest summed z-value
      // should be mentioned first (it's the top face)
      auto top_face = areas[0];
      size_t ix = 1;
      auto bot_face = areas[ix];

      // @@ the following check should never trigger for properly defined grids 
      while (same_face(top_face, bot_face, num_face_corners, face_corners)) {
          // the two faces are topologically the same.  Should generally not happen!
          cout << "Warning: duplicate face detected in cell " << cell << "!" << endl;
          if (ix == areas.size()-1)
            throw runtime_error("All faces in cell are identical!");
          bot_face = areas[++ix];
      }
      
      if (get<1>(top_face) > get<1>(bot_face))
        std::swap(top_face, bot_face);

      result.push_back({get<2>(top_face), get<2>(bot_face)});
    }
    
    return result;
  }
  
}; // end namespace
