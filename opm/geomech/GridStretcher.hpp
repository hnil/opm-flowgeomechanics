#pragma once

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/entityseed.hh>

#include <vector>
#include <algorithm>
#include <map>
#include <tuple>


namespace Opm {

struct BoundaryNormals {
  std::vector<Dune::FieldVector<double, 3>> bnode_normals;
  std::vector<Dune::FieldVector<double, 3>> bcell_normals;
};


  
class GridStretcher
{
public:

  using Grid = Dune::FoamGrid<2, 3>;
  using CellSeed = Grid::Codim<0>::EntitySeed;
  //using NodeSeed = Grid::Codim<2>::EntitySeed;
  // cell index -> {cell seed, boundary node ix 1, boundary node ix 2}
  using CellBnodeMap = std::map<size_t,std::tuple<CellSeed, size_t, size_t>>;
  using CoordType = Dune::FieldVector<double, 3>;
  
  GridStretcher(Grid& grid) :
    grid_(grid),
    nodecoords_(node_coordinates(grid_)), 
    bnindices_(boundary_node_indices(grid_)),
    iindices_(complement_of(bnindices_, grid_.leafGridView().size(2))), 
    iparam_(interior_parametrization(bnindices_, iindices_, nodecoords_)),
    c2bix_(compute_cell_2_bindices_mapping(grid_, bnindices_)),
    bcindices_(keyvec(c2bix_)),
    bcentroid_param_(bcentroid_param_mat(grid_, bnindices_, iindices_, iparam_, c2bix_)),
    boundary_normals_(boundary_normals(grid_, c2bix_, bcindices_, nodecoords_)) {}

  const std::vector<size_t>& boundaryNodeIndices() const { return bnindices_; }
  const std::vector<size_t>& boundaryCellIndices() const { return bcindices_; }

  // the vector should have one entry per boundary cell, expressing the distance
  // the boundary of that cell should be expanded outwards
  void expandBoundaryCells(const std::vector<double>& amounts); // will modify grid
                           

  // the vector should have two enties per boundary node, specifying its displacement
  // in the x and y direction
  void applyBoundaryNodeDisplacements(const std::vector<CoordType>& disp); // will modify grid
                                      
  std::vector<double>
  computeBoundaryNodeDisplacements(
         const std::vector<double>& amounts,
         const std::vector<CoordType> bnodenormals = std::vector<CoordType>()) const;
  
  std::vector<double> centroidEdgeDist() const;

  const std::vector<CoordType>& nodecoords() const {return nodecoords_;}
  const std::vector<CoordType>& bnodenormals() const {return boundary_normals_.bnode_normals;}
  const std::vector<CoordType>& bcellnormals() const {return boundary_normals_.bcell_normals;}

  // objective function and derivatives, when trying to stretch grid to a particular
  // target (in terms of distances between cell and edge centroids for boundary cells)
  double objective(const std::vector<double>& bndisp,
                   const std::vector<double>& dtarget,
                   std::vector<double>& grad,
                   bool fixed_cell_centroids = false);
  
private:

  // ----------------------- functions used by constructor -----------------------

  static BoundaryNormals boundary_normals(const Grid& grid,
                                          const CellBnodeMap& c2bix,
                                          const std::vector<size_t>& bcindices,
                                          const std::vector<CoordType>& nodecoords);
  static std::vector<CoordType> node_coordinates(const Grid& grid);
  static std::vector<size_t> boundary_node_indices(const Grid& grid);
  static std::vector<size_t> complement_of(const std::vector<size_t>& vec, const size_t N);
  static std::vector<double> interior_parametrization(const std::vector<size_t>& bix,
                                                      const std::vector<size_t>& iix,
                                                      const std::vector<CoordType>& coords);
  static CellBnodeMap compute_cell_2_bindices_mapping(const Grid& grid,
                                                      const std::vector<size_t>& bindices);
  static std::vector<double> bcentroid_param_mat(const Grid& grid,
                                                 const std::vector<size_t>& bnindices,
                                                 const std::vector<size_t>& iindices,
                                                 const std::vector<double>& iparam,
                                                 const CellBnodeMap& c2bix);
  template<typename Key, typename Value>
  static std::vector<Key> keyvec(std::map<Key, Value> map) {
    std::vector<Key> result;
    for (const auto& kv : map) result.push_back(kv.first);
    return result;
  }

    
  // ------------------------------- internal data -------------------------------

  Grid& grid_; // NB: mutable reference to externally owned grid!
  std::vector<CoordType> nodecoords_; // NB! should be updated when grid is updated.
    
  const std::vector<size_t> bnindices_; // indices of boundary gridnodes
  const std::vector<size_t> iindices_; // indices of internal gridnodes
  const std::vector<double> iparam_; // parametrization of internal nodes in
                                     // terms of boundary nodes
  const CellBnodeMap c2bix_; // map cell -> bindices
  const std::vector<size_t> bcindices_; // indices of boundary cells
  const std::vector<double> bcentroid_param_; // parametrization of boundary cell centroids

  BoundaryNormals boundary_normals_; // NB! should be updated when grid is updated.
    
};
  
  
  

}; // end namespace Opm
