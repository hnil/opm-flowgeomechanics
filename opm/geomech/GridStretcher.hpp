#pragma once

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/entityseed.hh>

#include <vector>
#include <algorithm>
#include <map>
#include <tuple>


namespace Opm {

class GridStretcher
{
public:

  using Grid = Dune::FoamGrid<2, 3>;
  using CellSeed = Grid::Codim<0>::EntitySeed;
  //using NodeSeed = Grid::Codim<2>::EntitySeed;
  using CellBnodeMap = std::map<size_t,std::tuple<CellSeed, size_t, size_t>>;
  using CoordType = Dune::FieldVector<double, 3>;
  
  GridStretcher(Grid& grid) :
    grid_(grid),
    bnindices_(boundary_node_indices(grid_)),
    iindices_(complement_of(bnindices_, grid_.leafGridView().size(2))), 
    iparam_(interior_parametrization(bnindices_, iindices_, nodecoords(grid_))),
    c2bix_(compute_cell_2_bindices_mapping(grid_, bnindices_)),
    bcindices_(keyvec(c2bix_)) {}

  const std::vector<size_t>& boundaryNodeIndices() const { return bnindices_; }
  const std::vector<size_t>& boundaryCellIndices() const { return bcindices_; }

  // the vector should have one entry per boundary cell, expressing the distance
  // the boundary of that cell should be expanded outwards
  void expandBoundaryCells(const std::vector<double>& amounts); // will modify grid

  // the vector should have two enties per boundary node, specifying its displacement
  // in the x and y direction
  void applyBoundaryNodeDisplacements(const std::vector<CoordType>& disp); // will modify grid

  std::vector<double> centroidEdgeDist() const;

  
private:

  // ----------------------- functions used by constructor -----------------------

  static std::vector<size_t> boundary_node_indices(const Grid& grid);
  static std::vector<size_t> complement_of(const std::vector<size_t>& vec, const size_t N);
  static std::vector<double> nodecoords(const Grid& grid);
  static std::vector<double> interior_parametrization(const std::vector<size_t>& bix,
                                                      const std::vector<size_t>& iix,
                                                      const std::vector<double>& coords);
  static CellBnodeMap compute_cell_2_bindices_mapping(const Grid& grid,
                                                      const std::vector<size_t>& bindices);
  template<typename Key, typename Value>
  static std::vector<Key> keyvec(std::map<Key, Value> map) {
    std::vector<Key> result;
    for (const auto& kv : map) result.push_back(kv.first);
    return result;
  }

  // ------------------------------- internal data -------------------------------

  Grid& grid_; // NB: mutable reference to externally owned grid!
  const std::vector<size_t> bnindices_; // indices of boundary gridnodes
  const std::vector<size_t> iindices_; // indices of internal gridnodes
  const std::vector<double> iparam_; // parametrization of internal nodes in
                                     // terms of boundary nodes
  const CellBnodeMap c2bix_; // map cell -> bindices
  const std::vector<size_t> bcindices_; // indices of boundary cells
};
  
  
  

}; // end namespace Opm
