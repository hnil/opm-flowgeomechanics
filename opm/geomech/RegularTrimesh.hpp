#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Opm {

typedef std::array<int, 3> CellRef; // (i, j, [0 | 1])
typedef std::array<int, 3> EdgeRef; // (i, j, [0 | 1 | 2])
typedef std::array<int, 2> NodeRef;  // (i, j)
typedef std::array<double, 3> Coord3D;
using Grid = Dune::FoamGrid<2, 3>;
  //typedef Dune::FoamGrid<2, 3> Grid;
//typedef int Grid; // @@ dummy


struct CellAttributes {
  CellAttributes() : flag(0) {};
  int flag;
}; // fill this with contents as needed


  
class RegularTrimesh
{
public:
  template<typename Iterator>
  RegularTrimesh(const Iterator cells_begin, const Iterator cells_end,
                 const std::array<double, 3>& origin = {0, 0, 0},
                 const std::array<double, 3>& axis1 = {1, 0, 0},
                 const std::array<double, 3>& axis2 = {0.5, std::sqrt(3)/2, 0},
                 const std::array<double, 2>& edgelen = {1, 1})
    : origin_(origin), axis1_(axis1), axis2_(axis2), edgelen_(edgelen)
  {
    for (auto it = cells_begin; it != cells_end; ++it) 
      cellinfo_[*it] = CellAttributes();
  }

  RegularTrimesh(const int layers = 1,
                 const std::array<double, 3>& origin = {0, 0, 0},
                 const std::array<double, 3>& axis1 = {1, 0, 0},
                 const std::array<double, 3>& axis2 = {0.5, std::sqrt(3)/2, 0},
                 const std::array<double, 2>& edgelen = {1, 1});
  
  // --------------------- Functions for inspecting the grid ---------------------
  std::vector<CellRef> cellIndices() const;
  std::vector<EdgeRef> edgeIndices() const;
  std::vector<NodeRef> nodeIndices() const;

  std::vector<EdgeRef> boundaryEdges() const;
  
  Coord3D cellCentroid(const CellRef& cell) const;
  Coord3D edgeCentroid(const EdgeRef& edge) const;
  Coord3D nodeCoord(const NodeRef& node) const;
  
  std::vector<Coord3D> cellCentroids() const;
  std::vector<Coord3D> edgeCentroids() const;
  std::vector<Coord3D> nodeCoords() const;

  std::vector<std::array<size_t, 3>> cellNodes() const;

  std::pair<NodeRef, NodeRef> edgeNodes(const EdgeRef& edge) const;
  std::vector<std::pair<size_t, size_t>> edgeNodeIndices(bool only_boundary=false) const;
  std::pair<Coord3D, Coord3D> edgeNodeCoords(const EdgeRef& edge) const;
  std::vector<std::pair<Coord3D, Coord3D>> edgeNodeCoords() const;

  bool isActive(const CellRef& cell) const;
  size_t numActive() const;
  size_t linearCellIndex(const CellRef& cell) const;
  CellRef cellIndex(const size_t idx) const;
  int getCellFlag(const CellRef& cell) const;
  std::vector<int> getCellFlags() const;

    // --------------------- Functions for modifying the grid ---------------------
  void setAllFlags(const int value);
  void setCellFlag(const CellRef& cell, const int value);
  bool setActive(const CellRef& cell);
  int expandGrid(const CellRef& cell);
  int expandGrid(const std::vector<CellRef>& cells);
  int expandGrid(); // uniform expansion all directions
  void removeSawtooths();
  
  // ---------------------- Functions for outputting grids ----------------------
  std::unique_ptr<Grid> createDuneGrid() const;
  void writeMatlabTriangulation(std::ostream& out) const;
  
private:
  // helper functions
  std::vector<EdgeRef> all_half_edges_() const ; // internal edges are duplicated

  // data members
  std::map<CellRef, CellAttributes> cellinfo_;
  const Coord3D origin_;
  const Coord3D axis1_;
  const Coord3D axis2_;
  const std::array<double, 2> edgelen_;
};

void writeMeshToVTK(const RegularTrimesh& mesh, const char* const filename);
void writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename);
  
} // namespace Opm
