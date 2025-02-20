#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <memory>

// // #include <dune/grid/common/mcmgmapper.hh> // mapper class
// // #include <dune/grid/io/file/vtk.hh>
// // #include <dune/grid/io/file/vtk/vtkwriter.hh>
// // #include <dune/grid/utility/persistentcontainer.hh>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>

// #include <dune/grid/yaspgrid.hh>
// #include <dune/grid/yaspgrid/partitioning.hh>
// #include <dune/istl/matrixmarket.hh>

namespace Opm {

typedef std::array<int, 3> CellRef; // (i, j, [0 | 1])
typedef std::array<int, 3> EdgeRef; // (i, j, [0 | 1 | 2])
typedef std::array<int, 2> NodeRef;  // (i, j)
typedef std::array<double, 3> Coord3D;
typedef Dune::FoamGrid<2, 3> Grid;
//typedef int Grid; // @@ dummy


struct CellAttributes {}; // fill this with contents as needed


  
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
  
  // functions for outputting grids
  std::unique_ptr<Grid> createDuneGrid() const;
  void writeMatlabTriangulation(std::ostream& out) const;
  
private:
  std::vector<EdgeRef> all_half_edges_() const ; // internal edges are duplicated
  
  std::map<CellRef, CellAttributes> cellinfo_;
  const Coord3D origin_;
  const Coord3D axis1_;
  const Coord3D axis2_;
  const std::array<double, 2> edgelen_;
};

  
} // namespace Opm
