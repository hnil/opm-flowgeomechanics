#include "RegularTrimesh.hpp"
#include "GeometryHelpers.hpp"
#include <opm/geomech/RegularTrimesh.hpp>
#include <set>
#include <algorithm>
#include <fstream>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Opm;
using namespace std;

namespace {

std::array<EdgeRef, 3> cell2edges(const CellRef& cell) {
  return (cell[2] == 0) ? std::array<EdgeRef, 3> {EdgeRef{cell[0], cell[1], 0},
                                                  EdgeRef{cell[0], cell[1], 1},
                                                  EdgeRef{cell[0], cell[1], 2}} :
                          std::array<EdgeRef, 3> {EdgeRef{cell[0], cell[1] + 1, 0},
                                                  EdgeRef{cell[0] + 1 , cell[1], 1},
                                                  EdgeRef{cell[0], cell[1], 2}};
}

std::array<CellRef, 2> edge2cells(const EdgeRef& edge) {
  return (edge[2] == 0) ?
        array<CellRef, 2> {CellRef{edge[0], edge[1], 0}, CellRef{edge[0], edge[1] - 1, 1}} :
    (edge[2] == 1) ?
        array<CellRef, 2> {CellRef{edge[0], edge[1], 0}, CellRef{edge[0]-1, edge[1], 1}} :
        array<CellRef, 2> {CellRef{edge[0], edge[1], 1}, CellRef{edge[0], edge[1], 0}};
}

std::array<CellRef, 3> cellNeighbors(const CellRef& cell) {
  return (cell[2] == 0) ?
    array<CellRef, 3> {CellRef{cell[0], cell[1], 1},
                       CellRef{cell[0] - 1, cell[1], 1},
                       CellRef{cell[0], cell[1] - 1, 1}} :
    array<CellRef, 3> {CellRef{cell[0], cell[1], 0},
                       CellRef{cell[0] + 1, cell[1], 0},
                       CellRef{cell[0], cell[1] + 1, 0}};
};
  
// += operator for Opm::Coord3D
Opm::Coord3D& operator+=(Opm::Coord3D& lhs, const Opm::Coord3D& rhs) {
  for (int i = 0; i != 3; ++i)
    lhs[i] += rhs[i];
  return lhs;
}

// *= operator for Opm::Coord3D
Opm::Coord3D& operator*=(Opm::Coord3D& lhs, const double rhs) {
  for (int i = 0; i != 3; ++i)
    lhs[i] *= rhs;
  return lhs;
}

// /= operator for Opm::Coord3D
Opm::Coord3D& operator/=(Opm::Coord3D& lhs, const double rhs) {
  for (int i = 0; i != 3; ++i)
    lhs[i] /= rhs;
  return lhs;
}  

// + operator for Opm::Coord3D
Opm::Coord3D operator+(const Opm::Coord3D& lhs, const Opm::Coord3D& rhs) {
  Opm::Coord3D result(lhs);
  result += rhs;
  return result;
}

// / operator for Opm::Coord3D
Opm::Coord3D operator/(const Opm::Coord3D& lhs, const double rhs) {
  Opm::Coord3D result(lhs);
  result /= rhs;
  return result;
}

// * operator for Opm::Coord3D
Opm::Coord3D operator*(const Opm::Coord3D& lhs, const double rhs) {
  Opm::Coord3D result(lhs);
  result *= rhs;
  return result;
}

// * operator for Opm::Coord3D with double on lhs
Opm::Coord3D operator*(const double lhs, const Opm::Coord3D& rhs) {
  return rhs * lhs;
};

// == operator for Opm::NodeRef
bool operator==(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs) {
  return lhs[0] == rhs[0] && lhs[1] == rhs[1];
}

// < operator for Opm::NodeRef
bool operator<(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs) {
  return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

// == operator for Opm::EdgeRef
bool operator==(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs) {
  return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
}

// < operator for Opm::EdgeRef
bool operator<(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs) {
  return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) ||
         (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
}
}; // namespace

namespace Opm
{
// ----------------------------------------------------------------------------
RegularTrimesh::RegularTrimesh(const int layers,
                               const std::array<double, 3>& origin,
                               const std::array<double, 3>& axis1,
                               const std::array<double, 3>& axis2, 
                               const std::array<double, 2>& edgelen)
// ----------------------------------------------------------------------------
  : origin_(origin),
    axis1_(RegularTrimesh::normalize(axis1)),
    axis2_(RegularTrimesh::normalize(axis2)),
    edgelen_(edgelen)
{
    cellinfo_[{0, 0, 0}] = CellAttributes(); // set a single seed cell
    for (int i = 0; i != layers; ++i)
      expandGrid();
}

// ----------------------------------------------------------------------------
std::vector<CellRef> RegularTrimesh::cellIndices() const
// ----------------------------------------------------------------------------
{
  std::vector<CellRef> indices;
  for (const auto& cell : cellinfo_)
    indices.push_back(cell.first);
  return indices;
}
  
// ----------------------------------------------------------------------------  
std::vector<EdgeRef> RegularTrimesh::edgeIndices() const
// ----------------------------------------------------------------------------
{
  std::vector<EdgeRef> all_edges = all_half_edges_();
  // keep only unique elements
  const auto last = std::unique(all_edges.begin(), all_edges.end());

  return std::vector<EdgeRef>(all_edges.begin(), last);
}
  
// ----------------------------------------------------------------------------
std::vector<NodeRef> RegularTrimesh::nodeIndices() const
// ----------------------------------------------------------------------------
{
  std::vector<NodeRef> indices;

  for (const auto& entry : cellinfo_) {
    const auto& cell = entry.first;
    for (int i = 0; i != 3; ++i)
      indices.push_back( i==0 ? (cell[2] == 0 ? NodeRef{cell[0],     cell[1]} :
                                                NodeRef{cell[0] + 1, cell[1] + 1}) :
                         i==1 ? NodeRef{cell[0] + 1, cell[1]} :
                                NodeRef{cell[0], cell[1] + 1} );
  }

  std::sort(indices.begin(), indices.end());
  const auto last = std::unique(indices.begin(), indices.end());

  return std::vector<NodeRef>(indices.begin(), last);
}

// ----------------------------------------------------------------------------  
std::vector<EdgeRef> RegularTrimesh::boundaryEdges() const
// ----------------------------------------------------------------------------
{
  // make a vector of all edges.  Count internal edges twice
  const std::vector<EdgeRef> all_edges = all_half_edges_();

  // boundary edges are those that are not duplicated
  std::vector<EdgeRef> result;
  for (auto it = all_edges.begin(); it != all_edges.end(); ++it)
    if (*it != *(it + 1))
      result.push_back(*it);
    else
      ++it; // skip the next one, which we already know is duplicated

  return result;
}

// ----------------------------------------------------------------------------
std::vector<CellRef> RegularTrimesh::boundaryCells() const
// ----------------------------------------------------------------------------
{
  std::vector<CellRef> result;
  for (const auto& edge : boundaryEdges())
    for (const auto& cell : edge2cells(edge))
      if (isActive(cell))
        result.push_back(cell);

  // remove any duplicates
  std::sort(result.begin(), result.end());
  const auto last = std::unique(result.begin(), result.end());
  return std::vector<CellRef>(result.begin(), last);
}

// ----------------------------------------------------------------------------
std::vector<EdgeRef> RegularTrimesh::all_half_edges_() const
// ----------------------------------------------------------------------------  
{
  // make a vector of all edges.  Count internal edges twice.  Sort the result
  std::vector<EdgeRef> all_edges;
  for (const auto& entry : cellinfo_) {
    const auto& cell = entry.first;
    if (cell[2] == 0) // "right way" triangle
      for (int i = 0; i != 3; ++i)
        all_edges.push_back( {cell[0], cell[1], i} );
    else
      for (int i = 0; i != 3; ++i)
        all_edges.push_back( i == 0 ? EdgeRef{cell[0],     cell[1] + 1, i} :
                             i == 1 ? EdgeRef{cell[0] + 1, cell[1],     i} :
                                      EdgeRef{cell[0],     cell[1],     i} );
  }
  std::sort(all_edges.begin(), all_edges.end());
  return all_edges;
}

// ----------------------------------------------------------------------------      
Coord3D RegularTrimesh::nodeCoord(const NodeRef& node) const
// ----------------------------------------------------------------------------
{
  Coord3D result(origin_);
  result += node[0] * edgelen_[0] * axis1_;
  result += node[1] * edgelen_[1] * axis2_;
  return result;
}
  
// ----------------------------------------------------------------------------    
Coord3D RegularTrimesh::cellCentroid(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
  Coord3D result = cell[2] == 0 ? nodeCoord({cell[0], cell[1]}) :
                                  nodeCoord({cell[0] + 1, cell[1] + 1});
  result += nodeCoord({cell[0] + 1, cell[1]});
  result += nodeCoord({cell[0], cell[1] + 1});

  result /= 3;
  return result;
}

// ----------------------------------------------------------------------------  
Coord3D RegularTrimesh::edgeCentroid(const EdgeRef& edge) const
// ----------------------------------------------------------------------------
{
  return
    (edge[2] == 0) ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0] + 1, edge[1]}))  / 2 :
    (edge[2] == 1) ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0], edge[1] + 1})) / 2 :
                     (nodeCoord({edge[0], edge[1] + 1}) + nodeCoord({edge[0] + 1, edge[1]})) / 2;
}

// ----------------------------------------------------------------------------  
std::vector<Coord3D> RegularTrimesh::cellCentroids() const
// ----------------------------------------------------------------------------  
{
  std::vector<Coord3D> result;
  for (const auto& entry : cellinfo_)
    result.push_back(cellCentroid(entry.first));
  return result;
}
  
// ----------------------------------------------------------------------------  
std::vector<Coord3D> RegularTrimesh::edgeCentroids() const
// ----------------------------------------------------------------------------  
{
  std::vector<Coord3D> result;
  for (const auto& edge : edgeIndices())
    result.push_back(edgeCentroid(edge));
  return result;
}

// ----------------------------------------------------------------------------  
std::vector<Coord3D> RegularTrimesh::nodeCoords() const
// ----------------------------------------------------------------------------  
{
  std::vector<Coord3D> result;
  for (const auto& node : nodeIndices())
    result.push_back(nodeCoord(node));
  return result;
}

//------------------------------------------------------------------------------  
vector<array<unsigned int, 3>> RegularTrimesh::getTriangles() const
//------------------------------------------------------------------------------  
{
  vector<array<unsigned int, 3>> result;
  for (const auto& cell : cellNodesLinear())
    result.push_back( {(unsigned int) cell[0],
                       (unsigned int) cell[1],
                       (unsigned int) cell[2]});
  return result;
}

//------------------------------------------------------------------------------  
vector<array<unsigned int, 3>> RegularTrimesh::getMultiresTriangles() const
//------------------------------------------------------------------------------  
{
  vector<array<unsigned int, 3>> result;

  // create mesh with boundary cells removed
  RegularTrimesh mesh(*this);
  for (int i = 0; i != 3; ++i)
    mesh.contractGrid();

  const RegularTrimesh coarsened = mesh.coarsen();

  // coarsened triangles

  const auto nodeindices = nodeIndices();
  // function to identify the index in nodeindices of a NodeRef
  const auto nodeix = [&nodeindices](const NodeRef& node) {
    return size_t(std::find(nodeindices.begin(), nodeindices.end(), node) - nodeindices.begin());
  };

  vector<CellRef> covered_finecells;
  for (const auto& cell : coarsened.cellIndices()) {

    for (const auto& cref : coarse_to_fine_(cell))
        covered_finecells.push_back(cref);

    const auto coarsenodes = cellNodes(cell);
    result.push_back( {(unsigned int) nodeix(coarse_to_fine_(coarsenodes[0])),
                       (unsigned int) nodeix(coarse_to_fine_(coarsenodes[1])),
                       (unsigned int) nodeix(coarse_to_fine_(coarsenodes[2]))});
  }
  sort(covered_finecells.begin(), covered_finecells.end());
  
  vector<CellRef> all_finecells = cellIndices();
  sort(all_finecells.begin(), all_finecells.end());

  vector<CellRef> uncovered_finecells;
  std::set_difference(all_finecells.begin(), all_finecells.end(),
                      covered_finecells.begin(), covered_finecells.end(),
                      std::back_inserter(uncovered_finecells));

  // cout << "Number of coarse cells: " << coarsened.numActive() << endl;
  // cout << "Number of covered finecells: " << covered_finecells.size() << endl;
  // cout << "Total number of finecells: " << all_finecells.size() << endl;
  // cout << "Number of uncovered finecells: " << uncovered_finecells.size() << endl;
  
  // insert fine triangles to result
  const auto cnodes = cellNodesLinear();
  for (const auto& cell : uncovered_finecells) {
    const auto nodes = cnodes[linearCellIndex(cell)];
    result.push_back({(unsigned int) nodes[0],
                      (unsigned int) nodes[1],
                      (unsigned int) nodes[2]});
  }
  
  return result; // @@ dummy, for now
}

//------------------------------------------------------------------------------
std::unique_ptr<Grid> RegularTrimesh::createDuneGrid(bool coarsen_interior) const 
//------------------------------------------------------------------------------  
{
  Dune::GridFactory<Grid> factory;

  // define points
  for (const auto& node : nodeCoords())
    factory.insertVertex( Dune::FieldVector<double, 3> {node[0], node[1], node[2]});

  // define triangles
  for (const auto& tri : (coarsen_interior? getMultiresTriangles() : getTriangles()))
        factory.insertElement(Dune::GeometryTypes::simplex(2),
                              vector<unsigned int> {tri[0], tri[1], tri[2]});
  
  return factory.createGrid();

}

// ----------------------------------------------------------------------------
std::pair<NodeRef, NodeRef> RegularTrimesh::edgeNodes(const EdgeRef& e) const
// ----------------------------------------------------------------------------
{
  return e[2] == 0 ? std::make_pair(NodeRef{e[0], e[1]}, NodeRef{e[0] + 1, e[1]}) :
         e[2] == 1 ? std::make_pair(NodeRef{e[0], e[1]}, NodeRef{e[0], e[1] + 1}) :
                     std::make_pair(NodeRef{e[0], e[1] + 1}, NodeRef{e[0] + 1, e[1]});
}

// ----------------------------------------------------------------------------
std::vector<std::pair<size_t, size_t>> RegularTrimesh::edgeNodeIndices(bool only_boundary) const
// ----------------------------------------------------------------------------  
{
  const auto nodeindices = nodeIndices();
  const auto edgeindices = only_boundary ? boundaryEdges() : edgeIndices();

  // make mapping from node to index
  std::map<NodeRef, size_t> node2index;
  for (size_t i = 0; i != nodeindices.size(); ++i)
    node2index[nodeindices[i]] = i;

  // map back to indices, for all edges in the mesh
  std::vector<std::pair<size_t, size_t>> result;  
  for (const auto& edge : edgeindices) {
    const auto nodes = edgeNodes(edge);
    result.push_back({node2index[nodes.first], node2index[nodes.second]});
  }

  return result;
}

// ----------------------------------------------------------------------------
std::pair<Coord3D, Coord3D> RegularTrimesh::edgeNodeCoords(const EdgeRef& edge) const
// ----------------------------------------------------------------------------  
{
  return std::make_pair(nodeCoord(edgeNodes(edge).first),
                        nodeCoord(edgeNodes(edge).second));
}

// ----------------------------------------------------------------------------
std::vector<std::pair<Coord3D, Coord3D>> RegularTrimesh::edgeNodeCoords() const
// ----------------------------------------------------------------------------
{
  const auto edgeindices = edgeIndices();
  std::vector<std::pair<Coord3D, Coord3D>> result;
  for (const auto& edge : edgeindices)
    result.push_back(edgeNodeCoords(edge));
  return result;
}

// ----------------------------------------------------------------------------
std::array<NodeRef, 3> RegularTrimesh::cellNodes(const CellRef& cell) 
// ----------------------------------------------------------------------------
{
        return cell[2] == 0 ?
          std::array<NodeRef, 3> {NodeRef{cell[0], cell[1]},
                                  NodeRef{cell[0] + 1, cell[1]},
                                  NodeRef{cell[0], cell[1] + 1}} :
          std::array<NodeRef, 3> {NodeRef{cell[0] + 1, cell[1]},
                                  NodeRef{cell[0], cell[1] + 1},
                                  NodeRef{cell[0] + 1, cell[1] + 1}};
}

// ----------------------------------------------------------------------------
std::vector<std::array<size_t, 3>> RegularTrimesh::cellNodesLinear() const
// ----------------------------------------------------------------------------  
{
  std::vector<std::array<size_t, 3>> result;
  const auto nodeindices = nodeIndices();
  const auto cellindices = cellIndices();

  const auto findnode = [&nodeindices](const NodeRef& node) {
    return size_t(std::find(nodeindices.begin(), nodeindices.end(), node) - nodeindices.begin());
  };
  
  for (const auto& cell : cellindices) {
    const auto noderefs = cellNodes(cell);
    result.push_back({findnode(noderefs[0]),
                      findnode(noderefs[1]),
                      findnode(noderefs[2])});
  }
  return result;
}
  
// ----------------------------------------------------------------------------
void RegularTrimesh::writeMatlabTriangulation(std::ostream& out) const
// ----------------------------------------------------------------------------  
{
  const auto nodeindices = nodeIndices();
  const auto nodecoords = nodeCoords();
  
  // define the vertices of the triangulation
  out << "vertices = [";
  for (const auto& node : nodecoords) {
    out << node[0] << " " << node[1] << " " << node[2] << "; ";
  }
  out << "];\n";
  
  // define the triangles of the triangulation
  const auto cellnodes = cellNodesLinear();
  out << "triangles = [";
  for (const auto& cell : cellnodes) {
    out << cell[0]+1 << " " << cell[1]+1 << " " << cell[2]+1 << "; ";
  }
  out << "];\n";

  // create a triangulation object
  out << "tri = triangulation(triangles, vertices);\n";

  // plot the triangulation, using triplot
  out << "figure;\n";
  out << "triplot(tri);\n";

  //set axis equal
  out << "axis equal;\n";
}


// ----------------------------------------------------------------------------
void writeMeshToVTK(const RegularTrimesh& mesh, const char* const filename,
                    bool coarsen_interior)
// ----------------------------------------------------------------------------
{
  const auto grid = mesh.createDuneGrid(coarsen_interior);

  // write grid to file
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(
                                           grid->leafGridView(),
                                           Dune::VTK::nonconforming);
  // write flag to file
  if (!coarsen_interior) {
    const auto flags = mesh.getCellFlags();
    vtkwriter->addCellData(flags, "flag");
  }
  
  vtkwriter->write(filename);
}

// ----------------------------------------------------------------------------  
void writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename)
// ----------------------------------------------------------------------------  
{
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  
  const std::vector<Coord3D> points = mesh.nodeCoords();
  const std::vector<std::pair<size_t, size_t>> edges = mesh.edgeNodeIndices(true);

  // header
  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  file << "  <UnstructuredGrid>\n";
  file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\""
       << edges.size() << "\">\n";
  // points
  file << "      <Points>\n";
  file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& point : points) 
    file << "          " << point[0] << " " << point[1] << " " << point[2] << "\n";
  file << "        </DataArray>\n";
  file << "      </Points>\n";

  // cells (edges)
  file << "      <Cells>\n";
  file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto& edge : edges) 
    file << "          " << edge.first << " " << edge.second << "\n";
  file << "        </DataArray>\n";
  file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i <= edges.size(); ++i) 
    file << "          " << i * 2 << "\n";
  file << "        </DataArray>\n";
  file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (size_t i = 0; i < edges.size(); ++i) 
        file << "          3\n";
  file << "        </DataArray>\n";
  file << "      </Cells>\n";

  // Footer
  file << "    </Piece>\n";
  file << "  </UnstructuredGrid>\n";
  file << "</VTKFile>\n";

  file.close();
}

// ----------------------------------------------------------------------------
size_t RegularTrimesh::numActive() const
// ----------------------------------------------------------------------------  
{
  return cellinfo_.size();
}
  
// ----------------------------------------------------------------------------  
bool RegularTrimesh::isActive(const CellRef& cell) const
// ----------------------------------------------------------------------------    
{
  return cellinfo_.find(cell) != cellinfo_.end();
}

// ----------------------------------------------------------------------------  
bool RegularTrimesh::setActive(const CellRef& cell)
// ----------------------------------------------------------------------------    
{
  if (isActive(cell))
    return false;
  cellinfo_.insert({cell, CellAttributes()});

  assert(cellinfo_[cell].flag == 0);
  return true;
}

// ----------------------------------------------------------------------------  
int RegularTrimesh::expandGrid(const CellRef& cell)
// ----------------------------------------------------------------------------    
{
  int result = 0;
  const auto edges = cell2edges(cell);
  for (const auto& e : edges) 
    for (const auto& c : edge2cells(e))
      result += setActive(c);

  return result;
}

// ----------------------------------------------------------------------------
int RegularTrimesh::contractGrid()
// ----------------------------------------------------------------------------  
{
  // identify boundary cells and remove them
  const vector<CellRef> bcells = boundaryCells();

  // identify all cells that are not boundary cells, using set_difference
  for (const auto& bc : bcells)
    cellinfo_.erase(bc);

  return (int) bcells.size();
}

  
// ----------------------------------------------------------------------------  
int RegularTrimesh::expandGrid(const std::vector<CellRef>& cells)
// ----------------------------------------------------------------------------    
{
  int result = 0;
  for (const auto& c : cells)
    result += expandGrid(c);
  return result;
}

// ----------------------------------------------------------------------------
int RegularTrimesh::expandGrid()
// ----------------------------------------------------------------------------
{
  // regular expansion in all directions
  return expandGrid(boundaryCells());
}

// ----------------------------------------------------------------------------
void RegularTrimesh::removeSawtooths()
// ----------------------------------------------------------------------------
{
  const auto bedges = boundaryEdges();
  std::vector<CellRef> candidates;
  for (const auto& edge : bedges) 
    for (const auto& cell : edge2cells(edge))
      if (!isActive(cell))
        candidates.push_back(cell);

  std::sort(candidates.begin(), candidates.end());

  // inactive cells adjacent to more than one boundary edge should be activated
  for (auto it = candidates.begin(); it != candidates.end(); ++it) {
    const auto range = std::equal_range(it, candidates.end(), *it);
    if (range.second - range.first > 1) {
        setActive(*it);
        it = range.second - 1;
    }
  }
}

// ----------------------------------------------------------------------------
size_t RegularTrimesh::linearCellIndex(const CellRef& cell) const
{
  // cellinfo is a map from CellRef to CellAttributes.
  assert(cellinfo_.find(cell) != cellinfo_.end());
  return std::distance(cellinfo_.begin(), cellinfo_.find(cell));
}

// ----------------------------------------------------------------------------
CellRef RegularTrimesh::cellIndex(const size_t index) const
{
  auto it = cellinfo_.begin();
  std::advance(it, index);
  return it->first;
}

// ----------------------------------------------------------------------------
void RegularTrimesh::setCellFlag(const CellRef& cell, const int value)
// ----------------------------------------------------------------------------
{
  cellinfo_[cell].flag = value;
}

// ----------------------------------------------------------------------------
void RegularTrimesh::setCellFlags(const std::vector<CellRef>& cells,
                                  const int value)
// ----------------------------------------------------------------------------
{
  for (const auto& cell : cells)
    setCellFlag(cell, value);
}

// ----------------------------------------------------------------------------
void RegularTrimesh::setAllFlags(const int value)
// ----------------------------------------------------------------------------  
{
  for (auto& it : cellinfo_)
    it.second.flag = value;
}

// ----------------------------------------------------------------------------
int RegularTrimesh::getCellFlag(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
  const auto it = cellinfo_.find(cell);
  return it->second.flag;
}

// ----------------------------------------------------------------------------
std::vector<int> RegularTrimesh::getCellFlags() const
// ----------------------------------------------------------------------------
{
  std::vector<int> result;
  for (const auto& el : cellinfo_)
    result.push_back(el.second.flag);
  return result;
}

// ----------------------------------------------------------------------------
RegularTrimesh RegularTrimesh::refine() const
// ----------------------------------------------------------------------------  
{
  map<CellRef, CellAttributes> new_cells;
  for (const auto& e : cellinfo_)
    for (const auto & c : coarse_to_fine_(e.first))
      new_cells[c] = e.second;

  return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0]/2, edgelen_[1]/2}}; 
}

// ----------------------------------------------------------------------------
RegularTrimesh RegularTrimesh::coarsen() const
// ----------------------------------------------------------------------------
{
  map<CellRef, CellAttributes> new_cells;
  for (const auto& e : cellinfo_) {
    const CellRef& c = e.first;
    const int k = (c[0]%2 + c[1]%2 + c[2]%2 > 1) ? 1 : 0;
    if ((k + c[2])% 2 == 1)
      new_cells[fine_to_coarse_(c)] = e.second;
  }
  return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0]*2, edgelen_[1]*2}}; 
}

// for (const auto& e : cellinfo_) {
//       const CellRef& c = e.first;
//       const CellAttributes& attr = e.second;
//       if (c[2] == 1 && c[0] % 2 == 0 && c[1] % 2 == 0) 
//         new_cells[{c[0]/2, c[1]/2, 0}] = attr;
//       else if (c[2] == 0 && c[0] % 2 == 1 && c[1] % 2 == 1) 
//         new_cells[{c[0]/2, c[1]/2, 1}] = attr;
// }
  
// ----------------------------------------------------------------------------
RegularTrimesh expand_to_criterion(const RegularTrimesh& mesh,
     function<std::vector<double>(const RegularTrimesh&)> score_function,
                                    double threshold)
// ----------------------------------------------------------------------------
{
  return RegularTrimesh(); // @@ dummy for now
}

// ----------------------------------------------------------------------------
std::array<CellRef, 4> RegularTrimesh::coarse_to_fine_(const CellRef& c)
// ----------------------------------------------------------------------------  
{
  return (c[2] == 0) ?
    array<CellRef, 4> { CellRef {2*c[0],   2*c[1],   0},
                        CellRef {2*c[0]+1, 2*c[1],   0},
                        CellRef {2*c[0],   2*c[1]+1, 0},
                        CellRef {2*c[0],   2*c[1],   1} } :
    array<CellRef, 4> { CellRef{2*c[0]+1, 2*c[1]+1, 0},
                        CellRef{2*c[0]+1, 2*c[1],   1}, 
                        CellRef{2*c[0],   2*c[1]+1, 1},
                        CellRef{2*c[0]+1, 2*c[1]+1, 1} } ;
}

// ----------------------------------------------------------------------------  
CellRef RegularTrimesh::fine_to_coarse_(const CellRef& cell)
// ----------------------------------------------------------------------------  
{
  const int k = cell[0]%2 + cell[1]%2 + cell[2]%2;
  return {cell[0]/2, cell[1]/2, k > 1 ? 1 : 0};
}

// ----------------------------------------------------------------------------
NodeRef RegularTrimesh::coarse_to_fine_(const NodeRef& node)
// ----------------------------------------------------------------------------
{
  return {2*node[0], 2*node[1]};
}

// ----------------------------------------------------------------------------
vector<CellRef> RegularTrimesh::interior_coarsegrid_() const
// ----------------------------------------------------------------------------  
{
  map<CellRef, int> cell_count;
  for (const auto& e : cellinfo_) {
    const CellRef coarse_cell = fine_to_coarse_(e.first);
    if (cell_count.find(coarse_cell) == cell_count.end())
      cell_count[coarse_cell] = 1;
    else
      cell_count[coarse_cell]++;
  }
  vector<CellRef> result;
  for (const auto& e : cell_count)
    if (e.second == 4) // this cell is fully covered by all its fine cells
      result.push_back(e.first);

  return result;
}


} // namespace

