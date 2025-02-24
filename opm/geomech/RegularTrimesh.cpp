#include "RegularTrimesh.hpp"
#include <opm/geomech/RegularTrimesh.hpp>
#include <set>
#include <algorithm>
#include <fstream>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

namespace {

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
std::unique_ptr<Grid> RegularTrimesh::createDuneGrid() const 
//------------------------------------------------------------------------------  
{
  Dune::GridFactory<Grid> factory;

  // define points
  for (const auto& node : nodeCoords())
    factory.insertVertex( Dune::FieldVector<double, 3> {node[0], node[1], node[2]});

  // define triangles 
  for (const auto& cell : cellNodes())
    factory.insertElement(Dune::GeometryTypes::simplex(2),
                          std::vector<unsigned int> {static_cast<unsigned int>(cell[0]),
                                                     static_cast<unsigned int>(cell[1]),
                                                     static_cast<unsigned int>(cell[2]) });
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
std::vector<std::array<size_t, 3>> RegularTrimesh::cellNodes() const
// ----------------------------------------------------------------------------  
{
  std::vector<std::array<size_t, 3>> result;
  const auto nodeindices = nodeIndices();
  const auto cellindices = cellIndices();

  const auto findnode = [&nodeindices](const NodeRef& node) {
    return size_t(std::find(nodeindices.begin(), nodeindices.end(), node) - nodeindices.begin());
  };
  
  for (const auto& cell : cellindices)
    if (cell[2] == 0)
      result.push_back({findnode({cell[0], cell[1]}),
                        findnode({cell[0] + 1, cell[1]}),
                        findnode({cell[0], cell[1] + 1})});
    else
      result.push_back({findnode({cell[0] + 1, cell[1]}),
                        findnode({cell[0], cell[1] + 1}),
                        findnode({cell[0] + 1, cell[1] + 1})});
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
  const auto cellnodes = cellNodes();
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
void writeMeshToVTK(const RegularTrimesh& mesh, const std::string& filename)
// ----------------------------------------------------------------------------
{
  const auto grid = mesh.createDuneGrid();

  // write grid to file
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(
                                           grid->leafGridView(),
                                           Dune::VTK::nonconforming);
  vtkwriter->write(filename);
}

// ----------------------------------------------------------------------------  
void writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const std::string& filename)
// ----------------------------------------------------------------------------  
{
  std::ofstream file(filename.c_str());
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


} // namespace
 
