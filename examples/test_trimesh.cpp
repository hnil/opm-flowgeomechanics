#include <iostream>
#include <vector>
#include <fstream>
#include "opm/geomech/RegularTrimesh.hpp"

#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace Opm;
using namespace std;

// // code to dump edges as vtu
void dumpEdges(const vector<pair<Coord3D, Coord3D>>& edges,
               const string& filename)
{
  ofstream file(filename);
  file << "# vtk DataFile Version 3.0\n";
  file << "Edges\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  file << "POINTS " << edges.size() * 2 << " float\n";
  for (const auto& edge : edges) {
    file << edge.first[0] << " " << edge.first[1] << " " << edge.first[2] << "\n";
    file << edge.second[0] << " " << edge.second[1] << " " << edge.second[2] << "\n";
  }
  file << "LINES " << edges.size() << " " << edges.size() * 3 << "\n";
  for (size_t i = 0; i < edges.size(); ++i) {
    file << "2 " << 2 * i << " " << 2 * i + 1 << "\n";
  }
  file.close();
}
  

// write a test progam
int main() {
  vector<array<int, 3>> cellspec = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 0, 1},
    {1, 1, 1},
    {1, 1, 0}
  };

  RegularTrimesh mesh(cellspec.begin(), cellspec.end());
  
  auto cells = mesh.cellIndices();
  cout << "Cells: \n";
  for (const auto& cell : cells)
    cout << cell[0] << " " << cell[1] << " " << cell[2] << "\n";
  cout << endl;
  
  auto edges = mesh.edgeIndices();
  cout << "Edges: \n";
  for (const auto& edge : edges)
    cout << edge[0] << " " << edge[1] << " " << edge[2] << "\n";
  cout << endl;
  
  auto nodes = mesh.nodeIndices();
  cout << "Nodes: \n";
  for (const auto& node : nodes)
    cout << node[0] << " " << node[1] << " \n";
  cout << endl;

  auto boundary_edges = mesh.boundaryEdges();
  cout << "Boundary edges: \n";
  for (const auto& edge : boundary_edges)
    cout << edge[0] << " " << edge[1] << " " << edge[2] << "\n";
  cout << endl;
  
  auto cell_centroids = mesh.cellCentroids();
  cout << "Cell centroids: \n";
  for (const auto& centroid : cell_centroids)
    cout << centroid[0] << " " << centroid[1] << " " << centroid[2] << "\n";
  cout << endl;
  
  auto edge_centroids = mesh.edgeCentroids();
  cout << "Edge centroids: \n";
  for (const auto& centroid : edge_centroids)
    cout << centroid[0] << " " << centroid[1] << " " << centroid[2] << "\n";
  cout << endl;
  
  auto node_coords = mesh.nodeCoords();
  cout << "Node coords: \n";
  for (const auto& coord : node_coords)
    cout << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
  cout << endl;

  // write to file as a matlab triangulation (open file as fstream)
  ofstream file("mesh.m");
  mesh.writeMatlabTriangulation(file);
  file.close();
  
  const auto grid = mesh.createDuneGrid();

  // write grid to file
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(
                                           grid->leafGridView(),
                                           Dune::VTK::nonconforming);
  vtkwriter->write("mesh");

  // write outer boundary edges to file
  vector<pair<Coord3D, Coord3D>> outer_boundary_edges;
  for (const auto& e : boundary_edges) 
    outer_boundary_edges.push_back(mesh.edgeNodeCoords(e));
    
  dumpEdges(outer_boundary_edges, "outer_boundary_edges.vtk");

}

