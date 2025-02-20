#include <iostream>
#include <vector>
#include <fstream>
#include "opm/geomech/RegularTrimesh.hpp"

using namespace Opm;
using namespace std;

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
  
  auto grid = mesh.createDuneGrid();
}

