#include <iostream>
#include <vector>
#include <fstream>
#include "opm/geomech/RegularTrimesh.hpp"

#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace Opm;
using namespace std;


int expand_grid_test(int turns)
{
  RegularTrimesh mesh;

  for (int i = 0; i != turns; ++i) {
    const auto cells = mesh.cellIndices();
    mesh.expandGrid(cells);
    mesh.removeSawtooths();
  }

  writeMeshToVTK(mesh, "expandedgrid");
  cout << "Number of cells: " << mesh.numActive() << endl;
  return 0;
}

int irregular_grid_test()
{
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

  //write grid to file
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(
                                           grid->leafGridView(),
                                           Dune::VTK::nonconforming);
  vtkwriter->write("mesh");

  // write outer boundary edges to file
  vector<pair<Coord3D, Coord3D>> outer_boundary_edges;
  for (const auto& e : boundary_edges) 
    outer_boundary_edges.push_back(mesh.edgeNodeCoords(e));
    
  //dumpEdges(outer_boundary_edges, "outer_boundary_edges.vtk");

  Opm::writeMeshBoundaryToVTK(mesh, "boundary.vtu");
  return 0;
}


// write a test progam
int main(int varnum, char** vararg) {

  if (varnum == 1)
    cout <<
      "Options are: \n" <<
      "1 - create irregular 5-cell grid\n" <<
      "2 - test grid expansion <n turns> \n" <<
      endl;
  else
    if (atoi(vararg[1]) == 1) return irregular_grid_test();
    else if (atoi(vararg[1]) == 2) return expand_grid_test(atoi(vararg[2]));
    else
      cout << "Invalid option given" << endl;
  return 1;
};
