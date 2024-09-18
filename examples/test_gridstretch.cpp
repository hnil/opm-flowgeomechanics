#include "opm/geomech/GridStretcher.hpp"

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <vector>
#include <string>
#include <iostream>

#include "opm/geomech/param_interior.hpp"
#include "opm/geomech/GridStretcher.hpp"

#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>



using namespace std;
using namespace Opm;

using Grid = Dune::FoamGrid<2, 3>;
using CoordType = GridStretcher::CoordType;
using Dune::VTKWriter;
using Dune::VTK::nonconforming;

namespace {

void translate(vector<CoordType>& disp, const CoordType& c) { for (auto& e : disp) e += c; }
void uniform_scale(vector<CoordType>& disp,
                   const vector<CoordType>& bcoords,
                   const double fac)
{
    for (size_t i = 0; i != disp.size(); ++i)
        disp[i] += bcoords[i] * (fac - 1);
}
  
}; // end anonymous namespace


int main(int varnum, char** vararg)
{
  const int choice = stoi(vararg[2]);
  string fname(vararg[1]);
  auto grid = Dune::GmshReader<Grid>::read(fname); // unique_ptr

  vector<CoordType> coords;
  for (const auto& vertex : vertices(grid->leafGridView()))
    coords.push_back({vertex.geometry().corner(0)[0],
                      vertex.geometry().corner(0)[1],
                      vertex.geometry().corner(0)[2]});
  
  GridStretcher gs(*grid);
  const vector<size_t>& bindices = gs.boundaryNodeIndices();

  vector<CoordType> bcoords;
  for (auto b : bindices) bcoords.push_back(coords[b]);

  gs.bcentroid_param_mat();
  return 0;
  
  switch (choice) {
  case 1: 
    // test boundary node displacements
    {
      cout << "Displacement test" << endl;
      vector<CoordType> disp(bindices.size(), {0, 0, 0});
      translate(disp, {4, 0, 0});
      uniform_scale(disp, bcoords, 1.5);
      gs.applyBoundaryNodeDisplacements(disp);
      break;
    }
  case 2:
    // test expansion
    {
      cout << "Expanding test" << endl;
      std::vector<double> amounts(bindices.size(), 0.05);
      amounts[4] = 0.4;
      amounts[5] = 0.3;
      amounts[10] = -0.1;
      gs.expandBoundaryCells(amounts);
    break;
    }
  default:
    cout << "Choice was not provided.  Must be 1 (displacement) or 2 (expansion).\n";
    return 1;
  }

  auto vtkwriter =
    make_unique<VTKWriter<Grid::LeafGridView>>(grid->leafGridView(), nonconforming);
  vtkwriter->write("transformed");

  return 0;
};
