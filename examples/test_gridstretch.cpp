#include "opm/geomech/GridStretcher.hpp"

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <vector>
#include <string>
#include <iostream>
#include <iterator>

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
  assert(varnum == 3);
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

  //gs.bcentroid_param_mat();
  //return 0;
  
  switch (choice) {
  case 1: 
    // test boundary node displacements
    {
      cout << "Displacement test" << endl;
      vector<CoordType> disp(bindices.size(), {0, 0, 0});
      translate(disp, {4, 0, 0});
      uniform_scale(disp, bcoords, 1.5);
      // vector<CoordType> disp(bindices.size(), {0, 0, 0});
      // disp[19] = {0.1, 0, 0};
      // disp[20] = {0.2, 0, 0};
      // disp[21] = {0.3, 0, 0};
      // disp[22] = {0.4, 0, 0};
      
      gs.applyBoundaryNodeDisplacements(disp);
      break;
    }
  case 2:
    // test expansion
    {
      cout << "Expanding test" << endl;
      // std::vector<double> amounts(bindices.size(), 0.05);
      // amounts[4] = 0.4;
      // amounts[5] = 0.3;
      // amounts[10] = -0.1;

      std::vector<double> amounts(bindices.size(), 0);
      amounts[0] = 0.1;
      amounts[1] = 0.2;
      //amounts[23] = 0.3;
      //amounts[20] = 0.4;
      
      gs.expandBoundaryCells(amounts);
    break;
    }
  case 3:
    // test gradient
  {
    vector<double> grad;
    //vector<double> disp(gs.boundaryNodeIndices().size(), 1);
    vector<double> disp(gs.boundaryNodeIndices().size(), 0);
    disp[0] = 1;
    disp[10] = 2;
    vector<double> target(gs.centroidEdgeDist());
    //target[0] *= 2;
    //for (int i = 0; i != target.size(); ++i) target[i] *= 2;
    const double obj = gs.objective(disp, target, grad);
    cout << "Objective value: " << obj << endl;
    cout << "Analytical gradient: ";
    copy(grad.begin(), grad.end(), ostream_iterator<double>(cout, " "));
    cout << endl;

    // computing numerical gradient
    double delta = 1e-5;
    for (int i = 0; i != disp.size(); ++i) {
      vector<double> dispDelta = disp;
      dispDelta[i] += delta;
      double obj2 = gs.objective(dispDelta, target, grad);
      cout << (obj2-obj) / delta << " " << grad[i] << endl;;
      //cout << i << " : " << (obj2-obj) / delta << endl;
    }
    cout << endl;
    
    return 0;
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
