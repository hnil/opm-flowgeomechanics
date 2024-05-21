//#include <dune/common/exceptions.hh> // We use exceptions
//#include <dune/common/version.hh>
//#include <dune/grid/utility/structuredgridfactory.hh>
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
//#include <dune/common/fmatrix.hh>
//#include <dune/common/version.hh>
//#include <dune/common/parallel/mpihelper.hh>
// #include <opm/grid/polyhedralgrid.hh>
// #include <opm/input/eclipse/Deck/Deck.hpp>
// #include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
// #include <opm/input/eclipse/EclipseState/EclipseState.hpp>
// #include <opm/input/eclipse/Parser/Parser.hpp>
// #include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
// #include <opm/common/utility/platform_dependent/reenable_warnings.h>
// #include <opm/input/eclipse/Python/Python.hpp>
// #include <opm/input/eclipse/Schedule/Schedule.hpp>
// #include <opm/input/eclipse/Schedule/Well/Well.hpp>
// #include <opm/simulators/utils/readDeck.hpp>
// #include <opm/input/eclipse/Parser/ParseContext.hpp>
// #if HAVE_ALUGRID
// #include <dune/alugrid/grid.hh>
// #include <dune/alugrid/common/fromtogridfactory.hh>
// #endif
// #include <opm/grid/utility/StopWatch.hpp>
// #include <opm/common/utility/parameters/ParameterGroup.hpp>

// #include <opm/simulators/linalg/PropertyTree.hpp>
// #include <opm/geomech/elasticity_solver.hpp>
// #include <opm/geomech/vem_elasticity_solver.hpp>
// #include <opm/elasticity/matrixops.hpp>

// #include <cstring>
// #include "vectorfunctions.hh"
// #include <unistd.h>
// #include <opm/geomech/boundaryutils.hh>

#include <iostream>

#include <dune/istl/matrixmarket.hh>
#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/gmsh4/gridcreators/lagrangegridcreator.hh>
#include <dune/grid/uggrid.hh> // debug
#include "opm/geomech/Fracture.hpp"

//#include <dune/common/filledarray.hh>
//#include <dune/common/parallel/mpihelper.hh>
//#include <dune/grid/yaspgrid.hh>

int main(int argc, char** argv) {

  auto refgrid = Dune::Gmsh4Reader<Dune::FoamGrid<2,3>>::createGridFromFile("Sphere.msh");
  
  // using Grid = Dune::FoamGrid<2, 3>;
  // Dune::GridFactory<Grid> factory;
  // Dune::Gmsh4::LagrangeGridCreator creator{factory};
  

  //using Grid = Dune::UGGrid<2>;
  //Dune::Gmsh4Reader<Grid>::createGridFromFile("filename.msh");

  
  using Point3D = Dune::FieldVector<double, 3>;

  // constructing and initializing fracture
  Opm::Fracture frac;
  
  int perf = 1;
  int wellcell = -1;
  Point3D origo {0, 0, 0};
  Point3D normal {0, 0, 1};
  std::cout << origo << std::endl;
  
  frac.init("testwell", perf, wellcell, origo, normal);

  // calling updateReservoirCells

  // @@ We might not have to do this, as it only creates a mapping between
  // global grid and fracture grid through filling in/updating
  // frac.reservoir_cells;
  
  //const Opm::EclipseGrid& eclgrid = eclState.getInputGrid();
  //external::cvf::ref<external::cvf::BoundingBoxTree> cellSearchTree;
  //external::buildBoundingBoxTree(cellSearchTree, eclgrid);
  //Grid3D& grid3D; 
  //frac.updateReservoirCells(cellSearchTree, grid3D);


  frac.updateReservoirProperties();
  frac.solve();

  frac.write();

  frac.printPressureMatrix();
  return 0;
}; // end main 
