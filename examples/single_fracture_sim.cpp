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
#include <fstream>

#include <dune/istl/matrixmarket.hh>

//#include "opm/geomech/Fracture.hpp"
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/istl/io.hh> // needed for printSparseMatrix??


#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>


#include "opm/geomech/DiscreteDisplacement.hpp"
//#include "GeometryHelpers.hpp" // ddm
//#include <dune/common/filledarray.hh>
//#include <dune/common/parallel/mpihelper.hh>
//#include <dune/grid/yaspgrid.hh>

int main(int argc, char** argv) {

  using Point3D = Dune::FieldVector<double, 3>;
  using Grid = Dune::FoamGrid<2, 3>;
  
  // read test grid from disk
  std::unique_ptr grid = Dune::GmshReader<Grid>::read("disk.msh");

  // write test grid to disk as VTK
  auto vtkwriter = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid->leafGridView(), Dune::VTK::nonconforming);
  vtkwriter->write("krull");

  // Assemble mechanic matrix
  int nc = grid->leafGridView().size(0);
  std::unique_ptr<Dune::DynamicMatrix<double>> frac_matrix(std::make_unique<Dune::DynamicMatrix<double>>());
  frac_matrix->resize(nc, nc);
  const double E = 1.0; //1e9;
  const double nu = 0.25;
  ddm::assembleMatrix(*frac_matrix,E, nu,*grid);

  std::ofstream os("dumped_matrix");
  Dune::printmatrix(os, *frac_matrix, "", "");
  os.close();
  
  // // create fracture
  // Opm::Fracture frac;
  // int perf = 1;
  // int wellcell = -1;
  // Point3D origo {0, 0, 0};
  // Point3D normal {0, 0, 1};
  
  // frac.init("testwell", perf, wellcell, origo, normal);
  // frac.setFractureGrid(grid);
  
  // frac.updateReservoirProperties();

  // // frac.solve();

  // // frac.write();

  // // frac.printPressureMatrix();
  // frac.printMechMatrix();
  
  return 0;
}




  
