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

#include "opm/geomech/Fracture.hpp"
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

  // // write test grid to disk as VTK
  // auto vtkwriter = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid->leafGridView(), Dune::VTK::nonconforming);
  // vtkwriter->write("krull");

  // // Assemble mechanic matrix
  // int nc = grid->leafGridView().size(0);
  // std::unique_ptr<Dune::DynamicMatrix<double>> frac_matrix(std::make_unique<Dune::DynamicMatrix<double>>());
  // frac_matrix->resize(nc, nc);
  // const double E = 1e9;
  // const double nu = 0.25;
  // ddm::assembleMatrix(*frac_matrix,E, nu,*grid);

  // std::ofstream os("dumped_matrix");
  // Dune::printmatrix(os, *frac_matrix, "", "");
  // os.close();
  
  // // create fracture
  Opm::Fracture frac;
  int perf = 1;
  int wellcell = -1;
  Point3D origo {0, 0, 0};
  Point3D normal {0, 0, 1};
  Opm::PropertyTree prm;
  prm.put("reservoir.perm", 1e-10);
  prm.put("reservoir.dist", 1.0);
  prm.put("reservoir.mobility", 1.0);
  prm.put("config.axis_scale", 1.0);
  prm.put("config.initial_fracture_width", 1e-6);
  prm.put("outputdir", std::string("."));
  prm.put("casename", std::string("testcase"));
  prm.put("solver.method", std::string("if"));

  //double perm = prm.get<double>("reservoir.perm");
  //  std::cout << "Perm: " << perm << std::endl;
  
  frac.init("testwell", perf, wellcell, origo, normal, prm);
  frac.setFractureGrid(grid);

  frac.solve();
  
  // frac.updateReservoirProperties();

  // // frac.solve();

  // // frac.write();

  // // frac.printPressureMatrix();
  // frac.printMechMatrix();
  
  return 0;
}




  
