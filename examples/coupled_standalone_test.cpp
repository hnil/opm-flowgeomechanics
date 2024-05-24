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

namespace{
  // definitions and constants
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  using Point3D = Dune::FieldVector<double, 3>;
  using Grid = Dune::FoamGrid<2, 3>;
  using Htrans = std::tuple<size_t,size_t, double, double>;
  const double E = 1e9;
  const double nu = 0.25;
  const double init_frac_press = 1e6;
}; // end anonymous namespace

// ============================================================================
std::tuple<std::unique_ptr<Matrix>, Htrans> pressureMatrixStructure()
// ============================================================================
{
  //@@ implement me
  return std::tuple<std::unique_ptr<Matrix>, Htrans>(std::unique_ptr<Matrix>(),
                                                     Htrans(0,0,0,0));
}

// ============================================================================
void updateTrans(std::tuple<std::unique_ptr<Matrix>, Htrans>& pmat,
                 const Dune::BlockVector<Dune::FieldVector<double, 1>>& aperture)
// ============================================================================  
{
  //@@ implement me
}

// ============================================================================
int main()
// ============================================================================
{


  // read test grid from disk
  std::unique_ptr grid = Dune::GmshReader<Grid>::read("disk.msh");
  
  // compute mech matrix
  int nc = grid->leafGridView().size(0);
  std::unique_ptr<Dune::DynamicMatrix<double>>
    frac_matrix(std::make_unique<Dune::DynamicMatrix<double>>());
  frac_matrix->resize(nc, nc);
  ddm::assembleMatrix(*frac_matrix, E, nu, *grid);
  *frac_matrix *= -1;
  
  // computing aperture
  Dune::BlockVector<Dune::FieldVector<double, 1>> frac_press(nc);
  Dune::BlockVector<Dune::FieldVector<double, 1>> frac_aperture(nc);
  frac_press = init_frac_press;
  frac_aperture = 0;
  
  frac_matrix->solve(frac_aperture, frac_press);
  
  std::cout << frac_aperture << std::endl;

  // compute flow matrix
  std::tuple<std::unique_ptr<Matrix>, Htrans> pmat = pressureMatrixStructure();

  updateTrans(pmat, frac_aperture);
  
  // std::ofstream os("dumped_matrix");
  // Dune::printmatrix(os, *frac_matrix, "", "");
  // os.close();
  
  return 0;
}; 

