#include <iostream>
#include <fstream>
#include <vector>

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
  using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
  using PressureMatrixInfo = std::tuple<std::unique_ptr<Matrix>, std::vector<Htrans>>;
  
  const double E = 1e9;
  const double nu = 0.25;
  const double init_frac_press = 1e6;
}; // end anonymous namespace

// ============================================================================
std::vector<Htrans> computeHtrans(const Grid& grid)
// ============================================================================
{
  std::vector<Htrans> result;
  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto geom = element.geometry();
    for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
      if (is.boundary())
        continue;
      const int nIdx = mapper.index(is.outside());
      const auto igeom = is.geometry();
      if( eIdx < nIdx) {
        const auto iscenter = igeom.center();
        const auto ncenter = is.outside().geometry().center() - iscenter;
        const auto ecenter = geom.center() - iscenter;
        const double area = igeom.volume();
        const double h1 = ecenter.two_norm() / area;
        const double h2 = ncenter.two_norm() / area;

        result.push_back({nIdx, eIdx, h1, h2});
      }
    }
  }
  return result;
}

// ============================================================================
PressureMatrixInfo pressureMatrixStructure(const std::unique_ptr<Grid>& grid)
// ============================================================================
{
  // computing pre-transmissibilities
  const std::vector<Htrans> htrans = computeHtrans(*grid);

  // setting up matrix and initializing its sparsity structure based on the
  // computed pre-transmissibilities
  std::unique_ptr<Matrix> pressure_matrix = std::make_unique<Matrix>();
  pressure_matrix->setBuildMode(Matrix::implicit);

  // map from dof=3*nodes at a cell (ca 3*3*3) to cell
  pressure_matrix->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
  const size_t nc = grid->leafGridView().size(0);
  pressure_matrix->setSize(nc, nc);
  for (const auto& elem : htrans) {
    const size_t ix(std::get<0>(elem));
    const size_t jx(std::get<1>(elem));
    pressure_matrix->entry(ix, jx) = 0.0;
    pressure_matrix->entry(jx, ix) = 0.0;
    pressure_matrix->entry(jx, jx) = 0.0;
    pressure_matrix->entry(ix, ix) = 0.0;
  }

  return PressureMatrixInfo(std::move(pressure_matrix), htrans);
}

// ============================================================================
void updateTrans(PressureMatrixInfo& pmat,
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
  PressureMatrixInfo pmat = pressureMatrixStructure(grid);

  updateTrans(pmat, frac_aperture);
  
  // std::ofstream os("dumped_matrix");
  // Dune::printmatrix(os, *frac_matrix, "", "");
  // os.close();
  
  return 0;
}; 

