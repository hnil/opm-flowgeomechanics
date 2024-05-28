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

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>



#include "opm/geomech/DiscreteDisplacement.hpp"

namespace{
  // definitions and constants
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
  using Point3D = Dune::FieldVector<double, 3>;
  using Grid = Dune::FoamGrid<2, 3>;
  using Htrans = std::tuple<size_t,size_t, double, double>;
  using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
  using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
  using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
  using PressureMatrixInfo = std::tuple<std::unique_ptr<Matrix>,std::vector<Htrans>>;
                                        
  using IntFloatPair = std::tuple<int, double>;


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
  pressure_matrix->compress();

  return PressureMatrixInfo(std::move(pressure_matrix), htrans);
}

// ============================================================================
template<typename T> inline T hmean(const T a, const T b) {return T(1) / ( T(1)/a + T(1)/b);};
// ============================================================================

// ============================================================================
template<typename T, typename Container>
inline bool inrange(T el, Container con) {return std::find(con.begin(), con.end(), el) != con.end();}
// ============================================================================

// ============================================================================
void updateTrans(PressureMatrixInfo& pmat,
                 const Vector& aperture,
                 const std::vector<size_t> imposed_vals_ixs = std::vector<size_t>())
// ============================================================================  
{
  auto& matrix = *std::get<0>(pmat);
  const auto& elems = std::get<1>(pmat);
  
  for (const auto& e : elems) {
    const size_t i = std::get<0>(e);
    const size_t j = std::get<1>(e);
    const double t1 = std::get<2>(e);
    const double t2 = std::get<3>(e);
    const double h1 = aperture[i];
    const double h2 = aperture[j];

    const double trans1 = h1 * h1 * t1 / 12.0;
    const double trans2 = h2 * h2 * t2 / 12.0;
    
    const double T = hmean(trans1, trans2);

    // imposed values will be associated with trivial equations 
    const bool i_imposed = inrange(i, imposed_vals_ixs);
    const bool j_imposed = inrange(j, imposed_vals_ixs);

    matrix[i][i] = i_imposed ? double(1) : double(matrix[i][i] + T);
    matrix[i][j] = i_imposed ? double(0) : double(matrix[i][j] - T);
    
    matrix[j][j] = j_imposed ? double(1) : double(matrix[j][j] + T);
    matrix[j][i] = j_imposed ? double(0) : double(matrix[j][i] - T);  
  }
}

// ============================================================================
Vector solvePressure(const Vector& aperture,
                     PressureMatrixInfo& pmat,
                     const std::vector<IntFloatPair> fixed_pvals = std::vector<IntFloatPair>())
// ============================================================================  
{
  // fill in pressure matrix with actual values, based on current aperture and
  // imposed pressure
  std::vector<size_t> imposed_ixs;
  for_each(fixed_pvals.begin(), fixed_pvals.end(), [&](const auto& el) {
    imposed_ixs.push_back(std::get<0>(el));});
    
  updateTrans(pmat, aperture, imposed_ixs);

  // prepare right-hand side
  Vector rhs(aperture.size()); rhs = 0;
  for (auto fp : fixed_pvals)
    rhs[std::get<0>(fp)] = std::get<1>(fp);

  // setup solver
  const auto prm = Opm::setupPropertyTree(Opm::FlowLinearSolverParameters(), true, true);
  auto op = PressureOperatorType(*std::get<0>(pmat));
  auto psolver = std::make_unique<FlexibleSolverType>(op, prm, std::function<Vector()>(), size_t(0));

  // solve pressure system
  Vector result(aperture.size());
  Dune::InverseOperatorResult r;
  psolver->apply(result, rhs, r);

  return result;  
}


// ============================================================================
template<int N>
std::array<size_t, N> n_closest(const Grid& grid)
// ============================================================================
{
  using Elem = std::tuple<size_t, double>;
  std::vector<Elem> distances;
  const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
  for (const auto& element : Dune::elements(grid.leafGridView())) {
    const int eIdx = mapper.index(element);
    const auto center = element.geometry().center();
    distances.push_back({distances.size(), center.two_norm()});
  }

  std::sort(distances.begin(), distances.end(), [](Elem& a, Elem& b) {return std::get<1>(a) < std::get<1>(b);});

  // std::for_each(distances.begin(), distances.end(), [&](const auto& t) {
  //   std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ")\n";
  // });

  std::array<size_t, N> result;
  for (int i = 0; i != N; ++i)
    result[i] = std::get<0>(distances[i]);
  
  return result;
}
  
// ============================================================================
int main()
// ============================================================================
{
  // read test grid from disk
  const auto grid = Dune::GmshReader<Grid>::read("disk.msh"); // unique_ptr
  
  // compute mech matrix
  const int nc = grid->leafGridView().size(0);
  std::unique_ptr<Dune::DynamicMatrix<double>>
    frac_matrix(std::make_unique<Dune::DynamicMatrix<double>>());
  frac_matrix->resize(nc, nc);
  ddm::assembleMatrix(*frac_matrix, E, nu, *grid);
  *frac_matrix *= -1;
  
  // computing aperture
  Vector frac_press(nc), frac_aperture(nc);
  frac_press = init_frac_press;
  frac_aperture = 0;
  
  frac_matrix->solve(frac_aperture, frac_press);
  
  //std::cout << frac_aperture << std::endl;

  // determining well cells and prescribing fixed pressure at these
  const auto wellcells = n_closest<2>(*grid);

  std::vector<IntFloatPair> fixed_pvals;
  const double pfixedval = 1.0;
  std::transform(wellcells.begin(), wellcells.end(), std::back_inserter(fixed_pvals),
                 [pfixedval] (const size_t ix) {return IntFloatPair(ix, pfixedval);});
  
  // prepare pressure matrix structure
  PressureMatrixInfo pmat = pressureMatrixStructure(grid);

  // solve pressure system
  frac_press = solvePressure(frac_aperture, pmat, fixed_pvals);

  std::cout << frac_press;
  
    // std::ofstream os("dumped_matrix");
  // Dune::printmatrix(os, *frac_matrix, "", "");
  // os.close();
  
  return 0;
}; 

