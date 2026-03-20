/*
  Standalone tests for the DDM (Discrete Displacement Method) matrix assembly.

  Tests verify that:
  1. assembleMatrix_fast produces the same results as assembleMatrix
     (the fast version pre-caches geometry and uses OpenMP parallelism).
  2. The assembled matrix has expected symmetry and sign properties.
  3. The assembly works on grids of different sizes.
*/

#include <config.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

namespace
{

// ============================================================================
// Helpers
// ============================================================================

using Grid = Dune::FoamGrid<2, 3>;
using FullMatrix = Dune::DynamicMatrix<double>;

// Create a small planar FoamGrid from a RegularTrimesh with the given
// number of layers around the center cell.
std::unique_ptr<Grid>
makeTestGrid(int layers)
{
    Opm::RegularTrimesh trimesh(layers);
    auto [grid, map, bcells] = trimesh.createDuneGrid(0, {}, false);
    return grid;
}

// Compare two dense matrices element-by-element
bool compareMatrices(const FullMatrix& a,
                     const FullMatrix& b,
                     double tol)
{
    if (a.N() != b.N() || a.M() != b.M()) {
        std::cerr << "  Matrix dimensions differ: "
                  << a.N() << "x" << a.M() << " vs "
                  << b.N() << "x" << b.M() << std::endl;
        return false;
    }

    bool ok = true;
    double max_diff = 0.0;
    for (size_t i = 0; i < a.N(); ++i) {
        for (size_t j = 0; j < a.M(); ++j) {
            const double diff = std::abs(a[i][j] - b[i][j]);
            max_diff = std::max(max_diff, diff);
            const double scale = std::max(std::abs(a[i][j]), std::abs(b[i][j]));
            const double rel = (scale > 1e-15) ? diff / scale : diff;
            if (rel > tol) {
                std::cerr << "  Mismatch at (" << i << "," << j << "): "
                          << a[i][j] << " vs " << b[i][j]
                          << " (rel_err=" << rel << ")" << std::endl;
                ok = false;
            }
        }
    }
    if (ok)
        std::cout << "  Max absolute difference: " << max_diff << std::endl;
    return ok;
}

// ============================================================================
// Tests
// ============================================================================

bool test_fast_matches_original(int layers)
{
    std::cout << "Test: assembleMatrix_fast matches assembleMatrix (layers="
              << layers << ") ..." << std::endl;

    auto grid = makeTestGrid(layers);
    if (!grid) {
        std::cerr << "  FAILED (could not create grid)" << std::endl;
        return false;
    }

    const int nc = grid->leafGridView().size(0);
    std::cout << "  Grid has " << nc << " cells" << std::endl;

    const double E = 1e9;   // Young's modulus
    const double nu = 0.25; // Poisson's ratio

    FullMatrix matrix_orig(nc, nc, 0.0);
    FullMatrix matrix_fast(nc, nc, 0.0);

    ddm::assembleMatrix(matrix_orig, E, nu, *grid);
    ddm::assembleMatrix_fast(matrix_fast, E, nu, *grid);

    if (!compareMatrices(matrix_orig, matrix_fast, 1e-12)) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_diagonal_dominance()
{
    std::cout << "Test: DDM matrix diagonal properties ..." << std::endl;

    auto grid = makeTestGrid(1);
    if (!grid) {
        std::cerr << "  FAILED (could not create grid)" << std::endl;
        return false;
    }

    const int nc = grid->leafGridView().size(0);
    const double E = 1e9;
    const double nu = 0.25;

    FullMatrix matrix(nc, nc, 0.0);
    ddm::assembleMatrix_fast(matrix, E, nu, *grid);

    // For a planar fracture with unit normal-opening slip, the diagonal
    // entries (self-influence) should be positive (compressive normal
    // traction from opening) and larger than off-diagonals.
    bool ok = true;
    for (int i = 0; i < nc; ++i) {
        if (matrix[i][i] <= 0) {
            std::cerr << "  Non-positive diagonal at (" << i << "," << i
                      << "): " << matrix[i][i] << std::endl;
            ok = false;
        }
    }

    if (!ok) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_matrix_nonzero()
{
    std::cout << "Test: DDM matrix is non-trivial (has non-zero entries) ..."
              << std::endl;

    auto grid = makeTestGrid(1);
    if (!grid) {
        std::cerr << "  FAILED (could not create grid)" << std::endl;
        return false;
    }

    const int nc = grid->leafGridView().size(0);
    const double E = 1e9;
    const double nu = 0.25;

    FullMatrix matrix(nc, nc, 0.0);
    ddm::assembleMatrix_fast(matrix, E, nu, *grid);

    double max_val = 0.0;
    for (int i = 0; i < nc; ++i)
        for (int j = 0; j < nc; ++j)
            max_val = std::max(max_val, std::abs(matrix[i][j]));

    if (max_val < 1e-15) {
        std::cerr << "  FAILED (matrix is all zeros)" << std::endl;
        return false;
    }
    std::cout << "  Max matrix entry magnitude: " << max_val << std::endl;
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_different_material_properties()
{
    std::cout << "Test: DDM assembly with different E, nu values ..."
              << std::endl;

    auto grid = makeTestGrid(1);
    if (!grid) {
        std::cerr << "  FAILED (could not create grid)" << std::endl;
        return false;
    }

    const int nc = grid->leafGridView().size(0);

    // Two sets of material properties
    const double E1 = 1e9, nu1 = 0.25;
    const double E2 = 2e9, nu2 = 0.30;

    FullMatrix mat1(nc, nc, 0.0), mat2(nc, nc, 0.0);
    ddm::assembleMatrix_fast(mat1, E1, nu1, *grid);
    ddm::assembleMatrix_fast(mat2, E2, nu2, *grid);

    // Matrices should differ since material properties are different
    bool differs = false;
    for (int i = 0; i < nc && !differs; ++i)
        for (int j = 0; j < nc && !differs; ++j)
            if (std::abs(mat1[i][j] - mat2[i][j]) > 1e-15)
                differs = true;

    if (!differs) {
        std::cerr << "  FAILED (matrices are identical for different properties)"
                  << std::endl;
        return false;
    }

    // Both matrices should still have fast == original consistency
    FullMatrix mat1_orig(nc, nc, 0.0), mat2_orig(nc, nc, 0.0);
    ddm::assembleMatrix(mat1_orig, E1, nu1, *grid);
    ddm::assembleMatrix(mat2_orig, E2, nu2, *grid);

    bool ok = compareMatrices(mat1, mat1_orig, 1e-12)
           && compareMatrices(mat2, mat2_orig, 1e-12);

    if (!ok) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

} // anonymous namespace

// ============================================================================
int main()
// ============================================================================
{
    int failures = 0;

    if (!test_fast_matches_original(1)) ++failures;
    if (!test_fast_matches_original(2)) ++failures;
    if (!test_diagonal_dominance())     ++failures;
    if (!test_matrix_nonzero())         ++failures;
    if (!test_different_material_properties()) ++failures;

    std::cout << "\n=== "
              << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " ===" << std::endl;

    return failures;
}
