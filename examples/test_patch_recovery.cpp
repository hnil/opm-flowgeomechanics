#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/geomech/vem/vemutils.hpp>
#include <opm/geomech/vem/vem.hpp>
#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/materials.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Minimal Eclipse deck describing a 3x3x3 Cartesian grid with given cell sizes.
static std::string
makeDeckString(double dx, double dy, double dz)
{
    std::ostringstream os;
    os << R"(
RUNSPEC
TITLE
PATCH_RECOVERY_TEST
DIMENS
3 3 3 /
OIL
WATER
FIELD
START
1 'JAN' 2000 /
TABDIMS
/
EQLDIMS
/
GRID
DX
27*)" << dx << R"( /
DY
27*)" << dy << R"( /
DZ
27*)" << dz << R"( /
TOPS
9*1000 /
PORO
27*0.3 /
PROPS
DENSITY
50 62.4 0.06 /
PVTW
4000 1.0 3.0E-6 0.5 0.0 /
ROCK
4000 3E-6 /
SWOF
0.2 0 1 0
1.0 1 0 0 /
SOLUTION
EQUIL
1150 4000 1200 0 1000 0 1 0 0 /
SCHEDULE
TSTEP
1 /
END
)";
    return os.str();
}

// --------------------------------------------------------------------------
// Helper: create a CpGrid from a deck string.
// --------------------------------------------------------------------------
static Dune::CpGrid
makeCpGrid(const std::string& deck_str)
{
    Opm::Parser parser;
    auto deck = parser.parseString(deck_str);
    Opm::EclipseState eclState(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclState.getInputGrid(), nullptr, false);
    return grid;
}

// --------------------------------------------------------------------------
// Test: apply patch recovery to a linear function f(x,y,z) = a + bx + cy + dz
// and verify that the recovered field reproduces the original cell values.
// --------------------------------------------------------------------------
static bool
testLinearPreservation(const Dune::CpGrid& grid, double tol = 1e-10)
{
    static constexpr int dim = 3;
    const auto& gv = grid.leafGridView();
    const int num_cells = grid.size(0);

    // Coefficients of the linear function.
    const double a0 = 1.5;
    const double a1 = 0.3;  // x-coefficient
    const double a2 = -0.7; // y-coefficient
    const double a3 = 0.2;  // z-coefficient

    // Evaluate linear function at each cell centroid.
    Dune::BlockVector<Dune::FieldVector<double,1>> cell_values(num_cells);
    for (const auto& elem : elements(gv)) {
        int idx = gv.indexSet().index(elem);
        auto c = elem.geometry().center();
        cell_values[idx][0] = a0 + a1 * c[0] + a2 * c[1] + a3 * c[2];
    }

    // Run patch recovery.
    auto recovered = vem::patchRecovery(grid, cell_values);

    // Check that recovered values match the original values.
    double maxErr = 0.0;
    for (int i = 0; i < num_cells; ++i) {
        double err = std::abs(recovered[i][0] - cell_values[i][0]);
        maxErr = std::max(maxErr, err);
    }

    std::cout << "  Linear preservation test: max error = " << maxErr << std::endl;
    return maxErr < tol;
}

// --------------------------------------------------------------------------
// Test: constant function should be preserved exactly.
// --------------------------------------------------------------------------
static bool
testConstantPreservation(const Dune::CpGrid& grid, double tol = 1e-12)
{
    const auto& gv = grid.leafGridView();
    const int num_cells = grid.size(0);

    const double constant_val = 42.0;
    Dune::BlockVector<Dune::FieldVector<double,1>> cell_values(num_cells);
    for (int i = 0; i < num_cells; ++i)
        cell_values[i][0] = constant_val;

    auto recovered = vem::patchRecovery(grid, cell_values);

    double maxErr = 0.0;
    for (int i = 0; i < num_cells; ++i) {
        double err = std::abs(recovered[i][0] - constant_val);
        maxErr = std::max(maxErr, err);
    }

    std::cout << "  Constant preservation test: max error = " << maxErr << std::endl;
    return maxErr < tol;
}

// --------------------------------------------------------------------------
// Tests of the physical-length Helmholtz smoothing (diffuseCellVector):
// a constant field must be preserved, the volume integral must be conserved,
// and a spike must be spread out (peak reduced, neighbors raised).
// --------------------------------------------------------------------------
static bool
testDiffuseCellVector(const Dune::CpGrid& grid, double length)
{
    const auto& gv = grid.leafGridView();
    const int num_cells = grid.size(0);
    bool ok = true;

    std::vector<double> vol(num_cells);
    for (const auto& elem : elements(gv))
        vol[gv.indexSet().index(elem)] = elem.geometry().volume();

    // constant preservation (relative tolerance: on extreme-aspect grids the
    // Helmholtz system has condition ~ length^2 * trans / volume, so the
    // attainable accuracy is a relative, not absolute, quantity)
    const double cval = 7.0;
    Dune::BlockVector<Dune::FieldVector<double,1>> cvec(num_cells);
    for (int i = 0; i < num_cells; ++i)
        cvec[i][0] = cval;
    auto csmooth = vem::diffuseCellVector(grid, cvec, length);
    double cerr = 0.0;
    for (int i = 0; i < num_cells; ++i)
        cerr = std::max(cerr, std::abs(csmooth[i][0] - cval));
    std::cout << "  diffuse constant preservation: max error = " << cerr << std::endl;
    if (cerr > 1e-6 * cval) {
        std::cerr << "FAILED: diffuse does not preserve constants" << std::endl;
        ok = false;
    }

    // spike: conservation and spreading.  Put a unit spike in the center cell (1,1,1)
    Dune::BlockVector<Dune::FieldVector<double,1>> spike(num_cells);
    spike = 0.0;
    const int center = 1 + 3 * (1 + 3 * 1);
    spike[center][0] = 1.0;
    auto ssmooth = vem::diffuseCellVector(grid, spike, length);

    double int_before = vol[center] * 1.0, int_after = 0.0, minval = 1e30;
    for (int i = 0; i < num_cells; ++i) {
        int_after += vol[i] * ssmooth[i][0];
        minval = std::min(minval, ssmooth[i][0]);
    }
    std::cout << "  diffuse spike: center " << ssmooth[center][0] << ", integral "
              << int_after / int_before << " of original, min " << minval << std::endl;
    if (std::abs(int_after - int_before) > 1e-6 * int_before) {
        std::cerr << "FAILED: diffuse does not conserve the volume integral" << std::endl;
        ok = false;
    }
    if (!(ssmooth[center][0] < 1.0 && ssmooth[center][0] > 0.0)) {
        std::cerr << "FAILED: spike not reduced (or negative)" << std::endl;
        ok = false;
    }
    if (minval < -1e-10) {
        std::cerr << "FAILED: diffuse produced negative values from a positive field"
                  << std::endl;
        ok = false;
    }
    return ok;
}

// --------------------------------------------------------------------------
// Test: 6-component nodal patch recovery (used for stress interpolation onto
// fracture surfaces).  Each component is a different linear function; the
// recovered nodal values must equal the functions at the node coordinates.
// --------------------------------------------------------------------------
static bool
testNodal6LinearExactness(const Dune::CpGrid& grid, double tol = 1e-9)
{
    static constexpr int dim = 3;
    const auto& gv = grid.leafGridView();
    const int num_cells = grid.size(0);

    const double coef[6][4] = {{1.5, 0.3, -0.7, 0.2},
                               {-2.0, 0.1, 0.4, -0.3},
                               {0.5, -0.2, 0.6, 0.1},
                               {3.0, 0.05, -0.15, 0.25},
                               {-1.0, 0.45, 0.2, -0.1},
                               {2.5, -0.35, 0.05, 0.4}};
    auto fval = [&](int c, const double* x) {
        return coef[c][0] + coef[c][1] * x[0] + coef[c][2] * x[1] + coef[c][3] * x[2];
    };

    Dune::BlockVector<Dune::FieldVector<double,6>> cell_values(num_cells);
    for (const auto& elem : elements(gv)) {
        const int idx = gv.indexSet().index(elem);
        const auto c = elem.geometry().center();
        const double x[3] = {c[0], c[1], c[2]};
        for (int comp = 0; comp < 6; ++comp)
            cell_values[idx][comp] = fval(comp, x);
    }

    const auto nodal = vem::patchRecoveryNodal6(grid, cell_values);

    double maxErr = 0.0, scale = 0.0;
    for (const auto& v : vertices(gv)) {
        const int nidx = gv.indexSet().index(v);
        const auto p = v.geometry().corner(0);
        const double x[3] = {p[0], p[1], p[2]};
        for (int comp = 0; comp < 6; ++comp) {
            const double ex = fval(comp, x);
            maxErr = std::max(maxErr, std::abs(nodal[nidx][comp] - ex));
            scale = std::max(scale, std::abs(ex));
        }
    }
    std::cout << "  Nodal 6-component linear exactness: max rel error = " << maxErr / scale
              << std::endl;
    return maxErr / scale < tol;
}

// --------------------------------------------------------------------------
// Test: the standalone Q1 hexahedron stiffness (vem::stiffness_matrix_fem_hex_3D,
// used by the array-based FEM/FEM_BBAR path) must reproduce the element stiffness
// of the legacy opm-upscaling Q1 implementation (Opm::Elasticity, used by the
// Dune/CpGrid FEM assembly path) entry by entry.
// --------------------------------------------------------------------------
static bool
testQ1MatchesOpmElasticity(const Dune::CpGrid& grid)
{
    constexpr int dim = 3;
    constexpr int comp = 6;
    constexpr int bfunc = 8;
    const double E_mod = 1e9, nu = 0.3;

    Opm::Elasticity::Elasticity upscaling_fem(grid);
    Dune::FieldMatrix<double, comp, comp> C;
    Opm::Elasticity::Isotropic material(1, E_mod, nu, 0.0);
    if (!material.getConstitutiveMatrix(C)) {
        std::cerr << "FAILED: could not build constitutive matrix" << std::endl;
        return false;
    }

    const auto& gv = grid.leafGridView();
    const auto elem = *gv.begin<0>();

    // legacy opm-upscaling Q1: 2x2x2 Gauss with Dune geometry
    Dune::FieldMatrix<double, dim * bfunc, dim * bfunc> K_legacy(0.0);
    const auto& rule = Dune::QuadratureRules<double, dim>::rule(elem.type(), 2);
    for (auto r = rule.begin(); r != rule.end(); ++r) {
        const auto jacInvTra = elem.geometry().jacobianInverseTransposed(r->position());
        const double detJ = elem.geometry().integrationElement(r->position());
        Dune::FieldMatrix<double, comp, dim * bfunc> lB;
        upscaling_fem.getBmatrix(lB, r->position(), jacInvTra);
        Dune::FieldMatrix<double, dim * bfunc, dim * bfunc> Aq;
        upscaling_fem.getStiffnessMatrix(Aq, lB, C, detJ * r->weight());
        K_legacy += Aq;
    }

    // standalone vem Q1 on the same corner coordinates (Dune reference ordering)
    std::array<double, 24> hc;
    for (int k = 0; k < 8; ++k) {
        const auto c = elem.geometry().corner(k);
        for (int d = 0; d < 3; ++d)
            hc[3 * k + d] = c[d];
    }
    std::vector<double> K_vem(24 * 24);
    vem::stiffness_matrix_fem_hex_3D(&hc[0], E_mod, nu, /*bbar=*/false, &K_vem[0]);

    double maxdiff = 0.0, scale = 0.0;
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j) {
            maxdiff = std::max(maxdiff, std::abs(K_legacy[i][j] - K_vem[i * 24 + j]));
            scale = std::max(scale, std::abs(K_legacy[i][j]));
        }
    std::cout << "  Q1 vs opm-upscaling element stiffness: max rel diff = " << maxdiff / scale
              << std::endl;
    return maxdiff / scale < 1e-12;
}

// --------------------------------------------------------------------------
static bool
runRecoveryTests(const Dune::CpGrid& grid, const std::string& label, double tol)
{
    bool ok = true;
    std::cout << "Running constant preservation test (" << label << ") ..." << std::endl;
    if (!testConstantPreservation(grid)) {
        std::cerr << "FAILED: constant function not preserved!" << std::endl;
        ok = false;
    } else {
        std::cout << "  PASSED" << std::endl;
    }

    std::cout << "Running linear preservation test (" << label << ") ..." << std::endl;
    if (!testLinearPreservation(grid, tol)) {
        std::cerr << "FAILED: linear function not preserved!" << std::endl;
        ok = false;
    } else {
        std::cout << "  PASSED" << std::endl;
    }

    std::cout << "Running nodal 6-component recovery test (" << label << ") ..." << std::endl;
    if (!testNodal6LinearExactness(grid)) {
        std::cerr << "FAILED: nodal 6-component recovery not linear-exact!" << std::endl;
        ok = false;
    } else {
        std::cout << "  PASSED" << std::endl;
    }
    return ok;
}

// --------------------------------------------------------------------------
int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    bool ok = true;

    std::cout << "Creating isotropic CpGrid (100x100x100 cells) ..." << std::endl;
    Dune::CpGrid grid = makeCpGrid(makeDeckString(100, 100, 100));
    std::cout << "  Grid has " << grid.size(0) << " cells, " << grid.size(3) << " vertices"
              << std::endl;
    ok = runRecoveryTests(grid, "isotropic", 1e-10) && ok;

    // high-aspect-ratio cells (100 x 100 x 0.1): stresses the conditioning of the
    // least-squares fit (the offsets are normalized per direction in tryFitAtNode)
    std::cout << "\nCreating anisotropic CpGrid (100x100x0.1 cells, aspect 1000) ..." << std::endl;
    Dune::CpGrid aniso_grid = makeCpGrid(makeDeckString(100, 100, 0.1));
    ok = runRecoveryTests(aniso_grid, "aspect 1000", 1e-9) && ok;

    std::cout << "\nRunning Helmholtz smoothing (diffuseCellVector) tests ..." << std::endl;
    ok = testDiffuseCellVector(grid, 150.0) && ok; // 1.5 cell sizes
    ok = testDiffuseCellVector(aniso_grid, 150.0) && ok;

    std::cout << "\nRunning Q1 vs opm-upscaling element stiffness comparison ..." << std::endl;
    if (!testQ1MatchesOpmElasticity(grid) || !testQ1MatchesOpmElasticity(aniso_grid)) {
        std::cerr << "FAILED: standalone Q1 does not match the opm-upscaling element!"
                  << std::endl;
        ok = false;
    } else {
        std::cout << "  PASSED" << std::endl;
    }

    if (ok) {
        std::cout << "\nAll patch recovery tests PASSED." << std::endl;
        return 0;
    } else {
        std::cerr << "\nSome patch recovery tests FAILED." << std::endl;
        return 1;
    }
}
