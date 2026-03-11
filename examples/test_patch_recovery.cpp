#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/geomech/vem/vemutils.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// Minimal Eclipse deck describing a 3x3x3 Cartesian grid.
static const char* deckString = R"(
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
27*100 /
DY
27*100 /
DZ
27*100 /
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
int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "Creating CpGrid from inline deck ..." << std::endl;
    Dune::CpGrid grid = makeCpGrid(deckString);
    std::cout << "  Grid has " << grid.size(0) << " cells, "
              << grid.size(3) << " vertices" << std::endl;

    bool ok = true;

    std::cout << "Running constant preservation test ..." << std::endl;
    if (!testConstantPreservation(grid)) {
        std::cerr << "FAILED: constant function not preserved!" << std::endl;
        ok = false;
    } else {
        std::cout << "  PASSED" << std::endl;
    }

    std::cout << "Running linear preservation test ..." << std::endl;
    if (!testLinearPreservation(grid)) {
        std::cerr << "FAILED: linear function not preserved!" << std::endl;
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
