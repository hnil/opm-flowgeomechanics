/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/geomech/vem/vemutils.hpp>
#include <opm/geomech/vem/vem.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Minimal deck describing a 3x3x3 Cartesian grid with given cell sizes.
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
// Test: Q1 and B-bar (vem::stiffness_matrix_fem_hex_3D) integrate the energy of a
// linear displacement exactly on every cell, and rigid motions cost nothing.
// --------------------------------------------------------------------------
static bool
testQ1LinearEnergy(const Dune::CpGrid& grid)
{
    const double E_mod = 1e9, nu = 0.3;
    const double lambda = E_mod * nu / ((1 + nu) * (1 - 2 * nu));
    const double mu = E_mod / (2 * (1 + nu));
    // u = G x with a symmetric part (strain) and a skew part (rotation)
    const double G[3][3] = {{1e-3, 4e-4, -2e-4}, {1e-4, -5e-4, 3e-4}, {2e-4, -1e-4, 7e-4}};
    double eps[3][3], tr = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            eps[i][j] = 0.5 * (G[i][j] + G[j][i]);
    for (int i = 0; i < 3; ++i)
        tr += eps[i][i];
    double density = lambda * tr * tr;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            density += 2 * mu * eps[i][j] * eps[i][j];

    const auto& gv = grid.leafGridView();
    double worst_e = 0.0, worst_r = 0.0;
    for (const auto& elem : elements(gv)) {
        std::array<double, 24> hc, u, r;
        for (int k = 0; k < 8; ++k) {
            const auto c = elem.geometry().corner(k);
            for (int d = 0; d < 3; ++d)
                hc[3 * k + d] = c[d];
            for (int d = 0; d < 3; ++d) {
                // relative to corner 0: a translation is in K's kernel but costs digits in u.K.u
                const auto c0 = elem.geometry().corner(0);
                u[3 * k + d] = G[d][0] * (c[0] - c0[0]) + G[d][1] * (c[1] - c0[1]) + G[d][2] * (c[2] - c0[2]);
                // rigid: translation plus a rotation about z
                r[3 * k + d] = (d == 0 ? 1.0 - c[1] : d == 1 ? 2.0 + c[0] : 3.0) * 1e-3;
            }
        }
        for (const bool bbar : {false, true}) {
            std::vector<double> K(24 * 24);
            vem::stiffness_matrix_fem_hex_3D(&hc[0], E_mod, nu, bbar, &K[0]);
            double e_lin = 0.0, rigid = 0.0, kmax = 0.0, rmax = 0.0;
            for (double v : r)
                rmax = std::max(rmax, std::abs(v));
            for (int i = 0; i < 24; ++i) {
                double Kr = 0.0;
                for (int j = 0; j < 24; ++j) {
                    e_lin += u[i] * K[i * 24 + j] * u[j];
                    Kr += K[i * 24 + j] * r[j];
                    kmax = std::max(kmax, std::abs(K[i * 24 + j]));
                }
                rigid = std::max(rigid, std::abs(Kr));
            }
            const double exact = density * elem.geometry().volume();
            worst_e = std::max(worst_e, std::abs(e_lin - exact) / exact);
            worst_r = std::max(worst_r, rigid / (kmax * rmax));
        }
    }
    std::cout << "  Q1 / B-bar: linear-field energy rel error = " << worst_e
              << ", rigid-mode residual = " << worst_r << std::endl;
    return worst_e < 1e-10 && worst_r < 1e-10;
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

    std::cout << "\nRunning Q1 linear-energy test ..." << std::endl;
    if (!testQ1LinearEnergy(grid) || !testQ1LinearEnergy(aniso_grid)) {
        std::cerr << "FAILED: Q1 element does not integrate linear fields exactly!" << std::endl;
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
