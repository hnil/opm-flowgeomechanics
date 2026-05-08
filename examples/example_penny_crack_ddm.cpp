/*
  Penny-shaped crack DDM example and validation test.

  This program sets up a circular (penny-shaped) crack using the Displacement
  Discontinuity Method (DDM) and compares the numerical crack opening profile
  and stress intensity factor K_I against analytic solutions from linear
  elastic fracture mechanics.

  Analytic solutions (Sneddon 1946, penny-shaped crack under uniform net
  pressure p, crack radius a):

    Total crack opening displacement (COD) at radial distance r from the center:
        w(r) = 8 * p * (1 - nu^2) / (pi * E) * sqrt(a^2 - r^2)

    Maximum opening at the center (r = 0):
        w_max = 8 * p * (1 - nu^2) * a / (pi * E)

    Mode-I stress intensity factor at the crack rim:
        K_I = 2 * p * sqrt(a) / sqrt(pi)

  The near-tip LEFM relation used to extract K_I from the DDM solution:
        K_I ~= E / (8*(1-nu^2)) * sqrt(2*pi/s) * w(s)
  where s is the distance from the boundary-element centre to the crack front.

  References:
    Sneddon, I. N. (1946). The distribution of stress in the neighbourhood of
    a crack in an elastic solid. Proc. R. Soc. London A, 187, 229-260.
    https://en.wikipedia.org/wiki/Stress_intensity_factor
*/

#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace
{

using Grid     = Dune::FoamGrid<2, 3>;
using GridView = Grid::LeafGridView;
using Mapper   = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

// ============================================================================
// Analytic solutions
// ============================================================================

/// Total crack opening displacement at radial distance r (Sneddon 1946).
/// w(r) = 8*p*(1-nu^2)/(pi*E) * sqrt(a^2 - r^2)
double
analytic_cod(const double r, const double a, const double p, const double E, const double nu)
{
    if (r >= a)
        return 0.0;
    return 8.0 * p * (1.0 - nu * nu) / (M_PI * E) * std::sqrt(a * a - r * r);
}

/// Mode-I stress intensity factor at the crack rim.
/// K_I = 2*p*sqrt(a) / sqrt(pi)
double
analytic_K1(const double a, const double p)
{
    return 2.0 * p * std::sqrt(a) / std::sqrt(M_PI);
}

/// Near-tip LEFM estimate of K_I from the total COD at distance s from tip.
/// From: w = 8*K_I*(1-nu^2)/E * sqrt(s/(2*pi))
/// =>  K_I = w * E / (8*(1-nu^2)) * sqrt(2*pi/s)
double
K1_from_cod(const double w, const double s, const double E, const double nu)
{
    if (s <= 0.0 || w <= 0.0)
        return 0.0;
    return w * E / (8.0 * (1.0 - nu * nu)) * std::sqrt(2.0 * M_PI / s);
}

// ============================================================================
// Single-resolution test
// ============================================================================

/// Run the penny-crack DDM test at a given number of refinement layers.
/// Returns true if all checks pass.
bool
run_test(const int layers,
         const double E,
         const double nu,
         const double p,
         const double cod_rel_tol,
         const double K1_rel_tol)
{
    // ------------------------------------------------------------------
    // Build the penny-shaped crack mesh.
    // Use the layers constructor with default edge length 1, so the
    // approximate crack radius is `layers` units.
    // ------------------------------------------------------------------
    const double h = 1.0; // element edge length

    Opm::RegularTrimesh trimesh(layers,
                                {0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.5, std::sqrt(3.0) / 2.0, 0.0},
                                {h, h});

    auto [grid, cellmap, bcells] = trimesh.createDuneGrid(0, {}, false);

    const int nc = grid->leafGridView().size(0);
    std::cout << "  layers=" << layers
              << "  cells=" << nc << std::endl;

    // ------------------------------------------------------------------
    // Determine the effective crack radius from the mesh.
    // Use the maximum element-centroid radius rounded up by h/2 so that
    // the analytic formula covers the full discrete crack.
    // ------------------------------------------------------------------
    Mapper mapper(grid->leafGridView(), Dune::mcmgElementLayout());

    double r_max = 0.0;
    for (auto elem : elements(grid->leafGridView())) {
        auto center = elem.geometry().center();
        const double r = std::sqrt(center[0] * center[0] + center[1] * center[1]);
        r_max = std::max(r_max, r);
    }
    // Effective crack radius: outermost element centroid + half edge length
    const double a_eff = r_max + 0.5 * h;

    std::cout << "  Effective crack radius a_eff = " << a_eff << " m" << std::endl;

    // ------------------------------------------------------------------
    // Assemble the DDM influence matrix.
    // A[i][j] = normal traction at centre of element i due to unit
    //           normal displacement discontinuity (opening) at element j.
    // ------------------------------------------------------------------
    Dune::DynamicMatrix<double> A(nc, nc, 0.0);
    ddm::assembleMatrix_fast(A, E, nu, *grid);

    // ------------------------------------------------------------------
    // Set up RHS: uniform net pressure p on all elements.
    // Solve:  A * w = p
    // where w is the vector of total crack openings per element.
    // ------------------------------------------------------------------
    Dune::DynamicVector<double> rhs(nc, p);
    Dune::DynamicVector<double> opening(nc, 0.0);

    // Note: DynamicMatrix::solve() performs LU factorisation in-place.
    A.solve(opening, rhs);

    // ------------------------------------------------------------------
    // Compare DDM openings with the analytic COD profile.
    // ------------------------------------------------------------------
    const double w_max_analytic = analytic_cod(0.0, a_eff, p, E, nu);
    std::cout << "  Analytic max COD (r=0): " << w_max_analytic << " m" << std::endl;

    double sum_rel_err   = 0.0;
    double max_rel_err   = 0.0;
    double sum_abs_err   = 0.0;
    int    n_interior    = 0;

    for (auto elem : elements(grid->leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double rx   = center[0];
        const double ry   = center[1];
        const double r    = std::sqrt(rx * rx + ry * ry);

        const double w_ddm      = opening[idx];
        const double w_analytic = analytic_cod(r, a_eff, p, E, nu);

        // Only compare interior elements (well away from the discretised rim)
        if (r < a_eff - 1.5 * h && w_analytic > 1e-16) {
            const double abs_err = std::abs(w_ddm - w_analytic);
            const double rel_err = abs_err / w_analytic;
            sum_rel_err += rel_err;
            sum_abs_err += abs_err;
            max_rel_err  = std::max(max_rel_err, rel_err);
            ++n_interior;
        }
    }

    const double avg_rel_err = (n_interior > 0) ? sum_rel_err / n_interior : 0.0;

    std::cout << "  Interior cells compared: " << n_interior << std::endl;
    std::cout << "  Max relative COD error:  " << max_rel_err * 100.0 << " %" << std::endl;
    std::cout << "  Avg relative COD error:  " << avg_rel_err * 100.0 << " %" << std::endl;

    // ------------------------------------------------------------------
    // Estimate K_I from the near-tip COD of boundary elements.
    // Use the LEFM relation K_I = w*E/(8*(1-nu^2)) * sqrt(2*pi/s)
    // where s = distance from element centroid to crack front (= a_eff - r).
    // ------------------------------------------------------------------
    const double K1_analytic = analytic_K1(a_eff, p);
    std::cout << "  Analytic K_I = " << K1_analytic << " Pa*sqrt(m)" << std::endl;

    double K1_sum   = 0.0;
    int    K1_count = 0;

    for (auto elem : elements(grid->leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double rx   = center[0];
        const double ry   = center[1];
        const double r    = std::sqrt(rx * rx + ry * ry);

        // Only use elements in the outermost ring (tip-distance in [0, 1.5*h])
        const double s = a_eff - r;
        if (s > 0.0 && s < 1.5 * h) {
            const double w   = opening[idx];
            const double K1  = K1_from_cod(w, s, E, nu);
            if (K1 > 0.0) {
                K1_sum  += K1;
                ++K1_count;
            }
        }
    }

    bool ok = true;

    if (K1_count > 0) {
        const double K1_ddm     = K1_sum / K1_count;
        const double K1_rel_err = std::abs(K1_ddm - K1_analytic) / K1_analytic;
        std::cout << "  DDM K_I (avg over " << K1_count
                  << " boundary elements) = " << K1_ddm << " Pa*sqrt(m)" << std::endl;
        std::cout << "  Relative K_I error: " << K1_rel_err * 100.0 << " %" << std::endl;

        if (K1_rel_err > K1_rel_tol) {
            std::cerr << "  FAILED: K_I relative error " << K1_rel_err * 100.0
                      << " % exceeds tolerance " << K1_rel_tol * 100.0 << " %\n";
            ok = false;
        }
    } else {
        std::cerr << "  WARNING: no boundary elements found for K_I estimate\n";
    }

    if (max_rel_err > cod_rel_tol) {
        std::cerr << "  FAILED: max COD relative error " << max_rel_err * 100.0
                  << " % exceeds tolerance " << cod_rel_tol * 100.0 << " %\n";
        ok = false;
    }

    std::cout << "  => " << (ok ? "PASSED" : "FAILED") << "\n";
    return ok;
}

// ============================================================================
// Convergence study helper
// ============================================================================

/// Print the COD profile along x-axis for the given mesh resolution.
void
print_cod_profile(const int layers,
                  const double E,
                  const double nu,
                  const double p)
{
    const double h = 1.0;

    Opm::RegularTrimesh trimesh(layers,
                                {0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.5, std::sqrt(3.0) / 2.0, 0.0},
                                {h, h});

    auto [grid, cellmap, bcells] = trimesh.createDuneGrid(0, {}, false);
    const int nc = grid->leafGridView().size(0);

    Mapper mapper(grid->leafGridView(), Dune::mcmgElementLayout());

    double r_max = 0.0;
    for (auto elem : elements(grid->leafGridView())) {
        auto center = elem.geometry().center();
        const double r = std::sqrt(center[0] * center[0] + center[1] * center[1]);
        r_max = std::max(r_max, r);
    }
    const double a_eff = r_max + 0.5 * h;

    Dune::DynamicMatrix<double> A(nc, nc, 0.0);
    ddm::assembleMatrix_fast(A, E, nu, *grid);

    Dune::DynamicVector<double> rhs(nc, p);
    Dune::DynamicVector<double> opening(nc, 0.0);
    A.solve(opening, rhs);

    std::cout << "\n  COD profile (layers=" << layers
              << ", a_eff=" << a_eff << "):\n";
    std::cout << "    r [m]      DDM w [m]   Analytic w [m]  Rel.err [%]\n";

    // Collect (r, w) pairs, sort by r, print
    std::vector<std::pair<double, double>> rv;
    for (auto elem : elements(grid->leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double r    = std::sqrt(center[0] * center[0] + center[1] * center[1]);
        rv.push_back({r, opening[idx]});
    }
    std::sort(rv.begin(), rv.end());

    for (auto& [r, w_ddm] : rv) {
        const double w_an  = analytic_cod(r, a_eff, p, E, nu);
        const double rel   = (w_an > 1e-16) ? std::abs(w_ddm - w_an) / w_an * 100.0 : 0.0;
        std::cout << "    " << r << "  " << w_ddm << "  " << w_an
                  << "  " << rel << "\n";
    }
}

} // anonymous namespace

// ============================================================================
int
main()
// ============================================================================
{
    // ------------------------------------------------------------------
    // Material and loading parameters
    // ------------------------------------------------------------------
    const double E  = 1.0e9; // Young's modulus [Pa]
    const double nu = 0.25;  // Poisson's ratio  [-]
    const double p  = 1.0e6; // Net internal pressure [Pa]

    std::cout << "======================================================\n";
    std::cout << " Penny-shaped crack DDM vs analytic comparison\n";
    std::cout << "======================================================\n";
    std::cout << "  E  = " << E  << " Pa\n";
    std::cout << "  nu = " << nu << "\n";
    std::cout << "  p  = " << p  << " Pa\n\n";

    // ------------------------------------------------------------------
    // Print COD profile for one resolution to aid visual inspection
    // ------------------------------------------------------------------
    std::cout << "--- COD profile (layers=5) ---\n";
    print_cod_profile(5, E, nu, p);

    // ------------------------------------------------------------------
    // Validation tests at two mesh resolutions.
    //
    // Tolerances are set conservatively since the DDM on a coarse
    // triangular mesh introduces discretisation errors.  The COD interior
    // tolerance of 15 % is appropriate for the coarse mesh (layers=5).
    // At finer resolution (layers=8) the error should be comfortably
    // below the same threshold.
    //
    // The K_I estimate uses the LEFM near-tip formula and a single ring of
    // boundary elements; a 30 % tolerance accommodates both the geometric
    // discretisation error and the variability in tip-distance s.
    // ------------------------------------------------------------------
    const double cod_tol = 0.15; // 15 % relative COD tolerance
    const double K1_tol  = 0.30; // 30 % relative K_I tolerance

    int failures = 0;

    std::cout << "\n--- Test: layers=5 ---\n";
    if (!run_test(5, E, nu, p, cod_tol, K1_tol))
        ++failures;

    std::cout << "\n--- Test: layers=8 ---\n";
    if (!run_test(8, E, nu, p, cod_tol, K1_tol))
        ++failures;

    std::cout << "\n======================================================\n";
    std::cout << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << "\n";
    std::cout << "======================================================\n";

    return failures;
}
