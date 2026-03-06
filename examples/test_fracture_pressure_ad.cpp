/*
  Unit tests for the AD-based fracture pressure assembly.

  Tests verify that:
  1. The pressure Jacobian (dR/dp) from the AD assembler matches the
     original manually coded assembly.
  2. The aperture coupling matrix (dR/dw) is correct by comparison with
     finite differences.
  3. Edge cases (width below min_width, different control types) are
     handled correctly.
*/

#include <config.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <opm/geomech/FracturePressureAssemblerAD.hpp>

namespace
{

// ============================================================================
// Helpers
// ============================================================================

// Compare two BCRS matrices element-by-element (iterate over sparsity of 'a')
bool compareMatrices(const Opm::BCRSMatrix1x1& a,
                     const Opm::BCRSMatrix1x1& b,
                     double tol)
{
    if (a.N() != b.N() || a.M() != b.M()) {
        std::cerr << "  Matrix dimensions differ: "
                  << a.N() << "x" << a.M() << " vs "
                  << b.N() << "x" << b.M() << std::endl;
        return false;
    }

    bool ok = true;
    for (auto row_it = a.begin(); row_it != a.end(); ++row_it) {
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            const size_t i = row_it.index();
            const size_t j = col_it.index();
            const double val_a = (*col_it)[0][0];

            double val_b = 0.0;
            auto b_row = b[i];
            auto b_col = b_row.find(j);
            if (b_col != b_row.end())
                val_b = (*b_col)[0][0];

            if (std::abs(val_a - val_b) > tol) {
                std::cerr << "  Mismatch at (" << i << "," << j << "): "
                          << val_a << " vs " << val_b
                          << " (diff=" << std::abs(val_a - val_b) << ")"
                          << std::endl;
                ok = false;
            }
        }
    }
    return ok;
}

// Compute residual R = A*p  where A is the original pressure matrix
Opm::BlockVector1 computeResidual(const Opm::FracturePressureInput& input)
{
    auto mat = Opm::assemblePressureOriginal(input);
    const size_t n = input.num_cells + input.num_well_equations;
    Opm::BlockVector1 p(n), res(n);
    for (size_t i = 0; i < n; ++i)
        p[i] = input.fracture_pressure[i];
    res = 0;
    mat->umv(p, res); // res += A * p
    return res;
}

// Simple 4-cell linear chain for testing
Opm::FracturePressureInput makeSimpleTestInput()
{
    Opm::FracturePressureInput input;
    input.num_cells = 4;
    input.min_width = 1e-4;

    // Linear chain: 0 -- 1 -- 2 -- 3
    // Htrans tuples follow the original convention: (nIdx, eIdx, geo_trans_e, geo_trans_n)
    // where eIdx < nIdx to avoid double-counting faces.
    input.htrans = {
        Opm::Htrans{1, 0, 2.0, 2.0},
        Opm::Htrans{2, 1, 1.5, 1.5},
        Opm::Htrans{3, 2, 1.0, 1.0}
    };

    // Widths (some well above min_width)
    input.fracture_width = {1e-3, 2e-3, 5e-4, 3e-3};

    // Pressures
    input.fracture_pressure = {1e6, 8e5, 1.2e6, 9e5};

    // Density and viscosity per cell (with zero pressure derivatives for
    // backward-compatible tests; mobility = density/viscosity = 1000/0.001 = 1e6).
    // Using constant properties here: rho=1000 kg/m3, mu=1e-4 Pa.s
    // gives mobility = 1e4 (same as the old reservoir_mobility values).
    input.density   = {{1000.0, 0.0}, {1000.0, 0.0}, {1000.0, 0.0}, {1000.0, 0.0}};
    input.viscosity = {{0.1, 0.0}, {0.1, 0.0}, {0.1, 0.0}, {0.1, 0.0}};

    // Leakoff
    input.leakof = {10.0, 20.0, 15.0, 25.0};

    // Default: rate control, no well equations
    input.control_type = "rate";
    input.num_well_equations = 0;

    return input;
}

// ============================================================================
// Tests
// ============================================================================

bool test_pressure_matrix_matches_original()
{
    std::cout << "Test 1: Comparing pressure matrices (AD vs original) ..."
              << std::endl;
    auto input = makeSimpleTestInput();

    auto ad_result = Opm::assemblePressureAD(input);
    auto orig_matrix = Opm::assemblePressureOriginal(input);

    if (!compareMatrices(*ad_result.pressure_matrix, *orig_matrix, 1e-12)) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_coupling_matrix_finite_differences()
{
    std::cout << "Test 2: Verifying coupling matrix via finite differences ..."
              << std::endl;
    auto input = makeSimpleTestInput();
    auto ad_result = Opm::assemblePressureAD(input);

    // base residual
    auto res_base = computeResidual(input);

    const double eps = 1e-7;
    bool ok = true;

    for (size_t k = 0; k < input.num_cells; ++k) {
        // perturb width[k]
        auto input_pert = input;
        input_pert.fracture_width[k] += eps;

        auto res_pert = computeResidual(input_pert);

        // check column k of the coupling matrix
        for (size_t i = 0; i < input.num_cells; ++i) {
            const double fd = (res_pert[i][0] - res_base[i][0]) / eps;

            double cmat_val = 0.0;
            auto row = (*ad_result.coupling_matrix)[i];
            auto col_it = row.find(k);
            if (col_it != row.end())
                cmat_val = (*col_it)[0][0];

            double rel_err = 0.0;
            const double scale = std::max(std::abs(fd), std::abs(cmat_val));
            if (scale > 1e-10)
                rel_err = std::abs(fd - cmat_val) / scale;

            if (rel_err > 1e-4) {
                std::cerr << "  FD mismatch at dR[" << i << "]/dw[" << k
                          << "]: AD=" << cmat_val << " FD=" << fd
                          << " rel_err=" << rel_err << std::endl;
                ok = false;
            }
        }
    }

    if (!ok) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_coupling_zero_below_min_width()
{
    std::cout << "Test 3: Coupling matrix zero when all widths < min_width ..."
              << std::endl;
    auto input = makeSimpleTestInput();

    // set all widths well below min_width
    input.min_width = 1e-2;
    input.fracture_width = {1e-4, 1e-4, 1e-4, 1e-4};

    auto ad_result = Opm::assemblePressureAD(input);

    bool ok = true;
    auto& cmat = *ad_result.coupling_matrix;
    for (auto row_it = cmat.begin(); row_it != cmat.end(); ++row_it) {
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            if (std::abs((*col_it)[0][0]) > 1e-15) {
                std::cerr << "  Non-zero coupling at ("
                          << row_it.index() << "," << col_it.index()
                          << "): " << (*col_it)[0][0] << std::endl;
                ok = false;
            }
        }
    }

    if (!ok) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_pressure_control()
{
    std::cout << "Test 4: Pressure matrix with pressure control ..."
              << std::endl;
    auto input = makeSimpleTestInput();
    input.control_type = "pressure";
    input.perfinj = {std::tuple<int,double>{0, 100.0}};

    auto ad_result = Opm::assemblePressureAD(input);
    auto orig_matrix = Opm::assemblePressureOriginal(input);

    if (!compareMatrices(*ad_result.pressure_matrix, *orig_matrix, 1e-12)) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_rate_well_control()
{
    std::cout << "Test 5: Pressure matrix with rate_well control ..."
              << std::endl;
    auto input = makeSimpleTestInput();
    input.control_type = "rate_well";
    input.num_well_equations = 1;
    input.perfinj = {std::tuple<int,double>{0, 50.0}};
    input.total_wellindex = 200.0;
    input.mobility_water_perf = 1000.0;

    // extra entry for well equation pressure
    input.fracture_pressure.push_back(1e6);

    auto ad_result = Opm::assemblePressureAD(input);
    auto orig_matrix = Opm::assemblePressureOriginal(input);

    if (!compareMatrices(*ad_result.pressure_matrix, *orig_matrix, 1e-12)) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_residual_consistency()
{
    std::cout << "Test 6: Residual R = A*p consistency ..."
              << std::endl;
    auto input = makeSimpleTestInput();

    auto ad_result = Opm::assemblePressureAD(input);

    // compute A*p using the pressure matrix from the AD result
    const size_t n = input.num_cells;
    Opm::BlockVector1 p(n), Ap(n);
    for (size_t i = 0; i < n; ++i)
        p[i] = input.fracture_pressure[i];
    Ap = 0;
    ad_result.pressure_matrix->umv(p, Ap);

    // compare with residual from the AD assembly
    bool ok = true;
    for (size_t i = 0; i < n; ++i) {
        const double diff = std::abs(Ap[i][0] - ad_result.residual[i][0]);
        if (diff > 1e-10) {
            std::cerr << "  Residual mismatch at " << i << ": "
                      << Ap[i][0] << " vs " << ad_result.residual[i][0]
                      << " (diff=" << diff << ")" << std::endl;
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

bool test_mobility_pressure_derivatives()
{
    std::cout << "Test 7: Mobility pressure derivatives via finite differences ..."
              << std::endl;

    // Use values where flow terms are dominant and well-resolved.
    // Large widths make transmissibility O(1); moderate pressures keep
    // the residual at a scale where FD perturbations are detectable.
    Opm::FracturePressureInput input;
    input.num_cells = 4;
    input.min_width = 1e-6;

    input.htrans = {
        Opm::Htrans{1, 0, 2.0, 2.0},
        Opm::Htrans{2, 1, 1.5, 1.5},
        Opm::Htrans{3, 2, 1.0, 1.0}
    };

    input.fracture_width = {0.05, 0.08, 0.04, 0.06};
    input.fracture_pressure = {2.0, 1.5, 2.5, 1.8};
    input.leakof = {0.0, 0.0, 0.0, 0.0}; // no leakoff to focus on flow terms
    input.control_type = "rate";
    input.num_well_equations = 0;

    // Density and viscosity with non-zero pressure derivatives
    input.density   = {{1000.0, 0.5}, {1010.0, 0.6}, {1020.0, 0.4}, {1030.0, 0.7}};
    input.viscosity = {{1e-3, -1e-7}, {1.1e-3, -1.2e-7}, {0.9e-3, -0.8e-7}, {1.2e-3, -1.1e-7}};

    auto ad_result = Opm::assemblePressureAD(input);
    auto res_base = computeResidual(input);

    const double eps = 1e-7;
    bool ok = true;

    for (size_t k = 0; k < input.num_cells; ++k) {
        auto input_pert = input;
        input_pert.fracture_pressure[k] += eps;
        input_pert.density[k].value   += input.density[k].dval_dp * eps;
        input_pert.viscosity[k].value += input.viscosity[k].dval_dp * eps;

        auto res_pert = computeResidual(input_pert);

        for (size_t i = 0; i < input.num_cells; ++i) {
            const double fd = (res_pert[i][0] - res_base[i][0]) / eps;

            double pmat_val = 0.0;
            auto row = (*ad_result.pressure_matrix)[i];
            auto col_it = row.find(k);
            if (col_it != row.end())
                pmat_val = (*col_it)[0][0];

            double rel_err = 0.0;
            const double scale = std::max(std::abs(fd), std::abs(pmat_val));
            if (scale > 1e-10)
                rel_err = std::abs(fd - pmat_val) / scale;

            if (rel_err > 1e-4) {
                std::cerr << "  FD mismatch at dR[" << i << "]/dp[" << k
                          << "]: AD=" << pmat_val << " FD=" << fd
                          << " rel_err=" << rel_err << std::endl;
                ok = false;
            }
        }
    }

    if (!ok) {
        std::cerr << "  FAILED" << std::endl;
        return false;
    }
    std::cout << "  PASSED" << std::endl;
    return true;
}

bool test_varying_fluid_properties()
{
    std::cout << "Test 8: Non-uniform density/viscosity across cells ..."
              << std::endl;
    auto input = makeSimpleTestInput();

    // Set different density/viscosity per cell (no pressure derivatives)
    input.density   = {{800.0, 0.0}, {900.0, 0.0}, {1000.0, 0.0}, {1100.0, 0.0}};
    input.viscosity = {{1e-1, 0.0}, {2e-1, 0.0}, {5e-2, 0.0}, {1e-1, 0.0}};

    auto ad_result = Opm::assemblePressureAD(input);
    auto orig_matrix = Opm::assemblePressureOriginal(input);

    if (!compareMatrices(*ad_result.pressure_matrix, *orig_matrix, 1e-12)) {
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

    if (!test_pressure_matrix_matches_original())  ++failures;
    if (!test_coupling_matrix_finite_differences()) ++failures;
    if (!test_coupling_zero_below_min_width())      ++failures;
    if (!test_pressure_control())                   ++failures;
    if (!test_rate_well_control())                  ++failures;
    if (!test_residual_consistency())               ++failures;
    if (!test_mobility_pressure_derivatives())      ++failures;
    if (!test_varying_fluid_properties())           ++failures;

    std::cout << "\n=== "
              << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << " ===" << std::endl;

    return failures;
}
