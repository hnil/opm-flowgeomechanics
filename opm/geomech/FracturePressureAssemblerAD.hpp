/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 NTNU

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

#ifndef OPM_FRACTURE_PRESSURE_ASSEMBLER_AD_HPP
#define OPM_FRACTURE_PRESSURE_ASSEMBLER_AD_HPP

#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

namespace Opm
{

// ============================================================================
// Simple forward-mode automatic differentiation scalar.
//
// This is a minimal AD type similar to, but simpler than, the DenseAd
// Evaluation type used in opm-simulators.  It tracks a scalar value and
// a fixed number of partial derivatives (determined at compile time).
// ============================================================================
template <int NumDerivs>
class LocalAD
{
public:
    double value;
    std::array<double, NumDerivs> derivatives;

    LocalAD() : value(0.0) { derivatives.fill(0.0); }
    explicit LocalAD(double v) : value(v) { derivatives.fill(0.0); }

    static LocalAD variable(double v, int deriv_index)
    {
        LocalAD result(v);
        result.derivatives[deriv_index] = 1.0;
        return result;
    }

    static LocalAD constant(double v)
    {
        return LocalAD(v);
    }

    // --- arithmetic operators -----------------------------------------------

    LocalAD operator+(const LocalAD& rhs) const
    {
        LocalAD r;
        r.value = value + rhs.value;
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = derivatives[i] + rhs.derivatives[i];
        return r;
    }

    LocalAD operator-(const LocalAD& rhs) const
    {
        LocalAD r;
        r.value = value - rhs.value;
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = derivatives[i] - rhs.derivatives[i];
        return r;
    }

    LocalAD operator*(const LocalAD& rhs) const
    {
        LocalAD r;
        r.value = value * rhs.value;
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = derivatives[i] * rhs.value + value * rhs.derivatives[i];
        return r;
    }

    LocalAD operator/(const LocalAD& rhs) const
    {
        assert(rhs.value != 0.0);
        LocalAD r;
        r.value = value / rhs.value;
        const double inv2 = 1.0 / (rhs.value * rhs.value);
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = (derivatives[i] * rhs.value - value * rhs.derivatives[i]) * inv2;
        return r;
    }

    LocalAD& operator+=(const LocalAD& rhs)
    {
        value += rhs.value;
        for (int i = 0; i < NumDerivs; ++i)
            derivatives[i] += rhs.derivatives[i];
        return *this;
    }

    LocalAD& operator-=(const LocalAD& rhs)
    {
        value -= rhs.value;
        for (int i = 0; i < NumDerivs; ++i)
            derivatives[i] -= rhs.derivatives[i];
        return *this;
    }

    LocalAD operator-() const
    {
        LocalAD r;
        r.value = -value;
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = -derivatives[i];
        return r;
    }

    // scalar * AD
    friend LocalAD operator*(double s, const LocalAD& a)
    {
        LocalAD r;
        r.value = s * a.value;
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = s * a.derivatives[i];
        return r;
    }

    // scalar / AD
    friend LocalAD operator/(double s, const LocalAD& a)
    {
        assert(a.value != 0.0);
        LocalAD r;
        r.value = s / a.value;
        const double neg_s_over_a2 = -s / (a.value * a.value);
        for (int i = 0; i < NumDerivs; ++i)
            r.derivatives[i] = neg_s_over_a2 * a.derivatives[i];
        return r;
    }
};

// max(ad, constant_threshold):  passes through derivatives when ad > threshold
template <int N>
LocalAD<N> adMax(const LocalAD<N>& a, double threshold)
{
    if (a.value >= threshold)
        return a;
    else
        return LocalAD<N>::constant(threshold);
}

// ============================================================================
// Data types used by the assembler
// ============================================================================

using Htrans = std::tuple<size_t, size_t, double, double>;
using BCRSMatrix1x1 = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using BlockVector1 = Dune::BlockVector<Dune::FieldVector<double, 1>>;

// Precomputed fluid property for a single fracture cell, carrying the
// value and its derivative with respect to the fracture cell pressure.
struct CellFluidProperty
{
    double value = 0.0;
    double dval_dp = 0.0;  // d(value)/d(fracture_pressure)
};

// Input data for the pressure assembly (standalone, independent of Fracture class)
struct FracturePressureInput
{
    std::vector<Htrans> htrans;
    std::vector<double> fracture_width;
    std::vector<double> fracture_pressure;
    std::vector<CellFluidProperty> density;    // per-cell density with pressure derivative
    std::vector<CellFluidProperty> viscosity;  // per-cell viscosity with pressure derivative
    std::vector<double> leakof;
    double min_width = 1e-6;
    size_t num_cells = 0;

    // Well data
    std::vector<std::tuple<int, double>> perfinj;
    std::string control_type = "rate";
    double total_wellindex = 0.0;
    double mobility_water_perf = 1.0;
    size_t num_well_equations = 0;
};

// Output from the AD-based assembly
struct PressureAssemblyADResult
{
    std::unique_ptr<BCRSMatrix1x1> pressure_matrix;  // dR/dp
    std::unique_ptr<BCRSMatrix1x1> coupling_matrix;  // dR/dw (aperture coupling)
    BlockVector1 residual;                            // R(p, w)
};

// ============================================================================
// Build the sparse matrix sparsity structure from connectivity information
// ============================================================================
inline std::unique_ptr<BCRSMatrix1x1>
buildMatrixStructure(const std::vector<Htrans>& htrans,
                     size_t nc,
                     const std::vector<std::tuple<int, double>>& perfinj = {},
                     size_t num_well_equations = 0)
{
    const size_t total_size = nc + num_well_equations;
    auto mat = std::make_unique<BCRSMatrix1x1>(total_size, total_size, 4, 0.4,
                                                BCRSMatrix1x1::implicit);

    for (const auto& ht : htrans) {
        const size_t i = std::get<0>(ht);
        const size_t j = std::get<1>(ht);
        mat->entry(i, j) = 0.0;
        mat->entry(j, i) = 0.0;
        mat->entry(i, i) = 0.0;
        mat->entry(j, j) = 0.0;
    }

    // ensure every diagonal exists
    for (size_t i = 0; i < total_size; ++i)
        mat->entry(i, i) = 0.0;

    if (num_well_equations > 0) {
        const size_t weqix = total_size - 1;
        for (const auto& pi : perfinj) {
            const int cell = std::get<0>(pi);
            mat->entry(cell, weqix) = 0.0;
            mat->entry(weqix, cell) = 0.0;
            mat->entry(cell, cell) = 0.0;
        }
        mat->entry(weqix, weqix) = 0.0;
    }

    mat->compress();
    return mat;
}

// ============================================================================
// Assemble the pressure equations using forward AD.
//
// For each internal face the flux is computed with AD scalars that track
// derivatives with respect to the two local pressures and the two local
// apertures.  The derivatives are then scattered into global matrices
// for the pressure Jacobian (dR/dp) and the aperture coupling (dR/dw).
//
// The residual is: R_i = sum_j T_ij(w) * (p_i - p_j) + leakof_i * p_i
//                        + [well coupling terms]
// ============================================================================
inline PressureAssemblyADResult
assemblePressureAD(const FracturePressureInput& input)
{
    const size_t nc = input.num_cells;
    const size_t total_size = nc + input.num_well_equations;

    PressureAssemblyADResult result;
    result.pressure_matrix = buildMatrixStructure(
        input.htrans, nc, input.perfinj, input.num_well_equations);
    result.coupling_matrix = buildMatrixStructure(
        input.htrans, nc, input.perfinj, input.num_well_equations);
    result.residual.resize(total_size);
    result.residual = 0;

    auto& pmat = *result.pressure_matrix;
    auto& cmat = *result.coupling_matrix;
    pmat = 0;
    cmat = 0;

    // Local AD type with 4 derivatives: p_i, p_j, w_i, w_j
    using AD4 = LocalAD<4>;
    constexpr int P_I = 0, P_J = 1, W_I = 2, W_J = 3;

    // ----- flow between fracture cells (face loop) -----
    for (const auto& ht : input.htrans) {
        const size_t i = std::get<0>(ht);
        const size_t j = std::get<1>(ht);
        const double t1 = std::get<2>(ht);
        const double t2 = std::get<3>(ht);

        // create AD variables for the local face stencil
        AD4 p_i_ad = AD4::variable(input.fracture_pressure[i], P_I);
        AD4 p_j_ad = AD4::variable(input.fracture_pressure[j], P_J);
        AD4 w_i_ad = AD4::variable(input.fracture_width[i], W_I);
        AD4 w_j_ad = AD4::variable(input.fracture_width[j], W_J);

        // apply minimum width (matches original std::max(width, min_width))
        AD4 h1 = adMax(w_i_ad, input.min_width);
        AD4 h2 = adMax(w_j_ad, input.min_width);

        // Poiseuille (cubic law) transmissibility – same formula as Fracture::assemblePressure
        AD4 h1_cubed = h1 * h1 * h1;
        AD4 h2_cubed = h2 * h2 * h2;
        AD4 inv_trans = AD4::constant(12.0) / (h1_cubed * AD4::constant(t1))
                      + AD4::constant(12.0) / (h2_cubed * AD4::constant(t2));

        // Compute per-cell mobility = density / viscosity as AD quantities.
        // Pressure derivatives of density and viscosity are mapped to the
        // local derivative slots P_I (for cell i) and P_J (for cell j).
        AD4 rho_i(input.density[i].value);
        rho_i.derivatives[P_I] = input.density[i].dval_dp;
        AD4 mu_i(input.viscosity[i].value);
        mu_i.derivatives[P_I] = input.viscosity[i].dval_dp;
        AD4 mob_i = rho_i / mu_i;

        AD4 rho_j(input.density[j].value);
        rho_j.derivatives[P_J] = input.density[j].dval_dp;
        AD4 mu_j(input.viscosity[j].value);
        mu_j.derivatives[P_J] = input.viscosity[j].dval_dp;
        AD4 mob_j = rho_j / mu_j;

        AD4 mobility = AD4::constant(0.5) * (mob_i + mob_j);
        AD4 trans = mobility / inv_trans;

        // flux from cell i to cell j
        AD4 flux = trans * (p_i_ad - p_j_ad);

        // scatter value into residual
        result.residual[i] += flux.value;
        result.residual[j] -= flux.value;

        // scatter pressure derivatives (dR/dp)
        pmat[i][i] += flux.derivatives[P_I];
        pmat[i][j] += flux.derivatives[P_J];
        pmat[j][i] -= flux.derivatives[P_I];
        pmat[j][j] -= flux.derivatives[P_J];

        // scatter width/aperture coupling derivatives (dR/dw)
        cmat[i][i] += flux.derivatives[W_I];
        cmat[i][j] += flux.derivatives[W_J];
        cmat[j][i] -= flux.derivatives[W_I];
        cmat[j][j] -= flux.derivatives[W_J];
    }

    // ----- leakoff (diagonal, no width dependence) -----
    for (size_t i = 0; i < input.leakof.size() && i < nc; ++i) {
        pmat[i][i] += input.leakof[i];
        result.residual[i] += input.leakof[i] * input.fracture_pressure[i];
    }

    // ----- well coupling -----
    if (input.control_type == "pressure" || input.control_type == "perf_pressure") {
        for (const auto& perfinj : input.perfinj) {
            const int cell = std::get<0>(perfinj);
            const double value = std::get<1>(perfinj);
            pmat[cell][cell] += value;
            result.residual[cell] += value * input.fracture_pressure[cell];
        }
    } else if (input.control_type == "rate_well") {
        assert(input.num_well_equations == 1);
        const double lambda = 1.0;
        const double WI_lambda = input.total_wellindex * lambda;
        pmat[total_size - 1][total_size - 1] = WI_lambda;

        for (const auto& pi : input.perfinj) {
            const int cell = std::get<0>(pi);
            const double value = std::get<1>(pi) * input.mobility_water_perf;
            pmat[total_size - 1][cell] = -value;
            pmat[total_size - 1][total_size - 1] += value;
            pmat[cell][total_size - 1] = -value;
            pmat[cell][cell] += value;
        }
    }

    return result;
}

// ============================================================================
// Assemble the pressure matrix using the original (non-AD) method.
// This is a standalone version of Fracture::assemblePressure(), suitable
// for comparison in unit tests.
// ============================================================================
inline std::unique_ptr<BCRSMatrix1x1>
assemblePressureOriginal(const FracturePressureInput& input)
{
    const size_t nc = input.num_cells;
    const size_t total_size = nc + input.num_well_equations;

    auto mat = buildMatrixStructure(
        input.htrans, nc, input.perfinj, input.num_well_equations);
    auto& matrix = *mat;
    matrix = 0;

    // flow terms (same as Fracture::assemblePressure)
    for (const auto& ht : input.htrans) {
        const size_t i = std::get<0>(ht);
        const size_t j = std::get<1>(ht);
        const double t1 = std::get<2>(ht);
        const double t2 = std::get<3>(ht);

        const double h1 = std::max(input.fracture_width[i], input.min_width);
        const double h2 = std::max(input.fracture_width[j], input.min_width);

        double value = 12.0 / (h1 * h1 * h1 * t1) + 12.0 / (h2 * h2 * h2 * t2);
        const double mob_i = input.density[i].value / input.viscosity[i].value;
        const double mob_j = input.density[j].value / input.viscosity[j].value;
        const double mobility = 0.5 * (mob_i + mob_j);
        value = 1.0 / value;
        value *= mobility;

        matrix[i][j] -= value;
        matrix[j][i] -= value;
        matrix[i][i] += value;
        matrix[j][j] += value;
    }

    // leakoff
    for (size_t i = 0; i < input.leakof.size() && i < nc; ++i) {
        matrix[i][i] += input.leakof[i];
    }

    // well coupling
    if (input.control_type == "pressure" || input.control_type == "perf_pressure") {
        for (const auto& perfinj : input.perfinj) {
            const int cell = std::get<0>(perfinj);
            const double value = std::get<1>(perfinj);
            matrix[cell][cell] += value;
        }
    } else if (input.control_type == "rate_well") {
        assert(input.num_well_equations == 1);
        const double lambda = 1.0;
        const double WI_lambda = input.total_wellindex * lambda;
        matrix[total_size - 1][total_size - 1] = WI_lambda;

        for (const auto& pi : input.perfinj) {
            const int cell = std::get<0>(pi);
            const double value = std::get<1>(pi) * input.mobility_water_perf;
            matrix[total_size - 1][cell] = -value;
            matrix[total_size - 1][total_size - 1] += value;
            matrix[cell][total_size - 1] = -value;
            matrix[cell][cell] += value;
        }
    }

    return mat;
}

} // namespace Opm

#endif // OPM_FRACTURE_PRESSURE_ASSEMBLER_AD_HPP
