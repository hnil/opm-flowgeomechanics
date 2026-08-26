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
        const double inv_rhs = 1.0 / rhs.value;
        r.value = value * inv_rhs;
        const double inv2 = inv_rhs * inv_rhs;
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
        const double inv_a = 1.0 / a.value;
        r.value = s * inv_a;
        const double neg_s_over_a2 = -(s * inv_a) * inv_a;
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
    std::vector<double> face_mobility;
    std::vector<CellFluidProperty> density;    // per-cell density with pressure derivative
    std::vector<CellFluidProperty> viscosity;  // per-cell viscosity with pressure derivative
    std::vector<double> leakof;
    double min_width = 1e-6;
    size_t num_cells = 0;

    // Well data
    std::vector<std::tuple<int, double>> perfinj;
    std::string control_type = "rate";
    double total_wellindex = 0.0;
    double well_target_pressure = 0.0; // bhp_well: p_w target at the perf-ref datum
    double mobility_water_perf = 1.0;
    size_t num_well_equations = 0;
    bool use_fluid_mobility_for_internal_faces = false;
    bool use_fluid_mobility_for_well_connections = false;
};

inline double
cellMobilityValue(const CellFluidProperty& density,
                  const CellFluidProperty& viscosity)
{
    assert(viscosity.value != 0.0);
    return density.value / viscosity.value;
}

inline double
constantFaceMobilityValue(const FracturePressureInput& input,
                          const size_t cell)
{
    if (cell < input.face_mobility.size())
        return input.face_mobility[cell];
    return cellMobilityValue(input.density[cell], input.viscosity[cell]);
}

template <int NumDerivs>
LocalAD<NumDerivs>
cellMobilityAD(const CellFluidProperty& density,
               const CellFluidProperty& viscosity,
               const int pressure_deriv_index)
{
    LocalAD<NumDerivs> rho(density.value);
    rho.derivatives[pressure_deriv_index] = density.dval_dp;
    LocalAD<NumDerivs> mu(viscosity.value);
    mu.derivatives[pressure_deriv_index] = viscosity.dval_dp;
    return rho / mu;
}

inline double
faceMobilityValue(const FracturePressureInput& input,
                  const size_t i,
                  const size_t j)
{
    if (input.use_fluid_mobility_for_internal_faces) {
        const double mob_i = cellMobilityValue(input.density[i], input.viscosity[i]);
        const double mob_j = cellMobilityValue(input.density[j], input.viscosity[j]);
        return 0.5 * (mob_i + mob_j);
    }

    const double mob_i = constantFaceMobilityValue(input, i);
    const double mob_j = constantFaceMobilityValue(input, j);
    return 0.5 * (mob_i + mob_j);
}

template <int NumDerivs>
LocalAD<NumDerivs>
faceMobilityAD(const FracturePressureInput& input,
               const size_t i,
               const size_t j,
               const int deriv_i,
               const int deriv_j)
{
    using AD = LocalAD<NumDerivs>;
    if (!input.use_fluid_mobility_for_internal_faces)
        return AD::constant(faceMobilityValue(input, i, j));

    return AD::constant(0.5)
        * (cellMobilityAD<NumDerivs>(input.density[i], input.viscosity[i], deriv_i)
         + cellMobilityAD<NumDerivs>(input.density[j], input.viscosity[j], deriv_j));
}

inline double
wellConnectionMobilityValue(const FracturePressureInput& input,
                            const int cell)
{
    if (input.use_fluid_mobility_for_well_connections)
        return cellMobilityValue(input.density[cell], input.viscosity[cell]);
    return input.mobility_water_perf;
}

template <int NumDerivs>
LocalAD<NumDerivs>
wellConnectionMobilityAD(const FracturePressureInput& input,
                        const int cell,
                        const int pressure_deriv_index)
{
    using AD = LocalAD<NumDerivs>;
    if (!input.use_fluid_mobility_for_well_connections)
        return AD::constant(input.mobility_water_perf);
    return cellMobilityAD<NumDerivs>(input.density[cell], input.viscosity[cell], pressure_deriv_index);
}

template <int NumDerivs>
void
assembleRateWellControlAD(const FracturePressureInput& input,
                          BCRSMatrix1x1& pmat,
                          BlockVector1& residual)
{
    using AD = LocalAD<NumDerivs>;
    constexpr int CELL_P = 0;
    constexpr int WELL_P = 1;

    assert(input.num_well_equations == 1);
    const size_t total_size = input.num_cells + input.num_well_equations;
    const size_t weqix = total_size - 1;
    const double WI_lambda = input.total_wellindex;

    pmat[weqix][weqix] += WI_lambda;
    residual[weqix] += WI_lambda * input.fracture_pressure[weqix];

    for (const auto& pi : input.perfinj) {
        const int cell = std::get<0>(pi);
        AD cell_pressure = AD::variable(input.fracture_pressure[cell], CELL_P);
        AD well_pressure = AD::variable(input.fracture_pressure[weqix], WELL_P);
        AD mobility = wellConnectionMobilityAD<NumDerivs>(input, cell, CELL_P);
        AD flux = AD::constant(std::get<1>(pi)) * mobility * (well_pressure - cell_pressure);

        residual[weqix] += flux.value;
        residual[cell] -= flux.value;

        pmat[weqix][cell] += flux.derivatives[CELL_P];
        pmat[weqix][weqix] += flux.derivatives[WELL_P];
        pmat[cell][cell] -= flux.derivatives[CELL_P];
        pmat[cell][weqix] -= flux.derivatives[WELL_P];
    }
}

// Constraint-row scale for the bhp_well row: the natural magnitude is the
// matrix-perforation WI sum (same as the rate row's leading diagonal term).
inline double
bhpWellRowScale(const FracturePressureInput& input)
{
    return std::max(input.total_wellindex, 1e-14);
}

// bhp_well: the well DOF is kept, the cell<->well fWI couplings are identical
// to the rate row, but the WELL row is the constraint  scale*(p_w - target) = 0.
// A control switch (rate <-> bhp) therefore changes one row's content only.
template <int NumDerivs>
void
assembleBhpWellControlAD(const FracturePressureInput& input,
                         BCRSMatrix1x1& pmat,
                         BlockVector1& residual)
{
    using AD = LocalAD<NumDerivs>;
    constexpr int CELL_P = 0;
    constexpr int WELL_P = 1;

    assert(input.num_well_equations == 1);
    const size_t total_size = input.num_cells + input.num_well_equations;
    const size_t weqix = total_size - 1;

    const double scale = bhpWellRowScale(input);
    pmat[weqix][weqix] += scale;
    residual[weqix] += scale * (input.fracture_pressure[weqix] - input.well_target_pressure);

    for (const auto& pi : input.perfinj) {
        const int cell = std::get<0>(pi);
        AD cell_pressure = AD::variable(input.fracture_pressure[cell], CELL_P);
        AD well_pressure = AD::variable(input.fracture_pressure[weqix], WELL_P);
        AD mobility = wellConnectionMobilityAD<NumDerivs>(input, cell, CELL_P);
        AD flux = AD::constant(std::get<1>(pi)) * mobility * (well_pressure - cell_pressure);

        // cells receive the flux; the well row is the constraint, not a balance
        residual[cell] -= flux.value;
        pmat[cell][cell] -= flux.derivatives[CELL_P];
        pmat[cell][weqix] -= flux.derivatives[WELL_P];
    }
}

inline void
assembleBhpWellControlOriginal(const FracturePressureInput& input,
                               BCRSMatrix1x1& matrix)
{
    assert(input.num_well_equations == 1);
    const size_t total_size = input.num_cells + input.num_well_equations;
    const size_t weqix = total_size - 1;

    matrix[weqix][weqix] += bhpWellRowScale(input);

    for (const auto& pi : input.perfinj) {
        const int cell = std::get<0>(pi);
        const double value = std::get<1>(pi) * wellConnectionMobilityValue(input, cell);
        matrix[cell][weqix] -= value;
        matrix[cell][cell] += value;
    }
}

inline void
assembleRateWellControlOriginal(const FracturePressureInput& input,
                                BCRSMatrix1x1& matrix)
{
    assert(input.num_well_equations == 1);
    const size_t total_size = input.num_cells + input.num_well_equations;
    const size_t weqix = total_size - 1;
    const double WI_lambda = input.total_wellindex;

    matrix[weqix][weqix] += WI_lambda;

    for (const auto& pi : input.perfinj) {
        const int cell = std::get<0>(pi);
        const double value = std::get<1>(pi) * wellConnectionMobilityValue(input, cell);
        matrix[weqix][cell] -= value;
        matrix[weqix][weqix] += value;
        matrix[cell][weqix] -= value;
        matrix[cell][cell] += value;
    }
}

// Output from the AD-based assembly
struct PressureAssemblyADResult
{
    std::unique_ptr<BCRSMatrix1x1> pressure_matrix;  // dR/dp
    std::unique_ptr<BCRSMatrix1x1> coupling_matrix;  // dR/dw (aperture coupling)
    BlockVector1 residual;                            // R(p, w)
};

// dR/dw has one column per fracture cell, so it is (nc+nw) x nc - NOT square.
// A square version makes every matrix-vector product read one entry past the
// end of the width vector (garbage: 0*NaN = NaN in the Krylov solve).
inline std::unique_ptr<BCRSMatrix1x1>
buildCouplingMatrixStructure(const std::vector<Htrans>& htrans,
                             size_t nc,
                             size_t num_well_equations)
{
    auto mat = std::make_unique<BCRSMatrix1x1>(
        nc + num_well_equations, nc, 4, 0.4, BCRSMatrix1x1::implicit);

    for (const auto& ht : htrans) {
        const size_t i = std::get<0>(ht);
        const size_t j = std::get<1>(ht);
        mat->entry(i, j) = 0.0;
        mat->entry(j, i) = 0.0;
        mat->entry(i, i) = 0.0;
        mat->entry(j, j) = 0.0;
    }
    for (size_t i = 0; i < nc; ++i)
        mat->entry(i, i) = 0.0;

    mat->compress();
    return mat;
}

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
    // implicit build: avg 4 entries/row; the overflow buffer must also hold the
    // well row, which has one entry per well-source cell (many with
    // solver.well_source_all_perfs), so size it from perfinj
    const double overflow = 0.4 + (total_size > 0 ? 2.0 * perfinj.size() / static_cast<double>(total_size) : 0.0);
    auto mat = std::make_unique<BCRSMatrix1x1>(total_size, total_size, 4, overflow,
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
    result.coupling_matrix
        = buildCouplingMatrixStructure(input.htrans, nc, input.num_well_equations);
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

        AD4 mobility = faceMobilityAD<4>(input, i, j, P_I, P_J);
        AD4 trans = mobility * (AD4::constant(1.0) / inv_trans);

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
        assembleRateWellControlAD<2>(input, pmat, result.residual);
    } else if (input.control_type == "bhp_well") {
        assembleBhpWellControlAD<2>(input, pmat, result.residual);
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
        const double mobility = faceMobilityValue(input, i, j);
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
        assembleRateWellControlOriginal(input, matrix);
    } else if (input.control_type == "bhp_well") {
        assembleBhpWellControlOriginal(input, matrix);
    }

    return mat;
}

} // namespace Opm

#endif // OPM_FRACTURE_PRESSURE_ASSEMBLER_AD_HPP
