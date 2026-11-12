/*
  Copyright 2025 Equinor ASA.

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

#include "config.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <numeric>

#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>


#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/DiagonalScalar.hpp>
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

#include <fmt/format.h>

using namespace std;
namespace
{
static int DEBUG_COUNT = 0;
[[maybe_unused]] string
debug_filename(const string& prefix, const string& suffix = ".txt")
{
    std::ostringstream oss;
    oss << prefix << DEBUG_COUNT++ << suffix;
    return oss.str();
}

// ========================== Convenience definitions ==========================
using Dune::Indices::_0;
using Dune::Indices::_1;

using ResVector = Opm::Fracture::Vector;
using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
using FMatrix = Dune::DynamicMatrix<double>; // full matrix

using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;

using Htrans = std::tuple<size_t, size_t, double, double>;
using LinearPrecondType = Opm::FractureMechanicsPreconditioner;
using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;

// ----------------------------------------------------------------------------
// Block scaling for the coupled fracture linear system.
//
// The mechanics block A (stiffness/traction, O(1e6-1e10)) and the flow block M
// (transmissibility, O(1e-6)) differ by many orders of magnitude, which wrecks
// the conditioning of the coupled Krylov solve. A symmetric diagonal scaling
//   D = diag(s_m on mech rows, s_p on pressure rows),  S -> D S D,  b -> D b,
// solving for y and recovering x = D y, normalises the two blocks. This mirrors
// the validated --scaling=blockmax path in examples/replay_fracture_linear_system.cpp
// (see tests/replay_sweep.sh). s_m = 1/sqrt(max|A|), s_p = 1/sqrt(max|M|).
// ----------------------------------------------------------------------------
struct BlockScaling
{
    double mechanics = 1.0;
    double pressure = 1.0;
    bool active() const { return mechanics != 1.0 || pressure != 1.0; }
};

double
max_abs_entry(const FMatrix& m)
{
    double res = 0.0;
    for (size_t i = 0; i != m.N(); ++i)
        for (size_t j = 0; j != m.M(); ++j)
            res = std::max(res, std::abs(m[i][j]));
    return res;
}

double
max_abs_entry(const SMatrix& m)
{
    double res = 0.0;
    for (auto row = m.begin(); row != m.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            res = std::max(res, std::abs((*col)[0][0]));
    return res;
}

void
scale_sparse(SMatrix& m, const double factor)
{
    for (auto row = m.begin(); row != m.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            (*col)[0][0] *= factor;
}

void
scale_dense(FMatrix& m, const double factor)
{
    for (size_t i = 0; i != m.N(); ++i)
        for (size_t j = 0; j != m.M(); ++j)
            m[i][j] *= factor;
}

// Bit-exact FNV-1a hash of a double sequence; used only by the
// solver.linsolver.debug_checksum diagnostic to tell whether two runs assemble
// the same coupled system.
struct Hash64
{
    uint64_t h {1469598103934665603ull};
    void operator()(double v)
    {
        uint64_t b;
        std::memcpy(&b, &v, sizeof b);
        h = (h ^ b) * 1099511628211ull;
    }
};

template <class SM>
void
hash_sparse(Hash64& h, const SM& m)
{
    for (auto row = m.begin(); row != m.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col) {
            h(static_cast<double>(col.index()));
            h((*col)[0][0]);
        }
}

// Compute s_m, s_p from the current diagonal blocks for the requested mode.
// "none" (default) returns identity scaling.
BlockScaling
compute_block_scaling(const SystemMatrix& S, const std::string& mode)
{
    if (mode == "none")
        return {};
    if (mode != "blockmax")
        OPM_THROW(std::runtime_error, "Unsupported fracture linsolver scaling mode: " + mode);

    const double a_max = max_abs_entry(S[_0][_0]);
    const double m_max = max_abs_entry(S[_1][_1]);
    if (a_max <= 0.0 || m_max <= 0.0)
        return {}; // degenerate block; fall back to no scaling rather than throw
    return {1.0 / std::sqrt(a_max), 1.0 / std::sqrt(m_max)};
}

// Scale the four blocks of S by D _ D (matrix only). fac=+1 applies the scaling,
// fac=-1 (reciprocal) restores the original matrix.
void
scale_system_matrix(SystemMatrix& S, const BlockScaling& s, const int sign)
{
    const double mm = (sign > 0) ? s.mechanics * s.mechanics : 1.0 / (s.mechanics * s.mechanics);
    const double mp = (sign > 0) ? s.mechanics * s.pressure : 1.0 / (s.mechanics * s.pressure);
    const double pp = (sign > 0) ? s.pressure * s.pressure : 1.0 / (s.pressure * s.pressure);
    scale_dense(S[_0][_0], mm);
    scale_sparse(S[_0][_1], mp);
    scale_sparse(S[_1][_0], mp);
    scale_sparse(S[_1][_1], pp);
}

// Apply D S D to the system and D b to the rhs (in place).
void
apply_block_scaling(SystemMatrix& S, VectorHP& rhs, const BlockScaling& s)
{
    scale_system_matrix(S, s, +1);
    for (auto& v : rhs[_0]) v[0] *= s.mechanics;
    for (auto& v : rhs[_1]) v[0] *= s.pressure;
}

// Direct dense solve of the coupled system S dx = rhs. The mechanics block A is a
// dense BEM (displacement-discontinuity) influence matrix, so the coupled system
// is well posed but extremely ill-conditioned (width O(1) vs pressure O(1e8-1e11)),
// which defeats Krylov methods and their cheap preconditioners. For SMALL fractures
// a direct dense LU of the assembled 2N x 2N system is robust and affordable; an
// exact Schur-complement preconditioner would cost the same order (A^-1 is dense).
// Cost ~O((2N)^3) time, O((2N)^2) memory — guard on N before calling.
bool
solve_coupled_direct(const SystemMatrix& S, const VectorHP& rhs, VectorHP& dx)
{
    const size_t n = S[_0][_0].N();
    const size_t m = S[_1][_1].N(); // n + numWellEquations (rate_well appends a row)
    FMatrix K(n + m, n + m, 0.0);
    const auto& A = S[_0][_0];
    for (size_t i = 0; i != n; ++i)
        for (size_t j = 0; j != n; ++j)
            K[i][j] = A[i][j];
    auto fill = [&](const SMatrix& Sb, const size_t ro, const size_t co) {
        for (auto r = Sb.begin(); r != Sb.end(); ++r)
            for (auto c = r->begin(); c != r->end(); ++c)
                K[ro + r.index()][co + c.index()] = (*c)[0][0];
    };
    fill(S[_0][_1], 0, n);
    fill(S[_1][_0], n, 0);
    fill(S[_1][_1], n, n);

    ResVector b(n + m), y(n + m);
    for (size_t i = 0; i != n; ++i) b[i] = rhs[_0][i];
    for (size_t i = 0; i != m; ++i) b[n + i] = rhs[_1][i];
    y = 0;
    try {
        K.solve(y, b); // dense LU (factorises K in place)
    } catch (const Dune::Exception&) {
        return false; // singular
    }
    for (size_t i = 0; i != n + m; ++i)
        if (!std::isfinite(y[i][0]))
            return false;
    for (size_t i = 0; i != n; ++i) dx[_0][i] = y[i];
    for (size_t i = 0; i != m; ++i) dx[_1][i] = y[n + i];
    return true;
}

void
validate_state_vector(const ResVector& values, const char* const label, const double max_abs_value)
{
    for (size_t i = 0; i < values.size(); ++i) {
        const double value = values[i][0];
        if (!std::isfinite(value)) {
            OPM_THROW(std::runtime_error,
                      std::string("Fracture nonlinear state contains non-finite ") + label
                          + " at index " + std::to_string(i));
        }
        if (std::abs(value) > max_abs_value) {
            OPM_THROW(std::runtime_error,
                      std::string("Fracture nonlinear state contains excessive ") + label
                          + " at index " + std::to_string(i)
                          + ": " + std::to_string(value));
        }
    }
}

void
validate_fracture_state(const VectorHP& x, const Opm::PropertyTree& prm)
{
    if (!prm.get<bool>("solver.guard_state", true)) {
        return;
    }
    const double max_abs_width = prm.get<double>("solver.guard_max_abs_width", 1e30);
    const double max_abs_pressure = prm.get<double>("solver.guard_max_abs_pressure", 1e30);
    validate_state_vector(x[_0], "width", max_abs_width);
    validate_state_vector(x[_1], "pressure", max_abs_pressure);
}

size_t
count_toggled_cells(const std::vector<int>& previous, const std::vector<int>& current)
{
    if (previous.size() != current.size()) {
        return current.size();
    }
    size_t count = 0;
    for (size_t i = 0; i < current.size(); ++i) {
        count += (previous[i] != current[i]);
    }
    return count;
}

// Count only cells that re-open (closed -> open). Establishing contact
// (open -> closed) must never be throttled: the contact set legitimately needs
// to grow by hundreds of cells in a single iteration when a fracture starts
// fully compressed, and throttling that growth freezes the closed set at empty
// and prevents the mechanics residual from ever converging. Chattering, the
// behaviour the toggle guard is meant to suppress, shows up as re-openings, so
// the guard only counts those.
size_t
count_reopened_cells(const std::vector<int>& previous, const std::vector<int>& current)
{
    if (previous.size() != current.size()) {
        return 0;
    }
    size_t count = 0;
    for (size_t i = 0; i < current.size(); ++i) {
        count += (previous[i] != 0 && current[i] == 0);
    }
    return count;
}

// ============================ Debugging functions ============================
// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_matrix(const FMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    const size_t rows = m.N();
    const size_t cols = m.M();
    std::ofstream os(name);
    for (size_t i = 0; i != rows; ++i)
        for (size_t j = 0; j != cols; ++j)
            os << m[i][j] << ((j == cols - 1) ? "\n" : " ");
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_matrix(const SMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    std::ofstream os(name);
    for (auto rowIt = m.begin(); rowIt != m.end(); ++rowIt) {
        const size_t i = rowIt.index();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const size_t j = colIt.index();
            os << i + 1 << " " << j + 1 << " " << m[i][j] << "\n";
        }
    }
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const vector<int>& v, const char* const name, bool append = false)
// ----------------------------------------------------------------------------
{
    // open for append is requested
    std::ofstream os(name, append ? std::ios::app : std::ios::out);
    for (size_t i = 0; i != v.size(); ++i)
        os << v[i] << "\n";
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const ResVector& v, const char* const name, bool append = false)
// ----------------------------------------------------------------------------
{
    if (!name)
        for (size_t i = 0; i != v.size(); ++i)
            std::cout << v[i] << "\n";
    if (!name)
        return;

    // open for append is requested
    std::ofstream os(name, append ? std::ios::app : std::ios::out);
    for (size_t i = 0; i != v.size(); ++i)
        os << v[i] << "\n";
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const VectorHP& v, const char* const name1, const char* const name2, bool append = false)
// ----------------------------------------------------------------------------
{
    dump_vector(v[_0], name1, append);
    dump_vector(v[_1], name2, append);
}

// ----------------------------------------------------------------------------
double
sparse_value_or_zero(const SMatrix& matrix, const size_t row, const size_t col)
// ----------------------------------------------------------------------------
{
    auto entry = matrix[row].find(col);
    return (entry == matrix[row].end()) ? 0.0 : (*entry)[0][0];
}

// ----------------------------------------------------------------------------
double
max_abs_matrix_difference(const SMatrix& lhs, const SMatrix& rhs)
// ----------------------------------------------------------------------------
{
    assert(lhs.N() == rhs.N());
    assert(lhs.M() == rhs.M());

    double max_diff = 0.0;
    for (auto row_it = lhs.begin(); row_it != lhs.end(); ++row_it) {
        const size_t row = row_it.index();
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            const size_t col = col_it.index();
            const double diff = std::abs((*col_it)[0][0] - sparse_value_or_zero(rhs, row, col));
            max_diff = std::max(max_diff, diff);
        }
    }

    for (auto row_it = rhs.begin(); row_it != rhs.end(); ++row_it) {
        const size_t row = row_it.index();
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            const size_t col = col_it.index();
            const double diff = std::abs((*col_it)[0][0] - sparse_value_or_zero(lhs, row, col));
            max_diff = std::max(max_diff, diff);
        }
    }

    return max_diff;
}

// ----------------------------------------------------------------------------
std::string
linear_system_dump_stem(const Opm::PropertyTree& prm,
                        const int nlin_iteration,
                        const std::string& reason)
// ----------------------------------------------------------------------------
{
    const std::string prefix = prm.get<std::string>("solver.linear_system_dump_prefix",
                                                    "fracture_linear_system");
    std::ostringstream oss;
    oss << prefix << "_nlin_" << nlin_iteration << "_" << reason;
    return oss.str();
}

// ----------------------------------------------------------------------------
std::string
linear_system_info_filename(const std::string& stem)
// ----------------------------------------------------------------------------
{
    const auto info_filename = std::filesystem::path(stem + "_info.json");
    return std::filesystem::absolute(info_filename).string();
}

// ----------------------------------------------------------------------------
void
write_property_tree_json(const Opm::PropertyTree& tree, const std::string& filename)
// ----------------------------------------------------------------------------
{
    std::ofstream os(filename);
    tree.write_json(os, true);
}

// ============================= Helper functions =============================

// ----------------------------------------------------------------------------
SMatrix
makeIdentity(size_t num,
             double right_padding = 0,
             double fac = 1,
             std::vector<int> zero_rows = std::vector<int>())
// ----------------------------------------------------------------------------
{
    // build a sparse matrix to represent the identity matrix of a given size
    SMatrix M;
    M.setBuildMode(SMatrix::implicit);
    M.setImplicitBuildModeParameters(1, 0.1);
    M.setSize(num, num + right_padding);
    for (size_t i = 0; i != num; ++i)
        M.entry(i, i) = fac;

    // set zero rows, if applicable
    for (size_t i = 0; i != zero_rows.size(); ++i)
        if (zero_rows[i])
            M.entry(i, i) = 0;

    M.compress();
    return M;
}
void modifyIdentity(SMatrix& M,
                  double fac = 1,
                  std::vector<int> zero_rows = std::vector<int>()){
                  int num =   M.N();
    for (int i = 0; i != num; ++i){
        M[i][i] = fac;
    }
    // set zero rows, if applicable
    for (size_t i = 0; i != zero_rows.size(); ++i){
        if (zero_rows[i]){
            M[i][i] = 0;
        }
    }
}

// ----------------------------------------------------------------------------
template <class DstMatrix, class SrcMatrix>
void
copyMatrixValuesWithSameSparsity(DstMatrix& dst, const SrcMatrix& src)
// ----------------------------------------------------------------------------
{
    assert(dst.N() == src.N());
    assert(dst.M() == src.M());

    auto srcRowIt = src.begin();
    auto dstRowIt = dst.begin();
    for (; srcRowIt != src.end() && dstRowIt != dst.end(); ++srcRowIt, ++dstRowIt) {
        assert(srcRowIt.index() == dstRowIt.index());

        auto srcColIt = srcRowIt->begin();
        auto dstColIt = dstRowIt->begin();
        for (; srcColIt != srcRowIt->end() && dstColIt != dstRowIt->end(); ++srcColIt, ++dstColIt) {
            assert(srcColIt.index() == dstColIt.index());
            *dstColIt = *srcColIt;
        }

        assert(srcColIt == srcRowIt->end());
        assert(dstColIt == dstRowIt->end());
    }

    assert(srcRowIt == src.end());
    assert(dstRowIt == dst.end());
}

// ----------------------------------------------------------------------------
std::unique_ptr<LinearSolverBase>
setupLinearSolver(const std::string& solver_type,
                  const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>& linop,
                  LinearPrecondType& precond,
                  const Opm::PropertyTree& prm,
                  const double lintol,
                  const int max_iter,
                  const int verbosity)
// ----------------------------------------------------------------------------
{
    if (solver_type == "bicgstab") {
        return std::make_unique<Dune::BiCGSTABSolver<VectorHP>>(linop,
                                                                 precond,
                                                                 lintol,
                                                                 max_iter,
                                                                 verbosity);
    }

    if (solver_type == "gmres") {
        const int restart = prm.get<int>("solver.linsolver.restart", 100);
        return std::make_unique<Dune::RestartedGMResSolver<VectorHP>>(linop,
                                                                       precond,
                                                                       lintol,
                                                                       restart,
                                                                       max_iter,
                                                                       verbosity);
    }

    if (solver_type == "fgmres") {
        const int restart = prm.get<int>("solver.linsolver.restart", 100);
        return std::make_unique<Dune::RestartedFlexibleGMResSolver<VectorHP>>(linop,
                                                                               precond,
                                                                               lintol,
                                                                               restart,
                                                                               max_iter,
                                                                               verbosity);
    }

    OPM_THROW(std::runtime_error, "Unknown linear solver type: for fracture");
}

// ----------------------------------------------------------------------------
std::unique_ptr<Opm::Fracture::Matrix>
initCouplingMatrixSparsity(const std::vector<Htrans>& htrans, size_t num_cells, size_t num_wells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    auto C = std::make_unique<Opm::Fracture::Matrix>(
        num_cells + num_wells, num_cells, 4, 0.4, Opm::Fracture::Matrix::implicit);
    for (const auto& e : htrans) {
        const size_t i = std::get<0>(e);
        const size_t j = std::get<1>(e);
        assert(i != j);
        C->entry(i, j) = 0;
        C->entry(j, i) = 0;
        C->entry(i, i) = 0;
        C->entry(j, j) = 0;
    }
    C->compress();
    return C;
}

// ----------------------------------------------------------------------------
void
updateCouplingMatrix(std::unique_ptr<Opm::Fracture::Matrix>& Cptr,
                     const size_t num_cells,
                     const size_t num_wells,
                     const std::vector<Htrans>& htrans,
                     const std::vector<double>& face_dgh, // empty: no gravity
                     const ResVector& pressure,
                     const ResVector& aperture,
                     const std::vector<int>& closed_cells,
                     const std::vector<double>& reservoir_mobility,
                     double min_width)
{
    OPM_TIMEFUNCTION();
    // create C if not done already
    if (!Cptr)
        Cptr = initCouplingMatrixSparsity(htrans, num_cells, num_wells);

    // NB: If the implementation of the pressure matrix changes (i.e. `assemblePressure()`),
    //        the code below might need to be updated accordingly as well.
    auto& C = *Cptr;
    C = 0;
    for (size_t k = 0; k < htrans.size(); ++k) {
        const auto& e = htrans[k];
        const size_t i = std::get<0>(e);
        const size_t j = std::get<1>(e);
        assert(i != j);

        const double t1 = std::get<2>(e);
        const double t2 = std::get<3>(e);
        double  dh1 = 1.0;
        double  dh2 = 1.0;
        if(aperture[i][0] <= min_width) dh1 = 0.0;
        if(aperture[j][0] <= min_width) dh2 = 0.0;
        assert(aperture[i][0]>=0.0);
        assert(aperture[j][0]>=0.0);
        const double h1 = std::max(aperture[i][0],min_width);
        const double h2 = std::max(aperture[j][0],min_width);
        const double g_ij = face_dgh.empty() ? 0.0 : face_dgh[k];
        const double p1 = pressure[i] - g_ij;
        const double p2 = pressure[j];

        const double q = (h1 * h1 * h1) * (h2 * h2 * h2) * (t1 * t2); // numerator
        const double d1q = 3 * (h1 * h1) * (h2 * h2 * h2) * (t1 * t2)* dh1;
        const double d2q = 3 * (h1 * h1 * h1) * (h2 * h2) * (t1 * t2)* dh2;

        const double r = 12 * (h1 * h1 * h1 * t1 + h2 * h2* h2 * t2); // denominator
        const double d1r = 36 * (h1 * h1) * t1 *dh1;
        const double d2r = 36 * (h2 * h2) * t2 *dh2;
        assert(r >= 0.0);
        assert(q >= 0.0);
        const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
        const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);

        const double mobility = 0.5 * (reservoir_mobility[i] + reservoir_mobility[j]);
        // diagonal elements
        C[i][i] += dTdh1 * (p1 - p2) * mobility;
        C[j][j] += dTdh2 * (p2 - p1) * mobility;

        // off-diagonal elements
        C[i][j] += dTdh2 * (p1 - p2) * mobility;
        C[j][i] += dTdh1 * (p2 - p1) * mobility;
    }
    // zeroing out columns corresponding to closed cells
    if(false){
      for (size_t col = 0; col != closed_cells.size(); ++col)
        if (closed_cells[col])
          for (size_t row = 0; row != C.N(); ++row)
            if (C.exists(row, col))
              C[row][col] = 0;
    }
}

// ----------------------------------------------------------------------------
std::unique_ptr<SMatrix>
assemble_coupling_original(const Opm::FracturePressureInput& input,
                           const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    ResVector pressure(input.num_cells + input.num_well_equations);
    ResVector aperture(input.num_cells + input.num_well_equations);
    std::vector<double> mobility(input.num_cells, 0.0);

    for (size_t i = 0; i < pressure.size(); ++i) {
        pressure[i][0] = input.fracture_pressure[i];
        aperture[i][0] = input.fracture_width[i];
    }

    for (size_t i = 0; i < input.num_cells; ++i) {
        mobility[i] = input.density[i].value / input.viscosity[i].value;
    }

    auto coupling = initCouplingMatrixSparsity(input.htrans,
                                               input.num_cells,
                                               input.num_well_equations);
    updateCouplingMatrix(coupling,
                         input.num_cells,
                         input.num_well_equations,
                         input.htrans,
                         input.face_gravity,
                         pressure,
                         aperture,
                         closed_cells,
                         mobility,
                         input.min_width);
    return coupling;
}

// ----------------------------------------------------------------------------
void
dump_linear_system_snapshot(const std::string& stem,
                            const SystemMatrix& system,
                            const VectorHP& rhs,
                            const VectorHP& x,
                            const std::vector<int>& closed_cells,
                            const bool use_ad,
                            const Opm::PropertyTree& prm,
                            const Opm::FracturePressureInput& input,
                            const Dune::InverseOperatorResult* iores = nullptr)
// ----------------------------------------------------------------------------
{
    dump_matrix(system[_0][_0], (stem + "_A.txt").c_str());
    Dune::storeMatrixMarket(system[_0][_1], stem + "_I.mtx");
    Dune::storeMatrixMarket(system[_1][_0], stem + "_C_active.mtx");
    Dune::storeMatrixMarket(system[_1][_1], stem + "_M_active.mtx");
    Dune::storeMatrixMarket(rhs[_0], stem + "_rhs_w.mtx");
    Dune::storeMatrixMarket(rhs[_1], stem + "_rhs_p.mtx");
    Dune::storeMatrixMarket(x[_0], stem + "_x_w.mtx");
    Dune::storeMatrixMarket(x[_1], stem + "_x_p.mtx");
    dump_vector(closed_cells, (stem + "_closed_cells.txt").c_str());

    const auto ad_result = assemblePressureAD(input);
    const auto original_pressure = assemblePressureOriginal(input);
    const auto original_coupling = assemble_coupling_original(input, closed_cells);

    Dune::storeMatrixMarket(*ad_result.pressure_matrix, stem + "_M_ad.mtx");
    Dune::storeMatrixMarket(*ad_result.coupling_matrix, stem + "_C_ad.mtx");
    Dune::storeMatrixMarket(*original_pressure, stem + "_M_original.mtx");
    Dune::storeMatrixMarket(*original_coupling, stem + "_C_original.mtx");

    Opm::PropertyTree info;
    info.put("use_ad_pressure_assembly", use_ad);
    info.put("active_pressure_matrix", std::string(use_ad ? "ad" : "original"));
    info.put("mechanics_size", system[_0][_0].N());
    info.put("pressure_size", system[_1][_1].N());
    info.put("pressure_max_abs_diff_ad_vs_original",
             max_abs_matrix_difference(*ad_result.pressure_matrix, *original_pressure));
    info.put("coupling_max_abs_diff_ad_vs_original",
             max_abs_matrix_difference(*ad_result.coupling_matrix, *original_coupling));
    if (iores) {
        info.put("linear.converged", iores->converged);
        info.put("linear.iterations", iores->iterations);
    }

    write_property_tree_json(info, stem + "_info.json");
    write_property_tree_json(prm.get_child("solver.linsolver"), stem + "_linsolver.json");
    write_property_tree_json(prm.get_child("solver.linsolver.preconditioner.flow_solver"),
                             stem + "_flow_solver.json");
}

// ----------------------------------------------------------------------------
inline bool
convergence_test(const VectorHP& res, const double tol_flow, double tol_mech, int verbosity)
// ----------------------------------------------------------------------------
{
    bool converged = res[_0].infinity_norm() < tol_mech && res[_1].infinity_norm() < tol_flow;
    if(converged && verbosity > 0){
      Opm::OpmLog::debug(fmt::format("Converged with residual norm[0] is {}\nConverged with residual norm[1] is {}",
                                res[_0].infinity_norm(), res[_1].infinity_norm()));
    }
    return converged;
}

// ----------------------------------------------------------------------------
template <typename Mat>
ResVector
diagvec(const Mat& M)
// ----------------------------------------------------------------------------
{
    ResVector res(M.M());
    for (size_t i = 0; i != res.size(); ++i)
        res[i] = Opm::diagScalar(M[i][i]);
    return res;
}


// ----------------------------------------------------------------------------
class TailoredPrecondDiag : public Dune::Preconditioner<VectorHP, VectorHP>
// ----------------------------------------------------------------------------
{
public:
    TailoredPrecondDiag(const SystemMatrix& S)
        : A_diag_(diagvec(S[_0][_0]))
        , M_diag_(diagvec(S[_1][_1]))
    {
    }
    virtual void apply(VectorHP& v, const VectorHP& d)
    {
        for (size_t i = 0; i != A_diag_.size(); ++i)
            v[_0][i] = d[_0][i] / A_diag_[i];
        for (size_t i = 0; i != M_diag_.size(); ++i)
            v[_1][i] = d[_1][i] / M_diag_[i];
    };
    virtual void post(VectorHP& /*v*/) {};
    virtual void pre(VectorHP& /*x*/, VectorHP& /*b*/) {};
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const ResVector A_diag_;
    const ResVector M_diag_;
};


// ----------------------------------------------------------------------------
[[maybe_unused]] double
estimate_step_fac(const VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------
{
    // estimate what might be a safe step size to avoid exiting the convergence radius
    const double x0_norm = x[_0].infinity_norm();
    const double x1_norm = x[_1].infinity_norm();

    const double f1 = x0_norm == 0 ? 1 : dx[_0].infinity_norm() / x0_norm;
    const double f2 = x1_norm == 0 ? 1 : dx[_1].infinity_norm() / x1_norm;
    const double fmax = std::max(f1, f2);
    const double threshold = 0.95;
    const double fac_min = 1e-1;
    return (fmax < threshold) ? 1.0 : std::max(threshold / fmax, fac_min);
}

// ----------------------------------------------------------------------------
FMatrix
modified_fracture_matrix(const FMatrix& A, const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    Dune::DynamicMatrix result(A);

    for (size_t row = 0; row != A.N(); ++row)
        if (closed_cells[row])
            for (size_t col = 0; col != A.N(); ++col)
                result[row][col] = (row == col);
    return result;
}

void modified_fracture_matrix(FMatrix& Aold, const FMatrix& A, const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    assert(Aold.N() == A.N());
    assert(Aold.M() == A.M());
    for(size_t i=0; i < A.N(); ++i){
      for(size_t j=0; j < A.M(); ++j){
        Aold[i][j] = A[i][j];
      }
    }
    for (size_t row = 0; row != A.N(); ++row)
        if (closed_cells[row])
            for (size_t col = 0; col != A.N(); ++col)
                Aold[row][col] = (row == col);
}

}; // end anonymous namespace

namespace Opm
{
// ----------------------------------------------------------------------------
std::vector<int>
Fracture::identify_closed(const FMatrix& A, const VectorHP& x, const ResVector& rhs)
// ----------------------------------------------------------------------------
{

    OPM_TIMEFUNCTION();
    const InnerSolveOptions& opt = inner_;
    const std::string& closing_type = opt.closingType;
    const std::string& closed_cell_policy = opt.closedCellPolicy;
    const double close_force_tolerance = opt.closeForceTol;
    const double reopen_force_tolerance = opt.reopenForceTol;
    const double reopen_width_tolerance = opt.reopenWidthTol;
    // "none": never close any cell — keep the fracture open and allow negative
    // aperture (no contact constraint). This makes the mechanics block linear, so
    // the keep-C coupled Newton stays smooth and the open/close chatter cannot
    // occur. Physical only where genuine contact is negligible; the flow already
    // floors transmissibility at min_width, so w<0 is benign for the flow.
    // Fischer-Burmeister handles contact through smooth mech-row weights in
    // fullSystemIteration, not a binary active set. Both report no hard-closed
    // cells here (the binary set is only used for tracking/metrics and the C-drop check).
    if (closed_cell_policy == "fischer_burmeister" || closed_cell_policy == "none") {
        return std::vector<int>(A.N(), 0);
    }
    const bool legacy_closed_cells = (closed_cell_policy == "legacy");
    // Primal-dual active set (semismooth Newton): the contact set is chosen by a
    // single consistently-scaled criterion  active(closed) <=> tmp_i - c_i w_i > tol
    // with c_i = |A_ii| (BEM self-stiffness, converts width to force units), instead
    // of the separate force/width tests + sticky hysteresis. Folding the open/close
    // decision into one smooth complementarity function is what stops the active-set
    // chatter that prevents the keep-C true-Newton iteration from converging.
    const bool pdas_closed_cells = (closed_cell_policy == "pdas");
    const auto was_closed = [this](const size_t index) {
        return index < closed_cells_.size() && closed_cells_[index] != 0;
    };

    // Opt-in (default off): when testing whether a CLOSED cell should reopen, use
    // the pressure that would develop once it opens, not its pinned (closed)
    // pressure. The static reopen test ignores that fluid flows into a cell as it
    // opens and raises its pressure — so it underestimates the opening tendency and
    // can lock cells closed (the "half-open fracture hard to open" pathology). We
    // estimate the pressure flow could communicate as the max fracture pressure
    // over currently-open cells.
    const bool reopen_use_open_pressure = opt.reopenUseOpenPressure;

    // computing rhs - A x[0] - I x[1]
    std::vector<int> result;
    if(closing_type == "org"){
        const ResVector& h = x[_0];
        const ResVector& p = x[_1];
        double p_open_estimate = 0.0;
        if (reopen_use_open_pressure)
            for (size_t i = 0; i != A.N(); ++i)
                if (!was_closed(i))
                    p_open_estimate = std::max(p_open_estimate, p[i][0]);
        const ResVector tmp = contactForce(A, h, p, rhs);
        for (size_t i = 0; i != A.N(); ++i) {
            if (legacy_closed_cells) {
                result.push_back(tmp[i] >= 0.0 && h[i] <= 0.0);
                continue;
            }
            if (pdas_closed_cells) {
                // Semismooth-Newton / primal-dual active set: one criterion,
                // stateless (no sticky hysteresis), consistently scaled by the
                // mech self-stiffness c_i = |A_ii|.
                const double c_i = std::abs(A[i][i]);
                result.push_back((tmp[i] - c_i * h[i][0]) > close_force_tolerance);
                continue;
            }
            const bool close_candidate = tmp[i] >= close_force_tolerance && h[i] <= 0.0;
            // Re-opening is force-driven. A closed cell has its width pinned at 0,
            // so a width gate (h[i] >= reopen_width_tolerance) with a positive
            // tolerance can never be satisfied and would lock the cell closed
            // forever, suppressing fracture growth. reopen_width_tolerance is kept
            // for backward compatibility but no longer gates on the pinned width.
            (void)reopen_width_tolerance;
            // tmp[i] = (rhs - A h - I p)[i]; replacing the cell pressure p[i] by the
            // higher pressure-on-opening makes the reopen test easier (more tensile),
            // since I is the identity on cell rows.
            double reopen_force = tmp[i];
            if (reopen_use_open_pressure && was_closed(i))
                reopen_force = tmp[i] - (p_open_estimate - p[i][0]);
            const bool reopen_candidate = reopen_force <= -reopen_force_tolerance;
            result.push_back(was_closed(i) ? !reopen_candidate : close_candidate);
        }
    } else if(closing_type == "open") {
            result.resize(A.N(), 1);
    } else if(closing_type == "simple"){
     // sum forces over faces
        result.resize(A.N(), 0);
        for (auto element : Dune::elements(grid_->leafGridView())) {
            int i = grid_->leafGridView().indexSet().index(element);
            double pressure = fracture_pressure_[i][0];
            //NB maybe use other pressure if closed
            if (legacy_closed_cells) {
                result[i] = (this->fractureForce(i) - pressure) < 0;
                continue;
            }
            const double closure_balance = this->fractureForce(i) - pressure;
            const bool close_candidate = closure_balance < -close_force_tolerance;
            // Force-driven reopen (see "org" branch): width is pinned at 0 for
            // closed cells, so we do not gate reopening on it.
            (void)reopen_width_tolerance;
            const bool reopen_candidate = closure_balance > reopen_force_tolerance;
            result[i] = was_closed(i) ? !reopen_candidate : close_candidate;
        }
    } else {
        std::stringstream ss;
        ss << "Unknown closing type: " << closing_type;
        OPM_THROW(std::runtime_error, ss.str());
    }

    return result;
}

// ----------------------------------------------------------------------------
void
Fracture::assemblePressureAndCouplingAD(const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    updateLeakoff();
    const auto input = makePressureAssemblyInput();

    // Run the AD-based assembly (produces both pressure and coupling matrices)
    auto result = assemblePressureAD(input);

    // Store the pressure matrix
    pressure_matrix_ = std::move(result.pressure_matrix);

    // Store the coupling matrix
    coupling_matrix_ = std::move(result.coupling_matrix);

}

// ----------------------------------------------------------------------------
FracturePressureInput
Fracture::makePressureAssemblyInput() const
// ----------------------------------------------------------------------------
{
    const size_t nc = numFractureCells();
    const size_t nw = numWellEquations();
    const bool have_water_phase_props = static_cast<bool>(fracture_water_property_evaluator_);
    const bool use_fluid_system_water_properties =
        have_water_phase_props
        && prm_.get<bool>("solver.use_fluid_system_water_properties", false);

    FracturePressureInput input;
    input.num_cells = nc;
    input.min_width = min_width_;
    input.htrans = htrans_;
    input.control_type = effectiveControlType();
    if (input.control_type == "bhp_well" && !perfinj_.empty()) {
        const double density = reservoir_.density[std::get<0>(perfinj_[0])];
        input.well_target_pressure = well_target_bhp_
            + density * gravity_ * (perf_ref_depth_ - well_ref_depth_);
    }
    input.num_well_equations = nw;
    input.perfinj.assign(perfinj_.begin(), perfinj_.end());
    input.total_wellindex = total_wellindex_;
    input.mobility_water_perf = mobility_water_perf_;
    input.face_mobility.resize(nc + nw, 0.0);
    for (size_t i = 0; i < reservoir_.mobility.size() && i < nc; ++i)
        input.face_mobility[i] = reservoir_.mobility[i];
    input.leakof = leakof_;
    input.use_fluid_mobility_for_internal_faces =
        use_fluid_system_water_properties
        && prm_.get<bool>("solver.use_water_phase_mobility_for_fracture_flow", false);
    input.use_fluid_mobility_for_well_connections =
        use_fluid_system_water_properties
        && prm_.get<bool>("solver.use_water_phase_mobility_for_well_fracture_flow",
                          input.use_fluid_mobility_for_internal_faces);

    input.fracture_width.resize(nc + nw);
    input.fracture_pressure.resize(nc + nw);
    input.density.resize(nc + nw);
    input.viscosity.resize(nc + nw);

    for (size_t i = 0; i < nc + nw; ++i) {
        // fracture_width_ has no entry for well equations (rate_well)
        input.fracture_width[i] = (i < fracture_width_.size()) ? fracture_width_[i][0] : 0.0;
        input.fracture_pressure[i] = fracture_pressure_[i][0];
    }
    input.face_gravity = faceGravityHeads();

    for (size_t i = 0; i < nc; ++i) {
        if (use_fluid_system_water_properties) {
            const auto fluid_properties = fracture_water_property_evaluator_(i, fracture_pressure_[i][0]);
            input.density[i] = fluid_properties.first;
            input.viscosity[i] = fluid_properties.second;
        }
        else {
            input.density[i] = {reservoir_.density[i], 0.0};
            input.viscosity[i] = {1.0, 0.0};
        }
    }
    for (size_t i = nc; i < nc + nw; ++i) {
        input.density[i] = {1.0, 0.0};
        input.viscosity[i] = {1.0, 0.0};
    }

    return input;
}

// ----------------------------------------------------------------------------
Fracture::InnerSolveOptions
Fracture::InnerSolveOptions::read(const PropertyTree& prm)
// ----------------------------------------------------------------------------
{
    InnerSolveOptions o;
    o.useAd = prm.get<bool>("solver.use_ad_pressure_assembly", false);
    o.verbosity = prm.get<int>("solver.verbosity", 0);
    o.closingType = prm.get<std::string>("solver.closing_type", "org");
    o.closedCellPolicy = prm.get<std::string>("solver.closed_cell_policy", "sticky");
    const std::string& policy = o.closedCellPolicy;
    if (policy != "sticky" && policy != "legacy" && policy != "pdas"
        && policy != "fischer_burmeister" && policy != "none") {
        OPM_THROW(std::runtime_error, "Unknown closed cell policy: " + policy);
    }
    o.isFb = (policy == "fischer_burmeister");
    o.closeForceTol = prm.get<double>("solver.close_force_tolerance", 0.0);
    o.reopenForceTol = prm.get<double>("solver.reopen_force_tolerance", o.closeForceTol);
    o.reopenWidthTol = prm.get<double>("solver.reopen_width_tolerance", 0.0);
    o.reopenUseOpenPressure = prm.get<bool>("solver.reopen_use_open_pressure", false);
    o.freezeClosedAfter = prm.get<int>("solver.freeze_closed_after", -1);
    o.chatterGuard = (prm.get<std::string>("solver.toggle_guard_mode", "count") == "chatter");
    if (o.chatterGuard) {
        o.chatterFlipLimit = prm.get<int>("solver.chatter_flip_limit", 3);
    } else {
        o.maxToggleCount = prm.get<int>("solver.max_closed_cell_toggle_count", 1000000000);
        o.maxToggleFraction = prm.get<double>("solver.max_closed_cell_toggle_fraction", 1.0);
    }
    if (o.isFb) {
        o.fbMu = prm.get<double>("solver.fb_mu", 0.0);
        o.fbScaledConvergence = prm.get<bool>("solver.fb_scaled_convergence", false);
        o.fbLineSearch = prm.get<bool>("solver.fb_line_search", false);
        if (o.fbLineSearch) {
            o.fbLineSearchMin = prm.get<double>("solver.fb_line_search_min", 0.0625);
        }
        o.fbClampWidth = prm.get<bool>("solver.fb_clamp_width", true);
    }
    o.dropCoupling = prm.get<bool>("solver.drop_fluid_mech_linearization", true);
    if (!o.dropCoupling) {
        o.dropTolH = prm.get<double>("solver.drop_tol_h", 1e1);
        o.dropTolP = prm.get<double>("solver.drop_tol_p", 1.0);
    }
    o.storageTerm = prm.get<bool>("solver.fracture_storage_term", false);
    o.scaling = prm.get<std::string>("solver.linsolver.scaling", "none");
    o.linearSolver = prm.get<std::string>("solver.linsolver.solver");
    if (o.linearSolver == "direct") {
        o.directMaxCells = prm.get<int>("solver.linsolver.direct_max_cells", 4000);
        o.directFallbackSolver = prm.get<std::string>("solver.linsolver.direct_fallback_solver", "fgmres");
    }
    o.precond = prm.get_child("solver.linsolver.preconditioner");
    o.updateMechOnReuse = o.precond.get<bool>("update_mech_on_reuse", false);
    o.linTol = prm.get<double>("solver.linsolver.tol");
    o.linAtol = prm.get<double>("solver.linsolver.atol");
    o.linMaxIter = static_cast<int>(prm.get<double>("solver.linsolver.max_iter"));
    o.linVerbosity = static_cast<int>(prm.get<double>("solver.linsolver.verbosity", 0));
    o.dumpMinNonlinIter = prm.get<int>("solver.dump_case_min_nonlin_iter", 0);
    o.writeCoupledSystem = prm.get<bool>("solver.write_coupled_linear_system", false);
    o.dumpOnFailure = prm.get<bool>("solver.dump_linear_system_on_failure", true);
    o.debugChecksum = prm.get<bool>("solver.linsolver.debug_checksum", false);
    o.failurePolicy = prm.get<std::string>("solver.linsolver.failure_policy", "throw");
    o.gmresFallback = prm.get<bool>("solver.linsolver.fallback_gmres_on_breakdown", true);
    o.damping = prm.get<double>("solver.damping");
    o.coupledStepLimit = prm.get<double>("solver.coupled_step_limit", 0.0);
    o.maxDwidth = prm.get<double>("solver.max_dwidth");
    o.maxDp = prm.get<double>("solver.max_dp");
    return o;
}

// ----------------------------------------------------------------------------
Fracture::ResVector
Fracture::contactForce(const FMatrix& A, const ResVector& w, const ResVector& p, const ResVector& traction)
// ----------------------------------------------------------------------------
{
    // traction - A w - p on the fracture cells; p may carry trailing well rows
    ResVector force(traction);
    A.mmv(w, force);
    for (size_t i = 0; i != A.N(); ++i) {
        force[i] -= p[i];
    }
    return force;
}

// ----------------------------------------------------------------------------
void
Fracture::updateFlowCoefficients_()
// ----------------------------------------------------------------------------
{
    updateLeakoff();
    if (!inner_.useAd) {
        assemblePressure();
    }
    addSource();
}

// ----------------------------------------------------------------------------
Fracture::ContactUpdate
Fracture::updateContactSet_(const VectorHP& x, const ResVector& traction, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    const int nlin_verbosity = opt.verbosity;
    // freeze_closed_after >= 0 reuses the previous contact set after that iteration,
    // separating active-set chatter from a bad Newton step for a fixed set.
    std::vector<int> closed_cells;
    if (opt.freezeClosedAfter >= 0 && nlin_iteration > opt.freezeClosedAfter
        && closed_cells_.size() == fracture_width_.size()) {
        closed_cells = closed_cells_; // frozen set
    } else {
        closed_cells = identify_closed(fractureMatrix(), x, traction);
    }
    if (opt.chatterGuard) {
        // Per-cell anti-chatter guard: chatter is the same cell flipping repeatedly
        // within one solve, while en-masse first-time re-opening is legitimate. A
        // global re-opening throttle suppresses an opening fracture, so pin only
        // cells that exceeded chatter_flip_limit changes this solve.
        if (nlin_iteration == 0 || cell_flip_counts_.size() != closed_cells.size()) {
            cell_flip_counts_.assign(closed_cells.size(), 0);
        }
        if (closed_cells_.size() == closed_cells.size()) {
            const int flip_limit = opt.chatterFlipLimit;
            size_t pinned = 0;
            for (size_t i = 0; i < closed_cells.size(); ++i) {
                if (closed_cells[i] != closed_cells_[i]) {
                    if (cell_flip_counts_[i] >= flip_limit) {
                        closed_cells[i] = closed_cells_[i]; // pin chattering cell
                        ++pinned;
                    } else {
                        ++cell_flip_counts_[i];
                    }
                }
            }
            if (pinned > 0 && nlin_verbosity > 0) {
                OpmLog::debug(fmt::format("Chatter guard pinned {} cells that changed contact state more than {} times this solve",
                                          pinned, flip_limit));
            }
        }
    } else if (!closed_cells_.empty() && closed_cells_.size() == closed_cells.size()) {
        // Only re-openings are throttled; establishing contact is always allowed.
        // A finite limit also rejects legitimate en-masse re-opening.
        const size_t toggled_cells = count_reopened_cells(closed_cells_, closed_cells);
        const size_t max_toggle_count = static_cast<size_t>(std::max(0, opt.maxToggleCount));
        const size_t fraction_limit = static_cast<size_t>(opt.maxToggleFraction * closed_cells.size());
        const size_t allowed_toggles = std::min(max_toggle_count, std::max<size_t>(fraction_limit, 0));
        if (toggled_cells > allowed_toggles) {
            if (nlin_verbosity > 0) {
                OpmLog::debug(fmt::format("Keeping previous closed-cell state after {} re-openings exceeded guard limit {}",
                                          toggled_cells, allowed_toggles));
            }
            closed_cells = closed_cells_;
        }
    }
    bool closed_cells_changed = false;
    if(closed_cells_.size() != closed_cells.size()){
      closed_cells_changed = true;
    }else{
      for(size_t i=0; i!=closed_cells.size(); ++i){
        if(closed_cells_[i] != closed_cells[i]){
          closed_cells_changed = true;
          break;
        }
      }
    }
    // Contact-chatter metric: cells that flipped since the previous iteration
    // (same grid only; a size mismatch is a grid change).
    if (closed_cells_.size() == closed_cells.size()) {
        last_solve_stats_.closed_cell_toggles +=
            static_cast<int>(count_toggled_cells(closed_cells_, closed_cells));
        // When each cell last toggled, for the propagation veto: sticky chatters by
        // nature, so what matters is whether the cell settled, not how often.
        if (cell_last_toggle_iter_.size() != closed_cells.size())
            cell_last_toggle_iter_.assign(closed_cells.size(), -1000);
        for (size_t i = 0; i < closed_cells.size(); ++i)
            if (closed_cells_[i] != closed_cells[i])
                cell_last_toggle_iter_[i] = nlin_iteration;
    } else {
        cell_last_toggle_iter_.assign(closed_cells.size(), -1000);
    }
    current_solve_iter_ = nlin_iteration;
    closed_cells_ = closed_cells;
    return {std::move(closed_cells), closed_cells_changed};
}

// ----------------------------------------------------------------------------
Fracture::FbRows
Fracture::assembleMechanicsRows_(const VectorHP& x,
                                 ResVector& mech_rhs,
                                 const std::vector<int>& closed_cells,
                                 const bool first_iteration)
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    FbRows fb;
    // Fischer-Burmeister contact replaces the binary open/closed mech rows by the
    // smooth complementarity row phi_FB(c_i w_i, lambda_i) = 0 with lambda the
    // contact force and c_i = |A_ii|. Jacobian row: A_FB = da c_i e_i - db A_i,
    // I_FB_ii = -db, with da = 1 - a/r, db = 1 - b/r, r = sqrt(a^2+b^2).
    if (opt.isFb) {
        const FMatrix& Aorig = fractureMatrix();
        const size_t n = Aorig.N();
        const ResVector& w = x[_0];
        const ResVector& p = x[_1];
        fb.phi.assign(n, 0.0);
        fb.scale.assign(n, 1.0);
        fb_cell_residual_.assign(n, 0.0);
        if (first_iteration) A_.reset(new FMatrix(Aorig));
        FMatrix& Afb = *A_;
        if (first_iteration) I_.reset(new SMatrix(makeIdentity(n, numWellEquations())));
        SMatrix& Ifb = *I_;
        const double eps = 1e-30;
        const ResVector force = contactForce(Aorig, w, p, mech_rhs);
        // fb_mu > 0 is Chen-Mangasarian smoothing of the a=b=0 kink; it relaxes
        // complementarity by ~mu, so it must be small against (c_i w)^2.
        const double mu = opt.fbMu;
        for (size_t i = 0; i < n; ++i) {
            const double c_i = std::max(std::abs(Aorig[i][i]), eps);
            const double a = c_i * w[i][0];
            const double b = force[i][0];
            const double r = std::sqrt(a * a + b * b + 2.0 * mu) + eps;
            const double da = 1.0 - a / r;
            const double db = 1.0 - b / r;
            fb.phi[i] = a + b - std::sqrt(a * a + b * b + 2.0 * mu);
            fb.scale[i] = 1.0 + std::abs(a) + std::abs(b);
            fb_cell_residual_[i] = std::abs(fb.phi[i]) / fb.scale[i];
            for (size_t j = 0; j < n; ++j) Afb[i][j] = -db * Aorig[i][j];
            Afb[i][i] += da * c_i;
            Ifb[i][i] = -db;
        }
    } else {
        if (first_iteration) {
            A_.reset(new FMatrix(modified_fracture_matrix(fractureMatrix(), closed_cells)));
        }else{
            modified_fracture_matrix(*A_, fractureMatrix(), closed_cells);
        }
        // Closed cells are held fixed in the mechanics residual and Jacobian block.
        for (size_t i = 0; i != closed_cells.size(); ++i)
            if (closed_cells[i])
                mech_rhs[i] = 0;
    }
    return fb;
}

// ----------------------------------------------------------------------------
void
Fracture::assembleFlowRows_(const ResVector& mech_rhs, const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    if (opt.useAd) {
        assemblePressureAndCouplingAD(closed_cells);
    } else {
        updateCouplingMatrix(coupling_matrix_,
                             pressure_matrix_->N() - numWellEquations(), // num cells
                             numWellEquations(), // num wells
                             htrans_,
                             faceGravityHeads(),
                             fracture_pressure_,
                             fracture_width_,
                             closed_cells,
                             reservoir_.mobility,
                             min_width_);
    }
    if (opt.dropCoupling) {
        *coupling_matrix_ = 0;
    } else if (mech_rhs.infinity_norm() > opt.dropTolH
               || rhs_pressure_.infinity_norm() > opt.dropTolP) {
        *coupling_matrix_ = 0;
    }
    // External-pressure mode: pressure rows become p = p_external (identity,
    // zero residual, no width coupling), so the Newton reduces to mechanics
    // and contact at the pressure the flow's aux cells hold.
    // Cells without an external value (no flow DOF yet) keep their own pressure
    // equation, so a front that grew inside the step can still be pressurised.
    if (external_pressure_) {
        const std::size_t nc = numFractureCells();
        const bool haveMask = (external_cell_mask_.size() == nc);
        auto& Mx = *pressure_matrix_;
        auto& Cx = *coupling_matrix_;
        for (auto row = Mx.begin(); row != Mx.end(); ++row) {
            const std::size_t r = row.index();
            const bool pinned = (r >= nc) || !haveMask || external_cell_mask_[r];
            if (!pinned) {
                continue;
            }
            for (auto col = row->begin(); col != row->end(); ++col) {
                *col = (col.index() == r) ? 1.0 : 0.0;
            }
            if (r < Cx.N()) {
                for (auto& entry : Cx[r]) {
                    entry = 0.0;
                }
            }
            rhs_pressure_[r] = fracture_pressure_[r];
        }
    }
}

// ----------------------------------------------------------------------------
void
Fracture::assembleSystemMatrix_(const std::vector<int>& closed_cells, const bool first_iteration)
// ----------------------------------------------------------------------------
{
    const auto& A = *A_;
    const auto& M = *pressure_matrix_;
    const auto& C = *coupling_matrix_;
    // Fischer-Burmeister already built I_ with its smooth -db_i weights.
    if (!inner_.isFb) {
        if (first_iteration) {
            I_.reset(new SMatrix(makeIdentity(A.N(), numWellEquations(), 1, closed_cells)));
        }else{
            modifyIdentity(*I_, 1, closed_cells);
        }
    }
    const auto& I = *I_;
    if(first_iteration){
        S_.reset(new SystemMatrix({{A, I},
                                   {C, M}}));
    } else{
        OPM_TIMEBLOCK(Modify_System_Matrix);
        // modify S in place to avoid copying the whole matrix
        auto& S = *S_;
        for(size_t i=0; i!=A.N(); ++i){
            for(size_t j=0; j!=A.M(); ++j){
                S[_0][_0][i][j] = A[i][j];
            }
        }
        copyMatrixValuesWithSameSparsity(S[_0][_1], I);
        copyMatrixValuesWithSameSparsity(S[_1][_0], C);
        copyMatrixValuesWithSameSparsity(S[_1][_1], M);
    }
}

// ----------------------------------------------------------------------------
Fracture::VectorHP
Fracture::assembleResidual_(const VectorHP& x, const ResVector& mech_rhs, const FbRows& fb) const
// ----------------------------------------------------------------------------
{
    const auto& S = *S_;
    const auto& M = *pressure_matrix_;
    VectorHP rhs {x};
    rhs[_0] = mech_rhs;
    rhs[_1] = rhs_pressure_;

    {
        // A block whose column count disagrees with the vector it multiplies
        // reads past the end of that vector: silent garbage in the residual and
        // in every Krylov matrix-vector product. Fail by name instead.
        const size_t nm = rhs[_0].size(), np = rhs[_1].size();
        if (x[_0].size() != nm || x[_1].size() != np)
            OPM_THROW(std::runtime_error,
                      "Fracture coupled solve: iterate/residual size mismatch: dx="
                          + std::to_string(x[_0].size()) + "/" + std::to_string(x[_1].size())
                          + " rhs=" + std::to_string(nm) + "/" + std::to_string(np));
        auto bad = [](const auto& B, size_t n, size_t m) { return B.N() != n || B.M() != m; };
        if (bad(S[_0][_0], nm, nm) || bad(S[_0][_1], nm, np) || bad(S[_1][_0], np, nm)
            || bad(S[_1][_1], np, np))
            OPM_THROW(std::runtime_error,
                      "Fracture coupled system block dimensions inconsistent: nm=" + std::to_string(nm)
                          + " np=" + std::to_string(np) + " A=" + std::to_string(S[_0][_0].N()) + "x"
                          + std::to_string(S[_0][_0].M()) + " I=" + std::to_string(S[_0][_1].N()) + "x"
                          + std::to_string(S[_0][_1].M()) + " C=" + std::to_string(S[_1][_0].N()) + "x"
                          + std::to_string(S[_1][_0].M()) + " M=" + std::to_string(S[_1][_1].N()) + "x"
                          + std::to_string(S[_1][_1].M()));
    }

    S[_1][_1].mmv(x[_1], rhs[_1]);
    S[_0][_0].mmv(x[_0], rhs[_0]);
    S[_0][_1].mmv(x[_1], rhs[_0]);

    // Fracture volume storage A_i (w_i - w_i^n)/dt: dominant while a cell opens
    // fast, where it makes opening cost pressure. Cells grown inside the step
    // store their full opening; assumes no renumbering within the step.
    if (!external_pressure_ && inner_.storageTerm && current_dt_ > 0.0) {
        const std::vector<double> areas = cellAreas();
        const size_t ncell = std::min(areas.size(), numFractureCells());
        for (size_t i = 0; i < ncell && i < rhs[_1].size() && i < x[_0].size(); ++i) {
            const double w_prev = (i < fracture_width_prev_.size())
                ? fracture_width_prev_[i][0] : 0.0;
            rhs[_1][i][0] -= areas[i] * (x[_0][i][0] - w_prev) / current_dt_;
        }
    }

    // Fischer-Burmeister: the mechanics residual is the complementarity function.
    if (inner_.isFb)
        for (size_t i = 0; i < fb.phi.size(); ++i)
            rhs[_0][i] = -fb.phi[i];

    // Fail recoverably (upstream catch -> timestep chop) instead of aborting GMRES
    // with an opaque "defect=nan"; a non-finite residual means S or rhs is bad.
    for (size_t i = 0; i < rhs[_0].size(); ++i)
        if (!std::isfinite(rhs[_0][i][0]))
            OPM_THROW(std::runtime_error,
                      "Fracture full-system residual non-finite in mechanics row "
                          + std::to_string(i));
    for (size_t i = 0; i < rhs[_1].size(); ++i)
        if (!std::isfinite(rhs[_1][i][0]))
            OPM_THROW(std::runtime_error,
                      "Fracture full-system residual non-finite in pressure row "
                          + std::to_string(i));
    if (numWellEquations() > 0) {
        // a zero well-equation pivot passes the finite checks but makes the
        // preconditioner emit NaN mid-Krylov ("defect=nan")
        const auto wix = M.N() - 1;
        const double wdiag = M[wix][wix][0][0];
        if (!std::isfinite(wdiag) || wdiag == 0.0)
            OPM_THROW(std::runtime_error,
                      "rate_well equation diagonal is " + std::to_string(wdiag)
                          + " (total_wellindex " + std::to_string(total_wellindex_)
                          + ", well_rate " + std::to_string(well_rate_) + ")");
    }
    return rhs;
}

// ----------------------------------------------------------------------------
bool
Fracture::isConverged_(const VectorHP& rhs, const FbRows& fb, const double tol) const
// ----------------------------------------------------------------------------
{
    // The same tolerance for both blocks, on the unscaled residual.
    const double tol_flow = tol;
    const double tol_mech = tol;
    // fb_scaled_convergence tests phi_i/(1+|a_i|+|b_i|), the same accuracy on every cell.
    if (inner_.isFb && inner_.fbScaledConvergence) {
        VectorHP rtest(rhs);
        for (size_t i = 0; i < fb.scale.size() && i < rtest[_0].size(); ++i)
            rtest[_0][i][0] /= fb.scale[i];
        return convergence_test(rtest, tol_mech, tol_flow, inner_.verbosity);
    }
    return convergence_test(rhs, tol_mech, tol_flow, inner_.verbosity);
}

// ----------------------------------------------------------------------------
bool
Fracture::fullSystemIteration(const double tol, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    ++last_solve_stats_.nonlinear_iterations;

    ++DEBUG_COUNT;

    const InnerSolveOptions& opt = inner_;
    const int nlin_verbosity = opt.verbosity;
    const bool use_ad = opt.useAd;
    updateFlowCoefficients_();

    VectorHP x {fracture_width_, fracture_pressure_};
    validate_fracture_state(x, prm_);

    ResVector mech_rhs(fracture_width_);
    normalFractureTraction(mech_rhs, false);

    const ContactUpdate contact = updateContactSet_(x, mech_rhs, nlin_iteration);
    const std::vector<int>& closed_cells = contact.closed;
    const bool closed_cells_changed = contact.changed;
    const bool first_iteration = (nlin_iteration == 0);
    const FbRows fb = assembleMechanicsRows_(x, mech_rhs, closed_cells, first_iteration);
    assembleFlowRows_(mech_rhs, closed_cells);
    assembleSystemMatrix_(closed_cells, first_iteration);
    auto& S = *S_;

    VectorHP dx = x;
    dx = 0; // Newton update for the current nonlinear iterate
    VectorHP rhs = assembleResidual_(x, mech_rhs, fb);
    auto rhs_org(rhs);

    if (isConverged_(rhs, fb, tol))
        return true;

    // Optional symmetric block scaling of the coupled system (opt-in; default
    // "none" reproduces the unscaled behaviour). Computed from the current
    // diagonal blocks and applied in place to S and rhs; the Newton update dx is
    // scaled back after the linear solve. The convergence test above is done on
    // the unscaled residual so tolerances keep their physical meaning.
    const BlockScaling block_scaling =
        compute_block_scaling(S, opt.scaling);
    if (block_scaling.active()) {
        apply_block_scaling(S, rhs, block_scaling);
    }

    const std::string& solver_type = opt.linearSolver;
    Dune::InverseOperatorResult iores;

    // Direct dense coupled solve (opt-in, solver.linsolver.solver = "direct"): for
    // small fractures the dense BEM-coupled system is best solved directly — Krylov
    // methods struggle with its O(1e11) conditioning and the only robust
    // preconditioner (exact Schur) costs the same as a direct solve anyway. Guarded
    // by direct_max_cells so large fractures fall back to the iterative path.
    bool solved_direct = false;
    if (solver_type == "direct") {
        const int direct_max_cells = opt.directMaxCells;
        const size_t ncoupled = fracture_width_.size();
        if (static_cast<int>(ncoupled) <= direct_max_cells) {
            solved_direct = solve_coupled_direct(S, rhs, dx);
            iores.converged = solved_direct;
            iores.iterations = 1;
            ++last_solve_stats_.linear_solves;
            if (!solved_direct && nlin_verbosity > 0)
                OpmLog::debug("Direct coupled solve failed (singular); will report non-convergence");
        } else if (nlin_verbosity > 0) {
            OpmLog::debug(fmt::format("Direct coupled solve skipped: {} cells > direct_max_cells {}; using iterative solver",
                                      ncoupled, direct_max_cells));
        }
    }

    if (!solved_direct) {
    // solve system equations
    if(first_iteration){
        S_linop_ = std::make_unique<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>>(S); // store for use in preconditioner
    }
    const Opm::PropertyTree& lprm = opt.precond;
    // When solver=="direct" but we fell back here (too large / singular), pick a
    // real Krylov type for the iterative path.
    const std::string iter_solver_type = (solver_type == "direct")
        ? opt.directFallbackSolver
        : solver_type;
    bool new_precond = first_iteration || closed_cells_changed;
    if(new_precond ){
        frac_flow_precond_.reset(new Opm::FractureMechanicsPreconditioner(S, lprm));
    }else{
        frac_flow_precond_->update(S, opt.updateMechOnReuse);
    }
    // make possible use abs tolerance for linear solver
    double res = rhs.two_norm2();
    const double linsolve_tol = opt.linTol;
    const double linsolve_atol = opt.linAtol;
    double lintol = std::max(linsolve_tol, linsolve_atol / res);

    const int max_iter = opt.linMaxIter;
    const int verbosity = opt.linVerbosity;
    if (nlin_verbosity > 1) {
      int num_closed_cells = std::accumulate(closed_cells.begin(), closed_cells.end(), 0);
      OpmLog::debug(fmt::format("Nonlinear iteration: {} rhs:  {} {} closed_cells: {} cells m: {} cells p: {}",
                                nlin_iteration, rhs[_0].infinity_norm(), rhs[_1].infinity_norm(),
                                num_closed_cells, rhs[_0].size(), rhs[_1].size()));
    }

    auto& precond = *frac_flow_precond_;
    if (new_precond) {
        psolver_ = setupLinearSolver(iter_solver_type, *S_linop_, precond, prm_, lintol, max_iter, verbosity);
    }

    // Dump a snapshot of the (physical, unscaled) coupled system either on explicit
    // request (write_coupled_linear_system) or when the solve is "struggling" — i.e.
    // the nonlinear iteration count reached dump_case_min_nonlin_iter (>0). The
    // latter captures a *relevant* hard case for later offline/standalone study
    // (replay tool, --contact-report) without dumping every system.
    const int dump_case_min_nlin = opt.dumpMinNonlinIter;
    const bool want_manual_dump = opt.writeCoupledSystem
        || (dump_case_min_nlin > 0 && nlin_iteration >= dump_case_min_nlin);
    if (want_manual_dump) {
        if (block_scaling.active()) scale_system_matrix(S, block_scaling, -1);
        dump_linear_system_snapshot(linear_system_dump_stem(prm_, nlin_iteration, "manual"),
                                    S,
                                    rhs_org,
                                    x,
                                    closed_cells,
                                    use_ad,
                                    prm_,
                                    makePressureAssemblyInput());
        if (block_scaling.active()) scale_system_matrix(S, block_scaling, +1);
    }

    if (opt.debugChecksum) {
        Hash64 hs, hr;
        for (size_t i = 0; i != S[_0][_0].N(); ++i)
            for (size_t j = 0; j != S[_0][_0].M(); ++j)
                hs(S[_0][_0][i][j]);
        hash_sparse(hs, S[_0][_1]);
        hash_sparse(hs, S[_1][_0]);
        hash_sparse(hs, S[_1][_1]);
        for (size_t i = 0; i != rhs[_0].size(); ++i) hr(rhs[_0][i][0]);
        for (size_t i = 0; i != rhs[_1].size(); ++i) hr(rhs[_1][i][0]);
        // Characterize the preconditioner as an operator: apply it to a fixed
        // probe and hash the result, so a stale/garbage internal state shows up
        // even though S and rhs are identical.
        Hash64 hp;
        {
            VectorHP probe(rhs), out(rhs);
            for (size_t i = 0; i != probe[_0].size(); ++i) probe[_0][i][0] = 1.0;
            for (size_t i = 0; i != probe[_1].size(); ++i) probe[_1][i][0] = 1.0;
            out = 0;
            precond.apply(out, probe);
            for (size_t i = 0; i != out[_0].size(); ++i) hp(out[_0][i][0]);
            for (size_t i = 0; i != out[_1].size(); ++i) hp(out[_1][i][0]);
        }
        OpmLog::debug(fmt::format("linsys checksum nlin={} S={:x} rhs={:x} prec={:x} nm={} np={}",
                                  nlin_iteration, hs.h, hr.h, hp.h, rhs[_0].size(), rhs[_1].size()));
    }

    // Pristine copy of the (possibly scaled) rhs fed to the solver, so each retry
    // restarts from the same system (apply() overwrites rhs with the residual).
    const VectorHP rhs_solve_org = rhs;
    const std::string& failure_policy = opt.failurePolicy;

    // Run a solver from a clean state and accumulate stats.
    auto run_solver = [&](LinearSolverBase& solver) {
        rhs = rhs_solve_org;
        dx = 0;
        solver.apply(dx, rhs, iores); // NB: will modify 'rhs'
        ++last_solve_stats_.linear_solves;
        last_solve_stats_.linear_iterations += iores.iterations;
    };

    {
        OPM_TIMEBLOCK(SolveCoupledSystem);
        try {
            run_solver(*psolver_);
        } catch (const std::exception& e1) {
            // Identify which preconditioner block generates non-finite values -
            // the Krylov abort ("defect=nan") is otherwise opaque.
            try {
                VectorHP vprobe(rhs_solve_org), dprobe(rhs_solve_org);
                vprobe = 0;
                precond.apply(vprobe, dprobe);
                size_t bad_m = 0, bad_f = 0;
                for (size_t i = 0; i < vprobe[_0].size(); ++i)
                    if (!std::isfinite(vprobe[_0][i][0])) ++bad_m;
                for (size_t i = 0; i < vprobe[_1].size(); ++i)
                    if (!std::isfinite(vprobe[_1][i][0])) ++bad_f;
                OpmLog::warning(fmt::format("Coupled solver abort ({}); precond probe non-finite: mech {}/{}, flow {}/{}",
                                            e1.what(), bad_m, vprobe[_0].size(), bad_f, vprobe[_1].size()));
            } catch (const std::exception& ep) {
                OpmLog::warning(fmt::format("Coupled solver abort ({}); precond probe itself threw: {}",
                                            e1.what(), ep.what()));
            }
            const bool allow_gmres_fallback = opt.gmresFallback;
            iores.converged = false;
            if (allow_gmres_fallback && iter_solver_type == "bicgstab") {
                psolver_ = setupLinearSolver("gmres", *S_linop_, precond, prm_, lintol, max_iter, verbosity);
                if (failure_policy == "ladder") {
                    // A Krylov abort (overflow -> defect=nan) is just another
                    // failed rung: fall through to the rescue ladder instead of
                    // aborting the nonlinear iteration.
                    try {
                        run_solver(*psolver_);
                    } catch (const std::exception&) {
                        iores.converged = false;
                    }
                } else {
                    run_solver(*psolver_);
                }
            } else if (failure_policy != "ladder") {
                throw;
            }
        }
    }

    if (!iores.converged) {
        ++last_solve_stats_.linear_solve_failures;
    }

    // Opt-in robustness ladder (solver.linsolver.failure_policy = "ladder"):
    // if the configured solve did not converge, try progressively more robust
    // rescue strategies before giving up. Each rung is transient — the persistent
    // solver/preconditioner are untouched, and S is rebuilt in place at the top of
    // the next nonlinear iteration, so no cleanup is required. Rungs operate on the
    // already-scaled S/rhs, so the dx unscaling below applies uniformly.
    if (!iores.converged && failure_policy == "ladder") {
        const auto fixed_stress_prm = [&]() {
            Opm::PropertyTree p = lprm;
            p.put("mode_policy", std::string("manual"));
            p.put("mech_first", std::string("true"));
            p.put("fixed_stress", std::string("true"));
            return p;
        };
        // rung 1: flexible GMRES with the existing preconditioner.
        try {
            auto s = setupLinearSolver("fgmres", *S_linop_, precond, prm_, lintol, max_iter, verbosity);
            run_solver(*s);
        } catch (const std::exception&) { /* fall through to next rung */ }

        // rung 2: fixed-stress (Schur) preconditioner + flexible GMRES.
        if (!iores.converged) {
            try {
                auto fsprm = fixed_stress_prm();
                Opm::FractureMechanicsPreconditioner fsprecond(S, fsprm);
                auto s = setupLinearSolver("fgmres", *S_linop_, fsprecond, prm_, lintol, max_iter, verbosity);
                run_solver(*s);
            } catch (const std::exception&) {}
        }

        // rung 3: drop the flow->mech coupling block C (decoupled / Picard) and
        // solve the block-triangular system with a fixed-stress preconditioner.
        if (!iores.converged) {
            try {
                S[_1][_0] = 0; // overwritten again at the top of the next iteration
                auto dcprm = fixed_stress_prm();
                Opm::FractureMechanicsPreconditioner dcprecond(S, dcprm);
                auto s = setupLinearSolver("fgmres", *S_linop_, dcprecond, prm_, lintol, max_iter, verbosity);
                run_solver(*s);
            } catch (const std::exception&) {}
        }

        if (iores.converged) {
            ++last_solve_stats_.ladder_rescues;
            if (nlin_verbosity > 0) {
                OpmLog::debug(fmt::format("Fracture coupled linear solve rescued by fallback ladder at nlin iter {} ({} iters)",
                                          nlin_iteration, iores.iterations));
            }
        }
    }
    } // end if(!solved_direct) — iterative coupled solve path

    if (!iores.converged) {
        // Restore the physical (unscaled) system so the dumped failure case is
        // directly replayable. The function throws right after, so no re-scale.
        if (block_scaling.active()) scale_system_matrix(S, block_scaling, -1);
        std::string dump_info_filename;
        if (opt.dumpOnFailure) {
            const auto dump_stem = linear_system_dump_stem(prm_, nlin_iteration, "linear_failure");
            dump_info_filename = linear_system_info_filename(dump_stem);
            dump_linear_system_snapshot(dump_stem,
                                        S,
                                        rhs_org,
                                        x,
                                        closed_cells,
                                        use_ad,
                                        prm_,
                                        makePressureAssemblyInput(),
                                        &iores);
            OpmLog::error(fmt::format("Fracture coupled linear solver failed. Snapshot info written to {}",
                                      dump_info_filename));
        }

        last_solve_stats_.converged = false;
        OPM_THROW(std::runtime_error,
                  "Fracture coupled linear solver failed to converge after "
                      + std::to_string(iores.iterations)
                      + " iterations."
                      + (dump_info_filename.empty()
                             ? std::string()
                             : " Snapshot info: " + dump_info_filename));
    }

    // Undo the block scaling: the solver returned y, the true Newton update is
    // dx = D y. All subsequent step control (damping, clamps) acts on dx.
    if (block_scaling.active()) {
        for (auto& v : dx[_0]) v[0] *= block_scaling.mechanics;
        for (auto& v : dx[_1]) v[0] *= block_scaling.pressure;
    }

    if (nlin_verbosity > 2) {
        OpmLog::debug(fmt::format("x:  {} {}\ndx: {} {}\nrhs after linsolve:  {} {}",
                                  x[_0].infinity_norm(), x[_1].infinity_norm(),
                                  dx[_0].infinity_norm(), dx[_1].infinity_norm(),
                                  rhs[_0].infinity_norm(), rhs[_1].infinity_norm()));
    }
    if (nlin_verbosity > 2 && !block_scaling.active()) {
      // Skipped when block scaling is active: S has been overwritten with D S D
      // in place, so this unscaled residual check would be meaningless.
      auto x_new(dx);
      auto rhs2_tmp(rhs_org);
      S.mmv(x_new, rhs2_tmp);
      OpmLog::debug(fmt::format("Check res : rhs after linsolve:  {} {}",
                                rhs2_tmp[_0].infinity_norm(), rhs2_tmp[_1].infinity_norm()));


    }


    limitStep_(x, dx, mech_rhs, fb);
    applyUpdate_(x, dx);
    return false;
}

// ----------------------------------------------------------------------------
void
Fracture::limitStep_(const VectorHP& x, VectorHP& dx, const ResVector& mech_rhs, const FbRows& fb) const
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    const int nlin_verbosity = opt.verbosity;
    double step_fac = opt.damping;
    // coupled_step_limit > 0 caps the per-block relative change ||dx||/||x||, a
    // cheap trust-region surrogate for the badly conditioned keep-C step.
    const double coupled_step_limit = opt.coupledStepLimit;
    if (coupled_step_limit > 0.0) {
        const double xw = x[_0].infinity_norm(), xp = x[_1].infinity_norm();
        const double dw = dx[_0].infinity_norm(), dp = dx[_1].infinity_norm();
        const double fw = (xw > 0.0) ? dw / xw : 0.0;
        const double fp = (xp > 0.0) ? dp / xp : 0.0;
        const double fmax = std::max(fw, fp);
        if (fmax > coupled_step_limit) {
            const double lim = std::max(coupled_step_limit / fmax, 1e-3);
            step_fac *= lim;
            if (nlin_verbosity > 1)
                OpmLog::debug(fmt::format("Coupled step limited: rel-change {} > {}, scaling step by {}",
                                          fmax, coupled_step_limit, lim));
        }
    }
    // FB backtracking: far from the solution the full semismooth Newton step
    // overshoots; damp it until the merit ||phi||^2 decreases.
    if (opt.isFb && opt.fbLineSearch) {
        const FMatrix& Aorig = fractureMatrix();
        const size_t n = Aorig.N();
        double merit0 = 0.0;
        for (size_t i = 0; i < n; ++i) merit0 += fb.phi[i] * fb.phi[i];
        const double s_min = opt.fbLineSearchMin;
        double sstep = 1.0, best_s = 1.0;
        double best_merit = std::numeric_limits<double>::max();
        ResVector w_t(n), p_t(n);
        while (sstep >= s_min) {
            for (size_t j = 0; j < n; ++j)
                w_t[j] = std::max(0.0, x[_0][j][0] + sstep * dx[_0][j][0]);
            for (size_t i = 0; i < n; ++i)
                p_t[i] = x[_1][i][0] + sstep * dx[_1][i][0];
            const ResVector lam = contactForce(Aorig, w_t, p_t, mech_rhs);
            double merit = 0.0;
            for (size_t i = 0; i < n; ++i) {
                const double c_i = std::max(std::abs(Aorig[i][i]), 1e-30);
                const double a = c_i * w_t[i][0];
                const double b = lam[i][0];
                const double phi = a + b - std::sqrt(a * a + b * b);
                merit += phi * phi;
            }
            if (merit < best_merit) { best_merit = merit; best_s = sstep; }
            if (merit <= (1.0 - 1e-4 * sstep) * merit0) break; // sufficient decrease
            sstep *= 0.5;
        }
        step_fac *= best_s;
        if (nlin_verbosity > 1 && best_s < 1.0)
            OpmLog::debug(fmt::format("FB line search: step {} merit {} -> {}", best_s, merit0, best_merit));
    }
    if (nlin_verbosity > 2) {
        OpmLog::debug(fmt::format("fac: {}", step_fac));
    }
    dx *= step_fac;
    auto& dx0 = dx[_0];
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        dx0[i][0] = std::max(-opt.maxDwidth, std::min(opt.maxDwidth, dx0[i][0]));
    }
    auto& dx1 = dx[_1];
    for (size_t i = 0; i != fracture_pressure_.size(); ++i) {
        dx1[i][0] = std::max(-opt.maxDp, std::min(opt.maxDp, dx1[i][0]));
    }
}

// ----------------------------------------------------------------------------
void
Fracture::applyUpdate_(VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------
{
    x += dx;
    if (inner_.verbosity > 3) {
        OpmLog::debug(fmt::format("after: dx: {} {}", dx[_0].infinity_norm(), dx[_1].infinity_norm()));
    }
    fracture_width_ = x[_0];
    // Policy "none" permits negative aperture by design, and FB satisfies w >= 0
    // at convergence (fb_clamp_width=false skips the projection on its iterates).
    const bool allow_negative_width = (inner_.closedCellPolicy == "none")
        || (inner_.isFb && !inner_.fbClampWidth);
    if (!allow_negative_width) {
        for (size_t i = 0; i != fracture_width_.size(); ++i) {
            fracture_width_[i][0] = std::max(0.0, fracture_width_[i][0]); // ensure non-negativity
        }
    }
    fracture_pressure_ = x[_1];
}

}; // end namespace Opm
