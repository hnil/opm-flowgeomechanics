#include "config.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <numeric>

#include <dune/common/indices.hh> // needed for _0, _1, etc.
//
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>


#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/DiagonalScalar.hpp>
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

using namespace std;
#define EXP_REF 0
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

// using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
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
// template<class MatrixType> void dump_matrix(const MatrixType& m, const char* const name)
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
// template<>
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
                     const ResVector& pressure,// is realy head
                     const ResVector& aperture,
                     const std::vector<int>& closed_cells,
                     const std::vector<double>& reservoir_mobility,
                     double min_width)
{
    OPM_TIMEFUNCTION();
    // create C if not done already
    if (!Cptr)
        Cptr = initCouplingMatrixSparsity(htrans, num_cells, num_wells);
    // if (!Cptr) Cptr = std::make_unique<Opm::Fracture::Matrix>(*M); // copy sparsity of M

    // @@ NB: If the implementation of the pressure matrix changes (i.e. `assemblePressure()`),
    //        the code below might need to be updated accordingly as well.
    auto& C = *Cptr;
    C = 0;
    for (const auto& e : htrans) {
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
        const double h1 = std::max(aperture[i][0],min_width);//aperture[i] + min_width;
        const double h2 = std::max(aperture[j][0],min_width);//aperture[j] + min_width;
        const double p1 = pressure[i];//-fracture_dgh[i];
        const double p2 = pressure[j];//-fracture_dgh[j];; pressure should already be head

        const double q = (h1 * h1 * h1) * (h2 * h2 * h2) * (t1 * t2); // numerator
        const double d1q = 3 * (h1 * h1) * (h2 * h2 * h2) * (t1 * t2)* dh1;
        const double d2q = 3 * (h1 * h1 * h1) * (h2 * h2) * (t1 * t2)* dh2;

        //const double r = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
        const double r = 12 * (h1 * h1 * h1 * t1 + h2 * h2* h2 * t2); // denominator
        const double d1r = 36 * (h1 * h1) * t1 *dh1;
        const double d2r = 36 * (h2 * h2) * t2 *dh2;
        //const double d1r = 24 * (h1) * t1;
        //const double d2r = 24 * (h2) * t2;
        assert(r >= 0.0);
        assert(q >= 0.0);
        const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
        const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);

        const double mobility = 0.5 * (reservoir_mobility[i] + reservoir_mobility[j]);
        //const double krull = 1; // 1e4; // @@ Not sure if this should be removed?
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
    // std::cout << "Residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
    // std::cout << "Residual norm[1] is " << res[_1].infinity_norm()<<std::endl;

    // std::cout << "tol mech is: " << tol_mech << std::endl;
    // std::cout << "tol flow is: " << tol_flow << std::endl;

    bool converged = res[_0].infinity_norm() < tol_mech && res[_1].infinity_norm() < tol_flow;
    if(converged && verbosity > 0){
      std::cout << "Converged with residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
      std::cout << "Converged with residual norm[1] is " << res[_1].infinity_norm()<<std::endl;
    }
    // else{
    //   std::cout << "Not converged with residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
    //   std::cout << "Not converged with residual norm[1] is " << res[_1].infinity_norm()<<std::endl;
    // }
    return converged;
    // return res.infinity_norm() < tol;
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
    const double fac_min = 1e-1; // @@ 1e-2;
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
    //Dune::DynamicMatrix result(A);
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
Fracture::identify_closed(const FMatrix& A, const VectorHP& x, const ResVector& rhs, const int nwells)
// ----------------------------------------------------------------------------
{
    
    OPM_TIMEFUNCTION();
    std::string closing_type = prm_.get<std::string>("solver.closing_type", "org");
    const std::string closed_cell_policy = prm_.get<std::string>("solver.closed_cell_policy", "sticky");
    const double close_force_tolerance = prm_.get<double>("solver.close_force_tolerance", 0.0);
    const double reopen_force_tolerance = prm_.get<double>("solver.reopen_force_tolerance", close_force_tolerance);
    const double reopen_width_tolerance = prm_.get<double>("solver.reopen_width_tolerance", 0.0);
    if (closed_cell_policy != "sticky" && closed_cell_policy != "legacy") {
        OPM_THROW(std::runtime_error,
                  "Unknown closed cell policy: " + closed_cell_policy);
    }
    const bool legacy_closed_cells = (closed_cell_policy == "legacy");
    ResVector tmp(rhs);
    const auto I = makeIdentity(A.N(), nwells);
    const auto was_closed = [this](const size_t index) {
        return index < closed_cells_.size() && closed_cells_[index] != 0;
    };

    // computing rhs - A x[0] - I x[1]
    std::vector<int> result;
    if(closing_type == "org"){
        const ResVector& h = x[_0];
        const ResVector& p = x[_1];
        // maybe use other pressure if closed
        A.mmv(h, tmp);
        I.mmv(p, tmp);
        for (size_t i = 0; i != A.N(); ++i) {
            if (legacy_closed_cells) {
                result.push_back(tmp[i] >= 0.0 && h[i] <= 0.0);
                continue;
            }
            const bool close_candidate = tmp[i] >= close_force_tolerance && h[i] <= 0.0;
            // Re-opening is force-driven. A closed cell has its width pinned at 0,
            // so a width gate (h[i] >= reopen_width_tolerance) with a positive
            // tolerance can never be satisfied and would lock the cell closed
            // forever, suppressing fracture growth. reopen_width_tolerance is kept
            // for backward compatibility but no longer gates on the pinned width.
            (void)reopen_width_tolerance;
            const bool reopen_candidate = tmp[i] <= -reopen_force_tolerance;
            result.push_back(was_closed(i) ? !reopen_candidate : close_candidate);
        }
    } else if(closing_type == "open") {
            result.resize(A.N(), 1);
    } else if(closing_type == "simple"){
     // sum forces over faces
        result.resize(A.N(), 0);
        //double sum_force = 0.0;
        for (auto element : Dune::elements(grid_->leafGridView())) {
           // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
            // Dune::mcmgElementLayout>::index(element);
            int i = grid_->leafGridView().indexSet().index(element);
            //double area = element.geometry().volume();
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
        
    // std::fill(result.begin(), result.end(), 1); // @@@
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

    // Zero out coupling matrix columns for closed cells
    // auto& C = *coupling_matrix_;
    // for (size_t col = 0; col < closed_cells.size(); ++col)
    //     if (closed_cells[col])
    //         for (size_t row = 0; row < C.N(); ++row)
    //             if (C.exists(row, col))
    //                 C[row][col] = 0;
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
    input.control_type = prm_.get_child("control").get<std::string>("type");
    input.num_well_equations = nw;
    input.perfinj.assign(perfinj_.begin(), perfinj_.end());
    input.total_wellindex = total_wellindex_;
    input.mobility_water_perf = mobility_water_perf_;
    input.face_mobility.resize(nc + nw, 0.0);
    for (size_t i = 0; i < reservoir_mobility_.size() && i < nc; ++i)
        input.face_mobility[i] = reservoir_mobility_[i];
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
        input.fracture_width[i] = fracture_width_[i][0];
        input.fracture_pressure[i] = fracture_pressure_[i][0];
        if (i < fracture_dgh_.size()) {
            input.fracture_pressure[i] -= fracture_dgh_[i];
        }
    }

    for (size_t i = 0; i < nc; ++i) {
        if (use_fluid_system_water_properties) {
            const auto fluid_properties = fracture_water_property_evaluator_(i, fracture_pressure_[i][0]);
            input.density[i] = fluid_properties.first;
            input.viscosity[i] = fluid_properties.second;
        }
        else {
            input.density[i] = {reservoir_density_[i], 0.0};
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
bool
Fracture::fullSystemIteration(const double tol, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    ++last_solve_stats_.nonlinear_iterations;

    ++DEBUG_COUNT;
    // min_width_ = prm_.get<double>("solver.min_width"); // min with only used for flow calculations
    //const double max_width = prm_.get<double>("solver.max_width");

    // update pressure matrix with the current values of `fracture_width_` and
    // `fracture_pressure_`
    // std::cout << "---- Assemble pressure" << std::endl;
    const bool use_ad = prm_.get<bool>("solver.use_ad_pressure_assembly", false);
    updateLeakoff();
    if (!use_ad) {
        assemblePressure(); // update pressure matrix (original method)
    }
    addSource(); // update right-hand side of pressure system;
    

    // std::cout << "---- Various " << std::endl;
    //  initialize vector of unknown, and vector represnting direction in tangent space
    VectorHP x {fracture_width_, fracture_pressure_};
    // dump_vector(x, "w", "p", true); // dump current state of fracture
    // dump_vector(x, debug_filename("w_").c_str(), debug_filename("p_").c_str());
    validate_fracture_state(x, prm_);
    const int nlin_verbosity = prm_.get<int>("solver.verbosity", 0);

    ResVector mech_rhs(fracture_width_);
    normalFractureTraction(mech_rhs, false);

    // make a version of the fracture matrix that has trivial equations for closed cells
    std::vector<int> closed_cells
        = identify_closed(fractureMatrix(), x, mech_rhs, numWellEquations());
    if (!closed_cells_.empty() && closed_cells_.size() == closed_cells.size()) {
        // Anti-chatter guard. Only re-openings (closed -> open) are throttled;
        // establishing contact (open -> closed) is always allowed so the contact
        // set can grow as fast as the physics demands. Throttling growth here was
        // observed to freeze the closed set at empty and prevent convergence on
        // fully-compressed fractures (regression decks). See count_reopened_cells.
        const size_t toggled_cells = count_reopened_cells(closed_cells_, closed_cells);
        const size_t max_toggle_count = static_cast<size_t>(std::max(
            0, prm_.get<int>("solver.max_closed_cell_toggle_count", 1000000000)));
        const double max_toggle_fraction = prm_.get<double>("solver.max_closed_cell_toggle_fraction", 1.0);
        const size_t fraction_limit = static_cast<size_t>(max_toggle_fraction * closed_cells.size());
        const size_t allowed_toggles = std::min(max_toggle_count, std::max<size_t>(fraction_limit, 0));
        if (toggled_cells > allowed_toggles) {
            if (nlin_verbosity > 0) {
                std::cout << "Keeping previous closed-cell state after " << toggled_cells
                          << " re-openings exceeded guard limit " << allowed_toggles << std::endl;
            }
            closed_cells = closed_cells_;
        }
    }
    // closed as changed
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
    //closed_cells_changed = true; // @@@
    // Contact-chatter metric: count cells that flipped state vs the previous
    // nonlinear iteration (same grid only; size mismatch = grid change, skip).
    if (closed_cells_.size() == closed_cells.size()) {
        last_solve_stats_.closed_cell_toggles +=
            static_cast<int>(count_toggled_cells(closed_cells_, closed_cells));
    }
    closed_cells_ = closed_cells;
    // dump_vector(closed_cells, debug_filename("closed_cells_").c_str());
    // dump_vector(closed_cells, "closed_cells", true);
    // const auto A = modified_fracture_matrix(fractureMatrix(), closed_cells);
    //bool update_mech_matrix = nlin_iteration > 0 || closed_cells_changed;
    bool first_iteration = (nlin_iteration == 0);
    if (first_iteration) {
        A_.reset(new FMatrix(modified_fracture_matrix(fractureMatrix(), closed_cells)));
    }else{
        modified_fracture_matrix(*A_, fractureMatrix(), closed_cells);
    }
    const auto& A = *A_;

    // Closed cells are held fixed in the mechanics residual and Jacobian block.
    for (size_t i = 0; i != closed_cells.size(); ++i)
        if (closed_cells[i])
            mech_rhs[i] = 0;

    if (use_ad) {
        // Use the standalone AD assembler for both pressure matrix and coupling matrix
        assemblePressureAndCouplingAD(closed_cells);
    } else {
        // update the coupling matrix (possibly create it if not already initialized)
        auto fracture_head(fracture_pressure_);
        assert(fracture_dgh_.size() == (fracture_pressure_.size()-numWellEquations()));
        for (size_t i = 0; i < fracture_dgh_.size(); ++i) {
            fracture_head[i] = fracture_pressure_[i] - fracture_dgh_[i];
        }

        updateCouplingMatrix(coupling_matrix_,
                             pressure_matrix_->N() - numWellEquations(), // num cells
                             numWellEquations(), // num wells
                             htrans_,
                             fracture_head,//Note use head here
                             fracture_width_,
                             closed_cells,
                             reservoir_mobility_,
                             min_width_);
    }
    // setup the full system
    if(prm_.get<bool>("solver.drop_fluid_mech_linearization",true)){
      *coupling_matrix_ = 0;
    }else{
            if(mech_rhs.infinity_norm() > prm_.get<double>("solver.drop_tol_h",1e1) ||
                 rhs_pressure_.infinity_norm() > prm_.get<double>("solver.drop_tol_p",1.0) ){
        *coupling_matrix_ = 0;
        // maybe change linear solver
      }

    }
    const auto& M = *pressure_matrix_;
    const auto& C = *coupling_matrix_;
    
    // const auto I = makeIdentity(A.N(), numWellEquations(), 1, closed_cells);
    if (first_iteration) {
        I_.reset(new SMatrix(makeIdentity(A.N(), numWellEquations(), 1, closed_cells)));
    }else{
        modifyIdentity(*I_, 1, closed_cells);
    }
    const auto& I = *I_;
#if EXP_REF    
    const auto& row1 = Dune::makeMultiTypeBlockVectorRef(A, I);
    const auto& row2 = Dune::makeMultiTypeBlockVectorRef(C, M);
    const auto& SS = Dune::makeMultiTypeBlockMatrixRef(row1, row2);
#endif
  
    //const auto SS_ref = std::make_unique<SystemMatrixRef>(SS); // copy to make sure we have a contiguous storage for the preconditioner
    // NB need to se how to modify this matrix as litle as possible and change only the part
    // of preconditioner need
    // system Jacobian (with cross term)  @@ should S be included as a member variable of Fracture?
    //S_.reset(
    //    new SystemMatrix({{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //                      {C, M}})); // flow system
    // this make a copy of the A,I,C,M

    if(first_iteration){
#if EXP_REF           
        S_ = std::make_unique<SystemMatrix>(SS); // copy to make sure we have a contiguous storage for the preconditioner
#else
        //S_ = std::make_unique<SystemMatrix>({{A, I}, // mechanics system (since A is negative, we leave I positive here)
        //                                      {C, M}}); // flow system
        S_.reset(new SystemMatrix({{A, I},                       
                                   {C, M}}));

#endif
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
    

    // SystemMatrix S = {{A, I},// mechanics system (since A is negative, we leave I positive here)
    //                   {C, M}};
    auto& S = *S_;

    // system equations
    //SystemMatrix S0 = S;
    //S0[_1][_0] = 0; // the equations themselves have no cross term
    // // dump_vector(rhs, debug_filename("rhs_w_").c_str(), debug_filename("rhs_p_").c_str());
    // // dump_vector(rhs, "rhs_w", "rhs_p", true);
    
    //S0.mmv(x, rhs); // rhs = rhs - S0 * x;   (we are working in the tanget plane)
    // write explicit using S

    VectorHP dx = x;
    dx = 0; // Newton update for the current nonlinear iterate

    // Assemble the Newton residual after the Jacobian blocks have been refreshed.
    VectorHP rhs {x};
    rhs[_0] = mech_rhs;
    rhs[_1] = rhs_pressure_;
    
       // S[_1][_0].mmv(x[_0], rhs[_1]);
        S[_1][_1].mmv(x[_1], rhs[_1]);
        S[_0][_0].mmv(x[_0], rhs[_0]);
        S[_0][_1].mmv(x[_1], rhs[_0]);   
    //}
   
    auto rhs_org(rhs);

    
    // Verify that equations have been chosen correctly
    // for (size_t i = 0; i != fracture_width_.size(); ++i)
    //   assert( !(rhs[_0][i] >= 0.0 && x[_0][i] <= 0.0));

    // check if system is already at a converged state (in which case we return immediately)
    //
    // for convergence, we use 'tol' directly for the mechanics system (where
    // residuals are expected to scale with pressure), and 'tol * ||M||' for the
    // flow system (where residuals scale with M*p)
    const double tol_flow = tol;//*std::max(A.infinity_norm(), M.infinity_norm()) * std::numeric_limits<double>::epsilon();
    const double tol_mech = tol;// * M.infinity_norm();
    if (convergence_test(rhs,
                         tol_mech,
                         tol_flow, nlin_verbosity))
        return true;

    // Optional symmetric block scaling of the coupled system (opt-in; default
    // "none" reproduces the unscaled behaviour). Computed from the current
    // diagonal blocks and applied in place to S and rhs; the Newton update dx is
    // scaled back after the linear solve. The convergence test above is done on
    // the unscaled residual so tolerances keep their physical meaning.
    const BlockScaling block_scaling =
        compute_block_scaling(S, prm_.get<std::string>("solver.linsolver.scaling", "none"));
    if (block_scaling.active()) {
        apply_block_scaling(S, rhs, block_scaling);
    }

    // solve system equations
    //const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S);
    if(first_iteration){
        S_linop_ = std::make_unique<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>>(S); // store for use in preconditioner
    }
    //TailoredPrecondDiag precond_old(S); // cannot be 'const' due to BiCGstabsolver interface
    const std::string solver_type = prm_.get<std::string>("solver.linsolver.solver");
    Opm::PropertyTree lprm = prm_.get_child("solver.linsolver.preconditioner");
    //Opm::FractureMechanicsPreconditioner precond(S, lprm);
    Dune::InverseOperatorResult iores;
    bool new_precond = first_iteration || closed_cells_changed;
    if(new_precond ){
        frac_flow_precond_.reset(new Opm::FractureMechanicsPreconditioner(S, lprm));
    }else{
        const bool update_mech_on_reuse = lprm.get<bool>("update_mech_on_reuse", false);
        frac_flow_precond_->update(S, update_mech_on_reuse);
    }
    // make possible use abs tolerance for linear solver
    double res = rhs.two_norm2();
    const double linsolve_tol = prm_.get<double>("solver.linsolver.tol");
    const double linsolve_atol = prm_.get<double>("solver.linsolver.atol");
    double lintol = std::max(linsolve_tol, linsolve_atol / res);
    
    const int max_iter = prm_.get<double>("solver.linsolver.max_iter");
    const int verbosity = prm_.get<double>("solver.linsolver.verbosity", 0);
    if (nlin_verbosity > 1) {
      int num_closed_cells = std::accumulate(closed_cells.begin(), closed_cells.end(), 0);
      std::cout << "Nonlinear iteration: " << nlin_iteration;// << std::endl; 
      std::cout << " rhs:  " << rhs[_0].infinity_norm() << " " << rhs[_1].infinity_norm();
      std::cout <<  " closed_cells: " << num_closed_cells;
      std::cout << " cells m: " << rhs[_0].size();
      std::cout << " cells p: " << rhs[_1].size();
      std::cout << std::endl;
    }

    //std::unique_ptr<LinearSolverBase> psolver;
    auto& precond = *frac_flow_precond_;
    //using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;
    if (new_precond) {
        psolver_ = setupLinearSolver(solver_type, *S_linop_, precond, prm_, lintol, max_iter, verbosity);
    }

    // Dump a snapshot of the (physical, unscaled) coupled system either on explicit
    // request (write_coupled_linear_system) or when the solve is "struggling" — i.e.
    // the nonlinear iteration count reached dump_case_min_nonlin_iter (>0). The
    // latter captures a *relevant* hard case for later offline/standalone study
    // (replay tool, --contact-report) without dumping every system.
    const int dump_case_min_nlin = prm_.get<int>("solver.dump_case_min_nonlin_iter", 0);
    const bool want_manual_dump = prm_.get<bool>("solver.write_coupled_linear_system", false)
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

    // Pristine copy of the (possibly scaled) rhs fed to the solver, so each retry
    // restarts from the same system (apply() overwrites rhs with the residual).
    const VectorHP rhs_solve_org = rhs;
    const std::string failure_policy =
        prm_.get<std::string>("solver.linsolver.failure_policy", "throw");

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
        } catch (const std::exception&) {
            const bool allow_gmres_fallback = prm_.get<bool>(
                "solver.linsolver.fallback_gmres_on_breakdown", true);
            if (allow_gmres_fallback && solver_type == "bicgstab") {
                psolver_ = setupLinearSolver("gmres", *S_linop_, precond, prm_, lintol, max_iter, verbosity);
                run_solver(*psolver_);
            } else {
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
                std::cout << "Fracture coupled linear solve rescued by fallback ladder at nlin iter "
                          << nlin_iteration << " (" << iores.iterations << " iters)" << std::endl;
            }
        }
    }

    if (!iores.converged) {
        // Restore the physical (unscaled) system so the dumped failure case is
        // directly replayable. The function throws right after, so no re-scale.
        if (block_scaling.active()) scale_system_matrix(S, block_scaling, -1);
        std::string dump_info_filename;
        if (prm_.get<bool>("solver.dump_linear_system_on_failure", true)) {
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
            std::cerr << "Fracture coupled linear solver failed. Snapshot info written to "
                      << dump_info_filename << std::endl;
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
        std::cout << "x:  " << x[_0].infinity_norm() << " " << x[_1].infinity_norm() << std::endl;
        std::cout << "dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm() << std::endl;
        std::cout << "rhs after linsolve:  " << rhs[_0].infinity_norm() <<
          " " << rhs[_1].infinity_norm()
                  << std::endl;
    }
    if (nlin_verbosity > 2 && !block_scaling.active()) {
      // Skipped when block scaling is active: S has been overwritten with D S D
      // in place, so this unscaled residual check would be meaningless.
      auto x_new(dx);
      //x_new += dx;
      auto rhs2_tmp(rhs_org);
      S.mmv(x_new, rhs2_tmp);
      std::cout << "Check res : rhs after linsolve:  " << rhs2_tmp[_0].infinity_norm() <<
        " " << rhs2_tmp[_1].infinity_norm() << std::endl;


    }
    
    // the following is a heuristic way to limit stepsize to stay within convergence radius
    const double damping = prm_.get<double>("solver.damping");
    const double step_fac = damping; // estimate_step_fac(x, dx) * damping;
    // std::cout << "fac: " << step_fac << std::endl;
    if (nlin_verbosity > 2) {
        std::cout << "fac: " << step_fac << std::endl;
    }
    dx *= step_fac;
    auto& dx0 = dx[_0];
    double max_dwidth = prm_.get<double>("solver.max_dwidth");
    // for(auto& dx0v: dx0){
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        dx0[i][0] = std::max(-max_dwidth, std::min(max_dwidth, dx0[i][0]));
    }
    auto& dx1 = dx[_1];
    double max_dp = prm_.get<double>("solver.max_dp");
    // for(auto& dx1v: dx1){
    for (size_t i = 0; i != fracture_pressure_.size(); ++i) {
        dx1[i][0] = std::max(-max_dp, std::min(max_dp, dx1[i][0]));
    }


    // dump_vector(dx, debug_filename("dx_w_").c_str(), debug_filename("dx_p_").c_str());
    // dump_vector(dx, "dx_w", "dx_p", true);
    x += dx;
    if (nlin_verbosity > 3) {
        std::cout << "after: dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm()
                  << std::endl;
    }
    // copying modified variables back to member variables
    fracture_width_ = x[_0];
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        fracture_width_[i][0] = std::max(0.0, fracture_width_[i][0]); // ensure non-negativity
        // fracture_width_[i][0] = std::min(max_width, fracture_width_[i][0]); // ensure non-negativity
    }
    fracture_pressure_ = x[_1];
    // if (nlin_verbosity > 2) {
    //   VectorHP rhs_tmp {x}; // same size as system, content set below
    //   normalFractureTraction(rhs_tmp[_0], false); // right-hand side equals the normal fracture traction
    //   rhs_tmp[_1] = rhs_pressure_;
    //   for (size_t i = 0; i != closed_cells.size(); ++i)
    //     if (closed_cells[i])
    //         rhs_tmp[_0][i] = 0;
    //   S0.mmv(x, rhs_tmp);
    //   std::cout << "New: rhs after linsolve:  " << rhs_tmp[_0].infinity_norm() << 
    //     " " << rhs_tmp[_1].infinity_norm() << std::endl;
    // }

    
    return false;
}

}; // end namespace Opm
