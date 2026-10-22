#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dune/istl/matrixmarket.hh>
#include <dune/istl/solvers.hh>

#include <opm/geomech/FractureMechanicsPreconditioner.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace
{
using Opm::FMatrix;
using Opm::FractureMechanicsPreconditioner;
using Opm::SMatrix;
using SystemMatrix = Opm::FractureSystemMatrix;
using Opm::Vector;
using Opm::VectorHP;
using CoupledOperator = Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>;
using CoupledSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;
using FlowOperator = Dune::MatrixAdapter<SMatrix, Vector, Vector>;
using FlowSolver = Dune::FlexibleSolver<FlowOperator>;

struct BlockScaling {
    double mechanics = 1.0;
    double pressure = 1.0;
};

void print_usage(const char* executable)
{
    std::cerr << "Usage: " << executable << " <dump-prefix> [--coupled]"
              << " [--pressure=active|ad|original] [--solver-json=<path>]"
              << " [--scaling=none|blockmax] [--drop-coupling] [--contact-report]" << std::endl;
}

FMatrix load_dense_matrix(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open dense matrix file: " + filename);
    }

    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream line_stream(line);
        std::vector<double> row;
        double value = 0.0;
        while (line_stream >> value) {
            row.push_back(value);
        }

        if (row.empty()) {
            continue;
        }

        if (!rows.empty() && row.size() != rows.front().size()) {
            throw std::runtime_error("Inconsistent row length in dense matrix file: " + filename);
        }

        rows.push_back(std::move(row));
    }

    FMatrix matrix;
    matrix.resize(rows.size(), rows.empty() ? 0 : rows.front().size());
    for (size_t row = 0; row < rows.size(); ++row) {
        for (size_t col = 0; col < rows[row].size(); ++col) {
            matrix[row][col] = rows[row][col];
        }
    }

    return matrix;
}

SMatrix load_sparse_matrix(const std::string& filename)
{
    SMatrix matrix;
    Dune::loadMatrixMarket(matrix, filename);
    return matrix;
}

Vector load_vector(const std::string& filename)
{
    Vector vector;
    Dune::loadMatrixMarket(vector, filename);
    return vector;
}

// Load a dumped closed-cell flag vector (one integer 0/1 per line, as written by
// dump_vector in Fracture_fullSystemIteration.cpp).
std::vector<int> load_int_vector(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open closed-cell file: " + filename);
    }
    std::vector<int> result;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        result.push_back(static_cast<int>(std::lround(std::stod(line))));
    }
    return result;
}

double max_abs_entry(const FMatrix& matrix)
{
    double result = 0.0;
    for (size_t row = 0; row < matrix.N(); ++row) {
        for (size_t col = 0; col < matrix.M(); ++col) {
            result = std::max(result, std::abs(matrix[row][col]));
        }
    }
    return result;
}

double max_abs_entry(const SMatrix& matrix)
{
    double result = 0.0;
    for (auto row_it = matrix.begin(); row_it != matrix.end(); ++row_it) {
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            result = std::max(result, std::abs((*col_it)[0][0]));
        }
    }
    return result;
}

void scale_matrix(FMatrix& matrix, const double factor)
{
    for (size_t row = 0; row < matrix.N(); ++row) {
        for (size_t col = 0; col < matrix.M(); ++col) {
            matrix[row][col] *= factor;
        }
    }
}

void scale_matrix(SMatrix& matrix, const double factor)
{
    for (auto row_it = matrix.begin(); row_it != matrix.end(); ++row_it) {
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
            (*col_it)[0][0] *= factor;
        }
    }
}

void scale_vector(Vector& vector, const double factor)
{
    for (auto& value : vector) {
        value[0] *= factor;
    }
}

BlockScaling compute_block_scaling(const FMatrix& mechanics,
                                   const SMatrix& pressure,
                                   const std::string& scaling_mode)
{
    if (scaling_mode == "none") {
        return {};
    }

    if (scaling_mode != "blockmax") {
        throw std::runtime_error("Unsupported scaling mode: " + scaling_mode);
    }

    const double mechanics_max = max_abs_entry(mechanics);
    const double pressure_max = max_abs_entry(pressure);
    if (mechanics_max <= 0.0 || pressure_max <= 0.0) {
        throw std::runtime_error("Cannot compute block scaling from zero-valued diagonal blocks");
    }

    BlockScaling scaling;
    scaling.mechanics = 1.0 / std::sqrt(mechanics_max);
    scaling.pressure = 1.0 / std::sqrt(pressure_max);
    return scaling;
}

void apply_block_scaling(SystemMatrix& system, VectorHP& rhs, const BlockScaling& scaling)
{
    scale_matrix(system[Dune::Indices::_0][Dune::Indices::_0], scaling.mechanics * scaling.mechanics);
    scale_matrix(system[Dune::Indices::_0][Dune::Indices::_1], scaling.mechanics * scaling.pressure);
    scale_matrix(system[Dune::Indices::_1][Dune::Indices::_0], scaling.mechanics * scaling.pressure);
    scale_matrix(system[Dune::Indices::_1][Dune::Indices::_1], scaling.pressure * scaling.pressure);
    scale_vector(rhs[Dune::Indices::_0], scaling.mechanics);
    scale_vector(rhs[Dune::Indices::_1], scaling.pressure);
}

std::unique_ptr<CoupledSolverBase>
make_coupled_solver(const std::string& solver_type,
                    const CoupledOperator& op,
                    FractureMechanicsPreconditioner& preconditioner,
                    const Opm::PropertyTree& linsolver_prm,
                    const double tolerance,
                    const int max_iter,
                    const int verbosity)
{
    if (solver_type == "bicgstab") {
        return std::make_unique<Dune::BiCGSTABSolver<VectorHP>>(op,
                                                                preconditioner,
                                                                tolerance,
                                                                max_iter,
                                                                verbosity);
    }

    if (solver_type == "gmres") {
        const int restart = linsolver_prm.get<int>("restart", 100);
        return std::make_unique<Dune::RestartedGMResSolver<VectorHP>>(op,
                                                                      preconditioner,
                                                                      tolerance,
                                                                      restart,
                                                                      max_iter,
                                                                      verbosity);
    }

    if (solver_type == "fgmres") {
        const int restart = linsolver_prm.get<int>("restart", 100);
        return std::make_unique<Dune::RestartedFlexibleGMResSolver<VectorHP>>(op,
                                                                              preconditioner,
                                                                              tolerance,
                                                                              restart,
                                                                              max_iter,
                                                                              verbosity);
    }

    throw std::runtime_error("Unsupported coupled solver type: " + solver_type);
}

int replay_coupled_system(const std::string& prefix,
                          const std::string& solver_json,
                          const std::string& scaling_mode,
                          const bool drop_coupling)
{
    const FMatrix aperture = load_dense_matrix(prefix + "_A.txt");
    const SMatrix identity = load_sparse_matrix(prefix + "_I.mtx");
    SMatrix coupling = load_sparse_matrix(prefix + "_C_active.mtx");
    const SMatrix pressure = load_sparse_matrix(prefix + "_M_active.mtx");
    // Mirror the production ladder's drop-C (decoupled / Picard) rescue rung: zero
    // the flow->mech coupling block so the system becomes block-triangular.
    if (drop_coupling) {
        scale_matrix(coupling, 0.0);
    }
    const Vector rhs_w = load_vector(prefix + "_rhs_w.mtx");
    const Vector rhs_p = load_vector(prefix + "_rhs_p.mtx");

    const Opm::PropertyTree linsolver_prm(solver_json.empty() ? prefix + "_linsolver.json" : solver_json);

    const SystemMatrix original_system({{aperture, identity},
                                        {coupling, pressure}});
    const VectorHP original_rhs{rhs_w, rhs_p};

    auto system = original_system;
    auto rhs = original_rhs;
    const auto scaling = compute_block_scaling(aperture, pressure, scaling_mode);
    apply_block_scaling(system, rhs, scaling);

    VectorHP dy{rhs_w, rhs_p};
    dy = 0;

    CoupledOperator op(system);
    FractureMechanicsPreconditioner preconditioner(system,
                                                   linsolver_prm.get_child("preconditioner"));

    const double rhs_norm_sq = rhs.two_norm2();
    const double tol = linsolver_prm.get<double>("tol");
    const double atol = linsolver_prm.get<double>("atol", 0.0);
    const double lintol = std::max(tol, atol / std::max(rhs_norm_sq, 1e-30));
    const int max_iter = linsolver_prm.get<int>("max_iter");
    const int verbosity = linsolver_prm.get<int>("verbosity", 0);

    auto solver = make_coupled_solver(linsolver_prm.get<std::string>("solver"),
                                      op,
                                      preconditioner,
                                      linsolver_prm,
                                      lintol,
                                      max_iter,
                                      verbosity);

    Dune::InverseOperatorResult result;
    auto rhs_work(rhs);
    solver->apply(dy, rhs_work, result);

    auto dx = dy;
    scale_vector(dx[Dune::Indices::_0], scaling.mechanics);
    scale_vector(dx[Dune::Indices::_1], scaling.pressure);

    auto residual = original_rhs;
    original_system.mmv(dx, residual);

    std::cout << "Coupled scaling mode: " << scaling_mode << std::endl;
    std::cout << "Scaling factors: " << scaling.mechanics << " " << scaling.pressure << std::endl;
    std::cout << "Coupled solve converged: " << result.converged << std::endl;
    std::cout << "Coupled iterations: " << result.iterations << std::endl;
    std::cout << "Residual infinity norms: "
              << residual[Dune::Indices::_0].infinity_norm() << " "
              << residual[Dune::Indices::_1].infinity_norm() << std::endl;
    std::cout << "Update infinity norms: "
              << dx[Dune::Indices::_0].infinity_norm() << " "
              << dx[Dune::Indices::_1].infinity_norm() << std::endl;

    return result.converged ? 0 : 2;
}

// Offline contact / fracture-opening report. Loads a dumped nonlinear fracture
// state and analyses the per-cell contact force balance the way Fracture::
// identify_closed does. The dumped rhs_w equals identify_closed's decision
// variable tmp = mech_rhs - A*w - I*p, so a cell's contact tendency is read
// directly from rhs_w[i] without any solve:
//   open  cell closes  when tmp >= close_force_tolerance and width <= 0
//   closed cell reopens when tmp <= -reopen_force_tolerance
// This exposes the "opened-too-large fracture is hard to re-open" problem:
// closed cells whose force balance sits just inside the reopen threshold are
// stuck closed even though they carry no compressive load.
int replay_contact_report(const std::string& prefix, const std::string& solver_json)
{
    const Vector force = load_vector(prefix + "_rhs_w.mtx"); // tmp = mech_rhs - A*w - I*p
    const Vector width = load_vector(prefix + "_x_w.mtx");
    const std::vector<int> closed = load_int_vector(prefix + "_closed_cells.txt");

    const std::size_t n = force.size();
    if (width.size() != n || closed.size() != n) {
        throw std::runtime_error("contact-report: mismatched sizes (force/width/closed)");
    }

    double close_tol = 0.0, reopen_tol = 0.0;
    if (!solver_json.empty()) {
        const Opm::PropertyTree prm(solver_json);
        close_tol = prm.get<double>("preconditioner.close_force_tolerance", 0.0);
        reopen_tol = prm.get<double>("preconditioner.reopen_force_tolerance", close_tol);
    }

    std::size_t n_closed = 0, n_open = 0;
    std::size_t open_barely = 0;     // open cells with near-zero width (half-open)
    std::size_t stuck_closed = 0;    // closed cells not tensile enough to reopen
    std::size_t open_wants_close = 0; // open cells with compressive force >= close_tol
    double max_abs_force = 0.0, max_width = 0.0;
    const double width_eps = 1e-4;   // "barely open" threshold (~initial_fracture_width)

    for (std::size_t i = 0; i < n; ++i) {
        const double f = force[i][0];
        const double w = width[i][0];
        max_abs_force = std::max(max_abs_force, std::abs(f));
        max_width = std::max(max_width, w);
        if (closed[i]) {
            ++n_closed;
            // would reopen only if f <= -reopen_tol; otherwise stuck closed
            if (f > -reopen_tol) ++stuck_closed;
        } else {
            ++n_open;
            if (std::abs(w) < width_eps) ++open_barely;
            if (f >= close_tol) ++open_wants_close;
        }
    }

    std::cout << "Contact / opening report for " << prefix << "\n"
              << "  cells: " << n << " (open " << n_open << ", closed " << n_closed << ")\n"
              << "  tolerances: close=" << close_tol << " reopen=" << reopen_tol << "\n"
              << "  max |force balance| = " << max_abs_force << ", max width = " << max_width << "\n"
              << "  open & barely-open (|w|<" << width_eps << "): " << open_barely << "\n"
              << "  open but force wants close (f>=close_tol): " << open_wants_close << "\n"
              << "  closed & stuck (f>-reopen_tol, cannot reopen): " << stuck_closed
              << " of " << n_closed << " closed\n";
    std::cout << "  => 'hard to open' signature = stuck-closed cells carrying little/no"
                 " compressive load (f near 0)." << std::endl;
    return 0;
}

// Direct dense solve of the full coupled system [[A,I],[C,M]] as a 2N x 2N matrix.
// Diagnostic: if the hard cases solve cleanly here, the coupled system is well
// posed and only the *iterative preconditioner* is at fault — which motivates an
// exact Schur-complement preconditioner using the dense BEM A's exact inverse
// (current fixed_stress uses only diag(A), a poor approximation for a dense A).
int replay_direct_solve(const std::string& prefix)
{
    const FMatrix A = load_dense_matrix(prefix + "_A.txt");
    const SMatrix I = load_sparse_matrix(prefix + "_I.mtx");
    const SMatrix C = load_sparse_matrix(prefix + "_C_active.mtx");
    const SMatrix M = load_sparse_matrix(prefix + "_M_active.mtx");
    const Vector bw = load_vector(prefix + "_rhs_w.mtx");
    const Vector bp = load_vector(prefix + "_rhs_p.mtx");
    const std::size_t n = A.N();

    // Assemble the dense 2n x 2n coupled matrix and rhs.
    FMatrix K(2 * n, 2 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            K[i][j] = A[i][j];
    auto add_sparse = [&](const SMatrix& S, std::size_t roff, std::size_t coff) {
        for (auto row = S.begin(); row != S.end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
                K[roff + row.index()][coff + col.index()] = (*col)[0][0];
    };
    add_sparse(I, 0, n);
    add_sparse(C, n, 0);
    add_sparse(M, n, n);

    Vector rhs(2 * n), x(2 * n);
    for (std::size_t i = 0; i < n; ++i) { rhs[i] = bw[i][0]; rhs[n + i] = bp[i][0]; }
    x = 0;

    FMatrix Kfac(K); // solve() factorizes in place
    Kfac.solve(x, rhs);

    // Residual r = rhs - K x (use the unfactored K).
    Vector r(rhs);
    K.mmv(x, r);
    double rw = 0.0, rp = 0.0, xw = 0.0, xp = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        rw = std::max(rw, std::abs(r[i][0]));     xw = std::max(xw, std::abs(x[i][0]));
        rp = std::max(rp, std::abs(r[n + i][0])); xp = std::max(xp, std::abs(x[n + i][0]));
    }
    std::cout << "Direct dense solve of [[A,I],[C,M]] for " << prefix << "\n"
              << "  size 2n = " << 2 * n << "\n"
              << "  residual inf-norm (mech / flow): " << rw << " " << rp << "\n"
              << "  solution inf-norm (w / p):       " << xw << " " << xp << std::endl;
    const double tol = 1e-6 * std::max(1.0, std::max(xw, xp));
    const bool ok = rw < tol && rp < tol;
    std::cout << "  => direct solve " << (ok ? "OK" : "FAILED") << std::endl;
    return ok ? 0 : 2;
}

int replay_pressure_system(const std::string& prefix,
                           const std::string& pressure_mode,
                           const std::string& solver_json)
{
    std::string matrix_file;
    if (pressure_mode == "active") {
        matrix_file = prefix + "_M_active.mtx";
    } else if (pressure_mode == "ad") {
        matrix_file = prefix + "_M_ad.mtx";
    } else if (pressure_mode == "original") {
        matrix_file = prefix + "_M_original.mtx";
    } else {
        throw std::runtime_error("Unsupported pressure replay mode: " + pressure_mode);
    }

    const SMatrix matrix = load_sparse_matrix(matrix_file);
    const Vector rhs = load_vector(prefix + "_rhs_p.mtx");
    const Opm::PropertyTree flow_solver_prm(solver_json.empty() ? prefix + "_flow_solver.json" : solver_json);

    FlowOperator op(matrix);
    FlowSolver solver(op, flow_solver_prm, std::function<Vector()>(), 1);

    Vector solution(rhs);
    solution = 0;

    Dune::InverseOperatorResult result;
    auto rhs_work(rhs);
    solver.apply(solution, rhs_work, result);

    auto residual(rhs);
    matrix.mmv(solution, residual);

    std::cout << "Pressure solve mode: " << pressure_mode << std::endl;
    std::cout << "Pressure solve converged: " << result.converged << std::endl;
    std::cout << "Pressure iterations: " << result.iterations << std::endl;
    std::cout << "Pressure residual infinity norm: " << residual.infinity_norm() << std::endl;
    std::cout << "Pressure solution infinity norm: " << solution.infinity_norm() << std::endl;

    return result.converged ? 0 : 2;
}
} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }

        const std::string prefix = argv[1];
        bool coupled_mode = true;
        std::string pressure_mode;
        std::string solver_json;
        std::string scaling_mode = "none";
        bool drop_coupling = false;
        bool contact_report = false;
        bool direct_solve = false;

        for (int arg = 2; arg < argc; ++arg) {
            const std::string option = argv[arg];
            if (option == "--coupled") {
                coupled_mode = true;
                pressure_mode.clear();
                continue;
            }

            const std::string pressure_prefix = "--pressure=";
            if (option.rfind(pressure_prefix, 0) == 0) {
                coupled_mode = false;
                pressure_mode = option.substr(pressure_prefix.size());
                continue;
            }

            const std::string solver_prefix = "--solver-json=";
            if (option.rfind(solver_prefix, 0) == 0) {
                solver_json = option.substr(solver_prefix.size());
                continue;
            }

            const std::string scaling_prefix = "--scaling=";
            if (option.rfind(scaling_prefix, 0) == 0) {
                scaling_mode = option.substr(scaling_prefix.size());
                continue;
            }

            if (option == "--drop-coupling") {
                drop_coupling = true;
                continue;
            }

            if (option == "--contact-report") {
                contact_report = true;
                continue;
            }

            if (option == "--direct") {
                direct_solve = true;
                continue;
            }

            throw std::runtime_error("Unknown option: " + option);
        }

        if (contact_report) {
            return replay_contact_report(prefix, solver_json);
        }

        if (direct_solve) {
            return replay_direct_solve(prefix);
        }

        if (coupled_mode) {
            return replay_coupled_system(prefix, solver_json, scaling_mode, drop_coupling);
        }

        return replay_pressure_system(prefix,
                                      pressure_mode.empty() ? "active" : pressure_mode,
                                      solver_json);
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}