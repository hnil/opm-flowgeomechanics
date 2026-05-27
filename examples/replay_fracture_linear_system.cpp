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
using Opm::SystemMatrix;
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
              << " [--scaling=none|blockmax]" << std::endl;
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
                          const std::string& scaling_mode)
{
    const FMatrix aperture = load_dense_matrix(prefix + "_A.txt");
    const SMatrix identity = load_sparse_matrix(prefix + "_I.mtx");
    const SMatrix coupling = load_sparse_matrix(prefix + "_C_active.mtx");
    const SMatrix pressure = load_sparse_matrix(prefix + "_M_active.mtx");
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

            throw std::runtime_error("Unknown option: " + option);
        }

        if (coupled_mode) {
            return replay_coupled_system(prefix, solver_json, scaling_mode);
        }

        return replay_pressure_system(prefix,
                                      pressure_mode.empty() ? "active" : pressure_mode,
                                      solver_json);
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}