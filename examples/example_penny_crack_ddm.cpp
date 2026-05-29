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

    Notes on method choice:
        The current standalone benchmark uses the CutDe full-space triangular
        dislocation kernel together with a cell-centered collocation DDM
        discretization. This is a reasonable baseline for coupling to cell-based
        flow models, where pressure is naturally represented per fracture cell.
        If higher mechanical accuracy is needed, the next exact step is usually a
        triangle-integrated or Galerkin DDM operator with dedicated self/near-cell
        treatment. If larger problems dominate runtime, a near/far split with FMM
        or hierarchical compression is the usual acceleration path while keeping
        exact near interactions.
*/

#include <config.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

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
using DenseVector = Dune::DynamicVector<double>;
using DenseElimVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using Clock = std::chrono::steady_clock;

enum class SolverBackend {
    Direct,
    LuReuse,
    MatrixFreeBiCGSTAB,
    MatrixFreeGMRes,
    All,
};

struct CommandLineOptions {
    double E = 1.0e9;
    double nu = 0.25;
    double p = 1.0e6;
    std::vector<int> validation_layers {5, 8};
    std::vector<int> benchmark_layers;
    int profile_layers = 5;
    bool print_profile = true;
    bool run_validation = true;
    SolverBackend backend = SolverBackend::Direct;
    int solve_repeats = 1;
    double cod_tol = 0.15;
    double K1_tol = 0.30;
    double linear_tol = 1e-8;
    int max_iter = 1000;
    int gmres_restart = 100;
    int solver_verbosity = 0;
};

struct SolveTimings {
    double dense_assembly_ms = 0.0;
    double cache_setup_ms = 0.0;
    double factorization_ms = 0.0;
    double preconditioner_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
};

struct SolveStats {
    int iterations = 0;
    int operator_applications = 0;
    bool converged = true;
};

struct SolveResult {
    std::unique_ptr<Grid> grid;
    DenseVector opening;
    double a_eff = 0.0;
    int nc = 0;
    SolveTimings timings;
    SolveStats stats;
};

class DenseLUSolver
{
public:
    class MatrixAccess : public Dune::DynamicMatrix<double>
    {
    public:
        static void factorize(Dune::DynamicMatrix<double>& matrix)
        {
            DenseElimVector tmp(matrix.N());
            typename Dune::DynamicMatrix<double>::template Elim<DenseElimVector> elim(tmp);
            Dune::Simd::Mask<double> nonsing(true);
            Dune::DynamicMatrix<double>::luDecomposition(matrix, elim, nonsing, false, false);
        }
    };

    explicit DenseLUSolver(const Dune::DynamicMatrix<double>& matrix)
        : lu_(matrix)
    {
    }

    void factorize()
    {
        MatrixAccess::factorize(lu_);
    }

    void solve(DenseVector& x, const DenseVector& rhs) const
    {
        x.resize(rhs.size());
        for (int i = 0; i < static_cast<int>(lu_.rows()); ++i) {
            x[i] = rhs[i];
            for (int j = 0; j < i; ++j) {
                x[i] -= lu_[i][j] * x[j];
            }
        }

        for (int i = static_cast<int>(lu_.rows()) - 1; i >= 0; --i) {
            for (size_t j = static_cast<size_t>(i) + 1; j < lu_.rows(); ++j) {
                x[i] -= lu_[i][j] * x[j];
            }
            x[i] /= lu_[i][i];
        }
    }

private:
    Dune::DynamicMatrix<double> lu_;
};

struct MatrixFreeGeometryCache {
    std::vector<Dune::FieldVector<double, 3>> centers;
    std::vector<Dune::FieldVector<double, 3>> normals;
    std::vector<std::array<ddm::Real3, 3>> tris;
};

MatrixFreeGeometryCache build_geometry_cache(const Grid& grid)
{
    MatrixFreeGeometryCache cache;
    Mapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
    const int nc = grid.leafGridView().size(0);

    cache.centers.resize(nc);
    cache.normals.resize(nc);
    cache.tris.resize(nc);

    for (auto elem : elements(grid.leafGridView())) {
        const int idx = mapper.index(elem);
        cache.centers[idx] = elem.geometry().center();
        cache.normals[idx] = ddm::normalOfElement(elem);
        cache.tris[idx] = ddm::getTri(elem);
    }

    return cache;
}

double matrix_free_normal_traction(const Dune::FieldVector<double, 3>& center,
                                   const Dune::FieldVector<double, 3>& normal,
                                   const std::array<ddm::Real3, 3>& tri,
                                   const double slip_value,
                                   const double E,
                                   const double nu)
{
    Dune::FieldVector<double, 3> slip(0.0);
    // The CutDe wrapper uses TDCS ordering with the normal-opening component first.
    slip[0] = slip_value;
    const auto strain = ddm::TDStrainFSTri(center, tri, slip, nu);
    const auto stress = ddm::strainToStress(E, nu, strain);
    return ddm::tractionSymTensor(stress, normal);
}

class MatrixFreeDDMOperator : public Dune::LinearOperator<DenseVector, DenseVector>
{
public:
    using field_type = typename DenseVector::field_type;

    MatrixFreeDDMOperator(const MatrixFreeGeometryCache& cache, const double E, const double nu)
        : cache_(cache)
        , E_(E)
        , nu_(nu)
    {
    }

    void apply(const DenseVector& x, DenseVector& y) const override
    {
        const int nc = static_cast<int>(cache_.centers.size());
        y.resize(nc);

#pragma omp parallel for
        for (int row = 0; row < nc; ++row) {
            double value = 0.0;
            for (int col = 0; col < nc; ++col) {
                value += matrix_free_normal_traction(cache_.centers[row],
                                                     cache_.normals[row],
                                                     cache_.tris[col],
                                                     x[col],
                                                     E_,
                                                     nu_);
            }
            y[row] = value;
        }

        ++apply_calls_;
    }

    void applyscaleadd(field_type alpha, const DenseVector& x, DenseVector& y) const override
    {
        DenseVector tmp;
        apply(x, tmp);
        if (y.size() != tmp.size()) {
            y.resize(tmp.size());
            y = 0.0;
        }
        for (size_t i = 0; i < tmp.size(); ++i)
            y[i] += alpha * tmp[i];
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    void resetApplyCount() const
    {
        apply_calls_ = 0;
    }

    int applyCount() const
    {
        return apply_calls_;
    }

private:
    const MatrixFreeGeometryCache& cache_;
    double E_;
    double nu_;
    mutable int apply_calls_ = 0;
};

class DiagonalPreconditioner : public Dune::Preconditioner<DenseVector, DenseVector>
{
public:
    explicit DiagonalPreconditioner(const DenseVector& diagonal)
        : inv_diag_(diagonal.size())
    {
        for (size_t i = 0; i < diagonal.size(); ++i) {
            const double value = diagonal[i];
            inv_diag_[i] = (std::abs(value) > 1e-30) ? 1.0 / value : 1.0;
        }
    }

    void apply(DenseVector& v, const DenseVector& d) override
    {
        v.resize(d.size());
        for (size_t i = 0; i < d.size(); ++i)
            v[i] = inv_diag_[i] * d[i];
    }

    void pre(DenseVector&, DenseVector&) override {}
    void post(DenseVector&) override {}

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    DenseVector inv_diag_;
};

DenseVector build_diagonal(const MatrixFreeGeometryCache& cache, const double E, const double nu)
{
    DenseVector diagonal(cache.centers.size());

#pragma omp parallel for
    for (int row = 0; row < static_cast<int>(cache.centers.size()); ++row) {
        diagonal[row] = matrix_free_normal_traction(cache.centers[row],
                                                    cache.normals[row],
                                                    cache.tris[row],
                                                    1.0,
                                                    E,
                                                    nu);
    }

    return diagonal;
}

double elapsed_ms(const Clock::time_point start, const Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

bool starts_with(const std::string& value, const std::string& prefix)
{
    return value.rfind(prefix, 0) == 0;
}

std::vector<int> parse_int_list(const std::string& text)
{
    std::vector<int> values;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            continue;
        }
        const int value = std::stoi(token);
        if (value <= 0) {
            throw std::invalid_argument("Layers must be positive integers");
        }
        values.push_back(value);
    }

    if (values.empty()) {
        throw std::invalid_argument("Expected at least one layer value");
    }

    return values;
}

SolverBackend parse_solver_backend(const std::string& text)
{
    if (text == "direct") {
        return SolverBackend::Direct;
    }
    if (text == "lu_reuse") {
        return SolverBackend::LuReuse;
    }
    if (text == "matrix_free_bicgstab") {
        return SolverBackend::MatrixFreeBiCGSTAB;
    }
    if (text == "matrix_free_gmres") {
        return SolverBackend::MatrixFreeGMRes;
    }
    if (text == "all") {
        return SolverBackend::All;
    }

    throw std::invalid_argument("Unknown solver backend: " + text);
}

const char* backend_name(const SolverBackend backend)
{
    switch (backend) {
    case SolverBackend::Direct:
        return "direct";
    case SolverBackend::LuReuse:
        return "lu_reuse";
    case SolverBackend::MatrixFreeBiCGSTAB:
        return "matrix_free_bicgstab";
    case SolverBackend::MatrixFreeGMRes:
        return "matrix_free_gmres";
    case SolverBackend::All:
        return "all";
    }
    return "unknown";
}

std::vector<SolverBackend> selected_backends(const SolverBackend backend)
{
    if (backend == SolverBackend::All) {
        return {
            SolverBackend::Direct,
            SolverBackend::LuReuse,
            SolverBackend::MatrixFreeBiCGSTAB,
            SolverBackend::MatrixFreeGMRes,
        };
    }
    return {backend};
}

bool is_matrix_free_backend(const SolverBackend backend)
{
    return backend == SolverBackend::MatrixFreeBiCGSTAB
        || backend == SolverBackend::MatrixFreeGMRes;
}

void print_usage(const char* executable)
{
    std::cout << "Usage: " << executable << " [options]\n"
              << "  --solver=direct|lu_reuse|matrix_free_bicgstab|matrix_free_gmres|all\n"
              << "                                   Solver backend to benchmark\n"
              << "  --layers=5,8                    Validation layers\n"
              << "  --benchmark-layers=18,24        Extra benchmark-only layer counts\n"
              << "  --profile-layers=5              Layer count for COD profile output\n"
              << "  --solve-repeats=1               Number of repeated solves for the same matrix\n"
              << "  --young=<value>                 Young's modulus [Pa]\n"
              << "  --nu=<value>                    Poisson ratio\n"
              << "  --pressure=<value>              Net pressure [Pa]\n"
              << "  --cod-tol=<value>               Relative COD tolerance\n"
              << "  --k1-tol=<value>                Relative K_I tolerance\n"
              << "  --tol=<value>                   Iterative solver residual reduction target\n"
              << "  --max-iter=<value>              Iterative solver iteration limit\n"
              << "  --gmres-restart=<value>         GMRES restart length\n"
              << "  --solver-verbosity=<value>      Iterative solver verbosity\n"
              << "  --skip-profile                  Skip COD profile printout\n"
              << "  --skip-validation               Skip analytic validation checks\n"
              << "  --help                          Show this help\n";
}

CommandLineOptions parse_options(int argc, char** argv)
{
    CommandLineOptions options;
    for (int arg = 1; arg < argc; ++arg) {
        const std::string option = argv[arg];
        if (option == "--help" || option == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
        if (option == "--skip-profile") {
            options.print_profile = false;
            continue;
        }
        if (option == "--skip-validation") {
            options.run_validation = false;
            continue;
        }
        if (starts_with(option, "--solver=")) {
            options.backend = parse_solver_backend(option.substr(9));
            continue;
        }
        if (starts_with(option, "--layers=")) {
            options.validation_layers = parse_int_list(option.substr(9));
            continue;
        }
        if (starts_with(option, "--benchmark-layers=")) {
            options.benchmark_layers = parse_int_list(option.substr(19));
            continue;
        }
        if (starts_with(option, "--profile-layers=")) {
            options.profile_layers = std::stoi(option.substr(17));
            if (options.profile_layers <= 0) {
                throw std::invalid_argument("Profile layers must be positive");
            }
            continue;
        }
        if (starts_with(option, "--solve-repeats=")) {
            options.solve_repeats = std::stoi(option.substr(16));
            if (options.solve_repeats <= 0) {
                throw std::invalid_argument("Solve repeats must be positive");
            }
            continue;
        }
        if (starts_with(option, "--young=")) {
            options.E = std::stod(option.substr(8));
            continue;
        }
        if (starts_with(option, "--nu=")) {
            options.nu = std::stod(option.substr(5));
            continue;
        }
        if (starts_with(option, "--pressure=")) {
            options.p = std::stod(option.substr(11));
            continue;
        }
        if (starts_with(option, "--cod-tol=")) {
            options.cod_tol = std::stod(option.substr(10));
            continue;
        }
        if (starts_with(option, "--k1-tol=")) {
            options.K1_tol = std::stod(option.substr(9));
            continue;
        }
        if (starts_with(option, "--tol=")) {
            options.linear_tol = std::stod(option.substr(6));
            continue;
        }
        if (starts_with(option, "--max-iter=")) {
            options.max_iter = std::stoi(option.substr(11));
            if (options.max_iter <= 0) {
                throw std::invalid_argument("Max iterations must be positive");
            }
            continue;
        }
        if (starts_with(option, "--gmres-restart=")) {
            options.gmres_restart = std::stoi(option.substr(16));
            if (options.gmres_restart <= 0) {
                throw std::invalid_argument("GMRES restart must be positive");
            }
            continue;
        }
        if (starts_with(option, "--solver-verbosity=")) {
            options.solver_verbosity = std::stoi(option.substr(19));
            continue;
        }

        throw std::invalid_argument("Unknown option: " + option);
    }

    return options;
}

SolveResult solve_penny_problem(const int layers,
                                const CommandLineOptions& options,
                                const SolverBackend backend)
{
    const auto total_start = Clock::now();
    const double h = 1.0;

    Opm::RegularTrimesh trimesh(layers,
                                {0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.5, std::sqrt(3.0) / 2.0, 0.0},
                                {h, h});

    auto [grid, cellmap, bcells] = trimesh.createDuneGrid(0, {}, false);
    const int nc = grid->leafGridView().size(0);

    double r_max = 0.0;
    for (auto elem : elements(grid->leafGridView())) {
        auto center = elem.geometry().center();
        const double r = std::sqrt(center[0] * center[0] + center[1] * center[1]);
        r_max = std::max(r_max, r);
    }
    const double a_eff = r_max + 0.5 * h;

    DenseVector rhs(nc, -options.p);
    DenseVector opening(nc, 0.0);

    SolveTimings timings;
    SolveStats stats;

    if (backend == SolverBackend::Direct) {
        Dune::DynamicMatrix<double> A(nc, nc, 0.0);
        const auto assembly_start = Clock::now();
        ddm::assembleMatrix_fast(A, options.E, options.nu, *grid);
        const auto assembly_end = Clock::now();
        timings.dense_assembly_ms = elapsed_ms(assembly_start, assembly_end);

        const auto solve_start = Clock::now();
        for (int repeat = 0; repeat < options.solve_repeats; ++repeat) {
            auto A_work = A;
            opening = 0.0;
            A_work.solve(opening, rhs);
        }
        const auto solve_end = Clock::now();
        timings.solve_ms = elapsed_ms(solve_start, solve_end);
        stats.iterations = options.solve_repeats;
    }
    else if (backend == SolverBackend::LuReuse) {
        Dune::DynamicMatrix<double> A(nc, nc, 0.0);
        const auto assembly_start = Clock::now();
        ddm::assembleMatrix_fast(A, options.E, options.nu, *grid);
        const auto assembly_end = Clock::now();
        timings.dense_assembly_ms = elapsed_ms(assembly_start, assembly_end);

        auto lu_solver = DenseLUSolver(A);

        const auto factor_start = Clock::now();
        lu_solver.factorize();
        const auto factor_end = Clock::now();
        timings.factorization_ms = elapsed_ms(factor_start, factor_end);

        const auto solve_start = Clock::now();
        for (int repeat = 0; repeat < options.solve_repeats; ++repeat) {
            opening = 0.0;
            lu_solver.solve(opening, rhs);
        }
        const auto solve_end = Clock::now();
        timings.solve_ms = elapsed_ms(solve_start, solve_end);
        stats.iterations = options.solve_repeats;
    }
    else if (is_matrix_free_backend(backend)) {
        const auto cache_start = Clock::now();
        const auto cache = build_geometry_cache(*grid);
        const auto cache_end = Clock::now();
        timings.cache_setup_ms = elapsed_ms(cache_start, cache_end);

        const auto preconditioner_start = Clock::now();
        const auto diagonal = build_diagonal(cache, options.E, options.nu);
        DiagonalPreconditioner preconditioner(diagonal);
        const auto preconditioner_end = Clock::now();
        timings.preconditioner_ms = elapsed_ms(preconditioner_start, preconditioner_end);

        MatrixFreeDDMOperator op(cache, options.E, options.nu);
        const auto solve_start = Clock::now();

        int total_iterations = 0;
        int total_apply_count = 0;
        bool converged = true;

        if (backend == SolverBackend::MatrixFreeBiCGSTAB) {
            Dune::BiCGSTABSolver<DenseVector> solver(op,
                                                     preconditioner,
                                                     options.linear_tol,
                                                     options.max_iter,
                                                     options.solver_verbosity);
            for (int repeat = 0; repeat < options.solve_repeats; ++repeat) {
                opening = 0.0;
                op.resetApplyCount();
                Dune::InverseOperatorResult result;
                auto rhs_work = rhs;
                solver.apply(opening, rhs_work, result);
                total_iterations += result.iterations;
                total_apply_count += op.applyCount();
                converged = converged && result.converged;
            }
        }
        else {
            Dune::RestartedGMResSolver<DenseVector> solver(op,
                                                           preconditioner,
                                                           options.linear_tol,
                                                           options.gmres_restart,
                                                           options.max_iter,
                                                           options.solver_verbosity);
            for (int repeat = 0; repeat < options.solve_repeats; ++repeat) {
                opening = 0.0;
                op.resetApplyCount();
                Dune::InverseOperatorResult result;
                auto rhs_work = rhs;
                solver.apply(opening, rhs_work, result);
                total_iterations += result.iterations;
                total_apply_count += op.applyCount();
                converged = converged && result.converged;
            }
        }

        const auto solve_end = Clock::now();
        timings.solve_ms = elapsed_ms(solve_start, solve_end);
        stats.iterations = total_iterations;
        stats.operator_applications = total_apply_count;
        stats.converged = converged;
    }

    timings.total_ms = elapsed_ms(total_start, Clock::now());

    return {std::move(grid), std::move(opening), a_eff, nc, timings, stats};
}

void print_benchmark_summary(const int layers,
                             const SolveResult& result,
                             const SolverBackend backend,
                             const int solve_repeats)
{
    std::cout << "  Benchmark backend=" << backend_name(backend)
              << " layers=" << layers
              << " cells=" << result.nc
              << " repeats=" << solve_repeats << "\n";
    if (backend == SolverBackend::LuReuse) {
        std::cout << "    dense_assembly_ms=" << result.timings.dense_assembly_ms << "\n"
                  << "    factorization_ms=" << result.timings.factorization_ms << "\n"
                  << "    triangular_solve_ms=" << result.timings.solve_ms;
        if (solve_repeats > 0) {
            std::cout << " (per_rhs=" << result.timings.solve_ms / solve_repeats << ")";
        }
        std::cout << "\n";
    }
    else if (backend == SolverBackend::Direct) {
        std::cout << "    dense_assembly_ms=" << result.timings.dense_assembly_ms << "\n"
                  << "    direct_solve_ms=" << result.timings.solve_ms;
        if (solve_repeats > 0) {
            std::cout << " (per_rhs=" << result.timings.solve_ms / solve_repeats << ")";
        }
        std::cout << "\n";
    }
    else {
        std::cout << "    cache_setup_ms=" << result.timings.cache_setup_ms << "\n"
                  << "    diagonal_preconditioner_ms=" << result.timings.preconditioner_ms << "\n"
                  << "    iterative_solve_ms=" << result.timings.solve_ms;
        if (solve_repeats > 0) {
            std::cout << " (per_rhs=" << result.timings.solve_ms / solve_repeats << ")";
        }
        std::cout << "\n"
                  << "    converged=" << (result.stats.converged ? "true" : "false") << "\n"
                  << "    iterations=" << result.stats.iterations;
        if (solve_repeats > 0) {
            std::cout << " (per_rhs=" << static_cast<double>(result.stats.iterations) / solve_repeats << ")";
        }
        std::cout << "\n"
                  << "    operator_applications=" << result.stats.operator_applications;
        if (solve_repeats > 0) {
            std::cout << " (per_rhs=" << static_cast<double>(result.stats.operator_applications) / solve_repeats << ")";
        }
        std::cout << "\n";
    }
    std::cout << "    total_ms=" << result.timings.total_ms << std::endl;
}

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
         const CommandLineOptions& options,
         const SolverBackend backend,
         const double cod_rel_tol,
         const double K1_rel_tol)
{
    const double h = 1.0; // element edge length
    auto result = solve_penny_problem(layers, options, backend);
    auto& grid = *result.grid;

    const int nc = result.nc;
    std::cout << "  layers=" << layers
              << "  cells=" << nc << std::endl;
    Mapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    const double a_eff = result.a_eff;

    std::cout << "  Effective crack radius a_eff = " << a_eff << " m" << std::endl;
    print_benchmark_summary(layers, result, backend, options.solve_repeats);

    if (is_matrix_free_backend(backend) && !result.stats.converged) {
        std::cerr << "  FAILED: iterative solver did not converge" << std::endl;
        return false;
    }

    // ------------------------------------------------------------------
    // Compare DDM openings with the analytic COD profile.
    // ------------------------------------------------------------------
    const double w_max_analytic = analytic_cod(0.0, a_eff, options.p, options.E, options.nu);
    std::cout << "  Analytic max COD (r=0): " << w_max_analytic << " m" << std::endl;

    double sum_rel_err   = 0.0;
    double max_rel_err   = 0.0;
    int    n_interior    = 0;

    for (auto elem : elements(grid.leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double rx   = center[0];
        const double ry   = center[1];
        const double r    = std::sqrt(rx * rx + ry * ry);

        const double w_ddm      = result.opening[idx];
        const double w_analytic = analytic_cod(r, a_eff, options.p, options.E, options.nu);

        // Only compare interior elements (well away from the discretised rim)
        if (r < a_eff - 1.5 * h && w_analytic > 1e-16) {
            const double abs_err = std::abs(w_ddm - w_analytic);
            const double rel_err = abs_err / w_analytic;
            sum_rel_err += rel_err;
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
    const double K1_analytic = analytic_K1(a_eff, options.p);
    std::cout << "  Analytic K_I = " << K1_analytic << " Pa*sqrt(m)" << std::endl;

    double K1_sum   = 0.0;
    int    K1_count = 0;

    for (auto elem : elements(grid.leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double rx   = center[0];
        const double ry   = center[1];
        const double r    = std::sqrt(rx * rx + ry * ry);

        // Only use elements in the outermost ring (tip-distance in [0, 1.5*h])
        const double s = a_eff - r;
        if (s > 0.0 && s < 1.5 * h) {
            const double w   = result.opening[idx];
            const double K1  = K1_from_cod(w, s, options.E, options.nu);
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
                  const CommandLineOptions& options,
                  const SolverBackend backend,
                  const int solve_repeats)
{
    auto result = solve_penny_problem(layers, options, backend);
    auto& grid = *result.grid;
    Mapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    print_benchmark_summary(layers, result, backend, solve_repeats);

    std::cout << "\n  COD profile (layers=" << layers
              << ", a_eff=" << result.a_eff << "):\n";
    std::cout << "    r [m]      DDM w [m]   Analytic w [m]  Rel.err [%]\n";

    // Collect (r, w) pairs, sort by r, print
    std::vector<std::pair<double, double>> rv;
    for (auto elem : elements(grid.leafGridView())) {
        const int  idx    = mapper.index(elem);
        auto       center = elem.geometry().center();
        const double r    = std::sqrt(center[0] * center[0] + center[1] * center[1]);
        rv.push_back({r, result.opening[idx]});
    }
    std::sort(rv.begin(), rv.end());

    for (auto& [r, w_ddm] : rv) {
        const double w_an  = analytic_cod(r, result.a_eff, options.p, options.E, options.nu);
        const double rel   = (w_an > 1e-16) ? std::abs(w_ddm - w_an) / w_an * 100.0 : 0.0;
        std::cout << "    " << r << "  " << w_ddm << "  " << w_an
                  << "  " << rel << "\n";
    }
}

} // anonymous namespace

// ============================================================================
int
main(int argc, char** argv)
// ============================================================================
{
    CommandLineOptions options;
    try {
        options = parse_options(argc, argv);
    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n\n";
        print_usage(argv[0]);
        return 2;
    }

    std::cout << "======================================================\n";
    std::cout << " Penny-shaped crack DDM vs analytic comparison\n";
    std::cout << "======================================================\n";
    std::cout << "  E  = " << options.E  << " Pa\n";
    std::cout << "  nu = " << options.nu << "\n";
    std::cout << "  p  = " << options.p  << " Pa\n";
    std::cout << "  solver = " << backend_name(options.backend) << "\n";
    std::cout << "  solve_repeats = " << options.solve_repeats << "\n\n";

    int failures = 0;
    for (const auto backend : selected_backends(options.backend)) {
        std::cout << "=== Backend: " << backend_name(backend) << " ===\n";

        if (options.print_profile) {
            std::cout << "--- COD profile (layers=" << options.profile_layers << ") ---\n";
            print_cod_profile(options.profile_layers, options, backend, options.solve_repeats);
        }

        if (!options.print_profile && !options.run_validation && options.benchmark_layers.empty()) {
            std::cout << "--- Benchmark only (layers=" << options.profile_layers << ") ---\n";
            const auto result = solve_penny_problem(options.profile_layers, options, backend);
            print_benchmark_summary(options.profile_layers,
                                    result,
                                    backend,
                                    options.solve_repeats);
        }

        for (const int layers : options.benchmark_layers) {
            std::cout << "--- Benchmark: layers=" << layers << " ---\n";
            const auto result = solve_penny_problem(layers, options, backend);
            print_benchmark_summary(layers, result, backend, options.solve_repeats);
        }

        if (options.run_validation) {
            for (const int layers : options.validation_layers) {
                std::cout << "\n--- Test: layers=" << layers << " ---\n";
                if (!run_test(layers, options, backend, options.cod_tol, options.K1_tol)) {
                    ++failures;
                }
            }
        }

        std::cout << std::endl;
    }

    std::cout << "\n======================================================\n";
    std::cout << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << "\n";
    std::cout << "======================================================\n";

    return failures;
}
