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

#include "config.h"

#include <opm/geomech/FractureCoupledLinearSolver.hpp>

#include <chrono>
#include <cmath>
#include <stdexcept>

#include <dune/common/indices.hh>
#include <dune/istl/solvers.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/Fracture.hpp>

#include <fmt/format.h>

namespace
{
using Dune::Indices::_0;
using Dune::Indices::_1;

using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using FMatrix = Dune::DynamicMatrix<double>;
using VectorHP = Opm::VectorHP;
using SystemMatrix = Opm::FractureSystemMatrix;

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

} // anonymous namespace

extern "C" void dgemv_(const char* trans, const int* m, const int* n, const double* alpha,
                       const double* a, const int* lda, const double* x, const int* incx,
                       const double* beta, double* y, const int* incy);

namespace
{
// S with the dense mechanics block copied to contiguous storage; refresh() after S changes.
class DenseBlasSystemOperator : public Dune::LinearOperator<Opm::VectorHP, Opm::VectorHP>
{
public:
    using SystemMatrix = Opm::FractureSystemMatrix;
    explicit DenseBlasSystemOperator(const SystemMatrix& S) : S_(S) {}

    void refresh()
    {
        using namespace Dune::Indices;
        const auto& A = S_[_0][_0];
        n_ = static_cast<int>(A.N());
        a_.resize(static_cast<std::size_t>(n_) * n_);
        for (int i = 0; i < n_; ++i)
            std::copy(A[i].begin(), A[i].end(), a_.begin() + static_cast<std::size_t>(i) * n_);
    }

    void apply(const Opm::VectorHP& x, Opm::VectorHP& y) const override
    {
        y = 0.0;
        applyscaleadd(1.0, x, y);
    }

    void applyscaleadd(field_type alpha, const Opm::VectorHP& x, Opm::VectorHP& y) const override
    {
        using namespace Dune::Indices;
        if (n_ > 0) {
            // row-major A is column-major A^T
            const char trans = 'T';
            const int one = 1;
            const double beta = 1.0;
            dgemv_(&trans, &n_, &n_, &alpha, a_.data(), &n_, &x[_0][0][0], &one, &beta, &y[_0][0][0], &one);
        }
        S_[_0][_1].usmv(alpha, x[_1], y[_0]);
        S_[_1][_0].usmv(alpha, x[_0], y[_1]);
        S_[_1][_1].usmv(alpha, x[_1], y[_1]);
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const SystemMatrix& S_;
    int n_{0};
    std::vector<double> a_;
};
} // namespace

namespace Opm
{

FractureCoupledLinearSolver::Options
FractureCoupledLinearSolver::Options::read(const PropertyTree& prm)
{
    Options o;
    o.scaling = prm.get<std::string>("solver.linsolver.scaling", "none");
    if (o.scaling != "none" && o.scaling != "blockmax") {
        OPM_THROW(std::runtime_error, "Unsupported fracture linsolver scaling mode: " + o.scaling);
    }
    o.solver = prm.get<std::string>("solver.linsolver.solver");
    o.precond = prm.get_child("solver.linsolver.preconditioner");
    o.updateMechOnReuse = o.precond.get<bool>("update_mech_on_reuse", false);
    o.tol = prm.get<double>("solver.linsolver.tol");
    o.atol = prm.get<double>("solver.linsolver.atol");
    o.maxIter = static_cast<int>(prm.get<double>("solver.linsolver.max_iter"));
    o.verbosity = static_cast<int>(prm.get<double>("solver.linsolver.verbosity", 0));
    o.restart = prm.get<int>("solver.linsolver.restart", 100);
    o.gmresFallback = prm.get<bool>("solver.linsolver.fallback_gmres_on_breakdown", true);
    o.blasOperator = prm.get<bool>("solver.linsolver.blas_operator", false);
    const auto reuse = o.precond.get<std::string>("mech_solver.reuse", "none");
    if (reuse != "none" && reuse != "rank_update") {
        OPM_THROW(std::runtime_error, "Unknown fracture mech_solver.reuse: " + reuse);
    }
    o.rankUpdate = (reuse == "rank_update") && !o.precond.get<bool>("diag_mech");
    o.mechFactor.inverse = o.precond.get<std::string>("mech_solver.type", "dune_lu") == "lapack_inv";
    o.mechFactor.tol = o.precond.get<double>("mech_solver.rank_update_tol", o.mechFactor.tol);
    o.mechFactor.maxRankFraction
        = o.precond.get<double>("mech_solver.rank_update_max_fraction", o.mechFactor.maxRankFraction);
    o.mechFactor.maxRankMin = o.precond.get<int>("mech_solver.rank_update_max_min", o.mechFactor.maxRankMin);
    o.mechFactor.refactorIterations
        = o.precond.get<int>("mech_solver.rank_update_refactor_iterations", o.mechFactor.refactorIterations);
    return o;
}

bool
FractureCoupledLinearSolver::scaled() const
{
    return scaleMech_ != 1.0 || scalePressure_ != 1.0;
}

// sign=+1 applies D S D, sign=-1 restores S.
void
FractureCoupledLinearSolver::scaleMatrix(const int sign)
{
    auto& S = *S_;
    const double sm = scaleMech_, sp = scalePressure_;
    const double mm = (sign > 0) ? sm * sm : 1.0 / (sm * sm);
    const double mp = (sign > 0) ? sm * sp : 1.0 / (sm * sp);
    const double pp = (sign > 0) ? sp * sp : 1.0 / (sp * sp);
    scale_dense(S[_0][_0], mm);
    scale_sparse(S[_0][_1], mp);
    scale_sparse(S[_1][_0], mp);
    scale_sparse(S[_1][_1], pp);
}

void
FractureCoupledLinearSolver::prepare(SystemMatrix& S, VectorHP& rhs)
{
    S_ = &S;
    scaleMech_ = 1.0;
    scalePressure_ = 1.0;
    // blockmax: D = diag(1/sqrt(max|A|), 1/sqrt(max|M|)) balances the
    // mechanics and flow blocks, which differ by many orders of magnitude.
    if (opt_.scaling == "blockmax") {
        const double a_max = max_abs_entry(S[_0][_0]);
        const double m_max = max_abs_entry(S[_1][_1]);
        if (a_max > 0.0 && m_max > 0.0) {
            scaleMech_ = 1.0 / std::sqrt(a_max);
            scalePressure_ = 1.0 / std::sqrt(m_max);
        }
    }
    if (scaled()) {
        scaleMatrix(+1);
        for (auto& v : rhs[_0]) v[0] *= scaleMech_;
        for (auto& v : rhs[_1]) v[0] *= scalePressure_;
    }
}

std::unique_ptr<FractureCoupledLinearSolver::Solver>
FractureCoupledLinearSolver::makeSolver(const std::string& type,
                                        FractureMechanicsPreconditioner& precond,
                                        const double tol) const
{
    if (type == "bicgstab") {
        return std::make_unique<Dune::BiCGSTABSolver<VectorHP>>(*operator_, precond, tol,
                                                                 opt_.maxIter, opt_.verbosity);
    }
    if (type == "gmres") {
        return std::make_unique<Dune::RestartedGMResSolver<VectorHP>>(*operator_, precond, tol,
                                                                       opt_.restart, opt_.maxIter,
                                                                       opt_.verbosity);
    }
    if (type == "fgmres") {
        return std::make_unique<Dune::RestartedFlexibleGMResSolver<VectorHP>>(*operator_, precond, tol,
                                                                               opt_.restart, opt_.maxIter,
                                                                               opt_.verbosity);
    }
    OPM_THROW(std::runtime_error, "Unknown linear solver type: for fracture");
}

void
FractureCoupledLinearSolver::solve(VectorHP& rhs,
                                   VectorHP& dx,
                                   const Iteration& it,
                                   FractureSolveStats& stats,
                                   const DumpHook& dump)
{
    auto& S = *S_;
    Dune::InverseOperatorResult iores;

    if (it.newSystem) {
        if (opt_.blasOperator)
            operator_ = std::make_unique<DenseBlasSystemOperator>(S);
        else
            operator_ = std::make_unique<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>>(S);
    }
    if (auto* blas = dynamic_cast<DenseBlasSystemOperator*>(operator_.get()))
        blas->refresh();
    using Clock = std::chrono::steady_clock;
    const auto since = [](const Clock::time_point& t0) {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };
    const DenseMechFactor* mech = nullptr;
    if (opt_.rankUpdate && it.mechRows) {
        if (!mechFactor_)
            mechFactor_ = std::make_unique<DenseMechFactor>(opt_.mechFactor);
        const auto t0 = Clock::now();
        const bool refactored = mechFactor_->update(S[_0][_0], scaleMech_ * scaleMech_, *it.mechRows);
        stats.mech_factor_seconds += since(t0);
        if (refactored) {
            ++stats.mech_factorizations;
        } else if (mechFactor_->rank() > 0) {
            ++stats.mech_rank_updates;
            stats.mech_rank_sum += mechFactor_->rank();
        }
        mech = mechFactor_.get();
    }
    const bool new_precond = it.newSystem || it.contactChanged;
    if (new_precond) {
        // free the old factor before building the new one
        precond_.reset();
        precond_.reset(new FractureMechanicsPreconditioner(S, opt_.precond, mech));
        stats.mech_factorizations += precond_->factorizations();
        stats.mech_factor_seconds += precond_->factorSeconds();
    } else {
        const int f0 = precond_->factorizations();
        const double s0 = precond_->factorSeconds();
        precond_->update(S, opt_.updateMechOnReuse);
        stats.mech_factorizations += precond_->factorizations() - f0;
        stats.mech_factor_seconds += precond_->factorSeconds() - s0;
    }
    const double apply0 = precond_->applySeconds();
    const auto krylov0 = Clock::now();
    // atol is an absolute tolerance on the squared residual norm
    const double res = rhs.two_norm2();
    const double lintol = std::max(opt_.tol, opt_.atol / res);
    if (it.verbosity > 1) {
        OpmLog::debug(fmt::format("Nonlinear iteration: {} rhs:  {} {} closed_cells: {} cells m: {} cells p: {}",
                                  it.index, rhs[_0].infinity_norm(), rhs[_1].infinity_norm(),
                                  it.numClosed, rhs[_0].size(), rhs[_1].size()));
    }

    auto& precond = *precond_;
    if (new_precond) {
        solver_ = makeSolver(opt_.solver, precond, lintol);
    }

    if (it.dump) {
        if (scaled()) scaleMatrix(-1);
        dump(S, "manual", nullptr);
        if (scaled()) scaleMatrix(+1);
    }

    // Each attempt restarts from the same rhs; apply() overwrites it.
    const VectorHP rhs_solve_org = rhs;
    auto run_solver = [&](Solver& solver) {
        rhs = rhs_solve_org;
        dx = 0;
        solver.apply(dx, rhs, iores);
        ++stats.linear_solves;
        stats.linear_iterations += iores.iterations;
    };

    {
        OPM_TIMEBLOCK(SolveCoupledSystem);
        try {
            run_solver(*solver_);
        } catch (const std::exception& e1) {
            // Name the preconditioner block that emits non-finite values;
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
            iores.converged = false;
            if (opt_.gmresFallback && opt_.solver == "bicgstab") {
                solver_ = makeSolver("gmres", precond, lintol);
                // a Krylov abort is one more failed rung
                try {
                    run_solver(*solver_);
                } catch (const std::exception&) {
                    iores.converged = false;
                }
            }
        }
    }

    if (!iores.converged) {
        ++stats.linear_solve_failures;
    }

    // Rescue ladder: transient attempts, each on the scaled S;
    // S is rebuilt at the next Newton iteration.
    if (!iores.converged) {
        const auto fixed_stress_prm = [&]() {
            PropertyTree p = opt_.precond;
            p.put("mode_policy", std::string("manual"));
            p.put("mech_first", std::string("true"));
            p.put("fixed_stress", std::string("true"));
            return p;
        };
        // rung 1: flexible GMRES with the existing preconditioner.
        try {
            auto s = makeSolver("fgmres", precond, lintol);
            run_solver(*s);
        } catch (const std::exception&) { /* fall through to next rung */ }

        // rung 2: fixed-stress (Schur) preconditioner + flexible GMRES.
        if (!iores.converged) {
            try {
                FractureMechanicsPreconditioner fsprecond(S, fixed_stress_prm(), mech);
                auto s = makeSolver("fgmres", fsprecond, lintol);
                run_solver(*s);
            } catch (const std::exception&) {}
        }

        // rung 3: drop the flow->mech coupling block C and solve the
        // block-triangular system with a fixed-stress preconditioner.
        if (!iores.converged) {
            try {
                S[_1][_0] = 0;
                FractureMechanicsPreconditioner dcprecond(S, fixed_stress_prm(), mech);
                auto s = makeSolver("fgmres", dcprecond, lintol);
                run_solver(*s);
            } catch (const std::exception&) {}
        }

        if (iores.converged) {
            ++stats.ladder_rescues;
            if (it.verbosity > 0) {
                OpmLog::debug(fmt::format("Fracture coupled linear solve rescued by fallback ladder at nlin iter {} ({} iters)",
                                          it.index, iores.iterations));
            }
        }
    }

    stats.linear_solve_seconds += since(krylov0);
    if (mechFactor_ && opt_.mechFactor.refactorIterations > 0
        && (!iores.converged || iores.iterations > opt_.mechFactor.refactorIterations))
        mechFactor_->invalidate();
    stats.precond_apply_seconds += precond_->applySeconds() - apply0;

    if (!iores.converged) {
        // the dump is the physical system; the throw makes a re-scale unnecessary
        if (scaled()) scaleMatrix(-1);
        const std::string dump_info_filename = dump(S, "linear_failure", &iores);
        stats.converged = false;
        OPM_THROW(std::runtime_error,
                  "Fracture coupled linear solver failed to converge after "
                      + std::to_string(iores.iterations)
                      + " iterations."
                      + (dump_info_filename.empty()
                             ? std::string()
                             : " Snapshot info: " + dump_info_filename));
    }

    // the solver returned y; the Newton update is D y
    if (scaled()) {
        for (auto& v : dx[_0]) v[0] *= scaleMech_;
        for (auto& v : dx[_1]) v[0] *= scalePressure_;
    }
}

} // namespace Opm
