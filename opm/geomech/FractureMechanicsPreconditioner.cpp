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
#include "FractureMechanicsPreconditioner.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <utility>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

#include <stdexcept>
#include <string>

extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a, const int* lda,
             const int* ipiv, double* b, const int* ldb, int* info);
void dgetri_(const int* n, double* a, const int* lda, const int* ipiv, double* work, const int* lwork, int* info);
void dgemv_(const char* trans, const int* m, const int* n, const double* alpha, const double* a, const int* lda,
            const double* x, const int* incx, const double* beta, double* y, const int* incy);
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha,
            const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c,
            const int* ldc);
}

namespace
{
double
secondsSince(const std::chrono::steady_clock::time_point& t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
} // namespace

namespace Opm
{
const char*
FractureMechanicsPreconditioner::modeName(ApplyMode mode)
{
    switch (mode) {
    case ApplyMode::MechLast:
        return "mech_last";
    case ApplyMode::MechFirst:
        return "mech_first";
    case ApplyMode::FixedStress:
        return "fixed_stress";
    }
    return "unknown";
}

double
FractureMechanicsPreconditioner::estimateCouplingIndicator(const FractureSystemMatrix& S) const
{
    const auto& C = S[_1][_0];
    const auto& I = S[_0][_1];
    const auto& M = S[_1][_1];
    constexpr double eps = 1e-30;

    double indicator = 0.0;
    for (auto rowIt = C.begin(); rowIt != C.end(); ++rowIt) {
        const size_t row = rowIt.index();
        double row_correction = 0.0;
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const size_t col = colIt.index();
            if (col >= A_diag_.size()) {
                continue;
            }
            const auto iit = I[col].find(col);
            if (iit == I[col].end()) {
                continue;
            }
            const double adiag = std::abs(A_diag_[col][0]);
            if (adiag <= eps) {
                continue;
            }
            row_correction += std::abs((*colIt)[0][0] * (*iit)[0][0] / adiag);
        }
        const double mdiag = std::abs(M[row][row][0][0]);
        indicator = std::max(indicator, row_correction / std::max(mdiag, eps));
    }
    return indicator;
}

FractureMechanicsPreconditioner::ApplyMode
FractureMechanicsPreconditioner::selectMode(const FractureSystemMatrix& S)
{
    if (mode_policy_ == "manual") {
        return fixed_stress_ ? ApplyMode::FixedStress
                             : (mech_first_ ? ApplyMode::MechFirst : ApplyMode::MechLast);
    }
    if (mode_policy_ == "coupling_auto") {
        last_coupling_indicator_ = estimateCouplingIndicator(S);
        if (last_coupling_indicator_ >= mode_switch_coupling_threshold_) {
            return ApplyMode::FixedStress;
        }
        return mech_first_ ? ApplyMode::MechFirst : ApplyMode::MechLast;
    }
    throw std::runtime_error("Unknown fracture preconditioner mode_policy: " + mode_policy_);
}

void
FractureMechanicsPreconditioner::rebuildFlowSolver(const Opm::FractureSystemMatrix& S)
{
    const int verbosity = prm_.get<int>("verbosity", 0);

    if (active_mode_ == ApplyMode::FixedStress) {
        updateFixedStressFlowMatrix(S);
        if (!diag_flow_) {
            fixed_stress_flowop_ = std::make_unique<FlowOperatorType>(*fixed_stress_matrix_);
            using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
            fixed_stress_flow_solver_ = std::make_unique<FlowSolverType>(
                *fixed_stress_flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);
        } else {
            setDiagvec(M_diag_, *fixed_stress_matrix_);
        }
    } else if (!diag_flow_) {
        flowop_ = std::make_unique<FlowOperatorType>(S[_1][_1]);
        using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
        flow_solver_ = std::make_unique<FlowSolverType>(
            *flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);
    } else {
        setDiagvec(M_diag_, S[_1][_1]);
    }

    if (verbosity > 0) {
        std::ostringstream msg;
        msg << "FractureMechanicsPreconditioner: active mode " << modeName(active_mode_);
        if (mode_policy_ == "coupling_auto") {
            msg << " coupling indicator=" << last_coupling_indicator_
                << " threshold=" << mode_switch_coupling_threshold_;
        }
        OpmLog::debug(msg.str());
    }
}

std::shared_ptr<const DdmInverse>
makeDdmInverse(const FMatrix& K, const std::uint64_t version)
{
    OPM_TIMEFUNCTION();
    auto inv = std::make_shared<DdmInverse>();
    const int n = K.N();
    inv->version = version;
    inv->K = &K;
    inv->n = n;
    inv->kinv.resize(std::size_t(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            inv->kinv[std::size_t(i) * n + j] = K[i][j];
    if (n == 0)
        return inv;
    // the row-major buffer is K^T for LAPACK; inverting it gives K^{-1} row-major
    std::vector<int> piv(n);
    int info = 0;
    dgetrf_(&n, &n, inv->kinv.data(), &n, piv.data(), &info);
    if (info != 0)
        OPM_THROW(std::runtime_error, "DDM inverse: dgetrf failed, info " + std::to_string(info));
    int lwork = -1;
    double wsize = 0.0;
    dgetri_(&n, inv->kinv.data(), &n, piv.data(), &wsize, &lwork, &info);
    lwork = std::max(n, static_cast<int>(wsize));
    std::vector<double> work(lwork);
    dgetri_(&n, inv->kinv.data(), &n, piv.data(), work.data(), &lwork, &info);
    if (info != 0)
        OPM_THROW(std::runtime_error, "DDM inverse: dgetri failed, info " + std::to_string(info));
    return inv;
}

FractureMechanicsPreconditioner::FractureMechanicsPreconditioner(const Opm::FractureSystemMatrix& S,
                                                                 Opm::PropertyTree prm,
                                                                 const MechBlockSolver* mech,
                                                                 std::shared_ptr<const DdmInverse> ddm)
    : A_(S)
    , mech_(mech)
    , A_diag_(diagvec(S[_0][_0]))
    , M_diag_(diagvec(S[_1][_1]))
    , prm_(prm)
{
    OPM_TIMEFUNCTION();
    diag_mech_ = prm.get<bool>("diag_mech");
    diag_flow_ = prm.get<bool>("diag_flow");
    fixed_stress_ = prm.get<bool>("fixed_stress", false);
    mode_policy_ = prm.get<std::string>("mode_policy", "manual");
    mode_switch_coupling_threshold_ = prm.get<double>("mode_switch_coupling_threshold", 1.0);
    int verbosity = prm.get<bool>("verbosity",0);
    mech_press_coupling_ = prm.get<bool>("mech_press_coupling", true);
    mech_first_ = prm_.get<bool>("mech_first");
    lu_pivoting_ = prm_.get<bool>("lu_pivoting", false);
    const auto mech_solver = prm_.get<std::string>("mech_solver.type", "dune_lu");
    if (mech_solver != "dune_lu" && mech_solver != "lapack_lu" && mech_solver != "lapack_inv"
        && mech_solver != "nearfield_ilu" && mech_solver != "fb_split") {
        OPM_THROW(std::invalid_argument, "Unknown fracture preconditioner mech_solver.type: " + mech_solver);
    }
    lapack_inv_ = (mech_solver == "lapack_inv");
    fb_split_ = (mech_solver == "fb_split");
    lapack_lu_ = (mech_solver == "lapack_lu") || lapack_inv_ || fb_split_;
    split_tol_ = prm_.get<double>("mech_solver.split_tol", 1e-2);
    split_max_fraction_ = prm_.get<double>("mech_solver.split_max_fraction", 0.5);
    ddm_ = std::move(ddm);
    verbosity_ = prm_.get<int>("mech_solver.verbosity", verbosity);
    if (!diag_mech_ && !mech_) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        if(verbosity>0){
        OpmLog::debug("FractureMechanicsPreconditioner: using full mechanics preconditioner, " + mech_solver);
        }
        const auto t0 = std::chrono::steady_clock::now();
        factorMechanics(S[_0][_0]);
        ++factorizations_;
        factor_seconds_ += secondsSince(t0);
    }
    active_mode_ = selectMode(S);
    if (verbosity > 0 && mode_policy_ == "manual") {
        OpmLog::debug(std::string("FractureMechanicsPreconditioner: configured mode ") + modeName(active_mode_));
    }
    rebuildFlowSolver(S);
}
void FractureMechanicsPreconditioner::update(const Opm::FractureSystemMatrix& S, bool new_lu_mech,
                                             std::shared_ptr<const DdmInverse> ddm){

    setDiagvec(A_diag_, S[_0][_0]);
    if (ddm)
        ddm_ = std::move(ddm);

    // fb_split is cheap to rebuild, so it always follows the current matrix
    if (!diag_mech_ && !mech_ && (new_lu_mech || fb_split_)) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        const auto t0 = std::chrono::steady_clock::now();
        factorMechanics(S[_0][_0]);
        ++factorizations_;
        factor_seconds_ += secondsSince(t0);
    }
    if(diag_mech_ && new_lu_mech){
        // if we are using a diagonal mechanics preconditioner, we need to update the diagonal entries
        setDiagvec(A_diag_,S[_0][_0]);
    }
    const ApplyMode new_mode = selectMode(S);
    const bool mode_changed = (new_mode != active_mode_);
    active_mode_ = new_mode;

    if (mode_changed) {
         OPM_TIMEBLOCK(SetupFlowPreconditioner);
         rebuildFlowSolver(S);
    } else if (active_mode_ == ApplyMode::FixedStress) {
         OPM_TIMEBLOCK(SetupFlowPreconditioner);
         updateFixedStressFlowMatrix(S);
         if (!diag_flow_) {
             fixed_stress_flow_solver_->preconditioner().update();
         } else {
             setDiagvec(M_diag_, *fixed_stress_matrix_);
         }
    } else if (!diag_flow_) {
         OPM_TIMEBLOCK(SetupFlowPreconditioner);
         flow_solver_->preconditioner().update();
    }else{
        // if we are using a diagonal flow preconditioner, we need to update the diagonal entries
        setDiagvec(M_diag_,S[_1][_1]);
    }
}

void
FractureMechanicsPreconditioner::apply(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    OPM_TIMEFUNCTION();
    const auto t0 = std::chrono::steady_clock::now();
    if (active_mode_ == ApplyMode::FixedStress) {
        applyfixed_stress(v, d);
    } else if (active_mode_ == ApplyMode::MechFirst){
        applymech_first(v, d);
    } else {
        applymech_last(v, d);
    }
    apply_seconds_ += secondsSince(t0);
}

void
FractureMechanicsPreconditioner::updateFixedStressFlowMatrix(const Opm::FractureSystemMatrix& S)
{
    const auto& C = S[_1][_0];
    const auto& I = S[_0][_1];
    const auto& M = S[_1][_1];

    if (!fixed_stress_matrix_) {
        fixed_stress_matrix_ = std::make_unique<SMatrix>();
        auto& fs = *fixed_stress_matrix_;
        fs.setBuildMode(SMatrix::implicit);
        size_t max_row_nnz = 1;
        for (size_t row = 0; row < M.N(); ++row) {
            std::vector<int> cols(M.M(), 0);
            size_t nnz = 0;
            for (auto colIt = M[row].begin(); colIt != M[row].end(); ++colIt) {
                cols[colIt.index()] = 1;
                ++nnz;
            }
            for (auto colIt = C[row].begin(); colIt != C[row].end(); ++colIt) {
                if (!cols[colIt.index()]) {
                    cols[colIt.index()] = 1;
                    ++nnz;
                }
            }
            max_row_nnz = std::max(max_row_nnz, nnz);
        }
        fs.setImplicitBuildModeParameters(max_row_nnz, 0.2);
        fs.setSize(M.N(), M.M());
        for (size_t row = 0; row < M.N(); ++row) {
            for (auto colIt = M[row].begin(); colIt != M[row].end(); ++colIt) {
                fs.entry(row, colIt.index()) = 0.0;
            }
            for (auto colIt = C[row].begin(); colIt != C[row].end(); ++colIt) {
                fs.entry(row, colIt.index()) = 0.0;
            }
        }
        fs.compress();
    }

    auto& fs = *fixed_stress_matrix_;
    fs = 0;
    // fs's pattern is M union C, so a lock-step copy (identical-sparsity) is not
    // safe when C carries an entry outside M's pattern - copy M by lookup, and
    // fail recoverably (upstream chop) on any inconsistency instead of UB.
    if (C.N() != M.N() || I.N() != C.M() || fs.N() != M.N() || fs.M() != M.M())
        OPM_THROW(std::runtime_error, "fixed-stress flow matrix: inconsistent block dimensions");
    for (auto rowIt = M.begin(); rowIt != M.end(); ++rowIt) {
        auto& fsrow = fs[rowIt.index()];
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            auto fit = fsrow.find(colIt.index());
            if (fit == fsrow.end())
                OPM_THROW(std::runtime_error,
                          "fixed-stress flow matrix: M entry outside stored pattern");
            *fit = *colIt;
        }
    }

    for (auto rowIt = C.begin(); rowIt != C.end(); ++rowIt) {
        const size_t row = rowIt.index();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const size_t col = colIt.index();
            if (col >= I.N() || col >= A_diag_.size())
                OPM_THROW(std::runtime_error,
                          "fixed-stress flow matrix: C column outside mechanics block");
            const auto iit = I[col].find(col);
            if (iit == I[col].end()) {
                continue;
            }
            const double coupling = (*colIt)[0][0];
            const double ij = (*iit)[0][0];
            const double adiag = A_diag_[col][0];
            if (adiag == 0.0) {
                continue;
            }
            auto fit = fs[row].find(col);
            if (fit == fs[row].end())
                OPM_THROW(std::runtime_error,
                          "fixed-stress flow matrix: C entry outside stored pattern");
            (*fit)[0][0] -= coupling * ij / adiag;
        }
    }
}

void
FractureMechanicsPreconditioner::solveMechanics(Opm::Vector& x, const Opm::Vector& rhs) const
{
    if (diag_mech_) {
        for (size_t i = 0; i != A_diag_.size(); ++i) {
            x[i] = rhs[i] / A_diag_[i];
        }
        return;
    }

    OPM_TIMEBLOCK(ApplyMechPreconditioner);
    if (mech_) {
        mech_->solve(x, rhs);
        return;
    }
    if (split_active_) {
        splitSolve(x, rhs);
        return;
    }
    if (lapack_lu_) {
        openBlockSolve(x, rhs);
        return;
    }
    auto tmp = rhs;
    backSolve(x, tmp);
}

void
FractureMechanicsPreconditioner::solveFlow(Opm::Vector& x, const Opm::Vector& rhs)
{
    if (diag_flow_) {
        for (size_t i = 0; i != M_diag_.size(); ++i) {
            x[i] = rhs[i] / M_diag_[i];
        }
        return;
    }

    OPM_TIMEBLOCK(ApplyFlowPreconditioner);
    Dune::InverseOperatorResult res;
    auto tmp = rhs;
    if (active_mode_ == ApplyMode::FixedStress) {
        fixed_stress_flow_solver_->apply(x, tmp, res);
    } else {
        flow_solver_->apply(x, tmp, res);
    }
}

void
FractureMechanicsPreconditioner::applymech_first(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    // FractureSystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //               {C, M}}; // flow system
    OPM_TIMEFUNCTION_LOCAL();
    solveMechanics(v[_0], d[_0]);
    auto rhs_flow = d[_1];
    if (mech_press_coupling_) {
        A_[_1][_0].mmv(v[_0], rhs_flow); // rhs_flow -= A_[_1][_0] * v[_0]
    }
    solveFlow(v[_1], rhs_flow);
};

void
FractureMechanicsPreconditioner::applyfixed_stress(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    OPM_TIMEFUNCTION_LOCAL();
    solveMechanics(v[_0], d[_0]);

    auto rhs_flow = d[_1];
    if (mech_press_coupling_) {
        A_[_1][_0].mmv(v[_0], rhs_flow);
    }
    solveFlow(v[_1], rhs_flow);
}

void
FractureMechanicsPreconditioner::applymech_last(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    // FractureSystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //               {C, M}}; // flow system
    OPM_TIMEFUNCTION_LOCAL();
    solveFlow(v[_1], d[_1]);

    auto rhs_mech = d[_0];
    if (mech_press_coupling_) {
        A_[_0][_1].mmv(v[_1], rhs_mech); // rhs_mech -= A_[_0][_1] * v[_1]
    }
    solveMechanics(v[_0], rhs_mech);



};

void
FractureMechanicsPreconditioner::factorMechanics(const FMatrix& A)
{
    split_active_ = fb_split_ && factorSplit(A);
    if (split_active_)
        return;
    if (lapack_lu_) {
        factorOpenBlock(A);
        return;
    }
    luM_ = A;
    if (lu_pivoting_) {
        MyDenseMatrix<double>::luDecomp(luM_, lu_pivot_);
    } else {
        lu_pivot_.clear();
        MyDenseMatrix<double>::luDecomp(luM_);
    }
    checkLuFactor();
}

void
FractureMechanicsPreconditioner::factorOpenBlock(const FMatrix& A)
{
    const int n = A.N();
    open_.clear();
    closed_.clear();
    for (int i = 0; i < n; ++i) {
        bool diagonal = true;
        for (int j = 0; j < n && diagonal; ++j) {
            diagonal = (j == i || A[i][j] == 0.0);
        }
        (diagonal ? closed_ : open_).push_back(i);
    }
    const int m = open_.size();
    const int k = closed_.size();

    closed_diag_.resize(k);
    for (int c = 0; c < k; ++c) {
        closed_diag_[c] = A[closed_[c]][closed_[c]];
        if (!std::isfinite(closed_diag_[c]) || closed_diag_[c] == 0.0)
            OPM_THROW(std::runtime_error,
                      "mechanics LAPACK LU: bad diagonal row " + std::to_string(closed_[c]));
    }
    open_lu_.resize(std::size_t(m) * m);
    open_closed_.resize(std::size_t(m) * k);
    for (int r = 0; r < m; ++r) {
        for (int q = 0; q < m; ++q)
            open_lu_[std::size_t(r) * m + q] = A[open_[r]][open_[q]];
        for (int c = 0; c < k; ++c)
            open_closed_[std::size_t(r) * k + c] = A[open_[r]][closed_[c]];
    }
    open_pivot_.resize(m);
    open_rhs_.resize(m);
    if (m == 0)
        return;
    int info = 0;
    dgetrf_(&m, &m, open_lu_.data(), &m, open_pivot_.data(), &info);
    bool finite = (info == 0);
    for (int r = 0; r < m && finite; ++r)
        finite = std::isfinite(open_lu_[std::size_t(r) * m + r]);
    if (!finite)
        OPM_THROW(std::runtime_error, "mechanics LAPACK LU: dgetrf failed, info " + std::to_string(info));
    if (lapack_inv_) {
        int lwork = -1;
        double wsize = 0.0;
        dgetri_(&m, open_lu_.data(), &m, open_pivot_.data(), &wsize, &lwork, &info);
        lwork = std::max(m, static_cast<int>(wsize));
        std::vector<double> work(lwork);
        dgetri_(&m, open_lu_.data(), &m, open_pivot_.data(), work.data(), &lwork, &info);
        if (info != 0)
            OPM_THROW(std::runtime_error, "mechanics LAPACK inverse: dgetri failed, info " + std::to_string(info));
    }
}

void
FractureMechanicsPreconditioner::openBlockSolve(Vector& x, const Vector& rhs) const
{
    const int m = open_.size();
    const int k = closed_.size();
    closed_x_.resize(k);
    for (int c = 0; c < k; ++c) {
        closed_x_[c] = rhs[closed_[c]][0] / closed_diag_[c];
        x[closed_[c]][0] = closed_x_[c];
    }
    if (m == 0)
        return;
    const double* xc = closed_x_.data();
    // four rows at a time; each row keeps its own summation order
    int r = 0;
    for (; r + 4 <= m; r += 4) {
        double v0 = rhs[open_[r]][0], v1 = rhs[open_[r + 1]][0];
        double v2 = rhs[open_[r + 2]][0], v3 = rhs[open_[r + 3]][0];
        const double* a0 = open_closed_.data() + std::size_t(r) * k;
        const double* a1 = a0 + k;
        const double* a2 = a1 + k;
        const double* a3 = a2 + k;
        for (int c = 0; c < k; ++c) {
            v0 = std::fma(-a0[c], xc[c], v0);
            v1 = std::fma(-a1[c], xc[c], v1);
            v2 = std::fma(-a2[c], xc[c], v2);
            v3 = std::fma(-a3[c], xc[c], v3);
        }
        open_rhs_[r] = v0;
        open_rhs_[r + 1] = v1;
        open_rhs_[r + 2] = v2;
        open_rhs_[r + 3] = v3;
    }
    for (; r < m; ++r) {
        double v = rhs[open_[r]][0];
        const double* a = open_closed_.data() + std::size_t(r) * k;
        for (int c = 0; c < k; ++c)
            v = std::fma(-a[c], xc[c], v);
        open_rhs_[r] = v;
    }
    // the row-major copy is the transpose in LAPACK's column-major layout
    const char trans = 'T';
    const int one = 1;
    if (lapack_inv_) {
        const double alpha = 1.0, beta = 0.0;
        open_x_.resize(m);
        dgemv_(&trans, &m, &m, &alpha, open_lu_.data(), &m, open_rhs_.data(), &one, &beta, open_x_.data(), &one);
        for (int r = 0; r < m; ++r)
            x[open_[r]][0] = open_x_[r];
        return;
    }
    int info = 0;
    dgetrs_(&trans, &m, &one, open_lu_.data(), &m, open_pivot_.data(), open_rhs_.data(), &m, &info);
    for (int r = 0; r < m; ++r)
        x[open_[r]][0] = open_rhs_[r];
}

bool
FractureMechanicsPreconditioner::factorSplit(const FMatrix& A)
{
    OPM_TIMEFUNCTION();
    const int n = A.N();
    if (!ddm_ || ddm_->n != n || ddm_->K == nullptr || n < 2) {
        if (verbosity_ > 0)
            OpmLog::debug("fb_split: no DDM inverse for this matrix, full factorisation");
        return false;
    }
    const FMatrix& K = *ddm_->K;
    split_rows_.clear();
    split_scale_.assign(n, 1.0);
    int n_shape = 0, n_diag = 0, n_zero = 0;
    double max_off = 0.0, max_g = 0.0;
    for (int i = 0; i < n; ++i) {
        int j0 = -1;
        double kmax = 0.0;
        for (int j = 0; j < n; ++j)
            if (j != i && std::abs(K[i][j]) > kmax) {
                kmax = std::abs(K[i][j]);
                j0 = j;
            }
        const double s = (j0 >= 0) ? -A[i][j0] / K[i][j0] : 0.0;
        bool reference = (s != 0.0) && std::isfinite(s);
        if (s == 0.0) {
            bool zero = true;
            for (int j = 0; j < n && zero; ++j)
                zero = (j == i || A[i][j] == 0.0);
            n_zero += zero;
        }
        // every off-diagonal must be -s K_ij, else the row is taken exactly
        double off = 0.0;
        for (int j = 0; j < n && reference; ++j)
            if (j != i)
                off = std::max(off, std::abs(A[i][j] + s * K[i][j]) / (std::abs(s) * kmax));
        if (reference && verbosity_ > 1) max_off = std::max(max_off, off);
        if (reference && off > 1e-10) { reference = false; ++n_shape; }
        if (reference) {
            const double g = std::abs(A[i][i] + s * K[i][i]) / std::abs(s * K[i][i]);
            if (verbosity_ > 1) max_g = std::max(max_g, g);
            if (g > split_tol_) { reference = false; ++n_diag; }
        }
        if (reference)
            split_scale_[i] = 1.0 / s;
        else
            split_rows_.push_back(i);
    }
    const int k = split_rows_.size();
    if (verbosity_ > 0)
        OpmLog::debug(fmt::format("fb_split: n {} exact rows {} ({:.3f}) diagonal-only {} shape {} diag {} max_off {:.2e} max_g {:.2e}",
                                  n, k, double(k) / n, n_zero, n_shape, n_diag, max_off, max_g));
    // the fallback LU factors only the m rows that are not diagonal-only: ~2/3 m^3 against ~2 n k^2
    const double m = n - n_zero;
    if (k > split_max_fraction_ * n || 3.0 * double(n) * k * k > m * m * m)
        return false;

    const double* kinv = ddm_->kinv.data();
    split_Ar_.resize(std::size_t(k) * n);
    split_KR_.resize(std::size_t(n) * k);
    for (int c = 0; c < k; ++c)
        for (int j = 0; j < n; ++j)
            split_Ar_[std::size_t(c) * n + j] = A[split_rows_[c]][j];
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < k; ++c)
            split_KR_[std::size_t(r) * k + c] = kinv[std::size_t(r) * n + split_rows_[c]];
    split_G_.resize(std::size_t(k) * k);
    split_piv_.resize(k);
    split_z_.resize(n);
    split_y_.resize(n);
    split_t_.resize(k);
    if (k == 0)
        return true;
    // G = A(R,:) K^{-1}(:,R), computed as G^T = KR^T Ar^T in column-major
    const double one = 1.0, zero = 0.0;
    const char nt = 'N';
    dgemm_(&nt, &nt, &k, &k, &n, &one, split_KR_.data(), &k, split_Ar_.data(), &n, &zero, split_G_.data(), &k);
    int info = 0;
    dgetrf_(&k, &k, split_G_.data(), &k, split_piv_.data(), &info);
    if (info != 0) {
        if (verbosity_ > 0)
            OpmLog::debug(fmt::format("fb_split: capacitance LU failed (info {}), full factorisation", info));
        return false;
    }
    return true;
}

// x = y - K^{-1}(:,R) G^{-1} (A(R,:) y - z_R) with z = scale .* rhs and y = -K^{-1} z
void
FractureMechanicsPreconditioner::splitSolve(Vector& x, const Vector& rhs) const
{
    const int n = split_scale_.size();
    const int k = split_rows_.size();
    for (int i = 0; i < n; ++i)
        split_z_[i] = rhs[i][0] * split_scale_[i];
    const char tr = 'T';
    const int inc = 1;
    const double one = 1.0, mone = -1.0, zero = 0.0;
    dgemv_(&tr, &n, &n, &mone, ddm_->kinv.data(), &n, split_z_.data(), &inc, &zero, split_y_.data(), &inc);
    if (k > 0) {
        dgemv_(&tr, &n, &k, &one, split_Ar_.data(), &n, split_y_.data(), &inc, &zero, split_t_.data(), &inc);
        for (int c = 0; c < k; ++c)
            split_t_[c] -= split_z_[split_rows_[c]];
        int info = 0;
        dgetrs_(&tr, &k, &inc, split_G_.data(), &k, split_piv_.data(), split_t_.data(), &k, &info);
        dgemv_(&tr, &k, &n, &mone, split_KR_.data(), &k, split_t_.data(), &inc, &one, split_y_.data(), &inc);
    }
    for (int i = 0; i < n; ++i)
        x[i][0] = split_y_[i];
}

void
FractureMechanicsPreconditioner::checkLuFactor() const
{
    // A non-finite or zero pivot poisons every later apply and surfaces only as
    // an opaque Krylov "defect=nan" - fail here, recoverably and by name.
    for (size_t i = 0; i < luM_.rows(); ++i) {
        const double piv = luM_[i][i];
        if (!std::isfinite(piv) || piv == 0.0)
            OPM_THROW(std::runtime_error,
                      "mechanics dense LU: bad pivot " + std::to_string(piv)
                          + " at row " + std::to_string(i));
    }
}

void
FractureMechanicsPreconditioner::backSolve(Opm::Vector& x, const Opm::Vector& rhs) const
{
    // apply the row permutation recorded by the pivoted decomposition
    auto prhs = rhs;
    for (size_t i = 0; i < lu_pivot_.size(); ++i)
        if (lu_pivot_[i] != i)
            std::swap(prhs[i], prhs[lu_pivot_[i]]);
    for (int i = 0; i < int(luM_.rows()); i++) {
        x[i] = prhs[i];
        for (int j = 0; j < i; j++) {
            x[i] -= luM_[i][j] * x[j];
        }
    }

    for (int i = luM_.rows() - 1; i >= 0; i--) {
        for (size_t j = i + 1; j < luM_.rows(); j++) {
            x[i] -= luM_[i][j] * x[j];
        }
        x[i] /= luM_[i][i];
    }
}
// ---------------------------------------------------------------------------
// DenseMechFactor

void
DenseMechFactor::factor(const FMatrix& B, const double scale, const MechRows& rows)
{
    const int n = B.N();
    n_ = n;
    generation_ = rows.generation;
    s0_ = scale;
    ratio_ = 1.0;
    alpha0_ = *rows.alpha;
    beta0_ = *rows.beta;
    f_.resize(std::size_t(n) * n);
    for (int i = 0; i < n; ++i)
        std::copy(B[i].begin(), B[i].end(), f_.begin() + std::size_t(i) * n);
    piv_.resize(n);
    int info = 0;
    if (n > 0)
        dgetrf_(&n, &n, f_.data(), &n, piv_.data(), &info);
    bool finite = (info == 0);
    for (int r = 0; r < n && finite; ++r)
        finite = std::isfinite(f_[std::size_t(r) * n + r]);
    if (!finite)
        OPM_THROW(std::runtime_error, "mechanics dense factor: dgetrf failed, info " + std::to_string(info));
    if (opt_.inverse && n > 0) {
        int lwork = -1;
        double wsize = 0.0;
        dgetri_(&n, f_.data(), &n, piv_.data(), &wsize, &lwork, &info);
        lwork = std::max(n, static_cast<int>(wsize));
        std::vector<double> work(lwork);
        dgetri_(&n, f_.data(), &n, piv_.data(), work.data(), &lwork, &info);
        if (info != 0)
            OPM_THROW(std::runtime_error, "mechanics dense inverse: dgetri failed, info " + std::to_string(info));
    }
    maxRank_ = std::min(n, std::max(opt_.maxRankMin, static_cast<int>(opt_.maxRankFraction * n)));
    slotOf_.assign(n, -1);
    slotRow_.clear();
    z_.resize(std::size_t(n) * maxRank_);
    ar_.resize(std::size_t(n) * maxRank_);
    g_.assign(std::size_t(maxRank_) * maxRank_, 0.0);
    cap_.clear();
    y_.resize(n);
    t_.resize(maxRank_);
}

bool
DenseMechFactor::addRows(const std::vector<int>& rows_in)
{
    const int n = n_;
    const int r0 = rank();
    const int nnew = rows_in.size();
    const int r1 = r0 + nnew;
    for (int q = 0; q < nnew; ++q) {
        const int j = rows_in[q];
        const int k = r0 + q;
        slotOf_[j] = k;
        slotRow_.push_back(j);
        double* z = z_.data() + std::size_t(k) * n;
        if (opt_.inverse) {
            // f_ is B0^{-1} row-major: gather column j
            for (int i = 0; i < n; ++i)
                z[i] = f_[std::size_t(i) * n + j];
        } else {
            std::fill(z, z + n, 0.0);
            z[j] = 1.0;
        }
    }
    if (!opt_.inverse) {
        const char trans = 'T';
        int info = 0;
        dgetrs_(&trans, &n, &nnew, f_.data(), &n, piv_.data(), z_.data() + std::size_t(r0) * n, &n, &info);
        if (info != 0)
            return false;
    }
    for (int k = r0; k < r1; ++k)
        for (int i = 0; i < n; ++i)
            if (!std::isfinite(z_[std::size_t(k) * n + i]))
                return false;
    return true;
}

bool
DenseMechFactor::factorCapacitance(const MechRows& rows)
{
    const int n = n_;
    const int r = rank();
    cap_.clear();
    if (r == 0)
        return true;
    const auto& alpha = *rows.alpha;
    const auto& beta = *rows.beta;
    va_.resize(r);
    vb_.resize(r);
    for (int k = 0; k < r; ++k) {
        const int j = slotRow_[k];
        va_[k] = s0_ * (alpha[j] - alpha0_[j]);
        vb_[k] = s0_ * (beta[j] - beta0_[j]);
    }
    // C = I + V^T Z, V^T(k,:) = va_k A(row_k,:) + vb_k e_{row_k}^T
    cap_.resize(std::size_t(r) * r);
    const int ld = maxRank_;
    for (int l = 0; l < r; ++l) {
        const double* z = z_.data() + std::size_t(l) * n;
        for (int k = 0; k < r; ++k)
            cap_[std::size_t(l) * r + k] = (k == l ? 1.0 : 0.0) + va_[k] * g_[std::size_t(l) * ld + k]
                + vb_[k] * z[slotRow_[k]];
    }
    capPiv_.resize(r);
    int info = 0;
    dgetrf_(&r, &r, cap_.data(), &r, capPiv_.data(), &info);
    if (info != 0)
        return false;
    for (int k = 0; k < r; ++k)
        if (!std::isfinite(cap_[std::size_t(k) * r + k]))
            return false;
    return true;
}

bool
DenseMechFactor::update(const FMatrix& B, const double scale, const MechRows& rows)
{
    const int n = B.N();
    if (!rows.ddm || !rows.alpha || !rows.beta || int(rows.ddm->N()) != n
        || int(rows.alpha->size()) != n || int(rows.beta->size()) != n)
        OPM_THROW(std::runtime_error, "mechanics dense factor: row coefficients do not match the block");
    if (f_.empty() || n != n_ || rows.generation != generation_ || stale_) {
        stale_ = false;
        factor(B, scale, rows);
        return true;
    }
    ratio_ = s0_ / scale;
    const auto& A = *rows.ddm;
    const auto& alpha = *rows.alpha;
    const auto& beta = *rows.beta;
    std::vector<std::pair<double, int>> moved; // relative row change, row
    for (int i = 0; i < n; ++i) {
        if (slotOf_[i] >= 0)
            continue;
        const double da = alpha[i] - alpha0_[i];
        const double db = beta[i] - beta0_[i];
        if (da == 0.0 && db == 0.0)
            continue;
        const double c = std::abs(A[i][i]);
        const double change = std::abs(da) * c + std::abs(db);
        const double size = std::abs(alpha[i]) * c + std::abs(beta[i]);
        if (change > opt_.tol * size)
            moved.emplace_back(size > 0.0 ? change / size : 1e300, i);
    }
    const int room = maxRank_ - rank();
    if (int(moved.size()) > room) {
        if (opt_.refactorIterations <= 0) {
            factor(B, scale, rows);
            return true;
        }
        // keep the largest changes; the rest only weaken the preconditioner
        std::nth_element(moved.begin(), moved.begin() + room, moved.end(),
                         [](const auto& a, const auto& b) { return a.first > b.first; });
        moved.resize(room);
    }
    std::vector<int> added;
    added.reserve(moved.size());
    for (const auto& m : moved)
        added.push_back(m.second);
    std::sort(added.begin(), added.end());
    if (!added.empty()) {
        const int r0 = rank();
        if (!addRows(added)) {
            factor(B, scale, rows);
            return true;
        }
        const int r1 = rank();
        for (int k = r0; k < r1; ++k) {
            const auto& row = A[slotRow_[k]];
            std::copy(row.begin(), row.end(), ar_.begin() + std::size_t(k) * n);
        }
        // new entries of G = Ar^T Z: columns r0..r1 for all rows, then rows r0..r1 of the old columns
        const char tr = 'T', nt = 'N';
        const double one = 1.0, zero = 0.0;
        const int ld = maxRank_;
        const int nnew = r1 - r0;
        dgemm_(&tr, &nt, &r1, &nnew, &n, &one, ar_.data(), &n, z_.data() + std::size_t(r0) * n, &n, &zero,
               g_.data() + std::size_t(r0) * ld, &ld);
        if (r0 > 0)
            dgemm_(&tr, &nt, &nnew, &r0, &n, &one, ar_.data() + std::size_t(r0) * n, &n, z_.data(), &n, &zero,
                   g_.data() + r0, &ld);
    }
    if (!factorCapacitance(rows)) {
        factor(B, scale, rows);
        return true;
    }
    return false;
}

void
DenseMechFactor::solve(Vector& x, const Vector& rhs) const
{
    const int n = n_;
    if (n == 0)
        return;
    const char trans = 'T', notrans = 'N';
    const int one = 1;
    if (opt_.inverse) {
        const double a = 1.0, b = 0.0;
        dgemv_(&trans, &n, &n, &a, f_.data(), &n, &rhs[0][0], &one, &b, y_.data(), &one);
    } else {
        for (int i = 0; i < n; ++i)
            y_[i] = rhs[i][0];
        int info = 0;
        dgetrs_(&trans, &n, &one, f_.data(), &n, piv_.data(), y_.data(), &n, &info);
    }
    const int r = rank();
    if (r > 0) {
        const double a = 1.0, b = 0.0, m1 = -1.0;
        dgemv_(&trans, &n, &r, &a, ar_.data(), &n, y_.data(), &one, &b, t_.data(), &one);
        for (int k = 0; k < r; ++k)
            t_[k] = va_[k] * t_[k] + vb_[k] * y_[slotRow_[k]];
        int info = 0;
        dgetrs_(&notrans, &r, &one, cap_.data(), &r, capPiv_.data(), t_.data(), &r, &info);
        dgemv_(&notrans, &n, &r, &m1, z_.data(), &n, t_.data(), &one, &a, y_.data(), &one);
    }
    for (int i = 0; i < n; ++i)
        x[i][0] = ratio_ * y_[i];
}

void
NearFieldMechILU::buildPattern(const std::vector<double>& xyz)
{
    const int n = n_;
    const auto d2 = [&](int i, int j) {
        double s = 0.0;
        for (int a = 0; a < 3; ++a) {
            const double d = xyz[3 * i + a] - xyz[3 * j + a];
            s += d * d;
        }
        return s;
    };
    // typical spacing: median nearest-neighbour distance over a sample of cells
    std::vector<double> nn;
    const int stride = std::max(1, n / 256);
    for (int i = 0; i < n; i += stride) {
        double best = std::numeric_limits<double>::max();
        for (int j = 0; j < n; ++j)
            if (j != i)
                best = std::min(best, d2(i, j));
        nn.push_back(best);
    }
    std::nth_element(nn.begin(), nn.begin() + nn.size() / 2, nn.end());
    const double h = nn.empty() ? 1.0 : std::sqrt(nn[nn.size() / 2]);
    const double R = radius_ * h;
    const double R2 = R * R;

    const auto bin = [&](int i, int a) { return static_cast<std::int64_t>(std::floor(xyz[3 * i + a] / R)); };
    const auto key = [](std::int64_t x, std::int64_t y, std::int64_t z) {
        return (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
    };
    std::unordered_map<std::int64_t, std::vector<int>> bins;
    for (int i = 0; i < n; ++i)
        bins[key(bin(i, 0), bin(i, 1), bin(i, 2))].push_back(i);

    ptr_.assign(1, 0);
    col_.clear();
    diag_.resize(n);
    std::vector<int> row;
    for (int i = 0; i < n; ++i) {
        row.clear();
        const std::int64_t bx = bin(i, 0), by = bin(i, 1), bz = bin(i, 2);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = bins.find(key(bx + dx, by + dy, bz + dz));
                    if (it == bins.end())
                        continue;
                    for (int j : it->second)
                        if (j == i || d2(i, j) <= R2)
                            row.push_back(j);
                }
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end()); // hash collisions
        for (int j : row) {
            if (j == i)
                diag_[i] = col_.size();
            col_.push_back(j);
        }
        ptr_.push_back(col_.size());
    }
    val_.resize(col_.size());
    havePattern_ = true;
}

void
NearFieldMechILU::update(const FMatrix& B, const MechRows& rows)
{
    const int n = B.N();
    if (!rows.centres || int(rows.centres->size()) != 3 * n)
        OPM_THROW(std::runtime_error, "near-field mechanics ILU: cell centres do not match the block");
    if (!havePattern_ || n != n_ || rows.generation != generation_) {
        n_ = n;
        generation_ = rows.generation;
        buildPattern(*rows.centres);
    }
    for (int i = 0; i < n; ++i)
        for (int k = ptr_[i]; k < ptr_[i + 1]; ++k)
            val_[k] = B[i][col_[k]];
    // ILU(0), IKJ order on the fixed pattern
    std::vector<int> pos(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int k = ptr_[i]; k < ptr_[i + 1]; ++k)
            pos[col_[k]] = k;
        for (int k = ptr_[i]; k < diag_[i]; ++k) {
            const int c = col_[k];
            const double l = (val_[k] /= val_[diag_[c]]);
            for (int q = diag_[c] + 1; q < ptr_[c + 1]; ++q)
                if (pos[col_[q]] >= 0)
                    val_[pos[col_[q]]] -= l * val_[q];
        }
        for (int k = ptr_[i]; k < ptr_[i + 1]; ++k)
            pos[col_[k]] = -1;
        double& p = val_[diag_[i]];
        if (!std::isfinite(p) || std::abs(p) <= 1e-12 * std::abs(B[i][i]))
            p = B[i][i]; // keep the factor usable on a breakdown
    }
}

void
NearFieldMechILU::solve(Vector& x, const Vector& rhs) const
{
    const int n = n_;
    x.resize(n);
    for (int i = 0; i < n; ++i) {
        double v = rhs[i][0];
        for (int k = ptr_[i]; k < diag_[i]; ++k)
            v -= val_[k] * x[col_[k]][0];
        x[i][0] = v;
    }
    for (int i = n - 1; i >= 0; --i) {
        double v = x[i][0];
        for (int k = diag_[i] + 1; k < ptr_[i + 1]; ++k)
            v -= val_[k] * x[col_[k]][0];
        x[i][0] = v / val_[diag_[i]];
    }
}
} // namespace Opm
