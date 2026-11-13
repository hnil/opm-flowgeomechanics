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
#include <cmath>
#include <sstream>
#include <utility>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <stdexcept>
#include <string>

extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a, const int* lda,
             const int* ipiv, double* b, const int* ldb, int* info);
void dgetri_(const int* n, double* a, const int* lda, const int* ipiv, double* work, const int* lwork, int* info);
void dgemv_(const char* trans, const int* m, const int* n, const double* alpha, const double* a, const int* lda,
            const double* x, const int* incx, const double* beta, double* y, const int* incy);
}

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

FractureMechanicsPreconditioner::FractureMechanicsPreconditioner(const Opm::FractureSystemMatrix& S,
                                                                 Opm::PropertyTree prm)
    : A_(S)
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
    if (mech_solver != "dune_lu" && mech_solver != "lapack_lu" && mech_solver != "lapack_inv") {
        OPM_THROW(std::invalid_argument, "Unknown fracture preconditioner mech_solver.type: " + mech_solver);
    }
    lapack_inv_ = (mech_solver == "lapack_inv");
    lapack_lu_ = (mech_solver == "lapack_lu") || lapack_inv_;
    if (!diag_mech_) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        if(verbosity>0){
        OpmLog::debug("FractureMechanicsPreconditioner: using full mechanics preconditioner, " + mech_solver);
        }
        factorMechanics(S[_0][_0]);
    }
    active_mode_ = selectMode(S);
    if (verbosity > 0 && mode_policy_ == "manual") {
        OpmLog::debug(std::string("FractureMechanicsPreconditioner: configured mode ") + modeName(active_mode_));
    }
    rebuildFlowSolver(S);
}
void FractureMechanicsPreconditioner::update(const Opm::FractureSystemMatrix& S,bool new_lu_mech){

    setDiagvec(A_diag_, S[_0][_0]);

    if (!diag_mech_ && new_lu_mech) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        factorMechanics(S[_0][_0]);
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
    if (active_mode_ == ApplyMode::FixedStress) {
        applyfixed_stress(v, d);
        return;
    }
    if (active_mode_ == ApplyMode::MechFirst){
        applymech_first(v, d);
    } else {
        applymech_last(v, d);
    }
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
    for (int c = 0; c < k; ++c)
        x[closed_[c]][0] = rhs[closed_[c]][0] / closed_diag_[c];
    if (m == 0)
        return;
    for (int r = 0; r < m; ++r) {
        double v = rhs[open_[r]][0];
        const double* a = open_closed_.data() + std::size_t(r) * k;
        for (int c = 0; c < k; ++c)
            v -= a[c] * x[closed_[c]][0];
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
} // namespace Opm
