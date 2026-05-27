#include "config.h"
#include "FractureMechanicsPreconditioner.hpp"
//#include <StrumpackSparseSolver.hpp>

namespace Opm
{
FractureMechanicsPreconditioner::FractureMechanicsPreconditioner(const Opm::SystemMatrix& S,
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
    int verbosity = prm.get<bool>("verbosity",0);
    mech_press_coupling_ = prm.get<bool>("mech_press_coupling", true);
    mech_first_ = prm_.get<bool>("mech_first");
    if (!diag_mech_) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        if(verbosity>0){
        std::cout << "FractureMechanicsPreconditioner: using full mechanics preconditioner" << std::endl;
        }
        // luM_ = S[_0][_0];
        luM_ = S[_0][_0];
        MyDenseMatrix<double>::luDecomp(luM_);
    }
        if (fixed_stress_) {
            updateFixedStressFlowMatrix(S);
            if(verbosity>0){
                std::cout << "FractureMechanicsPreconditioner: using fixed-stress flow preconditioner" << std::endl;
            }
            if (!diag_flow_) {
                fixed_stress_flowop_ = std::make_unique<Dune::MatrixAdapter<Opm::SMatrix, Opm::Vector, Opm::Vector>>(*fixed_stress_matrix_);
                using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
                fixed_stress_flow_solver_ = std::make_unique<FlowSolverType>(
                            *fixed_stress_flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);
            } else {
                setDiagvec(M_diag_, *fixed_stress_matrix_);
            }
        } else if (!diag_flow_) {
      if(verbosity>0){
        std::cout << "FractureMechanicsPreconditioner: using full flow preconditioner" << std::endl;
      }
      flowop_ = std::make_unique<Dune::MatrixAdapter<Opm::SMatrix, Opm::Vector, Opm::Vector>>(S[_1][_1]);
       using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
       flow_solver_ = std::make_unique<FlowSolverType>(
             *flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);
    }
}
void FractureMechanicsPreconditioner::update(const Opm::SystemMatrix& S,bool new_lu_mech){

    if (!diag_mech_ && new_lu_mech) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        //luM_ = S[_0][_0];
        //copy to dense matrix for lu factorization
        for(size_t i = 0; i != S[_0][_0].M(); ++i){
            for(size_t j = 0; j != S[_0][_0].N(); ++j){
                luM_[i][j] = S[_0][_0][i][j];
            }
        }
        MyDenseMatrix<double>::luDecomp(luM_);
    }
    if(diag_mech_ && new_lu_mech){
        // if we are using a diagonal mechanics preconditioner, we need to update the diagonal entries
        setDiagvec(A_diag_,S[_0][_0]);
    }
    if (fixed_stress_) {
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
        //  using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
        //  flow_solver_ = std::make_unique<FlowSolverType>(
        //     *flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);    
        //flow_solver_->update(*flowop_);
    }else{
        // if we are using a diagonal flow preconditioner, we need to update the diagonal entries
        setDiagvec(M_diag_,S[_1][_1]);
    }
}

void
FractureMechanicsPreconditioner::apply(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    //OPM_TIMEFUNCTION_LOCAL();
    OPM_TIMEFUNCTION();
    if (fixed_stress_) {
        applyfixed_stress(v, d);
        return;
    }
    if (mech_first_){
        applymech_first(v, d);
    } else {
        applymech_last(v, d);
    }
}

void
FractureMechanicsPreconditioner::updateFixedStressFlowMatrix(const Opm::SystemMatrix& S)
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
    copyMatrixValuesWithSameSparsity(fs, M);

    for (auto rowIt = C.begin(); rowIt != C.end(); ++rowIt) {
        const size_t row = rowIt.index();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const size_t col = colIt.index();
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
            fs[row][col][0][0] -= coupling * ij / adiag;
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
    if (fixed_stress_) {
        fixed_stress_flow_solver_->apply(x, tmp, res);
    } else {
        flow_solver_->apply(x, tmp, res);
    }
}

void
FractureMechanicsPreconditioner::applymech_first(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    // SystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //               {C, M}}; // flow system
    OPM_TIMEFUNCTION_LOCAL();
    solveMechanics(v[_0], d[_0]);
    auto rhs_flow = d[_1];
    if (mech_press_coupling_) {
        A_[_1][_0].mmv(v[_0], rhs_flow); // -1.0); // rhs_flow -= A_[_1][_0] * v[_0]
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
    // SystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //               {C, M}}; // flow system
    OPM_TIMEFUNCTION_LOCAL();
    solveFlow(v[_1], d[_1]);

    auto rhs_mech = d[_0];
    if (mech_press_coupling_) {
        A_[_0][_1].mmv(v[_1], rhs_mech); // -1.0); // rhs_flow -= A_[_1][_0] * v[_0]
    }
    solveMechanics(v[_0], rhs_mech);
    

    
};

void
FractureMechanicsPreconditioner::backSolve(Opm::Vector& x, const Opm::Vector& rhs) const
{
    // Vector& rhs = x;
    // rhs = rhs_in;
    for (int i = 0; i < int(luM_.rows()); i++) {
        x[i] = rhs[i];
        for (int j = 0; j < i; j++) {
            x[i] -= luM_[i][j] * x[j];
        }
        // rhs[i] = rhs[i]/luM_[i][i];
    }

    for (int i = luM_.rows() - 1; i >= 0; i--) {
        for (size_t j = i + 1; j < luM_.rows(); j++) {
            x[i] -= luM_[i][j] * x[j];
        }
        x[i] /= luM_[i][i];
    }
}
} // namespace Opm
