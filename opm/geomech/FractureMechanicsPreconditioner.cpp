#include "config.h"
#include "FractureMechanicsPreconditioner.hpp"
#include <StrumpackSparseSolver.hpp>

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
    mech_press_coupling_ = prm.get<bool>("mech_press_coupling", true);
    if (!diag_mech_) {
        OPM_TIMEBLOCK(SetupLuFactorization);
        std::cout << "FractureMechanicsPreconditioner: using full mechanics preconditioner" << std::endl;
        luM_ = S[_0][_0];
        MyDenseMatrix<double>::luDecomp(luM_);
    }
    if (!diag_flow_) {
        std::cout << "FractureMechanicsPreconditioner: using full flow preconditioner" << std::endl;
        flowop_
            = std::make_unique<Dune::MatrixAdapter<Opm::SMatrix, Opm::Vector, Opm::Vector>>(S[_1][_1]);
        using FlowSolverType = Dune::FlexibleSolver<FlowOperatorType>;
        flow_solver_ = std::make_unique<FlowSolverType>(
            *flowop_, prm_.get_child("flow_solver"), std::function<Vector()>(), 1);
    }
}

void
FractureMechanicsPreconditioner::apply(Opm::VectorHP& v, const Opm::VectorHP& d)
{
    // SystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //               {C, M}}; // flow system
    OPM_TIMEFUNCTION_LOCAL();
    if (diag_mech_) {
        for (size_t i = 0; i != A_diag_.size(); ++i) {
            v[_0][i] = d[_0][i] / A_diag_[i];
        }
    } else {
        // solve full mechanics system
        if (true) {
            auto tmp = d[_0];
            // A_[_0][_0].solve(v[_0],tmp);
            this->backSolve(v[_0], tmp);
        } else {
            // DenseMatrix<double> A(n, n);
            // for(int i=0; i< n; ++i){
            //   for(int j=0; i< n; ++j){
            //     A(i,j) = A[_0][_0][i][j];
            //   }
            // }
            //  structured::StructuredOptions<double> options;
            //  structured::ClusterTree tree(n);
            //  tree.refine(options.leaf_size());
            //  auto H = structured::construct_from_dense(A, options);
        }
        // auto diff = tmp;
        // diff -= v[_0];
        // auto err = diff.two_norm();///v[_0].two_norm();
        // assert(err<1e-8);
    }
    auto rhs_flow = d[_1];
    if (mech_press_coupling_) {
        A_[_1][_0].mmv(v[_0], rhs_flow); // -1.0); // rhs_flow -= A_[_1][_0] * v[_0]
    }

    if (diag_flow_) {
        for (size_t i = 0; i != M_diag_.size(); ++i) {
            v[_1][i] = d[_1][i] / M_diag_[i];
        }
    } else {
        Dune::InverseOperatorResult res;
        flow_solver_->apply(v[_1], rhs_flow, res);
        // throw std::runtime_error("FractureMechanicsPreconditioner: full flow preconditioner not
        // implemented");
    }
};

void
FractureMechanicsPreconditioner::backSolve(Opm::Vector& x, const Opm::Vector& rhs)
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
