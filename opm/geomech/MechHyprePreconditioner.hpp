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

#ifndef OPM_GEOMECH_MECH_HYPRE_PRECONDITIONER_HPP
#define OPM_GEOMECH_MECH_HYPRE_PRECONDITIONER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreInterface.hpp>

#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvercategory.hh>

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace Opm
{

// BoomerAMG for the mechanics system with the rigid rotations as interpolation vectors
// (interp_vec_variant > 0, needs nodal > 0). modes are in the solver's unknowns.
template <class Matrix, class Vector, class Comm>
class MechHyprePreconditioner : public Dune::PreconditionerWithUpdate<Vector, Vector>
{
public:
    MechHyprePreconditioner(const Matrix& A,
                            const PropertyTree& prm,
                            const std::vector<Vector>& modes,
                            const Comm& comm)
        : A_(A)
        , comm_(comm)
    {
        namespace HI = linalg::HypreInterface;
        HI::initialize(false);
        solver_ = HI::createAMGSolver();
        HI::setSolverParameters(solver_, prm, false);
        HYPRE_BoomerAMGSetNumFunctions(solver_, 3);
        if (prm.get<int>("p_max", -1) >= 0) {
            HYPRE_BoomerAMGSetPMaxElmts(solver_, prm.get<int>("p_max"));
        }
        if (prm.get<int>("keep_same_sign", 0) != 0) {
            HYPRE_BoomerAMGSetKeepSameSign(solver_, 1);
        }

        par_info_ = HI::setupHypreParallelInfo(comm_, A_);
        pattern_ = HI::setupSparsityPattern(A_, par_info_, par_info_.owner_first);
        host_.row_indexes
            = HI::computeRowIndexes(A_, pattern_.ncols, par_info_.local_dune_to_local_hypre, par_info_.owner_first);
        host_.indices.resize(par_info_.N_owned);
        std::iota(host_.indices.begin(), host_.indices.end(), par_info_.dof_offset);
        if (!par_info_.owner_first) {
            host_.continuous_vector_values.resize(par_info_.N_owned);
        }
        A_hypre_ = HI::createMatrix(par_info_.N_owned, par_info_.dof_offset, comm_);
        x_hypre_ = HI::createVector(par_info_.N_owned, par_info_.dof_offset, comm_);
        b_hypre_ = HI::createVector(par_info_.N_owned, par_info_.dof_offset, comm_);

        const int variant = prm.get<int>("interp_vec_variant", 0);
        if (variant > 0) {
            if (prm.get<int>("nodal", 0) <= 0) {
                OPM_THROW(std::invalid_argument, "mech hypre: interp_vec_variant needs nodal > 0");
            }
            for (const auto& m : modes) {
                if (m.size() != A_.N()) {
                    OPM_THROW(std::invalid_argument, "mech hypre: rotation modes do not match the matrix");
                }
                modes_hypre_.push_back(HI::createVector(par_info_.N_owned, par_info_.dof_offset, comm_));
                HI::transferVectorToHypre(m, modes_hypre_.back(), host_, device_, par_info_, false);
                HYPRE_ParVector pv;
                HYPRE_IJVectorGetObject(modes_hypre_.back(), reinterpret_cast<void**>(&pv));
                modes_par_.push_back(pv);
            }
            HYPRE_BoomerAMGSetInterpVectors(solver_, HYPRE_Int(modes_par_.size()), modes_par_.data());
            HYPRE_BoomerAMGSetInterpVecVariant(solver_, variant);
            HYPRE_BoomerAMGSetInterpVecQMax(solver_, prm.get<int>("interp_vec_qmax", 4));
        }
        update();
    }

    ~MechHyprePreconditioner() override
    {
        namespace HI = linalg::HypreInterface;
        HI::destroySolver(solver_);
        for (auto v : modes_hypre_) {
            HI::destroyVector(v);
        }
        HI::destroyVector(x_hypre_);
        HI::destroyVector(b_hypre_);
        HI::destroyMatrix(A_hypre_);
    }

    void update() override
    {
        namespace HI = linalg::HypreInterface;
        HI::updateMatrixValues(A_, A_hypre_, pattern_, host_, device_, false);
        HYPRE_ParCSRMatrix A;
        HYPRE_ParVector x, b;
        HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&A));
        HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&x));
        HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&b));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetup(solver_, A, b, x));
    }

    void pre(Vector& v, Vector&) override { comm_.copyOwnerToAll(v, v); }

    void apply(Vector& v, const Vector& d) override
    {
        namespace HI = linalg::HypreInterface;
        HI::transferVectorToHypre(v, x_hypre_, host_, device_, par_info_, false);
        HI::transferVectorToHypre(d, b_hypre_, host_, device_, par_info_, false);
        HYPRE_ParCSRMatrix A;
        HYPRE_ParVector x, b;
        HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&A));
        HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&x));
        HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&b));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSolve(solver_, A, b, x));
        HI::transferVectorFromHypre(x_hypre_, v, host_, device_, par_info_, false);
        comm_.copyOwnerToAll(v, v);
    }

    void post(Vector&) override {}

    Dune::SolverCategory::Category category() const override
    {
        return std::is_same_v<Comm, Dune::Amg::SequentialInformation> ? Dune::SolverCategory::sequential
                                                                       : Dune::SolverCategory::overlapping;
    }

    bool hasPerfectUpdate() const override { return true; }

private:
    const Matrix& A_;
    const Comm& comm_;
    linalg::HypreInterface::ParallelInfo par_info_;
    linalg::HypreInterface::SparsityPattern pattern_;
    linalg::HypreInterface::HostDataArrays host_;
    linalg::HypreInterface::DeviceDataArrays device_;
    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;
    std::vector<HYPRE_IJVector> modes_hypre_;
    std::vector<HYPRE_ParVector> modes_par_;
};

} // namespace Opm

#endif
