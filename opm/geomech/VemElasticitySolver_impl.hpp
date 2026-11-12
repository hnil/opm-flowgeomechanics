/*
  Copyright 2023 Equinor ASA.

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

#ifndef OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP
#define OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <tuple>
#include <vector>

namespace Opm {
namespace Elasticity {

namespace detail {
// VEM shear order xy, yz, zx to Voigt yz, zx, xy.
inline void vemToVoigt(VoigtField& voigt, const Vector& vem)
{
    for (std::size_t i = 0; i < voigt.size(); ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            voigt[i][k] = vem[i * 6 + k];
        }
        voigt[i][3] = vem[i * 6 + 4];
        voigt[i][5] = vem[i * 6 + 3];
        voigt[i][4] = vem[i * 6 + 5];
    }
}
} // namespace detail

template<class GridType>
void VemDiscretisation<GridType>::setOptions(int stability_choice_int,
                                             bool vem_source,
                                             bool stab_on_stress,
                                             bool vem_stress,
                                             bool diagonal_scaling)
{
    vem_source_ = vem_source;
    stab_on_stress_ = stab_on_stress;
    vem_stress_ = vem_stress;
    diagonal_scaling_ = diagonal_scaling;
    std::stringstream os;
    os << "Setting stability choise to " << stability_choice_int
       << (diagonal_scaling_ ? " with diagonal scaling" : "");
    OpmLog::info(os.str());
    stability_choice_ = vem::check_stability_choice(stability_choice_int);
}

template<class GridType>
void VemDiscretisation<GridType>::fixNodes(const std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes)
{
    // bc_nodes are sorted and unique (nodesAtBoundary/nodesAtOuterBoundary), so the dofs are too
    std::vector<int> fixed_dof_ixs;
    std::vector<double> fixed_dof_values;
    for (const auto& [node, value] : bc_nodes) {
        for (int m = 0; m < 3; ++m) {
            if (value.fixeddir[m]) {
                fixed_dof_ixs.push_back(3 * node + m);
                fixed_dof_values.push_back(value.disp[m]);
            }
        }
    }
    dirichlet_ = {int(fixed_dof_ixs.size()), fixed_dof_ixs, fixed_dof_values};
}

template<class GridType>
void VemDiscretisation<GridType>::setBodyForce(double gravity, const std::vector<double>& density)
{
    const int num_cells = grid_.leafGridView().size(0);
    body_force_.assign(3 * num_cells, 0.0);
    for (int i = 0; i < num_cells; ++i) {
        body_force_[3 * i + 2] = density[i] * gravity;
    }
}

template<class GridType>
void VemDiscretisation<GridType>::buildGeometry()
{
    coords_.clear();
    num_cell_faces_.clear();
    num_face_corners_.clear();
    face_corners_.clear();
    corner_mask_.clear();
    vem::getGridVectors(grid_, coords_, num_cell_faces_, num_face_corners_, face_corners_,
                        &corner_mask_);
    num_cells_ = grid_.leafGridView().size(0);
}

template<class GridType>
void VemDiscretisation<GridType>::assembleMatrix()
{
    OPM_TIMEBLOCK(assembleVEMSystem);
    if (coords_.empty()) {
        buildGeometry();
    }
    std::vector<std::tuple<int, int, double>> A_entries;
    rhs_force_.clear();
    vem::assemble_mech_system_3D_dune(grid_, &coords_[0], num_cells_, &num_cell_faces_[0],
                                      &num_face_corners_[0], &face_corners_[0], &ymodule_[0],
                                      &pratio_[0], &body_force_[0],
                                      /*num_neumann_faces*/ 0, nullptr, nullptr, A_entries,
                                      rhs_force_, stability_choice_, corner_mask_.data());
    const int num_dofs = int(rhs_force_.size());
    matrix_ = Matrix();
    matrix_.setBuildMode(Matrix::random);
    matrix_.setSize(num_dofs, num_dofs);
    makeDuneMatrixCompressed(A_entries, matrix_);
    applyDirichlet();
    scale_.resize(0);
    if (diagonal_scaling_) {
        applyDiagonalScaling();
    }
}

template<class GridType>
void VemDiscretisation<GridType>::buildOperators()
{
    if (coords_.empty()) {
        buildGeometry();
    }
    const auto& fixed_dof_ixs = std::get<1>(dirichlet_);
    const int num_dofs = grid_.leafGridView().size(3) * 3;
    {
        OPM_TIMEBLOCK(setUpExtraStructuresForVEM);
        std::vector<double> rhs_tmp(num_cells_, 0);
        std::vector<double> pressure_tmp(num_cells_, 0);
        std::vector<std::tuple<int, int, double>> divmat;
        if (vem_source_) {
            vem::potential_gradient_force_3D(&coords_[0], num_cells_, &num_cell_faces_[0],
                                             &num_face_corners_[0], &face_corners_[0],
                                             &pressure_tmp[0], rhs_tmp, divmat, true,
                                             corner_mask_.data());
        } else {
            vem::potential_gradient_force_3D_dune(grid_, &coords_[0], num_cells_,
                                                  &num_cell_faces_[0], &num_face_corners_[0],
                                                  &face_corners_[0], &pressure_tmp[0], rhs_tmp,
                                                  divmat, true, corner_mask_.data());
        }
        std::vector<char> fixed(num_dofs, 0);
        for (const int dof : fixed_dof_ixs) {
            fixed[dof] = 1;
        }
        // fixed dofs get no load
        std::vector<std::tuple<int, int, double>> divmatdof;
        for (const auto& [row, col, val] : divmat) {
            if (!fixed[row]) {
                divmatdof.emplace_back(row, col, val);
            }
        }
        divmat_ = Matrix();
        divmat_.setBuildMode(Matrix::random);
        divmat_.setSize(num_dofs, num_cells_);
        makeDuneMatrixCompressed(divmatdof, divmat_);
    }
    std::vector<double> dispall(num_dofs);
    std::vector<std::array<double, 6>> stresstmp(num_cells_);
    for (const bool is_stress : {true, false}) {
        OPM_TIMEBLOCK(setUpStressStrainMatrix);
        std::vector<std::tuple<int, int, double>> entries;
        vem::compute_stress_3D(&coords_[0], num_cells_, &num_cell_faces_[0],
                               &num_face_corners_[0], &face_corners_[0], &ymodule_[0],
                               &pratio_[0], dispall, stresstmp, entries, true, is_stress,
                               stability_choice_, stab_on_stress_, corner_mask_.data());
        Matrix& mat = is_stress ? stressmat_ : strainmat_;
        mat = Matrix();
        mat.setBuildMode(Matrix::random);
        mat.setSize(num_cells_ * 6, dispall.size());
        makeDuneMatrixCompressed(entries, mat);
    }
}

template<class GridType>
void VemDiscretisation<GridType>::rhs(const Vector& load, Vector& b) const
{
    OPM_TIMEBLOCK(updateRhsWithGrad);
    b.resize(rhs_force_.size());
    divmat_.mv(load, b);
    for (std::size_t i = 0; i < rhs_force_.size(); ++i) {
        b[i] += rhs_force_[i];
    }
}

// Symmetric elimination: fixed columns move into the rhs; a fixed row becomes
// d u_j = d g_j with d the mean diagonal, so the scaling of A is unchanged.
template<class GridType>
void VemDiscretisation<GridType>::applyDirichlet()
{
    const auto& fixed_dof_ixs = std::get<1>(dirichlet_);
    const auto& fixed_dof_values = std::get<2>(dirichlet_);
    double diag_sum = 0.0;
    for (std::size_t i = 0; i < matrix_.N(); ++i) {
        diag_sum += matrix_[i][i];
    }
    const double d = matrix_.N() > 0 ? diag_sum / matrix_.N() : 1.0;
    std::vector<double> g(matrix_.N(), 0.0);
    std::vector<char> fixed(matrix_.N(), 0);
    for (std::size_t k = 0; k < fixed_dof_ixs.size(); ++k) {
        fixed[fixed_dof_ixs[k]] = 1;
        g[fixed_dof_ixs[k]] = fixed_dof_values[k];
    }
    // the pattern is symmetric, so row j lists the rows that hold column j
    for (const int j : fixed_dof_ixs) {
        for (auto col = matrix_[j].begin(); col != matrix_[j].end(); ++col) {
            const auto i = col.index();
            if (!fixed[i]) {
                rhs_force_[i] -= matrix_[i][j] * g[j];
                matrix_[i][j] = 0.0;
            }
        }
    }
    for (const int j : fixed_dof_ixs) {
        matrix_[j] = 0.0;
        matrix_[j][j] = d;
        rhs_force_[j] = d * g[j];
    }
}

template<class GridType>
void VemDiscretisation<GridType>::stress(const Vector& dispall, VoigtField& out) const
{
    OPM_TIMEBLOCK(calculateStressPrecomputed);
    if (vem_stress_) {
        Vector s(6 * num_cells_);
        stressmat_.mv(dispall, s);
        out.resize(num_cells_);
        detail::vemToVoigt(out, s);
    } else {
        out = vem::computeStressFem(grid_, dispall, ymodule_, pratio_);
    }
}

template<class GridType>
void VemDiscretisation<GridType>::strain(const Vector& dispall, VoigtField& out) const
{
    OPM_TIMEBLOCK(calculateStrainPrecomputed);
    Vector s(6 * num_cells_);
    strainmat_.mv(dispall, s);
    out.resize(num_cells_);
    detail::vemToVoigt(out, s);
}

// A -> S A S with S = diag(1/sqrt(A_ii)): equilibrates the diagonal of high-aspect
// cells; the solver recovers x = S y.
template<class GridType>
void VemDiscretisation<GridType>::applyDiagonalScaling()
{
    auto& mat = matrix_;
    scale_.resize(mat.N());
    for (auto row = mat.begin(), rowEnd = mat.end(); row != rowEnd; ++row) {
        const auto i = row.index();
        double d = 0.0;
        for (auto col = row->begin(), colEnd = row->end(); col != colEnd; ++col) {
            if (col.index() == i) {
                d = (*col)[0][0];
            }
        }
        scale_[i] = (d > 0.0) ? 1.0 / std::sqrt(d) : 1.0;
    }
    for (auto row = mat.begin(), rowEnd = mat.end(); row != rowEnd; ++row) {
        const auto i = row.index();
        for (auto col = row->begin(), colEnd = row->end(); col != colEnd; ++col) {
            (*col) *= scale_[i][0] * scale_[col.index()][0];
        }
    }
}

template<class GridType>
void VemDiscretisation<GridType>::makeDuneMatrixCompressed(
    const std::vector<std::tuple<int, int, double>>& entries, Matrix& mat)
{
    OPM_TIMEBLOCK(makeDuneMatrixCompressed);
    {
        OPM_TIMEBLOCK(buildMatrixPattern);
        std::vector<std::pair<int, int>> pattern;
        pattern.reserve(entries.size());
        for (const auto& e : entries) {
            pattern.emplace_back(std::get<0>(e), std::get<1>(e));
        }
        std::sort(pattern.begin(), pattern.end());
        pattern.erase(std::unique(pattern.begin(), pattern.end()), pattern.end());
        std::vector<int> row_size(mat.N(), 0);
        for (const auto& [i, j] : pattern) {
            ++row_size[i];
        }
        for (std::size_t i = 0; i < mat.N(); ++i) {
            mat.setrowsize(i, row_size[i]);
        }
        mat.endrowsizes();
        for (const auto& [i, j] : pattern) {
            mat.addindex(i, j);
        }
        mat.endindices();
    }
    OPM_TIMEBLOCK(setMatrixValues);
    mat = 0;
    for (const auto& [i, j, val] : entries) {
        mat[i][j] += val;
    }
}

template<class GridType>
void MechLinearSolver<GridType>::setCopyRowsToZero(Matrix& A)
{
    assert(A.N() == A.M());
    assert(A.N() == dofParallelIndexSet_->size());
    for (const auto& index : *dofParallelIndexSet_) {
        const auto lind = index.local().local();
        if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::copy) {
            A[lind] = 0.0;
            A[lind][lind] = 1.0;
        }
    }
}

template<class GridType>
void MechLinearSolver<GridType>::setup(const PropertyTree& prm, Matrix& A,
                                       const std::vector<double>& coords)
{
    OPM_TIMEBLOCK(setupLinearSolver);
    tsolver_.reset();
    A_ = &A;
    addMechLinePlanePreconditioner(coords);
    // Opt-in: stop the warm-started solve when the residual of the full
    // load is below tol, not when the correction's own residual is.
    tol_relative_to_load_ = prm.get<bool>("tol_relative_to_load", true);
    warm_start_ = prm.get<bool>("warm_start", true);
    linear_tol_ = prm.get<double>("tol", 1e-5);
#if HAVE_MPI
    comm_.reset(new MechCommunication(grid_.comm(), Dune::SolverCategory::overlapping, /*free*/ false));
#endif
    const std::function<Vector()> weightsCalculator;
    const std::size_t pressureIndex = 0;
    if (comm_->communicator().size() > 1) {
        vertexParallelIndexSet_.reset(new MechCommunication::ParallelIndexSet(
            Opm::makeEntityEntityCommunication<3>(grid_)));
        assert(grid_.size(3) == vertexParallelIndexSet_->size());
        dofParallelIndexSet_.reset(new MechCommunication::ParallelIndexSet(
            Opm::entityToDofIndexSet(*vertexParallelIndexSet_, 3)));
        comm_->indexSet() = *dofParallelIndexSet_;
        assert(grid_.size(3) * 3 == dofParallelIndexSet_->size());
        comm_->remoteIndices().rebuild<false>();
        Opm::Elasticity::validateMechanicalPartition(grid_, *comm_);
        setCopyRowsToZero(A);
        comm_->communicator().barrier();
        auto pop = std::make_unique<ParOperatorType>(A, *comm_);
        using FlexibleSolverType = Dune::FlexibleSolver<ParOperatorType>;
        tsolver_ = std::make_unique<FlexibleSolverType>(*pop, *comm_, prm, weightsCalculator,
                                                        pressureIndex);
        pop_ = std::move(pop);
    } else {
        auto sop = std::make_unique<SeqOperatorType>(A);
        using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
        tsolver_ = std::make_unique<FlexibleSolverType>(*sop, prm, weightsCalculator,
                                                        pressureIndex);
        sop_ = std::move(sop);
    }
}

template<class GridType>
MechSolveResult MechLinearSolver<GridType>::solve(const Vector& rhs, const Vector& scale, Vector& x, bool cold)
{
    const bool scaled = scale.size() > 0;
    MechSolveResult result;
    // Copy rows hold an identity and a partial load, so only owner rows count.
    auto globalNorm2 = [this](const Vector& v) {
        double n2 = v.two_norm2();
#if HAVE_MPI
        if (grid_.comm().size() > 1)
            comm_->dot(v, v, n2);
#endif
        return n2;
    };
    try {
        Dune::InverseOperatorResult r;
        // with diagonal scaling the operator is S A S: solve (S A S) y = S b, x = S y
        Vector rhs_eff = rhs;
        if (scaled) {
            for (std::size_t i = 0; i < rhs_eff.size(); ++i)
                rhs_eff[i] *= scale[i][0];
        }
        const double b_norm = std::sqrt(globalNorm2(rhs_eff));
        const double target = linear_tol_ * b_norm;
        Vector y(rhs.size());
        if (x.size() != rhs.size() || !warm_start_ || cold) {
            y = 0;
            tsolver_->apply(y, rhs_eff, r);
        } else {
            y = x;
            if (scaled)
                for (std::size_t i = 0; i < y.size(); ++i)
                    y[i] = x[i][0] / scale[i][0];
            auto rhs_tmp = rhs_eff;
            A_->mmv(y, rhs_tmp);
            auto dy = y;
            dy = 0;
            if (tol_relative_to_load_) {
                // ||b - A(y+dy)|| <= tol ||b||: the correction only needs the reduction
                // that gets the full residual there, and none at all if it is already.
                const double have = std::sqrt(globalNorm2(rhs_tmp));
                if (have <= target) {
                    r.iterations = 0;
                    r.converged = true;
                    r.reduction = 1.0;
                } else {
                    tsolver_->apply(dy, rhs_tmp, std::min(0.5, target / have), r);
                }
            } else {
                tsolver_->apply(dy, rhs_tmp, r);
            }
            y += dy;
        }
        result.iterations = r.iterations;
        result.converged = r.converged;
        result.reduction = r.reduction;
        if (tol_relative_to_load_ && result.converged && b_norm > 0.0) {
            // GMRES stops on the preconditioned residual; accept only on the true one
            auto res = rhs_eff;
            A_->mmv(y, res);
            double res_norm = std::sqrt(globalNorm2(res));
            if (res_norm > target) {
                auto dy = y;
                dy = 0;
                Dune::InverseOperatorResult r2;
                auto rhs_tmp = res;
                tsolver_->apply(dy, rhs_tmp, std::min(0.5, target / res_norm), r2);
                y += dy;
                result.iterations += r2.iterations;
                res = rhs_eff;
                A_->mmv(y, res);
                res_norm = std::sqrt(globalNorm2(res));
            }
            if (res_norm > target) {
                result.converged = false;
                OpmLog::warning(fmt::format("Mechanical solve: true residual {:.3e} above "
                                            "tolerance {:.3e}", res_norm, target));
            }
        }
        x = y;
        if (scaled)
            for (std::size_t i = 0; i < x.size(); ++i)
                x[i] = scale[i][0] * y[i][0];
    } catch (Dune::ISTLError& e) {
        result.converged = false;
        result.failed = true;
        OpmLog::warning(fmt::format("exception thrown {}", e.what()));
    }
    return result;
}

} // namespace Elasticity
} // namespace Opm

#endif
