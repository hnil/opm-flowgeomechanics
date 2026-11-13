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
                                             bool stab_on_stress,
                                             bool diagonal_scaling)
{
    stab_on_stress_ = stab_on_stress;
    diagonal_scaling_ = diagonal_scaling;
    std::stringstream os;
    os << "Setting stability choise to " << stability_choice_int
       << (diagonal_scaling_ ? " with diagonal scaling" : "");
    OpmLog::info(os.str());
    stability_choice_ = vem::check_stability_choice(stability_choice_int);
}

template<class GridType>
void VemDiscretisation<GridType>::fixNodes(const std::vector<NodeConstraint>& constraints)
{
    // constraints are sorted and unique by node, so the dofs are too
    std::vector<int> fixed_dof_ixs;
    std::vector<double> fixed_dof_values;
    frame_of_node_.assign(grid_.leafGridView().size(3), -1);
    frames_.clear();
    for (const auto& c : constraints) {
        if (!c.axisAligned()) {
            frame_of_node_[c.node] = int(frames_.size());
            frames_.push_back(c.frame);
        }
        for (int m = 0; m < 3; ++m) {
            if (c.fixed[m]) {
                fixed_dof_ixs.push_back(3 * c.node + m);
                fixed_dof_values.push_back(c.value[m]);
            }
        }
    }
    dirichlet_ = {int(fixed_dof_ixs.size()), fixed_dof_ixs, fixed_dof_values};
}

// Rows (and columns) of rotated nodes become F A (F^T): u = F^T u~ with F the
// node's frame. Entries of other nodes are passed through unchanged.
template<class GridType>
void VemDiscretisation<GridType>::rotateRows(std::vector<std::tuple<int, int, double>>& entries,
                                              bool columns_too) const
{
    if (frames_.empty()) {
        return;
    }
    std::vector<std::tuple<int, int, double>> out;
    out.reserve(entries.size());
    const auto expand = [this](int dof, double val, auto&& emit) {
        const int f = frame_of_node_[dof / 3];
        if (f < 0) {
            emit(dof, val);
            return;
        }
        const int a = dof % 3;
        for (int i = 0; i < 3; ++i) {
            emit(3 * (dof / 3) + i, frames_[f][i][a] * val);
        }
    };
    for (const auto& [r, c, v] : entries) {
        expand(r, v, [&](int ri, double vr) {
            if (columns_too) {
                expand(c, vr, [&](int ci, double vc) { out.emplace_back(ri, ci, vc); });
            } else {
                out.emplace_back(ri, c, vr);
            }
        });
    }
    entries.swap(out);
}

// B -> F_a B F_b^T for the 3 x 3 block coupling nodes a and b.
template<class GridType>
void VemDiscretisation<GridType>::rotateBlock(int node_a, int node_b, double (&blk)[3][3]) const
{
    if (frames_.empty()) {
        return;
    }
    const int fa = frame_of_node_[node_a];
    const int fb = frame_of_node_[node_b];
    if (fa < 0 && fb < 0) {
        return;
    }
    double t[3][3];
    for (int i = 0; i != 3; ++i)
        for (int l = 0; l != 3; ++l) {
            double v = 0.0;
            for (int k = 0; k != 3; ++k)
                v += (fa < 0 ? double(i == k) : frames_[fa][i][k]) * blk[k][l];
            t[i][l] = v;
        }
    for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 3; ++j) {
            double v = 0.0;
            for (int l = 0; l != 3; ++l)
                v += t[i][l] * (fb < 0 ? double(j == l) : frames_[fb][j][l]);
            blk[i][j] = v;
        }
}

// Node-level pattern from the cells' face corners, three dofs per node.
template<class GridType>
void VemDiscretisation<GridType>::buildMatrixPattern(int num_dofs)
{
    OPM_TIMEBLOCK(buildMatrixPattern);
    const int num_nodes = num_dofs / 3;
    std::vector<std::vector<int>> adj(num_nodes);
    std::vector<int> cell_nodes;
    int fc = 0;
    int cf = 0;
    for (int c = 0; c < num_cells_; ++c) {
        cell_nodes.clear();
        for (int f = 0; f < num_cell_faces_[c]; ++f, ++cf) {
            cell_nodes.insert(cell_nodes.end(), &face_corners_[fc], &face_corners_[fc] + num_face_corners_[cf]);
            fc += num_face_corners_[cf];
        }
        std::sort(cell_nodes.begin(), cell_nodes.end());
        cell_nodes.erase(std::unique(cell_nodes.begin(), cell_nodes.end()), cell_nodes.end());
        for (const int a : cell_nodes) {
            adj[a].insert(adj[a].end(), cell_nodes.begin(), cell_nodes.end());
        }
    }
    std::size_t nnz = 0;
    for (auto& row : adj) {
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
        row.shrink_to_fit();
        nnz += 9 * row.size();
    }
    auto& A = *matrix_;
    A.setBuildMode(Matrix::random);
    A.setSize(num_dofs, num_dofs, nnz);
    for (int a = 0; a < num_nodes; ++a)
        for (int i = 0; i != 3; ++i)
            A.setrowsize(3 * a + i, 3 * adj[a].size());
    A.endrowsizes();
    for (int a = 0; a < num_nodes; ++a)
        for (int i = 0; i != 3; ++i)
            for (const int b : adj[a])
                for (int j = 0; j != 3; ++j)
                    A.addindex(3 * a + i, 3 * b + j);
    A.endindices();
}

template<class GridType>
template<class V>
void VemDiscretisation<GridType>::toLocal(V& u) const
{
    if (u.size() == 0) {
        return; // no start vector yet: zero in any frame
    }
    for (std::size_t k = 0; k < frame_of_node_.size(); ++k) {
        if (const int f = frame_of_node_[k]; f >= 0) {
            const double g[3] = {double(u[3 * k]), double(u[3 * k + 1]), double(u[3 * k + 2])};
            for (int i = 0; i < 3; ++i) {
                const auto& d = frames_[f][i];
                u[3 * k + i] = d[0] * g[0] + d[1] * g[1] + d[2] * g[2];
            }
        }
    }
}

template<class GridType>
template<class V>
void VemDiscretisation<GridType>::toGlobal(V& u) const
{
    if (u.size() == 0) {
        return;
    }
    for (std::size_t k = 0; k < frame_of_node_.size(); ++k) {
        if (const int f = frame_of_node_[k]; f >= 0) {
            const double l[3] = {double(u[3 * k]), double(u[3 * k + 1]), double(u[3 * k + 2])};
            const auto& F = frames_[f];
            for (int a = 0; a < 3; ++a) {
                u[3 * k + a] = F[0][a] * l[0] + F[1][a] * l[1] + F[2][a] * l[2];
            }
        }
    }
}

template<class GridType>
const std::vector<Vector>& VemDiscretisation<GridType>::rotationModes()
{
    const std::size_t nn = coords_.size() / 3;
    // Common centre on all ranks, so each mode is one rigid rotation globally.
    std::array<double, 4> sum = {0.0, 0.0, 0.0, double(nn)};
    for (std::size_t k = 0; k < nn; ++k) {
        for (int d = 0; d < 3; ++d) {
            sum[d] += coords_[3 * k + d];
        }
    }
    grid_.comm().sum(sum.data(), 4);
    const double c[3] = {sum[0] / sum[3], sum[1] / sum[3], sum[2] / sum[3]};
    rotation_modes_.assign(3, Vector(3 * nn));
    for (std::size_t k = 0; k < nn; ++k) {
        const double x = coords_[3 * k] - c[0];
        const double y = coords_[3 * k + 1] - c[1];
        const double z = coords_[3 * k + 2] - c[2];
        const double r[3][3] = {{-y, x, 0.0}, {0.0, -z, y}, {z, 0.0, -x}};
        for (int m = 0; m < 3; ++m) {
            for (int d = 0; d < 3; ++d) {
                rotation_modes_[m][3 * k + d] = r[m][d];
            }
        }
    }
    for (auto& mode : rotation_modes_) {
        toLocal(mode);
        for (const int dof : std::get<1>(dirichlet_)) {
            mode[dof] = 0.0;
        }
        if (scale_.size() == mode.size()) {
            for (std::size_t i = 0; i < mode.size(); ++i) {
                mode[i] /= scale_[i];
            }
        }
    }
    return rotation_modes_;
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
    const int num_dofs = 3 * grid_.leafGridView().size(3);
    matrix_ = std::make_unique<Matrix>(); // a built BCRSMatrix cannot be rebuilt in place
    buildMatrixPattern(num_dofs);
    auto& A = *matrix_;
    A = 0.0;
    rhs_force_.clear();
    const auto scatter = [this, &A](const int* nodes, const int n, const double* K) {
        const int m = 3 * n;
        for (int a = 0; a != n; ++a) {
            for (int b = 0; b != n; ++b) {
                double blk[3][3];
                for (int i = 0; i != 3; ++i)
                    for (int j = 0; j != 3; ++j)
                        blk[i][j] = K[(3 * a + i) * m + 3 * b + j];
                rotateBlock(nodes[a], nodes[b], blk);
                for (int i = 0; i != 3; ++i)
                    for (int j = 0; j != 3; ++j)
                        A[3 * nodes[a] + i][3 * nodes[b] + j] += blk[i][j];
            }
        }
    };
    vem::assemble_mech_system_3D_dune(grid_, &coords_[0], num_cells_, &num_cell_faces_[0],
                                      &num_face_corners_[0], &face_corners_[0], &ymodule_[0],
                                      &pratio_[0], &body_force_[0],
                                      /*num_neumann_faces*/ 0, nullptr, nullptr, scatter,
                                      rhs_force_, stability_choice_, corner_mask_.data(),
                                      hasOperators() ? nullptr : &pending_ops_);
    assert(int(rhs_force_.size()) == num_dofs);
    toLocal(rhs_force_);
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
        if (!pending_ops_.div.empty()) {
            divmat.swap(pending_ops_.div);
        } else {
            vem::potential_gradient_force_3D(&coords_[0], num_cells_, &num_cell_faces_[0],
                                             &num_face_corners_[0], &face_corners_[0],
                                             &pressure_tmp[0], rhs_tmp, divmat, true,
                                             corner_mask_.data());
        }
        rotateRows(divmat, /*columns_too*/ false);
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
        divmat_.setBuildMode(Matrix::random);
        divmat_.setSize(num_dofs, num_cells_);
        makeDuneMatrixCompressed(divmatdof, divmat_);
    }
    std::vector<double> dispall(num_dofs);
    std::vector<std::array<double, 6>> stresstmp(num_cells_);
    for (const bool is_stress : {true, false}) {
        OPM_TIMEBLOCK(setUpStressStrainMatrix);
        std::vector<std::tuple<int, int, double>> entries;
        auto& pending = is_stress ? pending_ops_.stress : pending_ops_.strain;
        if (!pending.empty()) {
            entries.swap(pending);
        } else {
            vem::compute_stress_3D(&coords_[0], num_cells_, &num_cell_faces_[0],
                                   &num_face_corners_[0], &face_corners_[0], &ymodule_[0],
                                   &pratio_[0], dispall, stresstmp, entries, true, is_stress,
                                   stability_choice_, stab_on_stress_, corner_mask_.data());
        }
        Matrix& mat = is_stress ? stressmat_ : strainmat_;
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
    for (std::size_t i = 0; i < matrix_->N(); ++i) {
        diag_sum += (*matrix_)[i][i];
    }
    const double d = matrix_->N() > 0 ? diag_sum / matrix_->N() : 1.0;
    std::vector<double> g(matrix_->N(), 0.0);
    std::vector<char> fixed(matrix_->N(), 0);
    for (std::size_t k = 0; k < fixed_dof_ixs.size(); ++k) {
        fixed[fixed_dof_ixs[k]] = 1;
        g[fixed_dof_ixs[k]] = fixed_dof_values[k];
    }
    // the pattern is symmetric, so row j lists the rows that hold column j
    for (const int j : fixed_dof_ixs) {
        for (auto col = (*matrix_)[j].begin(); col != (*matrix_)[j].end(); ++col) {
            const auto i = col.index();
            if (!fixed[i]) {
                rhs_force_[i] -= (*matrix_)[i][j] * g[j];
                (*matrix_)[i][j] = 0.0;
            }
        }
    }
    for (const int j : fixed_dof_ixs) {
        (*matrix_)[j] = 0.0;
        (*matrix_)[j][j] = d;
        rhs_force_[j] = d * g[j];
    }
}

template<class GridType>
void VemDiscretisation<GridType>::stress(const Vector& dispall, VoigtField& out) const
{
    OPM_TIMEBLOCK(calculateStressPrecomputed);
    Vector s(6 * num_cells_);
    stressmat_.mv(dispall, s);
    out.resize(num_cells_);
    detail::vemToVoigt(out, s);
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
    auto& mat = *matrix_;
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
void MechLinearSolver<GridType>::setup(const PropertyTree& prm_in, Matrix& A,
                                       const std::vector<double>& coords,
                                       const std::vector<Vector>& modes)
{
    OPM_TIMEBLOCK(setupLinearSolver);
    tsolver_.reset();
    A_ = &A;
    addMechPreconditioners(coords, modes);
    // hypre with interp_vec_variant needs the rotations, which only mech_hypre passes;
    // hypre uses them with nodal coarsening only.
    PropertyTree prm = prm_in;
    for (const std::string path : {"preconditioner", "preconditioner.inner"}) {
        if (prm.get<std::string>(path + ".type", "") != "hypre"
            || prm.get<int>(path + ".interp_vec_variant", 0) <= 0) {
            continue;
        }
        if (prm.get<int>(path + ".nodal", 0) > 0) {
            prm.put(path + ".type", std::string("mech_hypre"));
        } else if (grid_.comm().rank() == 0) {
            OpmLog::warning("mech solver: interp_vec_variant ignored, it needs nodal > 0");
        }
    }
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
