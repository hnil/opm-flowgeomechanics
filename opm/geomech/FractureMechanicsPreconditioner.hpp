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

#pragma once

#include <cstdint>
#include <vector>
#include <memory>
#include <string>
#include <opm/common/TimingMacros.hpp>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <opm/geomech/DiagonalScalar.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

namespace Opm{
  // ----------------------------------------------------------------------------
template <class DstMatrix, class SrcMatrix>
void
copyMatrixValuesWithSameSparsity(DstMatrix& dst, const SrcMatrix& src)
// ----------------------------------------------------------------------------
{
    assert(dst.N() == src.N());
    assert(dst.M() == src.M());

    auto srcRowIt = src.begin();
    auto dstRowIt = dst.begin();
    for (; srcRowIt != src.end() && dstRowIt != dst.end(); ++srcRowIt, ++dstRowIt) {
        assert(srcRowIt.index() == dstRowIt.index());

        auto srcColIt = srcRowIt->begin();
        auto dstColIt = dstRowIt->begin();
        for (; srcColIt != srcRowIt->end() && dstColIt != dstRowIt->end(); ++srcColIt, ++dstColIt) {
            assert(srcColIt.index() == dstColIt.index());
            *dstColIt = *srcColIt;
        }

        assert(srcColIt == srcRowIt->end());
        assert(dstColIt == dstRowIt->end());
    }

    assert(srcRowIt == src.end());
    assert(dstRowIt == dst.end());
}
}
namespace Opm{
  using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using VectorHP = Dune::MultiTypeBlockVector<Vector, Vector>;
    using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
  template <class value_t>
  class MyDenseMatrix : public Dune::DynamicMatrix<value_t>{
  public:
    static void luDecomp(Dune::DynamicMatrix<value_t>& mat){
      Vector tmp(mat.N());
      typename Dune::DynamicMatrix<value_t>::template Elim<Vector> elim(tmp);
      Dune::Simd::Mask<double> nonsing(true);
      Dune::DynamicMatrix<value_t>::luDecomposition(mat,elim, nonsing, false, false);
    }
    // Partially pivoted variant. The unpivoted decomposition overflows to
    // Inf/NaN when O(1) pinned closed-cell rows mix with O(1e9) BEM rows.
    using pivot_index = std::size_t; // == DenseMatrix::simd_index_type for scalar entries
    static void luDecomp(Dune::DynamicMatrix<value_t>& mat,
                         std::vector<pivot_index>& pivot){
      pivot.resize(mat.N());
      typename Dune::DynamicMatrix<value_t>::ElimPivot elim(pivot);
      Dune::Simd::Mask<double> nonsing(true);
      Dune::DynamicMatrix<value_t>::luDecomposition(mat, elim, nonsing, false, true);
    }
  };


  using FMatrix = Dune::DynamicMatrix<double>; // full matrix

using FractureSystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
using Dune::Indices::_0;
using Dune::Indices::_1;

// Explicit inverse of the elasticity (DDM) matrix K of one fracture mesh; valid while K and version (the DDM generation) match.
struct DdmInverse
{
    std::uint64_t version{0};
    const FMatrix* K{nullptr};
    int n{0};
    std::vector<double> kinv; // row-major
};
std::shared_ptr<const DdmInverse> makeDdmInverse(const FMatrix& K, std::uint64_t version);


// Unscaled mechanics block as diag(alpha) * A + diag(beta), A the elasticity matrix.
// generation changes whenever A is reassembled.
struct MechRows
{
    const FMatrix* ddm{nullptr};
    const std::vector<double>* alpha{nullptr};
    const std::vector<double>* beta{nullptr};
    std::uint64_t generation{0};
    const std::vector<double>* centres{nullptr}; // cell centres, xyz per cell
};

// Approximate solve with the mechanics block, used by the preconditioner.
class MechBlockSolver
{
public:
    virtual ~MechBlockSolver() = default;
    virtual void solve(Vector& x, const Vector& rhs) const = 0;
};

// Dense LU (or inverse) of the mechanics block B = s * (diag(alpha) A + diag(beta)),
// kept across Newton iterations and fracture solves while A is unchanged. Rows whose
// coefficients moved since the factorisation are added as a Woodbury update; rows that
// moved by less than tol are ignored, so the result is a preconditioner, exact when tol=0.
class DenseMechFactor : public MechBlockSolver
{
public:
    struct Options
    {
        bool inverse{true};
        double tol{0.01};
        double maxRankFraction{0.05};
        int maxRankMin{32};
        // refactor once a Krylov solve needs more iterations; 0: refactor when the slots are full
        int refactorIterations{10};
    };
    explicit DenseMechFactor(const Options& opt) : opt_(opt) {}
    // B is the current (scaled) block; returns true if it was refactored.
    bool update(const FMatrix& B, double scale, const MechRows& rows);
    void solve(Vector& x, const Vector& rhs) const override;
    int rank() const { return static_cast<int>(slotRow_.size()); }
    // the next update refactors
    void invalidate() { stale_ = true; }

private:
    void factor(const FMatrix& B, double scale, const MechRows& rows);
    bool addRows(const std::vector<int>& rows);
    bool factorCapacitance(const MechRows& rows);

    Options opt_;
    int n_{0};
    std::uint64_t generation_{0};
    double s0_{1.0};
    double ratio_{1.0}; // s0 / s
    std::vector<double> alpha0_, beta0_;
    std::vector<double> f_; // row-major B0: LU of its transpose, or B0^{-1}
    std::vector<int> piv_;
    // Woodbury: one slot per updated row j, Z(:,k) = B0^{-1} e_j, Ar(:,k) = A(j,:)^T
    std::vector<int> slotOf_;
    std::vector<int> slotRow_;
    std::vector<double> z_, ar_;
    std::vector<double> g_; // G(k,l) = A(row_k,:) z_l, column-major, leading dimension maxRank_
    int maxRank_{0};
    bool stale_{false};
    std::vector<double> va_, vb_;
    std::vector<double> cap_;
    std::vector<int> capPiv_;
    mutable std::vector<double> y_, t_;
};

// ILU(0) of the mechanics block restricted to cell pairs closer than radius times the
// typical cell spacing. Cheap enough to refactor at every Newton iteration.
class NearFieldMechILU : public MechBlockSolver
{
public:
    explicit NearFieldMechILU(double radius) : radius_(radius) {}
    void update(const FMatrix& B, const MechRows& rows);
    void solve(Vector& x, const Vector& rhs) const override;
    double nnzPerRow() const { return n_ > 0 ? double(col_.size()) / n_ : 0.0; }

private:
    void buildPattern(const std::vector<double>& xyz);

    double radius_;
    int n_{0};
    std::uint64_t generation_{0};
    bool havePattern_{false};
    std::vector<int> ptr_, col_, diag_;
    std::vector<double> val_;
};

// Block preconditioner for the coupled fracture width / pressure linear system.
//
// The operator is assembled as
//   [ A  I ]
//   [ C  M ]
// where A is the dense mechanics block, M is the sparse flow block, and I/C are
// the cross-coupling terms. The preconditioner supports three application modes:
//
// - mech_last: solve the flow block first and then apply the mechanics solve to
//   the pressure-corrected mechanics residual.
// - mech_first: solve the mechanics block first and then apply the flow solve to
//   the mechanics-corrected flow residual.
// - fixed_stress: keep the exact mechanics solve, but replace the flow solve by
//   a fixed-stress Schur approximation built from M and the diagonal of A.
//
// The mechanics solve is either diagonal or dense-LU based. The flow solve is
// either diagonal or delegated to FlexibleSolver, which makes it cheap to test
// exact sparse-flow solves on difficult coupled cases.
class FractureMechanicsPreconditioner :
   public Dune::Preconditioner<VectorHP, VectorHP>
{

public:
  // prm expects the fracture linsolver preconditioner subtree. The relevant
  // switches are diag_mech, diag_flow, mech_first, fixed_stress,
  // mech_press_coupling, lu_pivoting, mech_solver.type and flow_solver.*.
  // A non-null mech replaces the internal mechanics factorisation; it must outlive this object.
  FractureMechanicsPreconditioner(const FractureSystemMatrix& S, Opm::PropertyTree prm,
                                  const MechBlockSolver* mech = nullptr,
                                  std::shared_ptr<const DdmInverse> ddm = {});
  virtual void apply(VectorHP& v, const VectorHP& d);
    virtual void post(VectorHP& /*v*/) { };
    virtual void pre(VectorHP& /*x*/, VectorHP& /*b*/) { };
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }
    // Refresh the internal factors after the outer nonlinear step assembles a
    // new coupled system. new_lu_mech controls whether the dense mechanics LU is
    // rebuilt or only the flow-side data is refreshed.
    void update(const Opm::FractureSystemMatrix& S, bool new_lu_mech,
                std::shared_ptr<const DdmInverse> ddm = {});
    int factorizations() const { return factorizations_; }
    double factorSeconds() const { return factor_seconds_; }
    double applySeconds() const { return apply_seconds_; }
private:
enum class ApplyMode {
  MechLast,
  MechFirst,
  FixedStress,
};

// Flow-first block triangular application.
void applymech_last(Opm::VectorHP& v, const Opm::VectorHP& d);
// Mechanics-first block triangular application.
void applymech_first(Opm::VectorHP& v, const Opm::VectorHP& d);
// Mechanics-first application with a fixed-stress Schur approximation on flow.
void applyfixed_stress(Opm::VectorHP& v, const Opm::VectorHP& d);
  ApplyMode selectMode(const FractureSystemMatrix& S);
  double estimateCouplingIndicator(const FractureSystemMatrix& S) const;
  void rebuildFlowSolver(const Opm::FractureSystemMatrix& S);
  static const char* modeName(ApplyMode mode);
template <typename Mat>
Vector diagvec(const Mat& M)
  {
        Vector res(M.M());
        for (size_t i = 0; i != res.size(); ++i)
      res[i] = Opm::diagScalar(M[i][i]);
        return res;
   }
   template <typename Mat>
  void setDiagvec(Vector& res, const Mat& M)
  {
        assert(res.size() == M.M() && res.size() == M.N());
        for (size_t i = 0; i != res.size(); ++i)
      res[i] = Opm::diagScalar(M[i][i]);
   }
  // Assemble the fixed-stress flow operator M - C * diag(A)^{-1} * I using the
  // sparsity of M augmented with active C entries.
  void updateFixedStressFlowMatrix(const Opm::FractureSystemMatrix& S);
  void solveMechanics(Vector& x, const Vector& rhs) const;
  void solveFlow(Vector& x, const Vector& rhs);
  void factorMechanics(const FMatrix& A);
  void factorOpenBlock(const FMatrix& A);
  bool factorSplit(const FMatrix& A);
  void splitSolve(Vector& x, const Vector& rhs) const;
  void backSolve(Vector& x,const Vector& rhs_in) const;
  void openBlockSolve(Vector& x, const Vector& rhs) const;
  void checkLuFactor() const;
    const FractureSystemMatrix& A_;
    const MechBlockSolver* mech_{nullptr};
    int factorizations_{0};
    double factor_seconds_{0.0};
    double apply_seconds_{0.0};
    mutable FMatrix luM_;
    std::vector<MyDenseMatrix<double>::pivot_index> lu_pivot_;
    bool lu_pivoting_{false};
    // mech_solver.type: "dune_lu" (default), "lapack_lu" or "lapack_inv"
    bool lapack_lu_{false};
    // lapack_inv: open_lu_ holds the explicit inverse, applied with dgemv
    bool lapack_inv_{false};
    // lapack_lu: rows of A with only a diagonal entry are solved directly, the rest by dgetrf
    std::vector<int> open_, closed_;
    std::vector<double> open_lu_;   // row-major open block, factored as its transpose
    std::vector<int> open_pivot_;
    std::vector<double> closed_diag_;
    std::vector<double> open_closed_; // A(open, closed), row-major
    mutable std::vector<double> open_rhs_;
    mutable std::vector<double> open_x_;
    mutable std::vector<double> closed_x_;
    // fb_split: rows of A of the form -s_i K_i + g_i e_i with small g_i are taken as -s_i K_i
    // and solved with K^{-1}; the other rows R enter exactly through a |R|x|R| Woodbury correction.
    bool fb_split_{false};
    bool split_active_{false};
    double split_tol_{1e-2};
    double split_max_fraction_{0.5};
    std::shared_ptr<const DdmInverse> ddm_;
    std::vector<int> split_rows_;
    std::vector<double> split_scale_; // 1/s_i, or 1 on rows in R
    std::vector<double> split_Ar_;    // A(R, :), row-major
    std::vector<double> split_KR_;    // K^{-1}(:, R), row-major
    std::vector<double> split_G_;     // LU of A(R,:) K^{-1}(:,R)
    std::vector<int> split_piv_;
    mutable std::vector<double> split_z_, split_y_, split_t_;
    int verbosity_{0};
    mutable Vector A_diag_;
    mutable Vector M_diag_;
    Opm::PropertyTree prm_;
    using FlowOperatorType = Dune::MatrixAdapter<SMatrix, Vector, Vector>;
    std::unique_ptr< FlowOperatorType> flowop_;
    std::unique_ptr< Dune::FlexibleSolver<FlowOperatorType> > flow_solver_;
    std::unique_ptr<SMatrix> fixed_stress_matrix_;
    std::unique_ptr<FlowOperatorType> fixed_stress_flowop_;
    std::unique_ptr<Dune::FlexibleSolver<FlowOperatorType>> fixed_stress_flow_solver_;
  bool diag_mech_{true};
  bool diag_flow_{true};
  bool mech_press_coupling_{false};
  bool press_mech_coupling_{false};
  bool mech_first_{true};
  bool fixed_stress_{false};
  ApplyMode active_mode_{ApplyMode::MechFirst};
  std::string mode_policy_{"manual"};
  double mode_switch_coupling_threshold_{1.0};
  double last_coupling_indicator_{0.0};
};

} // namespace Opm
