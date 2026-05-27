#pragma once

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
  using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;//Opm::Fracture::Vector
using VectorHP = Dune::MultiTypeBlockVector<Vector, Vector>;
    using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
  // class MySwap{
  // public:
  //   void swap(int i, int j){}
  // };
  template <class value_t>
  class MyDenseMatrix : public Dune::DynamicMatrix<value_t>{
  public:
    static void luDecomp(Dune::DynamicMatrix<value_t>& mat){
      Vector tmp(mat.N());
      typename Dune::DynamicMatrix<value_t>::template Elim<Vector> elim(tmp);
      //AutonomousValue<MAT> A(asImp());
      Dune::Simd::Mask<double> nonsing(true);
      Dune::DynamicMatrix<value_t>::luDecomposition(mat,elim, nonsing, false, false);
    } 
  };


  using FMatrix = Dune::DynamicMatrix<double>; // full matrix

using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;
using Dune::Indices::_0;
using Dune::Indices::_1;

  
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
   //public Dune::PreconditionerWithUpdate<VectorHP, VectorHP>
{
    
public:
  // prm expects the fracture linsolver preconditioner subtree. The relevant
  // switches are diag_mech, diag_flow, mech_first, fixed_stress,
  // mech_press_coupling, and flow_solver.*.
  FractureMechanicsPreconditioner(const SystemMatrix& S, Opm::PropertyTree prm);
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
    void update(const Opm::SystemMatrix& S,bool new_lu_mech);
    // void update(){
    //     // default to recomputing the preconditioner, but in some cases (e.g. when only a few cells are closed) we can just update the existing preconditioner
    //     assert(false);
    //     update(A_);
    // };
    //bool hasPerfectUpdate() const override { return false; }
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
  ApplyMode selectMode(const SystemMatrix& S);
  double estimateCouplingIndicator(const SystemMatrix& S) const;
  void rebuildFlowSolver(const Opm::SystemMatrix& S);
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
  void updateFixedStressFlowMatrix(const Opm::SystemMatrix& S);
  void solveMechanics(Vector& x, const Vector& rhs) const;
  void solveFlow(Vector& x, const Vector& rhs);
  void backSolve(Vector& x,const Vector& rhs_in) const;
    const SystemMatrix& A_;
    mutable FMatrix luM_;
    mutable Vector A_diag_;
    mutable Vector M_diag_;
    Opm::PropertyTree prm_;
    using FlowOperatorType = Dune::MatrixAdapter<SMatrix, Vector, Vector>;
    std::unique_ptr< FlowOperatorType> flowop_;
    std::unique_ptr< Dune::FlexibleSolver<FlowOperatorType> > flow_solver_;
    std::unique_ptr<SMatrix> fixed_stress_matrix_;
    std::unique_ptr<FlowOperatorType> fixed_stress_flowop_;
    std::unique_ptr<Dune::FlexibleSolver<FlowOperatorType>> fixed_stress_flow_solver_;
  //
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

} // namespace Opm::Geomech
