#pragma once

#include <vector>
#include <memory>
#include <opm/common/TimingMacros.hpp>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
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

  
class FractureMechanicsPreconditioner : 
   public Dune::Preconditioner<VectorHP, VectorHP>
{
    
public:
  FractureMechanicsPreconditioner(const SystemMatrix& S, Opm::PropertyTree prm);
  virtual void apply(VectorHP& v, const VectorHP& d);
    virtual void post(VectorHP& /*v*/) { };
    virtual void pre(VectorHP& /*x*/, VectorHP& /*b*/) { };
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
template <typename Mat>
Vector diagvec(const Mat& M)
  {
        Vector res(M.M());
        for (size_t i = 0; i != res.size(); ++i)
            res[i] = M[i][i];
        return res;
   }
  void backSolve(Vector& x,const Vector& rhs_in);
    const SystemMatrix& A_;
    mutable FMatrix luM_;
    const Vector A_diag_;
    const Vector M_diag_;
    Opm::PropertyTree prm_;
    using FlowOperatorType = Dune::MatrixAdapter<SMatrix, Vector, Vector>;
    std::unique_ptr< FlowOperatorType> flowop_;
    std::unique_ptr< Dune::FlexibleSolver<FlowOperatorType> > flow_solver_;
  //
  bool diag_mech_{true};
  bool diag_flow_{true};
  bool mech_press_coupling_{false};
  bool press_mech_coupling_{false};
  bool mech_first_{true};
};

} // namespace Opm::Geomech
