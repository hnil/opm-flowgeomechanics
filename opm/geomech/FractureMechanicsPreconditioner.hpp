#pragma once

#include <vector>
#include <memory>
#include <dune/istl/preconditioners.hh>
#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <opm/simulators/linalg/PropertyTree.hpp>
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
    FractureMechanicsPreconditioner(const SystemMatrix& S, Opm::PropertyTree prm)
        : A_(S)
        , A_diag_(diagvec(S[_0][_0]))
        , M_diag_(diagvec(S[_1][_1]))
        , prm_(prm)          
    {
      OPM_TIMEFUNCTION();
      diag_mech_ = prm.get<bool>("diag_mech");
      diag_flow_ = prm.get<bool>("diag_flow");
      luM_ = S[_0][_0];
      MyDenseMatrix<double>::luDecomp(luM_);
    }
    virtual void apply(VectorHP& v, const VectorHP& d)
    {
      // SystemMatrix S {{A, I}, // mechanics system (since A is negative, we leave I positive here)
      //               {C, M}}; // flow system
      OPM_TIMEFUNCTION_LOCAL();
      if(diag_mech_){
        for (size_t i = 0; i != A_diag_.size(); ++i){
            v[_0][i] = d[_0][i] / A_diag_[i];
        }
      }else{
        // solve full mechanics system
        A_[_0][_0].solve(v[_0],d[_0]);
      }
      if(diag_flow_){
        for (size_t i = 0; i != M_diag_.size(); ++i){
            v[_1][i] = d[_1][i] / M_diag_[i];
        }
      }
    };
    virtual void post(VectorHP& /*v*/) { };
    virtual void pre(VectorHP& /*x*/, VectorHP& /*b*/) { };
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
template <typename Mat>
    Vector diagvec(const Mat& M){
        Vector res(M.M());
        for (size_t i = 0; i != res.size(); ++i)
            res[i] = M[i][i];
        return res;
    }
  void backSolve(Vector& x,const Vector& rhs_in){
    Vector& rhs = x;
    rhs = rhs_in;
    for(int i=0; i < luM_.rows(); i++) {
      for(int j=0; j<i-1; j++){
        rhs[i] -= luM_[i][j]*x[j];
      }
      rhs[i] = rhs[i]/luM_[i][i];
    }

    for(int i=luM_.rows()-1; i>=0; i--) {
        for (size_t j=i+1; j<luM_.rows(); j++){
          rhs[i] -= luM_[i][j]*x[j];
        }
        rhs[i] = rhs[i]/luM_[i][i];
      }
  }
    const SystemMatrix& A_;
    mutable FMatrix luM_;
    const Vector A_diag_;
    const Vector M_diag_;
    Opm::PropertyTree prm_;
  //
  bool diag_mech_;
  bool diag_flow_;
};

} // namespace Opm::Geomech
