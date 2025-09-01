#pragma once

#include <vector>
#include <memory>
#include <dune/istl/preconditioners.hh>
namespace Opm{
using Vector = Dune::BlockVector<double>;
using VectorHP = Dune::MultiTypeBlockVector<Vector, Vector>;
class FractureMechanicsPreconditioner : 
   public Dune::Preconditioner<VectorHP, VectorHP>
{
public:
    FractureMechanicsPreconditioner()(const SystemMatrix& S,Opm::PropertyTree prm)
        : A_(S[_0][_0])
        , A_diag_(diagvec(S[_0][_0]))
        , M_diag_(diagvec(S[_1][_1]))
        , prm_(prm)
    {
    }
    virtual void apply(VectorHP& v, const VectorHP& d)
    {
        for (size_t i = 0; i != A_diag_.size(); ++i){
            v[_0][i] = d[_0][i] / A_diag_[i];
        }
        for (size_t i = 0; i != M_diag_.size(); ++i){
            v[_1][i] = d[_1][i] / M_diag_[i];
        }
    };
    virtual void post(VectorHP& v) { };
    virtual void pre(VectorHP& x, VectorHP& b) { };
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const FullMatrix& A_;
    const Vector A_diag_;
    const Vector M_diag_;
    Opm::PropertyTree prm_;
};

} // namespace Opm::Geomech