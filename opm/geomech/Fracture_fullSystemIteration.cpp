#include "config.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <numeric>

#include <dune/common/indices.hh> // needed for _0, _1, etc.
//
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>


#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/DiagonalScalar.hpp>
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

using namespace std;
#define EXP_REF 0
namespace
{
static int DEBUG_COUNT = 0;
[[maybe_unused]] string
debug_filename(const string& prefix, const string& suffix = ".txt")
{
    std::ostringstream oss;
    oss << prefix << DEBUG_COUNT++ << suffix;
    return oss.str();
}

// ========================== Convenience definitions ==========================
using Dune::Indices::_0;
using Dune::Indices::_1;

// using ResVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using ResVector = Opm::Fracture::Vector;
using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
using FMatrix = Dune::DynamicMatrix<double>; // full matrix

using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;

using Htrans = std::tuple<size_t, size_t, double, double>;
using LinearPrecondType = Opm::FractureMechanicsPreconditioner;
using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;

// ============================ Debugging functions ============================
// ----------------------------------------------------------------------------
// template<class MatrixType> void dump_matrix(const MatrixType& m, const char* const name)
[[maybe_unused]] void
dump_matrix(const FMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    const size_t rows = m.N();
    const size_t cols = m.M();
    std::ofstream os(name);
    for (size_t i = 0; i != rows; ++i)
        for (size_t j = 0; j != cols; ++j)
            os << m[i][j] << ((j == cols - 1) ? "\n" : " ");
    os.close();
}

// ----------------------------------------------------------------------------
// template<>
[[maybe_unused]] void
dump_matrix(const SMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    std::ofstream os(name);
    for (auto rowIt = m.begin(); rowIt != m.end(); ++rowIt) {
        const size_t i = rowIt.index();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const size_t j = colIt.index();
            os << i + 1 << " " << j + 1 << " " << m[i][j] << "\n";
        }
    }
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const vector<int>& v, const char* const name, bool append = false)
// ----------------------------------------------------------------------------
{
    // open for append is requested
    std::ofstream os(name, append ? std::ios::app : std::ios::out);
    for (size_t i = 0; i != v.size(); ++i)
        os << v[i] << "\n";
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const ResVector& v, const char* const name, bool append = false)
// ----------------------------------------------------------------------------
{
    if (!name)
        for (size_t i = 0; i != v.size(); ++i)
            std::cout << v[i] << "\n";
    if (!name)
        return;

    // open for append is requested
    std::ofstream os(name, append ? std::ios::app : std::ios::out);
    for (size_t i = 0; i != v.size(); ++i)
        os << v[i] << "\n";
    os.close();
}

// ----------------------------------------------------------------------------
[[maybe_unused]] void
dump_vector(const VectorHP& v, const char* const name1, const char* const name2, bool append = false)
// ----------------------------------------------------------------------------
{
    dump_vector(v[_0], name1, append);
    dump_vector(v[_1], name2, append);
}

// ============================= Helper functions =============================

// ----------------------------------------------------------------------------
SMatrix
makeIdentity(size_t num,
             double right_padding = 0,
             double fac = 1,
             std::vector<int> zero_rows = std::vector<int>())
// ----------------------------------------------------------------------------
{
    // build a sparse matrix to represent the identity matrix of a given size
    SMatrix M;
    M.setBuildMode(SMatrix::implicit);
    M.setImplicitBuildModeParameters(1, 0.1);
    M.setSize(num, num + right_padding);
    for (size_t i = 0; i != num; ++i)
        M.entry(i, i) = fac;

    // set zero rows, if applicable
    for (size_t i = 0; i != zero_rows.size(); ++i)
        if (zero_rows[i])
            M.entry(i, i) = 0;

    M.compress();
    return M;
}
void modifyIdentity(SMatrix& M,
                  double fac = 1,
                  std::vector<int> zero_rows = std::vector<int>()){
                  int num =   M.N();
    for (int i = 0; i != num; ++i){
        M[i][i] = fac;
    }
    // set zero rows, if applicable
    for (size_t i = 0; i != zero_rows.size(); ++i){
        if (zero_rows[i]){
            M[i][i] = 0;
        }
    }
}

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

// ----------------------------------------------------------------------------
std::unique_ptr<LinearSolverBase>
setupLinearSolver(const std::string& solver_type,
                  const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>& linop,
                  LinearPrecondType& precond,
                  const Opm::PropertyTree& prm,
                  const double lintol,
                  const int max_iter,
                  const int verbosity)
// ----------------------------------------------------------------------------
{
    if (solver_type == "bicgstab") {
        return std::make_unique<Dune::BiCGSTABSolver<VectorHP>>(linop,
                                                                 precond,
                                                                 lintol,
                                                                 max_iter,
                                                                 verbosity);
    }

    if (solver_type == "gmres") {
        const int restart = prm.get<int>("solver.linsolver.restart", 100);
        return std::make_unique<Dune::RestartedGMResSolver<VectorHP>>(linop,
                                                                       precond,
                                                                       lintol,
                                                                       restart,
                                                                       max_iter,
                                                                       verbosity);
    }

    if (solver_type == "fgmres") {
        const int restart = prm.get<int>("solver.linsolver.restart", 100);
        return std::make_unique<Dune::RestartedFlexibleGMResSolver<VectorHP>>(linop,
                                                                               precond,
                                                                               lintol,
                                                                               restart,
                                                                               max_iter,
                                                                               verbosity);
    }

    OPM_THROW(std::runtime_error, "Unknown linear solver type: for fracture");
}

// ----------------------------------------------------------------------------
std::unique_ptr<Opm::Fracture::Matrix>
initCouplingMatrixSparsity(const std::vector<Htrans>& htrans, size_t num_cells, size_t num_wells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    auto C = std::make_unique<Opm::Fracture::Matrix>(
        num_cells + num_wells, num_cells, 4, 0.4, Opm::Fracture::Matrix::implicit);
    for (const auto& e : htrans) {
        const size_t i = std::get<0>(e);
        const size_t j = std::get<1>(e);
        assert(i != j);
        C->entry(i, j) = 0;
        C->entry(j, i) = 0;
        C->entry(i, i) = 0;
        C->entry(j, j) = 0;
    }
    C->compress();
    return C;
}

// ----------------------------------------------------------------------------
void
updateCouplingMatrix(std::unique_ptr<Opm::Fracture::Matrix>& Cptr,
                     const size_t num_cells,
                     const size_t num_wells,
                     const std::vector<Htrans>& htrans,
                     const ResVector& pressure,// is realy head
                     const ResVector& aperture,
                     const std::vector<int>& closed_cells,
                     const std::vector<double>& reservoir_mobility,
                     double min_width)
{
    OPM_TIMEFUNCTION();
    // create C if not done already
    if (!Cptr)
        Cptr = initCouplingMatrixSparsity(htrans, num_cells, num_wells);
    // if (!Cptr) Cptr = std::make_unique<Opm::Fracture::Matrix>(*M); // copy sparsity of M

    // @@ NB: If the implementation of the pressure matrix changes (i.e. `assemblePressure()`),
    //        the code below might need to be updated accordingly as well.
    auto& C = *Cptr;
    C = 0;
    for (const auto& e : htrans) {
        const size_t i = std::get<0>(e);
        const size_t j = std::get<1>(e);
        assert(i != j);

        const double t1 = std::get<2>(e);
        const double t2 = std::get<3>(e);
        double  dh1 = 1.0;
        double  dh2 = 1.0;
        if(aperture[i][0] <= min_width) dh1 = 0.0;
        if(aperture[j][0] <= min_width) dh2 = 0.0;
        assert(aperture[i][0]>=0.0);
        assert(aperture[j][0]>=0.0);        
        const double h1 = std::max(aperture[i][0],min_width);//aperture[i] + min_width;
        const double h2 = std::max(aperture[j][0],min_width);//aperture[j] + min_width;
        const double p1 = pressure[i];//-fracture_dgh[i];
        const double p2 = pressure[j];//-fracture_dgh[j];; pressure should already be head

        const double q = (h1 * h1 * h1) * (h2 * h2 * h2) * (t1 * t2); // numerator
        const double d1q = 3 * (h1 * h1) * (h2 * h2 * h2) * (t1 * t2)* dh1;
        const double d2q = 3 * (h1 * h1 * h1) * (h2 * h2) * (t1 * t2)* dh2;

        //const double r = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
        const double r = 12 * (h1 * h1 * h1 * t1 + h2 * h2* h2 * t2); // denominator
        const double d1r = 36 * (h1 * h1) * t1 *dh1;
        const double d2r = 36 * (h2 * h2) * t2 *dh2;
        //const double d1r = 24 * (h1) * t1;
        //const double d2r = 24 * (h2) * t2;
        assert(r >= 0.0);
        assert(q >= 0.0);
        const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
        const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);

        const double mobility = 0.5 * (reservoir_mobility[i] + reservoir_mobility[j]);
        //const double krull = 1; // 1e4; // @@ Not sure if this should be removed?
        // diagonal elements
        C[i][i] += dTdh1 * (p1 - p2) * mobility;
        C[j][j] += dTdh2 * (p2 - p1) * mobility;

        // off-diagonal elements
        C[i][j] += dTdh2 * (p1 - p2) * mobility;
        C[j][i] += dTdh1 * (p2 - p1) * mobility;
    }
    // zeroing out columns corresponding to closed cells
    if(true){
      for (size_t col = 0; col != closed_cells.size(); ++col)
        if (closed_cells[col])
          for (size_t row = 0; row != C.N(); ++row)
            if (C.exists(row, col))
              C[row][col] = 0;
    }
}

// ----------------------------------------------------------------------------
inline bool
convergence_test(const VectorHP& res, const double tol_flow, double tol_mech, int verbosity)
// ----------------------------------------------------------------------------
{
    // std::cout << "Residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
    // std::cout << "Residual norm[1] is " << res[_1].infinity_norm()<<std::endl;

    // std::cout << "tol mech is: " << tol_mech << std::endl;
    // std::cout << "tol flow is: " << tol_flow << std::endl;

    bool converged = res[_0].infinity_norm() < tol_mech && res[_1].infinity_norm() < tol_flow;
    if(converged && verbosity > 0){
      std::cout << "Converged with residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
      std::cout << "Converged with residual norm[1] is " << res[_1].infinity_norm()<<std::endl;
    }
    // else{
    //   std::cout << "Not converged with residual norm[0] is " << res[_0].infinity_norm()<<std::endl;
    //   std::cout << "Not converged with residual norm[1] is " << res[_1].infinity_norm()<<std::endl;
    // }
    return converged;
    // return res.infinity_norm() < tol;
}

// ----------------------------------------------------------------------------
template <typename Mat>
ResVector
diagvec(const Mat& M)
// ----------------------------------------------------------------------------
{
    ResVector res(M.M());
    for (size_t i = 0; i != res.size(); ++i)
        res[i] = Opm::diagScalar(M[i][i]);
    return res;
}


// ----------------------------------------------------------------------------
class TailoredPrecondDiag : public Dune::Preconditioner<VectorHP, VectorHP>
// ----------------------------------------------------------------------------
{
public:
    TailoredPrecondDiag(const SystemMatrix& S)
        : A_diag_(diagvec(S[_0][_0]))
        , M_diag_(diagvec(S[_1][_1]))
    {
    }
    virtual void apply(VectorHP& v, const VectorHP& d)
    {
        for (size_t i = 0; i != A_diag_.size(); ++i)
            v[_0][i] = d[_0][i] / A_diag_[i];
        for (size_t i = 0; i != M_diag_.size(); ++i)
            v[_1][i] = d[_1][i] / M_diag_[i];
    };
    virtual void post(VectorHP& /*v*/) {};
    virtual void pre(VectorHP& /*x*/, VectorHP& /*b*/) {};
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const ResVector A_diag_;
    const ResVector M_diag_;
};


// ----------------------------------------------------------------------------
[[maybe_unused]] double
estimate_step_fac(const VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------
{
    // estimate what might be a safe step size to avoid exiting the convergence radius
    const double x0_norm = x[_0].infinity_norm();
    const double x1_norm = x[_1].infinity_norm();

    const double f1 = x0_norm == 0 ? 1 : dx[_0].infinity_norm() / x0_norm;
    const double f2 = x1_norm == 0 ? 1 : dx[_1].infinity_norm() / x1_norm;
    const double fmax = std::max(f1, f2);
    const double threshold = 0.95;
    const double fac_min = 1e-1; // @@ 1e-2;
    return (fmax < threshold) ? 1.0 : std::max(threshold / fmax, fac_min);
}

// ----------------------------------------------------------------------------
FMatrix
modified_fracture_matrix(const FMatrix& A, const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    Dune::DynamicMatrix result(A);

    for (size_t row = 0; row != A.N(); ++row)
        if (closed_cells[row])
            for (size_t col = 0; col != A.N(); ++col)
                result[row][col] = (row == col);
    return result;
}

void modified_fracture_matrix(FMatrix& Aold, const FMatrix& A, const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    //Dune::DynamicMatrix result(A);
    assert(Aold.N() == A.N());
    assert(Aold.M() == A.M());
    for(size_t i=0; i < A.N(); ++i){
      for(size_t j=0; j < A.M(); ++j){
        Aold[i][j] = A[i][j];
      }
    }
    for (size_t row = 0; row != A.N(); ++row)
        if (closed_cells[row])
            for (size_t col = 0; col != A.N(); ++col)
                Aold[row][col] = (row == col);
}

}; // end anonymous namespace

namespace Opm
{
// ----------------------------------------------------------------------------
std::vector<int>
Fracture::identify_closed(const FMatrix& A, const VectorHP& x, const ResVector& rhs, const int nwells)
// ----------------------------------------------------------------------------
{
    
    OPM_TIMEFUNCTION();
    std::string closing_type = prm_.get<std::string>("solver.closing_type", "org");
    ResVector tmp(rhs);
    const auto I = makeIdentity(A.N(), nwells);

    // computing rhs - A x[0] - I x[1]
    std::vector<int> result;
    if(closing_type == "org"){
        const ResVector& h = x[_0];
        const ResVector& p = x[_1];
        // maybe use other pressure if closed
        A.mmv(h, tmp);
        I.mmv(p, tmp);
        for (size_t i = 0; i != A.N(); ++i)
            result.push_back(tmp[i] >= 0 && h[i] <= 0.0);
    } else if(closing_type == "open") {
            result.resize(A.N(), 1);
    } else if(closing_type == "simple"){
     // sum forces over faces
        result.resize(A.N(), 0);
        //double sum_force = 0.0;
        for (auto element : Dune::elements(grid_->leafGridView())) {
           // int i = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView,
            // Dune::mcmgElementLayout>::index(element);
            int i = grid_->leafGridView().indexSet().index(element);
            //double area = element.geometry().volume();
            double pressure = fracture_pressure_[i][0];
            //NB maybe use other pressure if closed
            if((this->fractureForce(i)-pressure) < 0){
                result[i] = 1;
            } else {
                result[i] = 0;
            }
        }
    } else {
        std::stringstream ss;
        ss << "Unknown closing type: " << closing_type;
        OPM_THROW(std::runtime_error, ss.str());
    }
        
    // std::fill(result.begin(), result.end(), 1); // @@@
    return result;
}

// ----------------------------------------------------------------------------
void
Fracture::assemblePressureAndCouplingAD(const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    updateLeakoff();

    const size_t nc = numFractureCells();
    const size_t nw = numWellEquations();
    const std::string control_type = prm_.get_child("control").get<std::string>("type");

    // Build standalone input from Fracture member data
    FracturePressureInput input;
    input.num_cells = nc;
    input.min_width = min_width_;
    input.htrans = htrans_;
    input.control_type = control_type;
    input.num_well_equations = nw;
    input.perfinj.assign(perfinj_.begin(), perfinj_.end());
    input.total_wellindex = total_wellindex_;
    input.mobility_water_perf = mobility_water_perf_;

    // Convert fracture_width_ (BlockVector<FieldVector<double,1>>) to std::vector<double>
    input.fracture_width.resize(nc + nw);
    for (size_t i = 0; i < nc + nw; ++i)
        input.fracture_width[i] = fracture_width_[i][0];

    // Convert fracture_pressure_ to std::vector<double>, using head (pressure - dgh)
    input.fracture_pressure.resize(nc + nw);
    for (size_t i = 0; i < nc + nw; ++i) {
        input.fracture_pressure[i] = fracture_pressure_[i][0];
        if (i < fracture_dgh_.size())
            input.fracture_pressure[i] -= fracture_dgh_[i];
    }

    // Map reservoir_mobility_ to density/viscosity: use mobility as density, viscosity = 1
    // This yields mobility = density/viscosity = reservoir_mobility_[i] / 1.0
    input.density.resize(nc + nw);
    input.viscosity.resize(nc + nw);
    for (size_t i = 0; i < nc; ++i) {
        input.density[i] = {reservoir_mobility_[i], 0.0};
        input.viscosity[i] = {1.0, 0.0};
    }
    for (size_t i = nc; i < nc + nw; ++i) {
        input.density[i] = {1.0, 0.0};
        input.viscosity[i] = {1.0, 0.0};
    }

    input.leakof = leakof_;

    // Run the AD-based assembly (produces both pressure and coupling matrices)
    auto result = assemblePressureAD(input);

    // Store the pressure matrix
    pressure_matrix_ = std::move(result.pressure_matrix);

    // Store the coupling matrix
    coupling_matrix_ = std::move(result.coupling_matrix);

    // Zero out coupling matrix columns for closed cells
    auto& C = *coupling_matrix_;
    for (size_t col = 0; col < closed_cells.size(); ++col)
        if (closed_cells[col])
            for (size_t row = 0; row < C.N(); ++row)
                if (C.exists(row, col))
                    C[row][col] = 0;
}

// ----------------------------------------------------------------------------
bool
Fracture::fullSystemIteration(const double tol, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();

    ++DEBUG_COUNT;
    // min_width_ = prm_.get<double>("solver.min_width"); // min with only used for flow calculations
    //const double max_width = prm_.get<double>("solver.max_width");

    // update pressure matrix with the current values of `fracture_width_` and
    // `fracture_pressure_`
    // std::cout << "---- Assemble pressure" << std::endl;
    const bool use_ad = prm_.get<bool>("solver.use_ad_pressure_assembly", false);
    if (!use_ad) {
        assemblePressure(); // update pressure matrix (original method)
    }
    addSource(); // update right-hand side of pressure system;

    // std::cout << "---- Various " << std::endl;
    //  initialize vector of unknown, and vector represnting direction in tangent space
    VectorHP x {fracture_width_, fracture_pressure_};
    // dump_vector(x, "w", "p", true); // dump current state of fracture
    // dump_vector(x, debug_filename("w_").c_str(), debug_filename("p_").c_str());

    VectorHP dx = x;
    dx = 0; // gradient of 'x' (which we aim to compute below)


    // set right hand side
    VectorHP rhs {x}; // same size as system, content set below
    normalFractureTraction(rhs[_0], false); // right-hand side equals the normal fracture traction
    rhs[_1] = rhs_pressure_; // should have been updated in call to `assemblePressure` above

    // make a version of the fracture matrix that has trivial equations for closed cells
    const std::vector<int> closed_cells
        = identify_closed(fractureMatrix(), x, rhs[_0], numWellEquations());
    // closed as changed
    bool closed_cells_changed = false;
    if(closed_cells_.size() != closed_cells.size()){
      closed_cells_changed = true;
    }else{
      for(size_t i=0; i!=closed_cells.size(); ++i){
        if(closed_cells_[i] != closed_cells[i]){
          closed_cells_changed = true;
          break;
        }
      }
    }
    //closed_cells_changed = true; // @@@
    closed_cells_ = closed_cells;
    // dump_vector(closed_cells, debug_filename("closed_cells_").c_str());
    // dump_vector(closed_cells, "closed_cells", true);
    // const auto A = modified_fracture_matrix(fractureMatrix(), closed_cells);
    //bool update_mech_matrix = nlin_iteration > 0 || closed_cells_changed;
    bool first_iteration = (nlin_iteration == 0);
    if (first_iteration) {
        A_.reset(new FMatrix(modified_fracture_matrix(fractureMatrix(), closed_cells)));
    }else{
        modified_fracture_matrix(*A_, fractureMatrix(), closed_cells);
    }
    const auto& A = *A_;

    // also modify right hand side for closed cells
    for (size_t i = 0; i != closed_cells.size(); ++i)
        if (closed_cells[i])
            rhs[_0][i] = 0;

    if (use_ad) {
        // Use the standalone AD assembler for both pressure matrix and coupling matrix
        assemblePressureAndCouplingAD(closed_cells);
    } else {
        // update the coupling matrix (possibly create it if not already initialized)
        auto fracture_head(fracture_pressure_);
        assert(fracture_dgh_.size() == (fracture_pressure_.size()-numWellEquations()));
        for (size_t i = 0; i < fracture_dgh_.size(); ++i) {
            fracture_head[i] = fracture_pressure_[i] - fracture_dgh_[i];
        }

        updateCouplingMatrix(coupling_matrix_,
                             pressure_matrix_->N() - numWellEquations(), // num cells
                             numWellEquations(), // num wells
                             htrans_,
                             fracture_head,//Note use head here
                             fracture_width_,
                             closed_cells,
                             reservoir_mobility_,
                             min_width_);
    }
    // setup the full system
    if(prm_.get<bool>("solver.drop_fluid_mech_linearization",true)){
      *coupling_matrix_ = 0;
    }else{
      if(rhs[_0].infinity_norm() > prm_.get<double>("solver.drop_tol_h",1e5) ||
         rhs[_1].infinity_norm() > prm_.get<double>("solver.drop_tol_p",1.0) ){
        *coupling_matrix_ = 0;
        // maybe change linear solver
      }

    }
    const auto& M = *pressure_matrix_;
    const auto& C = *coupling_matrix_;
    
    // const auto I = makeIdentity(A.N(), numWellEquations(), 1, closed_cells);
    if (first_iteration) {
        I_.reset(new SMatrix(makeIdentity(A.N(), numWellEquations(), 1, closed_cells)));
    }else{
        modifyIdentity(*I_, 1, closed_cells);
    }
    const auto& I = *I_;
#if EXP_REF    
    const auto& row1 = Dune::makeMultiTypeBlockVectorRef(A, I);
    const auto& row2 = Dune::makeMultiTypeBlockVectorRef(C, M);
    const auto& SS = Dune::makeMultiTypeBlockMatrixRef(row1, row2);
#endif
  
    //const auto SS_ref = std::make_unique<SystemMatrixRef>(SS); // copy to make sure we have a contiguous storage for the preconditioner
    // NB need to se how to modify this matrix as litle as possible and change only the part
    // of preconditioner need
    // system Jacobian (with cross term)  @@ should S be included as a member variable of Fracture?
    //S_.reset(
    //    new SystemMatrix({{A, I}, // mechanics system (since A is negative, we leave I positive here)
    //                      {C, M}})); // flow system
    // this make a copy of the A,I,C,M

    if(first_iteration){
#if EXP_REF           
        S_ = std::make_unique<SystemMatrix>(SS); // copy to make sure we have a contiguous storage for the preconditioner
#else
        //S_ = std::make_unique<SystemMatrix>({{A, I}, // mechanics system (since A is negative, we leave I positive here)
        //                                      {C, M}}); // flow system
        S_.reset(new SystemMatrix({{A, I},                       
                                   {C, M}}));

#endif
    } else{
        OPM_TIMEBLOCK(Modify_System_Matrix);
        // modify S in place to avoid copying the whole matrix
        auto& S = *S_;
        for(size_t i=0; i!=A.N(); ++i){
            for(size_t j=0; j!=A.M(); ++j){
                S[_0][_0][i][j] = A[i][j];
            }
        }
        copyMatrixValuesWithSameSparsity(S[_0][_1], I);
        copyMatrixValuesWithSameSparsity(S[_1][_0], C);
        copyMatrixValuesWithSameSparsity(S[_1][_1], M);
    }  
    

    // SystemMatrix S = {{A, I},// mechanics system (since A is negative, we leave I positive here)
    //                   {C, M}};
    auto& S = *S_;

    // system equations
    //SystemMatrix S0 = S;
    //S0[_1][_0] = 0; // the equations themselves have no cross term
    // // dump_vector(rhs, debug_filename("rhs_w_").c_str(), debug_filename("rhs_p_").c_str());
    // // dump_vector(rhs, "rhs_w", "rhs_p", true);
    
    //S0.mmv(x, rhs); // rhs = rhs - S0 * x;   (we are working in the tanget plane)
    // write explicit using S
    
       // S[_1][_0].mmv(x[_0], rhs[_1]);
        S[_1][_1].mmv(x[_1], rhs[_1]);
        S[_0][_0].mmv(x[_0], rhs[_0]);
        S[_0][_1].mmv(x[_1], rhs[_0]);   
    //}
   
    auto rhs_org(rhs);

    
    // Verify that equations have been chosen correctly
    // for (size_t i = 0; i != fracture_width_.size(); ++i)
    //   assert( !(rhs[_0][i] >= 0.0 && x[_0][i] <= 0.0));

    // check if system is already at a converged state (in which case we return immediately)
    //
    // for convergence, we use 'tol' directly for the mechanics system (where
    // residuals are expected to scale with pressure), and 'tol * ||M||' for the
    // flow system (where residuals scale with M*p)
    const double tol_flow = tol;//*std::max(A.infinity_norm(), M.infinity_norm()) * std::numeric_limits<double>::epsilon();
    const double tol_mech = tol;// * M.infinity_norm();
    const int nlin_verbosity = prm_.get<double>("solver.verbosity");
    if (convergence_test(rhs,
                         tol_mech,
                         tol_flow, nlin_verbosity))
        return true;


    // solve system equations
    //const Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S);
    if(first_iteration){
        S_linop_ = std::make_unique<Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>>(S); // store for use in preconditioner
    }
    //TailoredPrecondDiag precond_old(S); // cannot be 'const' due to BiCGstabsolver interface
    const std::string solver_type = prm_.get<std::string>("solver.linsolver.solver");
    Opm::PropertyTree lprm = prm_.get_child("solver.linsolver.preconditioner");
    //Opm::FractureMechanicsPreconditioner precond(S, lprm);
    Dune::InverseOperatorResult iores;
    bool new_precond = first_iteration || closed_cells_changed;
    if(new_precond ){
        frac_flow_precond_.reset(new Opm::FractureMechanicsPreconditioner(S, lprm));
    }else{
        // need to modify code to use original storage for
        frac_flow_precond_->update(S,/*new_lu_mech*/ false);
    }
    // make possible use abs tolerance for linear solver
    double res = rhs.two_norm2();
    const double linsolve_tol = prm_.get<double>("solver.linsolver.tol");
    const double linsolve_atol = prm_.get<double>("solver.linsolver.atol");
    double lintol = std::max(linsolve_tol, linsolve_atol / res);
    
    const int max_iter = prm_.get<double>("solver.linsolver.max_iter");
    const int verbosity = prm_.get<double>("solver.linsolver.verbosity");
    if (nlin_verbosity > 1) {
      int num_closed_cells = std::accumulate(closed_cells.begin(), closed_cells.end(), 0);
      std::cout << "Nonlinear iteration: " << nlin_iteration;// << std::endl; 
      std::cout << " rhs:  " << rhs[_0].infinity_norm() << " " << rhs[_1].infinity_norm();
      std::cout <<  " closed_cells: " << num_closed_cells;
      std::cout << " cells m: " << rhs[_0].size();
      std::cout << " cells p: " << rhs[_1].size();
      std::cout << std::endl;
    }

    //std::unique_ptr<LinearSolverBase> psolver;
    auto& precond = *frac_flow_precond_;
    //using LinearSolverBase = Dune::InverseOperator<VectorHP, VectorHP>;
    if (new_precond) {
        psolver_ = setupLinearSolver(solver_type, *S_linop_, precond, prm_, lintol, max_iter, verbosity);
    }

    {
        OPM_TIMEBLOCK(SolveCoupledSystem);
        psolver_->apply(dx, rhs, iores); // NB: will modify 'rhs'
    }

    if (nlin_verbosity > 2) {
        std::cout << "x:  " << x[_0].infinity_norm() << " " << x[_1].infinity_norm() << std::endl;
        std::cout << "dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm() << std::endl;
        std::cout << "rhs after linsolve:  " << rhs[_0].infinity_norm() <<
          " " << rhs[_1].infinity_norm()
                  << std::endl;
    }
    if (nlin_verbosity > 2) {
      auto x_new(dx);
      //x_new += dx;
      auto rhs2_tmp(rhs_org);
      S.mmv(x_new, rhs2_tmp);
      std::cout << "Check res : rhs after linsolve:  " << rhs2_tmp[_0].infinity_norm() << 
        " " << rhs2_tmp[_1].infinity_norm() << std::endl;

      
    }
    
    // the following is a heuristic way to limit stepsize to stay within convergence radius
    const double damping = prm_.get<double>("solver.damping");
    const double step_fac = damping; // estimate_step_fac(x, dx) * damping;
    // std::cout << "fac: " << step_fac << std::endl;
    if (nlin_verbosity > 2) {
        std::cout << "fac: " << step_fac << std::endl;
    }
    dx *= step_fac;
    auto& dx0 = dx[_0];
    double max_dwidth = prm_.get<double>("solver.max_dwidth");
    // for(auto& dx0v: dx0){
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        dx0[i][0] = std::max(-max_dwidth, std::min(max_dwidth, dx0[i][0]));
    }
    auto& dx1 = dx[_1];
    double max_dp = prm_.get<double>("solver.max_dp");
    // for(auto& dx1v: dx1){
    for (size_t i = 0; i != fracture_pressure_.size(); ++i) {
        dx1[i][0] = std::max(-max_dp, std::min(max_dp, dx1[i][0]));
    }


    // dump_vector(dx, debug_filename("dx_w_").c_str(), debug_filename("dx_p_").c_str());
    // dump_vector(dx, "dx_w", "dx_p", true);
    x += dx;
    if (nlin_verbosity > 3) {
        std::cout << "after: dx: " << dx[_0].infinity_norm() << " " << dx[_1].infinity_norm()
                  << std::endl;
    }
    // copying modified variables back to member variables
    fracture_width_ = x[_0];
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        fracture_width_[i][0] = std::max(0.0, fracture_width_[i][0]); // ensure non-negativity
        // fracture_width_[i][0] = std::min(max_width, fracture_width_[i][0]); // ensure non-negativity
    }
    fracture_pressure_ = x[_1];
    // if (nlin_verbosity > 2) {
    //   VectorHP rhs_tmp {x}; // same size as system, content set below
    //   normalFractureTraction(rhs_tmp[_0], false); // right-hand side equals the normal fracture traction
    //   rhs_tmp[_1] = rhs_pressure_;
    //   for (size_t i = 0; i != closed_cells.size(); ++i)
    //     if (closed_cells[i])
    //         rhs_tmp[_0][i] = 0;
    //   S0.mmv(x, rhs_tmp);
    //   std::cout << "New: rhs after linsolve:  " << rhs_tmp[_0].infinity_norm() << 
    //     " " << rhs_tmp[_1].infinity_norm() << std::endl;
    // }

    
    return false;
}

}; // end namespace Opm
