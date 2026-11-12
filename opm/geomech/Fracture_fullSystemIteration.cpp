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

#include "config.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <numeric>

#include <dune/common/indices.hh> // needed for _0, _1, etc.
#include <dune/common/fmatrix.hh> // Dune::FieldMatrix
#include <dune/istl/bcrsmatrix.hh> // Dune::BCRSMatrix
#include <dune/istl/bvector.hh> // Dune::BlockVector
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>


#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/FractureCoupledLinearSolver.hpp>
#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

#include <fmt/format.h>

using namespace std;
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

using ResVector = Opm::Fracture::Vector;
using VectorHP = Dune::MultiTypeBlockVector<ResVector, ResVector>;

using SMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>; // sparse matrix
using FMatrix = Dune::DynamicMatrix<double>; // full matrix

using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FMatrix, SMatrix>,
                                                Dune::MultiTypeBlockVector<SMatrix, SMatrix>>;


void
validate_state_vector(const ResVector& values, const char* const label, const double max_abs_value)
{
    for (size_t i = 0; i < values.size(); ++i) {
        const double value = values[i][0];
        if (!std::isfinite(value)) {
            OPM_THROW(std::runtime_error,
                      std::string("Fracture nonlinear state contains non-finite ") + label
                          + " at index " + std::to_string(i));
        }
        if (std::abs(value) > max_abs_value) {
            OPM_THROW(std::runtime_error,
                      std::string("Fracture nonlinear state contains excessive ") + label
                          + " at index " + std::to_string(i)
                          + ": " + std::to_string(value));
        }
    }
}

void
validate_fracture_state(const VectorHP& x, const Opm::PropertyTree& prm)
{
    if (!prm.get<bool>("solver.guard_state", true)) {
        return;
    }
    const double max_abs_width = prm.get<double>("solver.guard_max_abs_width", 1e30);
    const double max_abs_pressure = prm.get<double>("solver.guard_max_abs_pressure", 1e30);
    validate_state_vector(x[_0], "width", max_abs_width);
    validate_state_vector(x[_1], "pressure", max_abs_pressure);
}

size_t
count_toggled_cells(const std::vector<int>& previous, const std::vector<int>& current)
{
    if (previous.size() != current.size()) {
        return current.size();
    }
    size_t count = 0;
    for (size_t i = 0; i < current.size(); ++i) {
        count += (previous[i] != current[i]);
    }
    return count;
}

// ============================ Debugging functions ============================
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
std::string
linear_system_dump_stem(const Opm::PropertyTree& prm,
                        const int nlin_iteration,
                        const std::string& reason)
// ----------------------------------------------------------------------------
{
    const std::string prefix = prm.get<std::string>("solver.linear_system_dump_prefix",
                                                    "fracture_linear_system");
    std::ostringstream oss;
    oss << prefix << "_nlin_" << nlin_iteration << "_" << reason;
    return oss.str();
}

// ----------------------------------------------------------------------------
std::string
linear_system_info_filename(const std::string& stem)
// ----------------------------------------------------------------------------
{
    const auto info_filename = std::filesystem::path(stem + "_info.json");
    return std::filesystem::absolute(info_filename).string();
}

// ----------------------------------------------------------------------------
void
write_property_tree_json(const Opm::PropertyTree& tree, const std::string& filename)
// ----------------------------------------------------------------------------
{
    std::ofstream os(filename);
    tree.write_json(os, true);
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
void
dump_linear_system_snapshot(const std::string& stem,
                            const SystemMatrix& system,
                            const VectorHP& rhs,
                            const VectorHP& x,
                            const std::vector<int>& closed_cells,
                            const Opm::PropertyTree& prm,
                            const Dune::InverseOperatorResult* iores = nullptr)
// ----------------------------------------------------------------------------
{
    dump_matrix(system[_0][_0], (stem + "_A.txt").c_str());
    Dune::storeMatrixMarket(system[_0][_1], stem + "_I.mtx");
    Dune::storeMatrixMarket(system[_1][_0], stem + "_C_active.mtx");
    Dune::storeMatrixMarket(system[_1][_1], stem + "_M_active.mtx");
    Dune::storeMatrixMarket(rhs[_0], stem + "_rhs_w.mtx");
    Dune::storeMatrixMarket(rhs[_1], stem + "_rhs_p.mtx");
    Dune::storeMatrixMarket(x[_0], stem + "_x_w.mtx");
    Dune::storeMatrixMarket(x[_1], stem + "_x_p.mtx");
    dump_vector(closed_cells, (stem + "_closed_cells.txt").c_str());

    Opm::PropertyTree info;
    info.put("mechanics_size", system[_0][_0].N());
    info.put("pressure_size", system[_1][_1].N());
    if (iores) {
        info.put("linear.converged", iores->converged);
        info.put("linear.iterations", iores->iterations);
    }

    write_property_tree_json(info, stem + "_info.json");
    write_property_tree_json(prm.get_child("solver.linsolver"), stem + "_linsolver.json");
    write_property_tree_json(prm.get_child("solver.linsolver.preconditioner.flow_solver"),
                             stem + "_flow_solver.json");
}

// ----------------------------------------------------------------------------
inline bool
convergence_test(const VectorHP& res, const double tol_flow, double tol_mech, int verbosity)
// ----------------------------------------------------------------------------
{
    bool converged = res[_0].infinity_norm() < tol_mech && res[_1].infinity_norm() < tol_flow;
    if(converged && verbosity > 0){
      Opm::OpmLog::debug(fmt::format("Converged with residual norm[0] is {}\nConverged with residual norm[1] is {}",
                                res[_0].infinity_norm(), res[_1].infinity_norm()));
    }
    return converged;
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
Fracture::identify_closed(const FMatrix& A, const VectorHP& x, const ResVector& rhs)
// ----------------------------------------------------------------------------
{

    OPM_TIMEFUNCTION();
    const InnerSolveOptions& opt = inner_;
    // Fischer-Burmeister handles contact through smooth mech rows, not a binary set.
    if (opt.isFb) {
        return std::vector<int>(A.N(), 0);
    }
    const auto was_closed = [this](const size_t index) {
        return index < closed_cells_.size() && closed_cells_[index] != 0;
    };
    const ResVector& h = x[_0];
    const ResVector force = contactForce(A, h, x[_1], rhs);
    std::vector<int> result;
    for (size_t i = 0; i != A.N(); ++i) {
        // A closed cell has its width pinned at 0, so reopening is decided by force alone.
        const bool close_candidate = force[i] >= opt.closeForceTol && h[i] <= 0.0;
        const bool reopen_candidate = force[i] <= -opt.reopenForceTol;
        result.push_back(was_closed(i) ? !reopen_candidate : close_candidate);
    }
    return result;
}

// ----------------------------------------------------------------------------
void
Fracture::assemblePressureAndCouplingAD(const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    updateLeakoff();
    const auto input = makePressureAssemblyInput();

    // Run the AD-based assembly (produces both pressure and coupling matrices)
    auto result = assemblePressureAD(input);

    // Store the pressure matrix
    pressure_matrix_ = std::move(result.pressure_matrix);

    // Store the coupling matrix
    coupling_matrix_ = std::move(result.coupling_matrix);

}

// ----------------------------------------------------------------------------
FracturePressureInput
Fracture::makePressureAssemblyInput() const
// ----------------------------------------------------------------------------
{
    const size_t nc = numFractureCells();
    const size_t nw = numWellEquations();
    const bool have_water_phase_props = static_cast<bool>(fracture_water_property_evaluator_);
    const bool use_fluid_system_water_properties =
        have_water_phase_props
        && prm_.get<bool>("solver.use_fluid_system_water_properties", false);

    FracturePressureInput input;
    input.num_cells = nc;
    input.min_width = min_width_;
    input.htrans = htrans_;
    input.control_type = effectiveControlType();
    if (input.control_type == "bhp_well" && !perfinj_.empty()) {
        const double density = reservoir_.density[std::get<0>(perfinj_[0])];
        input.well_target_pressure = well_target_bhp_
            + density * gravity_ * (perf_ref_depth_ - well_ref_depth_);
    }
    input.num_well_equations = nw;
    input.perfinj.assign(perfinj_.begin(), perfinj_.end());
    input.total_wellindex = total_wellindex_;
    input.mobility_water_perf = mobility_water_perf_;
    input.face_mobility.resize(nc + nw, 0.0);
    for (size_t i = 0; i < reservoir_.mobility.size() && i < nc; ++i)
        input.face_mobility[i] = reservoir_.mobility[i];
    input.leakof = leakof_;
    input.use_fluid_mobility_for_internal_faces =
        use_fluid_system_water_properties
        && prm_.get<bool>("solver.use_water_phase_mobility_for_fracture_flow", false);
    input.use_fluid_mobility_for_well_connections =
        use_fluid_system_water_properties
        && prm_.get<bool>("solver.use_water_phase_mobility_for_well_fracture_flow",
                          input.use_fluid_mobility_for_internal_faces);

    input.fracture_width.resize(nc + nw);
    input.fracture_pressure.resize(nc + nw);
    input.density.resize(nc + nw);
    input.viscosity.resize(nc + nw);

    for (size_t i = 0; i < nc + nw; ++i) {
        // fracture_width_ has no entry for well equations (rate_well)
        input.fracture_width[i] = (i < fracture_width_.size()) ? fracture_width_[i][0] : 0.0;
        input.fracture_pressure[i] = fracture_pressure_[i][0];
    }
    input.face_gravity = faceGravityHeads();

    for (size_t i = 0; i < nc; ++i) {
        if (use_fluid_system_water_properties) {
            const auto fluid_properties = fracture_water_property_evaluator_(i, fracture_pressure_[i][0]);
            input.density[i] = fluid_properties.first;
            input.viscosity[i] = fluid_properties.second;
        }
        else {
            input.density[i] = {reservoir_.density[i], 0.0};
            input.viscosity[i] = {1.0, 0.0};
        }
    }
    for (size_t i = nc; i < nc + nw; ++i) {
        input.density[i] = {1.0, 0.0};
        input.viscosity[i] = {1.0, 0.0};
    }

    return input;
}

// ----------------------------------------------------------------------------
Fracture::InnerSolveOptions
Fracture::InnerSolveOptions::read(const PropertyTree& prm)
// ----------------------------------------------------------------------------
{
    InnerSolveOptions o;
    o.verbosity = prm.get<int>("solver.verbosity", 0);
    const auto policy = prm.get<std::string>("solver.closed_cell_policy", "sticky");
    if (policy != "sticky" && policy != "fischer_burmeister") {
        OPM_THROW(std::runtime_error, "Unknown closed cell policy: " + policy);
    }
    o.isFb = (policy == "fischer_burmeister");
    o.closeForceTol = prm.get<double>("solver.close_force_tolerance", 0.0);
    o.reopenForceTol = prm.get<double>("solver.reopen_force_tolerance", o.closeForceTol);
    if (o.isFb) {
        o.fbScaledConvergence = prm.get<bool>("solver.fb_scaled_convergence", false);
        o.fbLineSearch = prm.get<bool>("solver.fb_line_search", false);
        if (o.fbLineSearch) {
            o.fbLineSearchMin = prm.get<double>("solver.fb_line_search_min", 0.0625);
        }
    }
    o.dropCoupling = prm.get<bool>("solver.drop_fluid_mech_linearization", true);
    o.storageTerm = prm.get<bool>("solver.fracture_storage_term", false);
    o.linear = FractureCoupledLinearSolver::Options::read(prm);
    o.dumpMinNonlinIter = prm.get<int>("solver.dump_case_min_nonlin_iter", 0);
    o.writeCoupledSystem = prm.get<bool>("solver.write_coupled_linear_system", false);
    o.dumpOnFailure = prm.get<bool>("solver.dump_linear_system_on_failure", true);
    o.maxDwidth = prm.get<double>("solver.max_dwidth");
    return o;
}

// ----------------------------------------------------------------------------
Fracture::ResVector
Fracture::contactForce(const FMatrix& A, const ResVector& w, const ResVector& p, const ResVector& traction)
// ----------------------------------------------------------------------------
{
    // traction - A w - p on the fracture cells; p may carry trailing well rows
    ResVector force(traction);
    A.mmv(w, force);
    for (size_t i = 0; i != A.N(); ++i) {
        force[i] -= p[i];
    }
    return force;
}

// ----------------------------------------------------------------------------
void
Fracture::updateFlowCoefficients_()
// ----------------------------------------------------------------------------
{
    updateLeakoff();
    addSource();
}

// ----------------------------------------------------------------------------
Fracture::ContactUpdate
Fracture::updateContactSet_(const VectorHP& x, const ResVector& traction, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    std::vector<int> closed_cells = identify_closed(fractureMatrix(), x, traction);
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
    // Contact-chatter metric: cells that flipped since the previous iteration
    // (same grid only; a size mismatch is a grid change).
    if (closed_cells_.size() == closed_cells.size()) {
        last_solve_stats_.closed_cell_toggles +=
            static_cast<int>(count_toggled_cells(closed_cells_, closed_cells));
        // When each cell last toggled, for the propagation veto: sticky chatters by
        // nature, so what matters is whether the cell settled, not how often.
        if (cell_last_toggle_iter_.size() != closed_cells.size())
            cell_last_toggle_iter_.assign(closed_cells.size(), -1000);
        for (size_t i = 0; i < closed_cells.size(); ++i)
            if (closed_cells_[i] != closed_cells[i])
                cell_last_toggle_iter_[i] = nlin_iteration;
    } else {
        cell_last_toggle_iter_.assign(closed_cells.size(), -1000);
    }
    current_solve_iter_ = nlin_iteration;
    closed_cells_ = closed_cells;
    return {std::move(closed_cells), closed_cells_changed};
}

// ----------------------------------------------------------------------------
Fracture::FbRows
Fracture::assembleMechanicsRows_(const VectorHP& x,
                                 ResVector& mech_rhs,
                                 const std::vector<int>& closed_cells,
                                 const bool first_iteration)
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    FbRows fb;
    // Fischer-Burmeister contact replaces the binary open/closed mech rows by the
    // smooth complementarity row phi_FB(c_i w_i, lambda_i) = 0 with lambda the
    // contact force and c_i = |A_ii|. Jacobian row: A_FB = da c_i e_i - db A_i,
    // I_FB_ii = -db, with da = 1 - a/r, db = 1 - b/r, r = sqrt(a^2+b^2).
    if (opt.isFb) {
        const FMatrix& Aorig = fractureMatrix();
        const size_t n = Aorig.N();
        const ResVector& w = x[_0];
        const ResVector& p = x[_1];
        fb.phi.assign(n, 0.0);
        fb.scale.assign(n, 1.0);
        fb_cell_residual_.assign(n, 0.0);
        if (first_iteration) A_.reset(new FMatrix(Aorig));
        FMatrix& Afb = *A_;
        if (first_iteration) I_.reset(new SMatrix(makeIdentity(n, numWellEquations())));
        SMatrix& Ifb = *I_;
        const double eps = 1e-30;
        const ResVector force = contactForce(Aorig, w, p, mech_rhs);
        for (size_t i = 0; i < n; ++i) {
            const double c_i = std::max(std::abs(Aorig[i][i]), eps);
            const double a = c_i * w[i][0];
            const double b = force[i][0];
            const double r = std::sqrt(a * a + b * b) + eps;
            const double da = 1.0 - a / r;
            const double db = 1.0 - b / r;
            fb.phi[i] = a + b - std::sqrt(a * a + b * b);
            fb.scale[i] = 1.0 + std::abs(a) + std::abs(b);
            fb_cell_residual_[i] = std::abs(fb.phi[i]) / fb.scale[i];
            for (size_t j = 0; j < n; ++j) Afb[i][j] = -db * Aorig[i][j];
            Afb[i][i] += da * c_i;
            Ifb[i][i] = -db;
        }
    } else {
        if (first_iteration) {
            A_.reset(new FMatrix(modified_fracture_matrix(fractureMatrix(), closed_cells)));
        }else{
            modified_fracture_matrix(*A_, fractureMatrix(), closed_cells);
        }
        // Closed cells are held fixed in the mechanics residual and Jacobian block.
        for (size_t i = 0; i != closed_cells.size(); ++i)
            if (closed_cells[i])
                mech_rhs[i] = 0;
    }
    return fb;
}

// ----------------------------------------------------------------------------
void
Fracture::assembleFlowRows_(const std::vector<int>& closed_cells)
// ----------------------------------------------------------------------------
{
    assemblePressureAndCouplingAD(closed_cells);
    if (inner_.dropCoupling) {
        *coupling_matrix_ = 0;
    }
    // External-pressure mode: pressure rows become p = p_external (identity,
    // zero residual, no width coupling), so the Newton reduces to mechanics
    // and contact at the pressure the flow's aux cells hold.
    // Cells without an external value (no flow DOF yet) keep their own pressure
    // equation, so a front that grew inside the step can still be pressurised.
    if (external_pressure_) {
        const std::size_t nc = numFractureCells();
        const bool haveMask = (external_cell_mask_.size() == nc);
        auto& Mx = *pressure_matrix_;
        auto& Cx = *coupling_matrix_;
        for (auto row = Mx.begin(); row != Mx.end(); ++row) {
            const std::size_t r = row.index();
            const bool pinned = (r >= nc) || !haveMask || external_cell_mask_[r];
            if (!pinned) {
                continue;
            }
            for (auto col = row->begin(); col != row->end(); ++col) {
                *col = (col.index() == r) ? 1.0 : 0.0;
            }
            if (r < Cx.N()) {
                for (auto& entry : Cx[r]) {
                    entry = 0.0;
                }
            }
            rhs_pressure_[r] = fracture_pressure_[r];
        }
    }
}

// ----------------------------------------------------------------------------
void
Fracture::assembleSystemMatrix_(const std::vector<int>& closed_cells, const bool first_iteration)
// ----------------------------------------------------------------------------
{
    const auto& A = *A_;
    const auto& M = *pressure_matrix_;
    const auto& C = *coupling_matrix_;
    // Fischer-Burmeister already built I_ with its smooth -db_i weights.
    if (!inner_.isFb) {
        if (first_iteration) {
            I_.reset(new SMatrix(makeIdentity(A.N(), numWellEquations(), 1, closed_cells)));
        }else{
            modifyIdentity(*I_, 1, closed_cells);
        }
    }
    const auto& I = *I_;
    if(first_iteration){
        S_.reset(new SystemMatrix({{A, I},
                                   {C, M}}));
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
}

// ----------------------------------------------------------------------------
Fracture::VectorHP
Fracture::assembleResidual_(const VectorHP& x, const ResVector& mech_rhs, const FbRows& fb) const
// ----------------------------------------------------------------------------
{
    const auto& S = *S_;
    const auto& M = *pressure_matrix_;
    VectorHP rhs {x};
    rhs[_0] = mech_rhs;
    rhs[_1] = rhs_pressure_;

    {
        // A block whose column count disagrees with the vector it multiplies
        // reads past the end of that vector: silent garbage in the residual and
        // in every Krylov matrix-vector product. Fail by name instead.
        const size_t nm = rhs[_0].size(), np = rhs[_1].size();
        if (x[_0].size() != nm || x[_1].size() != np)
            OPM_THROW(std::runtime_error,
                      "Fracture coupled solve: iterate/residual size mismatch: dx="
                          + std::to_string(x[_0].size()) + "/" + std::to_string(x[_1].size())
                          + " rhs=" + std::to_string(nm) + "/" + std::to_string(np));
        auto bad = [](const auto& B, size_t n, size_t m) { return B.N() != n || B.M() != m; };
        if (bad(S[_0][_0], nm, nm) || bad(S[_0][_1], nm, np) || bad(S[_1][_0], np, nm)
            || bad(S[_1][_1], np, np))
            OPM_THROW(std::runtime_error,
                      "Fracture coupled system block dimensions inconsistent: nm=" + std::to_string(nm)
                          + " np=" + std::to_string(np) + " A=" + std::to_string(S[_0][_0].N()) + "x"
                          + std::to_string(S[_0][_0].M()) + " I=" + std::to_string(S[_0][_1].N()) + "x"
                          + std::to_string(S[_0][_1].M()) + " C=" + std::to_string(S[_1][_0].N()) + "x"
                          + std::to_string(S[_1][_0].M()) + " M=" + std::to_string(S[_1][_1].N()) + "x"
                          + std::to_string(S[_1][_1].M()));
    }

    S[_1][_1].mmv(x[_1], rhs[_1]);
    S[_0][_0].mmv(x[_0], rhs[_0]);
    S[_0][_1].mmv(x[_1], rhs[_0]);

    // Fracture volume storage A_i (w_i - w_i^n)/dt: dominant while a cell opens
    // fast, where it makes opening cost pressure. Cells grown inside the step
    // store their full opening; assumes no renumbering within the step.
    if (!external_pressure_ && inner_.storageTerm && current_dt_ > 0.0) {
        const std::vector<double> areas = cellAreas();
        const size_t ncell = std::min(areas.size(), numFractureCells());
        for (size_t i = 0; i < ncell && i < rhs[_1].size() && i < x[_0].size(); ++i) {
            const double w_prev = (i < fracture_width_prev_.size())
                ? fracture_width_prev_[i][0] : 0.0;
            rhs[_1][i][0] -= areas[i] * (x[_0][i][0] - w_prev) / current_dt_;
        }
    }

    // Fischer-Burmeister: the mechanics residual is the complementarity function.
    if (inner_.isFb)
        for (size_t i = 0; i < fb.phi.size(); ++i)
            rhs[_0][i] = -fb.phi[i];

    // Fail recoverably (upstream catch -> timestep chop) instead of aborting GMRES
    // with an opaque "defect=nan"; a non-finite residual means S or rhs is bad.
    for (size_t i = 0; i < rhs[_0].size(); ++i)
        if (!std::isfinite(rhs[_0][i][0]))
            OPM_THROW(std::runtime_error,
                      "Fracture full-system residual non-finite in mechanics row "
                          + std::to_string(i));
    for (size_t i = 0; i < rhs[_1].size(); ++i)
        if (!std::isfinite(rhs[_1][i][0]))
            OPM_THROW(std::runtime_error,
                      "Fracture full-system residual non-finite in pressure row "
                          + std::to_string(i));
    if (numWellEquations() > 0) {
        // a zero well-equation pivot passes the finite checks but makes the
        // preconditioner emit NaN mid-Krylov ("defect=nan")
        const auto wix = M.N() - 1;
        const double wdiag = M[wix][wix][0][0];
        if (!std::isfinite(wdiag) || wdiag == 0.0)
            OPM_THROW(std::runtime_error,
                      "rate_well equation diagonal is " + std::to_string(wdiag)
                          + " (total_wellindex " + std::to_string(total_wellindex_)
                          + ", well_rate " + std::to_string(well_rate_) + ")");
    }
    return rhs;
}

// ----------------------------------------------------------------------------
bool
Fracture::isConverged_(const VectorHP& rhs, const FbRows& fb, const double tol) const
// ----------------------------------------------------------------------------
{
    // The same tolerance for both blocks, on the unscaled residual.
    const double tol_flow = tol;
    const double tol_mech = tol;
    // fb_scaled_convergence tests phi_i/(1+|a_i|+|b_i|), the same accuracy on every cell.
    if (inner_.isFb && inner_.fbScaledConvergence) {
        VectorHP rtest(rhs);
        for (size_t i = 0; i < fb.scale.size() && i < rtest[_0].size(); ++i)
            rtest[_0][i][0] /= fb.scale[i];
        return convergence_test(rtest, tol_mech, tol_flow, inner_.verbosity);
    }
    return convergence_test(rhs, tol_mech, tol_flow, inner_.verbosity);
}

// ----------------------------------------------------------------------------
bool
Fracture::fullSystemIteration(const double tol, const int nlin_iteration)
// ----------------------------------------------------------------------------
{
    OPM_TIMEFUNCTION();
    ++last_solve_stats_.nonlinear_iterations;

    ++DEBUG_COUNT;

    const InnerSolveOptions& opt = inner_;
    const int nlin_verbosity = opt.verbosity;
    updateFlowCoefficients_();

    VectorHP x {fracture_width_, fracture_pressure_};
    validate_fracture_state(x, prm_);

    ResVector mech_rhs(fracture_width_);
    normalFractureTraction(mech_rhs, false);

    const ContactUpdate contact = updateContactSet_(x, mech_rhs, nlin_iteration);
    const std::vector<int>& closed_cells = contact.closed;
    const bool closed_cells_changed = contact.changed;
    const bool first_iteration = (nlin_iteration == 0);
    const FbRows fb = assembleMechanicsRows_(x, mech_rhs, closed_cells, first_iteration);
    assembleFlowRows_(closed_cells);
    assembleSystemMatrix_(closed_cells, first_iteration);
    auto& S = *S_;

    VectorHP dx = x;
    dx = 0; // Newton update for the current nonlinear iterate
    VectorHP rhs = assembleResidual_(x, mech_rhs, fb);
    auto rhs_org(rhs);

    if (isConverged_(rhs, fb, tol))
        return true;

    coupledSolver_.prepare(S, rhs);
    FractureCoupledLinearSolver::Iteration it;
    it.index = nlin_iteration;
    it.verbosity = nlin_verbosity;
    it.numClosed = std::accumulate(closed_cells.begin(), closed_cells.end(), 0);
    it.newSystem = first_iteration;
    it.contactChanged = closed_cells_changed;
    // dump_case_min_nonlin_iter > 0 captures the systems of struggling solves for replay
    it.dump = opt.writeCoupledSystem
        || (opt.dumpMinNonlinIter > 0 && nlin_iteration >= opt.dumpMinNonlinIter);
    const auto dump = [&](const SystemMatrix& system, const char* tag,
                          const Dune::InverseOperatorResult* result) -> std::string {
        const bool failure = (result != nullptr);
        if (failure && !opt.dumpOnFailure) {
            return {};
        }
        const auto stem = linear_system_dump_stem(prm_, nlin_iteration, tag);
        dump_linear_system_snapshot(stem, system, rhs_org, x, closed_cells, prm_, result);
        if (!failure) {
            return {};
        }
        const auto info = linear_system_info_filename(stem);
        OpmLog::error(fmt::format("Fracture coupled linear solver failed. Snapshot info written to {}", info));
        return info;
    };
    coupledSolver_.solve(rhs, dx, it, last_solve_stats_, dump);

    if (nlin_verbosity > 2) {
        OpmLog::debug(fmt::format("x:  {} {}\ndx: {} {}\nrhs after linsolve:  {} {}",
                                  x[_0].infinity_norm(), x[_1].infinity_norm(),
                                  dx[_0].infinity_norm(), dx[_1].infinity_norm(),
                                  rhs[_0].infinity_norm(), rhs[_1].infinity_norm()));
    }
    if (nlin_verbosity > 2 && !coupledSolver_.scaled()) {
      // S holds D S D when scaled, so this unscaled check would be meaningless.
      auto x_new(dx);
      auto rhs2_tmp(rhs_org);
      S.mmv(x_new, rhs2_tmp);
      OpmLog::debug(fmt::format("Check res : rhs after linsolve:  {} {}",
                                rhs2_tmp[_0].infinity_norm(), rhs2_tmp[_1].infinity_norm()));


    }


    limitStep_(x, dx, mech_rhs, fb);
    applyUpdate_(x, dx);
    return false;
}

// ----------------------------------------------------------------------------
void
Fracture::limitStep_(const VectorHP& x, VectorHP& dx, const ResVector& mech_rhs, const FbRows& fb) const
// ----------------------------------------------------------------------------
{
    const InnerSolveOptions& opt = inner_;
    const int nlin_verbosity = opt.verbosity;
    // FB backtracking: far from the solution the full semismooth Newton step
    // overshoots; damp it until the merit ||phi||^2 decreases.
    if (opt.isFb && opt.fbLineSearch) {
        const FMatrix& Aorig = fractureMatrix();
        const size_t n = Aorig.N();
        double merit0 = 0.0;
        for (size_t i = 0; i < n; ++i) merit0 += fb.phi[i] * fb.phi[i];
        const double s_min = opt.fbLineSearchMin;
        double sstep = 1.0, best_s = 1.0;
        double best_merit = std::numeric_limits<double>::max();
        ResVector w_t(n), p_t(n);
        while (sstep >= s_min) {
            for (size_t j = 0; j < n; ++j)
                w_t[j] = std::max(0.0, x[_0][j][0] + sstep * dx[_0][j][0]);
            for (size_t i = 0; i < n; ++i)
                p_t[i] = x[_1][i][0] + sstep * dx[_1][i][0];
            const ResVector lam = contactForce(Aorig, w_t, p_t, mech_rhs);
            double merit = 0.0;
            for (size_t i = 0; i < n; ++i) {
                const double c_i = std::max(std::abs(Aorig[i][i]), 1e-30);
                const double a = c_i * w_t[i][0];
                const double b = lam[i][0];
                const double phi = a + b - std::sqrt(a * a + b * b);
                merit += phi * phi;
            }
            if (merit < best_merit) { best_merit = merit; best_s = sstep; }
            if (merit <= (1.0 - 1e-4 * sstep) * merit0) break; // sufficient decrease
            sstep *= 0.5;
        }
        dx *= best_s;
        if (nlin_verbosity > 1 && best_s < 1.0)
            OpmLog::debug(fmt::format("FB line search: step {} merit {} -> {}", best_s, merit0, best_merit));
    }
    auto& dx0 = dx[_0];
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        dx0[i][0] = std::max(-opt.maxDwidth, std::min(opt.maxDwidth, dx0[i][0]));
    }
}

// ----------------------------------------------------------------------------
void
Fracture::applyUpdate_(VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------
{
    x += dx;
    if (inner_.verbosity > 3) {
        OpmLog::debug(fmt::format("after: dx: {} {}", dx[_0].infinity_norm(), dx[_1].infinity_norm()));
    }
    fracture_width_ = x[_0];
    for (size_t i = 0; i != fracture_width_.size(); ++i) {
        fracture_width_[i][0] = std::max(0.0, fracture_width_[i][0]);
    }
    fracture_pressure_ = x[_1];
}

}; // end namespace Opm
