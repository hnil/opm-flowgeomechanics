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

#ifndef OPM_VEM_ELASTICITY_SOLVER_HPP
#define OPM_VEM_ELASTICITY_SOLVER_HPP

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <opm/grid/CpGrid.hpp>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Schedule/BCState.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/geomech/MechLinePlanePreconditioner.hpp>
#include <opm/geomech/DuneCommunicationHelpers.hpp>
#include <opm/geomech/MechanicalPartitionValidation.hpp>
#include <opm/geomech/vem/vem.hpp>

#include <memory>
#include <tuple>
#include <vector>

namespace Opm {
namespace Elasticity {

using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
using VoigtField = Dune::BlockVector<Dune::FieldVector<double,6>>;

#if HAVE_MPI
using MechCommunication = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
using MechCommunication = Dune::Communication<int>;
#endif

//! Outcome of one mechanics linear solve.
struct MechSolveResult
{
    bool converged = true;
    int iterations = 0;
    double reduction = 1.0;
    bool failed = false; // the solver threw; the solution was not updated
};

//! Displacement (all node dofs, or the free ones with reduce_boundary) and cell fields.
struct MechState
{
    Vector displacement;
    VoigtField stress;
    VoigtField strain;
};

//! VEM operators on a fixed grid with fixed Dirichlet nodes: stiffness matrix,
//! load-to-rhs (divmat_), and displacement-to-stress/strain matrices.
template<class GridType>
class VemDiscretisation
{
public:
    explicit VemDiscretisation(const GridType& grid)
      : grid_(grid)
    {}

    void setMaterial(const std::vector<double>& ymodule, const std::vector<double>& pratio)
    {
        ymodule_ = ymodule;
        pratio_ = pratio;
    }

    void setOptions(int stability_choice_int,
                    bool vem_source,
                    bool stab_on_stress,
                    bool vem_stress,
                    bool diagonal_scaling,
                    bool reduce_boundary);

    void fixNodes(const std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes);
    void setBodyForce(double gravity, const std::vector<double>& density);

    //! Stiffness matrix and the fixed-load rhs (body force, Dirichlet values).
    void assembleMatrix();
    //! divmat_, stressmat_, strainmat_; independent of the body force.
    void buildOperators();
    bool hasOperators() const { return divmat_.N() > 0; }

    //! b = divmat * load + rhs_force
    void rhs(const Vector& load, Vector& b) const;
    void expandSolution(Vector& result, const Vector& u) const;
    void stress(const Vector& dispall, VoigtField& out) const;
    void strain(const Vector& dispall, VoigtField& out) const;

    Matrix& matrix() { return matrix_; }
    const Matrix& matrix() const { return matrix_; }
    //! Jacobi factors when the matrix is S A S (diagonal_scaling), else empty.
    const Vector& scale() const { return scale_; }
    const std::vector<double>& coords() const { return coords_; }
    bool reduceBoundary() const { return reduce_boundary_; }

private:
    void buildGeometry();
    void applyDiagonalScaling();
    static void makeDuneMatrixCompressed(const std::vector<std::tuple<int, int, double>>& entries,
                                         Matrix& mat);

    const GridType& grid_;
    int num_cells_ = 0;
    std::vector<double> coords_;
    std::vector<int> num_cell_faces_, num_face_corners_, face_corners_;
    std::vector<char> corner_mask_; // per face_corners_ entry: a grid corner of its face
    std::vector<int> idx_free_;
    std::vector<double> ymodule_;
    std::vector<double> pratio_;
    std::vector<double> body_force_;
    std::tuple<int, std::vector<int>, std::vector<double>> dirichlet_;
    std::vector<double> rhs_force_;

    Matrix matrix_;
    Matrix divmat_;    // from cell load to dofs
    Matrix stressmat_; // from all dofs (no elimination) to cell stress
    Matrix strainmat_;
    Vector scale_;

    vem::StabilityChoice stability_choice_ = vem::SIMPLE;
    bool vem_source_ = true;
    bool stab_on_stress_ = false;
    bool vem_stress_ = true;
    bool diagonal_scaling_ = false;
    bool reduce_boundary_ = false;
};

//! Krylov solve of the assembled mechanics system (FlexibleSolver), warm started
//! from the previous displacement.
template<class GridType>
class MechLinearSolver
{
public:
    explicit MechLinearSolver(const GridType& grid)
      : grid_(grid)
    {}

    //! Sets up the preconditioner on A; in parallel, A's copy rows become identity rows.
    void setup(const PropertyTree& prm, Matrix& A, const std::vector<double>& coords);
    //! x is the unscaled solution and the warm start; scale is empty or the Jacobi factors of A.
    MechSolveResult solve(const Vector& b, const Vector& scale, Vector& x);
    const MechCommunication* comm() const { return comm_.get(); }

private:
    void setCopyRowsToZero(Matrix& A);

    using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    using ParOperatorType = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, MechCommunication>;

    const GridType& grid_;
    const Matrix* A_ = nullptr;
    std::unique_ptr<SeqOperatorType> sop_;
    std::unique_ptr<ParOperatorType> pop_;
    std::unique_ptr<Dune::InverseOperator<Vector, Vector>> tsolver_;
    std::shared_ptr<MechCommunication> comm_;
    std::shared_ptr<MechCommunication::ParallelIndexSet> vertexParallelIndexSet_;
    std::shared_ptr<MechCommunication::ParallelIndexSet> dofParallelIndexSet_;
    bool tol_relative_to_load_ = false;
    bool warm_start_ = true;
    double linear_tol_ = 1e-5;
};

} // namespace Elasticity
} // namespace Opm

#include "VemElasticitySolver_impl.hpp"

#endif
