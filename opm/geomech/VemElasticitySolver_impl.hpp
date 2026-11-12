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

//==============================================================================
//!
//! \file elasticity_upscale_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class - template implementations
//!
//==============================================================================
#ifndef OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP
#define OPM_VEM_ELASTICITY_SOLVER_IMPL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <opm/common/TimingMacros.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>
namespace Opm {
namespace Elasticity {

#undef IMPL_FUNC
#define IMPL_FUNC(A,B) template<class GridType> \
                         A VemElasticitySolver<GridType>::B

  IMPL_FUNC(void,assignToVoigt(Dune::BlockVector< Dune::FieldVector<double,6> >& voigt_stress, const Dune::BlockVector< Dune::FieldVector<double,1> >& vemstress)){
    // This is used in flow!!!
    for (size_t i = 0; i < voigt_stress.size(); ++i) {
                for (size_t k = 0; k < 3; ++k) {
		  voigt_stress[i][k] = vemstress[i * 6 + k];
                }
		// from VEM shear order xy, yz, zx to Voigt yz, zx, xy
		voigt_stress[i][3] = vemstress[i * 6 + 4];
		voigt_stress[i][5] = vemstress[i * 6 + 3];
		voigt_stress[i][4] = vemstress[i * 6 + 5];
    }
  }
  IMPL_FUNC(void,assignToVoigtSymMat(Dune::BlockVector< Dune::FieldVector<double,6> >& voigt_stress, const std::vector< std::array<double,6> >& vemstress)){
    for (size_t i = 0; i < voigt_stress.size(); ++i) {
      for (size_t k = 0; k < 3; ++k) {
	voigt_stress[i][k] = vemstress[i][k];
      }
      // from VEM shear order xy, yz, zx to Voigt yz, zx, xy
      voigt_stress[i][3] = vemstress[i][4];
      voigt_stress[i][5] = vemstress[i][3];
      voigt_stress[i][4] = vemstress[i][5];
    }
  }

  IMPL_FUNC(void, calculateStressPrecomputed(const Vector& dispalldune))
    {

            OPM_TIMEBLOCK(calculateStressPrecomputed);
            //NB stressmat is defined in linear indices not block linear indices
            if(vem_stress_){
              Vector stress(6 * grid_.leafGridView().size(0));
              stressmat_.mv(dispalldune,stress);
              stress_.resize(num_cells_);

              assignToVoigt(stress_,stress);
            }else{
              stress_ = vem::computeStressFem(grid_,
                                dispalldune,
                                ymodule_,
                                pratio_);
            }

    }
    IMPL_FUNC(void, calculateStress())
    {
        bool expand = true;
            OPM_TIMEBLOCK(calculateStressFull);
            // assumes the grid structure is made
            num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
            // assemble the mechanical system
            stress_.resize(num_cells_);
            std::vector<std::array<double, 6>> stress;
            stress.resize(num_cells_);
            std::vector<double> dispall;
            this->expandDisp(dispall,expand);
            std::vector<std::tuple<int, int, double>> stressmat;
            vem::compute_stress_3D(&coords_[0],
                                   num_cells_,
                                   &num_cell_faces_[0],
                                   &num_face_corners_[0],
                                   &face_corners_[0],
                                   &ymodule_[0],
                                   &pratio_[0],
                                   dispall,
                                   stress,
                                   stressmat,
                                   false,
                                   true,
                                   stability_choice_,
                                   stab_on_stress_, corner_mask_.data()
                );
            // copy to dune definitions
            stress_.resize(num_cells_);
	    assignToVoigtSymMat(stress_,stress);

    }

    IMPL_FUNC(void, calculateStrainPrecomputed(const Vector& dispalldune))
    {

            OPM_TIMEBLOCK(calculateStrainPrecomputed);
            //NB stressmat is defined in linear indices not block linear indices
            Vector strain(6 * grid_.leafGridView().size(0));
            strainmat_.mv(dispalldune,strain);
            strain_.resize(num_cells_);
	    assignToVoigt(strain_,strain);
    }
    IMPL_FUNC(void, calculateStrain())
    {
        bool expand = true;// assumes that matrices is caculated with reduced_boundary = true
            OPM_TIMEBLOCK(calculateStressFull);
            // assumes the grid structure is made
            num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
            // assemble the mechanical system
            strain_.resize(num_cells_);
            std::vector<std::array<double, 6>> strain;
            strain.resize(num_cells_);

            std::vector<double> dispall;
            this->expandDisp(dispall,expand);
            std::vector<std::tuple<int, int, double>> stressmat;
            vem::compute_stress_3D(&coords_[0],
                                   num_cells_,
                                   &num_cell_faces_[0],
                                   &num_face_corners_[0],
                                   &face_corners_[0],
                                   &ymodule_[0],
                                   &pratio_[0],
                                   dispall,
                                   strain,
                                   stressmat,
                                   false,
                                   false,
                                   stability_choice_,
                                   stab_on_stress_, corner_mask_.data()

                );
            // copy to dune definitions
            strain_.resize(num_cells_);
	        assignToVoigtSymMat(strain_,strain);

    }

    IMPL_FUNC(void, assemble(const Vector& pressure, bool reduce_boundary))
{
    OPM_TIMEBLOCK(assemble);
    using namespace std;
    this->getLoadVector() = 0;
    this->getOperator() = 0;
    {
        vem::getGridVectors(grid_,coords_,
                            num_cell_faces_,
                            num_face_corners_,
                            face_corners_,
                            &corner_mask_);

        const int num_fixed_dofs = get<0>(dirichlet_);
        const vector<int>& fixed_dof_ixs = std::get<1>(dirichlet_);
        const vector<double>& fixed_dof_values = std::get<2>(dirichlet_);

        // neumann boundary conditions
        const int num_neumann_faces = 0;
        num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
        // assemble the mechanical system
        {// scope for timing and removing A_entries util assembled directly in dune matrices
        OPM_TIMEBLOCK(assembleVEMSystem);
        vector<tuple<int, int, double>> A_entries;

        vem::assemble_mech_system_3D_dune(grid_, &coords_[0], num_cells_, &num_cell_faces_[0], &num_face_corners_[0],
                       &face_corners_[0], &ymodule_[0], &pratio_[0], &body_force_[0],
                       num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                       num_neumann_faces, nullptr, nullptr,
                       A_entries, rhs_force_, stability_choice_,reduce_boundary, corner_mask_.data());


        this->makeDuneSystemMatrix(A_entries);
        if(diagonal_scaling_){
            this->applyDiagonalScaling();
        }
        }
        // hack
        if(divmat_.N() > 0){// not rebuilding divmat_, stressmat_,strainmat_ for now
            return;
        }
        {
        OPM_TIMEBLOCK(setUpExtraStructuresForVEM);
        // make indexing for div operator i.e. all nodes to dofs
        std::vector<int> dof_idx(grid_.leafGridView().size(3)*3);
        std::iota(dof_idx.begin(), dof_idx.end(),0);

        std::set_difference(dof_idx.begin(), dof_idx.end(), fixed_dof_ixs.begin(), fixed_dof_ixs.end(),std::back_inserter(idx_free_));
        vector<double> rhs_tmp(pressure.size(),0);
        vector<double> pressure_tmp(pressure.size(),0);
        vector<tuple<int, int, double>> divmat;
        if(vem_source_){
          vem::potential_gradient_force_3D(&coords_[0],
                                         num_cells_,
                                         &num_cell_faces_[0],
                                         &num_face_corners_[0],
                                         &face_corners_[0],
                                         &pressure_tmp[0],
                                         rhs_tmp,
                                         divmat,
                                         true, corner_mask_.data());
        }else{
          vem::potential_gradient_force_3D_dune(grid_,&coords_[0],
                                             num_cells_,
                                             &num_cell_faces_[0],
                                             &num_face_corners_[0],
                                             &face_corners_[0],
                                             &pressure_tmp[0],
                                             rhs_tmp,
                                             divmat,
                                             true, corner_mask_.data());
        }
        std::vector<int> global_to_dof(grid_.leafGridView().size(3)*3,-1);
        for(size_t i=0; i< idx_free_.size(); ++i){
            if(reduce_boundary){
                //new numbering of dofs
                global_to_dof[idx_free_[i]] = i;
            }else{
                global_to_dof[idx_free_[i]] = idx_free_[i];
            }
        }
        //renumber and eliminate dof fized
        vector<tuple<int, int, double>> divmatdof;
        for(const auto& elem: divmat){
            int I = global_to_dof[std::get<0>(elem)];
            if(I>-1){
                int J = get<1>(elem);
                double val = std::get<2>(elem);
                divmatdof.push_back(std::tuple<int, int, double>(I,J,val));
            }
        }
        // finaly make dune matrix
        divmat_.setBuildMode(Matrix::random);
        // map from dof=3*nodes at a cell (ca 3*3*3) to cell
        {
            const int rows = reduce_boundary ? int(idx_free_.size())
                                             : grid_.leafGridView().size(3)*3;
            divmat_.setSize(rows, num_cells_);
        }
        makeDuneMatrixCompressed(divmatdof,divmat_);
        }
        std::vector<double> dispall(grid_.leafGridView().size(3)*3);
        std::vector<std::array<double,6>> stresstmp(grid_.leafGridView().size(0));
        {
        OPM_TIMEBLOCK(setUpStressStressMatrix);
        // also make stress matrix
        std::vector<std::tuple<int, int, double>> stressmat;
        vem::compute_stress_3D(&coords_[0],
                               num_cells_,
                               &num_cell_faces_[0],
                               &num_face_corners_[0],
                               &face_corners_[0],
                               &ymodule_[0], &pratio_[0],
                               dispall,
                               stresstmp,
                               stressmat,
                               true,
                               true,
                               stability_choice_,
                               stab_on_stress_, corner_mask_.data()
        );
        stressmat_.setBuildMode(Matrix::random);
        stressmat_.setSize(num_cells_*6, dispall.size());
        makeDuneMatrixCompressed(stressmat, stressmat_);
        }
        {
        OPM_TIMEBLOCK(setUpStrainMatrix);
        std::vector<std::tuple<int, int, double>> strainmat;
        vem::compute_stress_3D(&coords_[0],
                               num_cells_,
                               &num_cell_faces_[0],
                               &num_face_corners_[0],
                               &face_corners_[0],
                               &ymodule_[0], &pratio_[0],
                               dispall,
                               stresstmp,
                               strainmat,
                               true,
                               false,
                               stability_choice_,
                               stab_on_stress_, corner_mask_.data()
        );
        strainmat_.setBuildMode(Matrix::random);
        strainmat_.setSize(num_cells_*6, dispall.size());
        makeDuneMatrixCompressed(strainmat, strainmat_);
        }
    }
}

IMPL_FUNC(void, solve())
{
  try {
    Dune::InverseOperatorResult r;
    ++num_solves_;
    Vector& rhs = this->getLoadVector();
    // when diagonal scaling is active the assembled operator is S A S; we solve
    // (S A S) y = S b and recover x = S y.  sol_ always holds the unscaled solution.
    Vector rhs_eff = rhs;
    if (diagonal_scaling_) {
        for (std::size_t i = 0; i < rhs_eff.size(); ++i)
            rhs_eff[i] *= scale_[i][0];
    }
    if(sol_.size() != rhs.size() || !warm_start_){
        sol_.resize(rhs.size());
        sol_ = 0;
        Vector y = sol_;
        tsolver_->apply(y, rhs_eff, r);
        sol_ = y;
        if (diagonal_scaling_)
            for (std::size_t i = 0; i < sol_.size(); ++i)
                sol_[i] = scale_[i][0]*y[i][0];
    }else{
        const auto& mat= this->getOperator();
        // warm start from the previous (unscaled) solution
        Vector y = sol_;
        if (diagonal_scaling_)
            for (std::size_t i = 0; i < y.size(); ++i)
                y[i] = sol_[i][0]/scale_[i][0];
        auto rhs_tmp = rhs_eff;
        mat.mmv(y,rhs_tmp);
        auto dy = y;
        dy=0;
        if (tol_relative_to_load_) {
            // ||b - A(y+dy)|| <= tol ||b||: the correction only needs the reduction
            // that gets the full residual there, and none at all if it is already.
            // Copy rows hold an identity and a partial load, so only owner rows count.
            double b2 = rhs_eff.two_norm2();
            double r2 = rhs_tmp.two_norm2();
#if HAVE_MPI
            if (grid_.comm().size() > 1) {
                comm_->dot(rhs_eff, rhs_eff, b2);
                comm_->dot(rhs_tmp, rhs_tmp, r2);
            }
#endif
            const double target = linear_tol_ * std::sqrt(b2);
            const double have = std::sqrt(r2);
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
        sol_ = y;
        if (diagonal_scaling_)
            for (std::size_t i = 0; i < sol_.size(); ++i)
                sol_[i] = scale_[i][0]*y[i][0];
    }
    last_linear_iterations_ = r.iterations;
    total_linear_iterations_ += r.iterations;
    last_linear_solve_converged_ = r.converged;
        // MAYBe do other initialization or shift solution.

  } catch (Dune::ISTLError& e) {
    last_linear_solve_converged_ = false;
    OpmLog::warning(fmt::format("exception thrown {}", e.what()));
  }
}

}} // namespace Opm, Elasticity

#endif
