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

#include <iostream>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <tuple>
#include <vector>
#include <algorithm>
#include <opm/common/TimingMacros.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>
namespace Opm {
namespace Elasticity {

#undef IMPL_FUNC
#define IMPL_FUNC(A,B) template<class GridType> \
                         A VemElasticitySolver<GridType>::B


// IMPL_FUNC(void, fixNodes(const std::vector<size_t>& fixed_nodes))
// {
//     // makestructure for vem assembly
//     this->fixNodesVem(fixed_nodes);

//   typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
//   const VertexLeafIterator itend = grid_.leafGridView().template end<dim>();

//   // make a mapper for codim 0 entities in the leaf grid
//   using LeafGridView = typename GridType::LeafGridView;
//   Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(grid_.leafGridView(), Dune::mcmgVertexLayout());

//   NodeValue zerovec;
//   zerovec = 0.0;
//   // iterate over vertices
//   for (VertexLeafIterator it = grid_.leafGridView().template begin<dim>(); it != itend; ++it) {
//       int indexi = mapper.index(*it);
//       //assert(indexi == grid_.leafGridView().indexSet().index(it));
//       bool exist = std::find(fixed_nodes.begin(), fixed_nodes.end(), indexi)
//           !=
//           fixed_nodes.end();
//       if(exist){
//           A.updateFixedNode(indexi,std::make_pair(XYZ,zerovec));
//       }
//   }
// }

    IMPL_FUNC(void, calculateStress(bool precomputed))
    {
        if (precomputed) {
            OPM_TIMEBLOCK(calculateStressPrecomputed);
            //NB stressmat is defined in linear indices not block linear indices
            Vector dispalldune;
            dispalldune.resize(3 * grid_.leafGridView().size(3));
            this->expandSolution(dispalldune, this->u);
            // Dune::BlockVector< DuneFieldVector<double,1> >
            Vector stress(6 * grid_.leafGridView().size(0));
            stressmat_.mv(dispalldune,stress);
            stress_.resize(num_cells_);
            for (size_t i = 0; i < num_cells_; ++i) {
                for (size_t k = 0; k < 6; ++k) {
                    stress_[i][k] = stress[i * 6 + k];
                }
            }
        } else {
            OPM_TIMEBLOCK(calculateStressFull);
            // assumes the grid structure is made
            const int num_neumann_faces = 0;
            num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
            // assemble the mechanical system
            vem::StabilityChoice stability_choice = vem::D_RECIPE;
            // const int numdof =
            //  const int tot_num_faces = accumulate(num_cell_faces_, num_cell_faces_ + num_cells_, 0);
            //  const int tot_num_fcorners = accumulate(num_face_corners_, &num_face_corners_[0] + tot_num_faces, 0);
            //  const int tot_num_nodes = *max_element(face_corners_, face_corners + tot_num_fcorners) + 1;
            stress_.resize(num_cells_);
            std::vector<std::array<double, 6>> stress;
            stress.resize(num_cells_);

            std::vector<double> dispall;
            dispall.resize(3 * grid_.leafGridView().size(3));
            {
                Vector dispalldune;
                dispalldune.resize(3 * grid_.leafGridView().size(3));
                this->expandSolution(dispalldune, this->u);
                for (size_t i = 0; i < dispall.size(); ++i) {
                    dispall[i] = dispalldune[i]; // fieldvector<double,1> can be converted to double
                }
            }
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
                                   stability_choice,
                                   stressmat,
                                   false,
                                   true
                );
            // copy to dune definitions
            stress_.resize(num_cells_);
            for (size_t i = 0; i < num_cells_; ++i) {
                for (size_t k = 0; k < 6; ++k) {
                    stress_[i][k] = stress[i][k];
                }
            }
        }
    }

    IMPL_FUNC(void, calculateStrain(bool precomputed))
    {
        if (precomputed) {
            OPM_TIMEBLOCK(calculateStrainPrecomputed);
            //NB stressmat is defined in linear indices not block linear indices
            Vector dispalldune;
            dispalldune.resize(3 * grid_.leafGridView().size(3));
            this->expandSolution(dispalldune, this->u);
            // Dune::BlockVector< DuneFieldVector<double,1> >
            Vector strain(6 * grid_.leafGridView().size(0));
            strainmat_.mv(dispalldune,strain);
            strain_.resize(num_cells_);
            for (size_t i = 0; i < num_cells_; ++i) {
                for (size_t k = 0; k < 6; ++k) {
                    strain_[i][k] = strain[i * 6 + k];
                }
            }
        } else {
            OPM_TIMEBLOCK(calculateStressFull);
            // assumes the grid structure is made
            const int num_neumann_faces = 0;
            num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
            // assemble the mechanical system
            vem::StabilityChoice stability_choice = vem::D_RECIPE;
            // const int numdof =
            //  const int tot_num_faces = accumulate(num_cell_faces_, num_cell_faces_ + num_cells_, 0);
            //  const int tot_num_fcorners = accumulate(num_face_corners_, &num_face_corners_[0] + tot_num_faces, 0);
            //  const int tot_num_nodes = *max_element(face_corners_, face_corners + tot_num_fcorners) + 1;
            strain_.resize(num_cells_);
            std::vector<std::array<double, 6>> strain;
            strain.resize(num_cells_);

            std::vector<double> dispall;
            dispall.resize(3 * grid_.leafGridView().size(3));
            {
                Vector dispalldune;
                dispalldune.resize(3 * grid_.leafGridView().size(3));
                this->expandSolution(dispalldune, this->u);
                for (size_t i = 0; i < dispall.size(); ++i) {
                    dispall[i] = dispalldune[i]; // fieldvector<double,1> can be converted to double
                }
            }
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
                                   stability_choice,
                                   stressmat,
                                   false,
                                   false
                );
            // copy to dune definitions
            strain_.resize(num_cells_);
            for (size_t i = 0; i < num_cells_; ++i) {
                for (size_t k = 0; k < 6; ++k) {
                    strain_[i][k] = strain[i][k];
                }
            }
        }
    }

    IMPL_FUNC(void, assemble(const Vector& pressure, bool do_matrix, bool do_vector))
{
    OPM_TIMEBLOCK(assemble);
    using namespace std;
    Vector& b = A.getLoadVector();
    b = 0;
    A.getLoadVector() = 0;
    if (do_matrix)
        A.getOperator() = 0;

    if(do_matrix){
        vem::getGridVectors(grid_,coords_,
                            num_cell_faces_,
                            num_face_corners_,
                            face_corners_);

        const int num_fixed_dofs = get<0>(dirichlet_);
        const vector<int>& fixed_dof_ixs = std::get<1>(dirichlet_);
        const vector<double>& fixed_dof_values = std::get<2>(dirichlet_);

        // neumann boundary conditions
        const int num_neumann_faces = 0;
        num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
        // assemble the mechanical system
        vector<tuple<int, int, double>> A_entries;
        vem::StabilityChoice stability_choice = vem::D_RECIPE;
        {
        OPM_TIMEBLOCK(assembleVEM);
        const int numdof =    vem::assemble_mech_system_3D(&coords_[0], num_cells_, &num_cell_faces_[0], &num_face_corners_[0],
                                         &face_corners_[0], &ymodule_[0], &pratio_[0], &body_force_[0],
                                         num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                                         num_neumann_faces, nullptr, nullptr,
                                         A_entries, rhs_force_, stability_choice);
        }


        this->makeDuneSystemMatrix(A_entries);

        {
        OPM_TIMEBLOCK(setUpExtraStructuresForVEM);
        // make indexing for div operator i.e. all nodes to dofs
        std::vector<int> dof_idx(grid_.leafGridView().size(3)*3);
        std::iota(dof_idx.begin(), dof_idx.end(),0);

        std::set_difference(dof_idx.begin(), dof_idx.end(), fixed_dof_ixs.begin(), fixed_dof_ixs.end(),std::back_inserter(idx_free_));
        //
        vector<double> rhs_tmp(pressure.size(),0);
        vector<double> pressure_tmp(pressure.size(),0);
        vector<tuple<int, int, double>> divmat;
        vem::potential_gradient_force_3D(&coords_[0],
                                         num_cells_,
                                         &num_cell_faces_[0],
                                         &num_face_corners_[0],
                                         &face_corners_[0],
                                         &pressure_tmp[0],
                                         rhs_tmp,
                                         divmat,
                                         true);
        // sort(divmat_.begin(),
        //      divmat_.end(),
        //      [](const auto& aa, const auto& bb) { return std::get<1>(aa) < std::get<1>(bb); });
        std::vector<int> global_to_dof(grid_.leafGridView().size(3)*3,-1);
        for(size_t i=0; i< idx_free_.size(); ++i){
            global_to_dof[idx_free_[i]] = i;
        }
        //renumber and eliminate dof fized
        vector<tuple<int, int, double>> divmatdof;
        for(const auto& elem: divmat){
            int I = global_to_dof[std::get<0>(elem)];
            if(I>-1){
                // renumber
                int J = get<1>(elem);
                double val = std::get<2>(elem);
                divmatdof.push_back(std::tuple<int, int, double>(I,J,val));
            }
        }
        // finaly make dune matrix
        divmat_.setBuildMode(Matrix::implicit);
        // map from dof=3*nodes at a cell (ca 3*3*3) to cell
        divmat_.setImplicitBuildModeParameters (3*3*3, 0.4);
        divmat_.setSize(idx_free_.size(), num_cells_);
        makeDuneMatrixCompressed(divmatdof,divmat_);
        // also make stress matrix
        std::vector<std::tuple<int, int, double>> stressmat;
        std::vector<double> dispall(grid_.leafGridView().size(3)*3);
        std::vector<std::array<double,6>> stresstmp(grid_.leafGridView().size(0));
        {
        OPM_TIMEBLOCK(setUpStressStressMatrix);
        vem::compute_stress_3D(&coords_[0],
                               num_cells_,
                               &num_cell_faces_[0],
                               &num_face_corners_[0],
                               &face_corners_[0],
                               &ymodule_[0], &pratio_[0],
                               dispall,
                               stresstmp,
                               stability_choice,
                               stressmat,
                               true,
                               true
        );
        }
        stressmat_.setBuildMode(Matrix::implicit);
        stressmat_.setImplicitBuildModeParameters (3*3*3, 0.4);
        stressmat_.setSize(num_cells_*6, dispall.size());
        makeDuneMatrixCompressed(stressmat, stressmat_);
        std::vector<std::tuple<int, int, double>> strainmat;
        {
        OPM_TIMEBLOCK(setUpStressStrainMatrix);
        vem::compute_stress_3D(&coords_[0],
                               num_cells_,
                               &num_cell_faces_[0],
                               &num_face_corners_[0],
                               &face_corners_[0],
                               &ymodule_[0], &pratio_[0],
                               dispall,
                               stresstmp,
                               stability_choice,
                               strainmat,
                               true,
                               false
        );
        }
        strainmat_.setBuildMode(Matrix::implicit);
        strainmat_.setImplicitBuildModeParameters (3*3*3, 0.4);
        strainmat_.setSize(num_cells_*6, dispall.size());
        makeDuneMatrixCompressed(strainmat, strainmat_);
        }
    }
    if(do_vector){
        OPM_TIMEBLOCK(calculateRHS);
        //NB rhs_force_ is calculated by matrix call
        vector<double> rhs_pressure;
        vector<double> std_pressure(pressure.size(), 0);
        for(size_t i = 0; i < pressure.size(); ++i){
            std_pressure[i] = pressure[i][0];
        }
        std::vector<std::tuple<int, int, double>> divmat;
        vem::potential_gradient_force_3D(&coords_[0],
                                         num_cells_,
                                         &num_cell_faces_[0],
                                         &num_face_corners_[0],
                                         &face_corners_[0],
                                         &std_pressure[0],
                                         rhs_pressure,
                                         divmat,
                                         false);

        //Sign is added here  i.e \div \sigma =
        vector<double> rhs(rhs_force_);
        assert(rhs_force_.size() == idx_free_.size());
        for(size_t i=0; i< idx_free_.size(); ++i){
            rhs[i] += rhs_pressure[idx_free_[i]];
        }
        b.resize(rhs.size());
        // end initialization
        for(int i=0; i < rhs.size(); ++i){
            b[i] = rhs[i];
        }
    }
}


IMPL_FUNC(void, solve())
{
  try {
    Dune::InverseOperatorResult r;
    Vector& rhs = A.getLoadVector();
    u.resize(rhs.size());
    u = 0;
    tsolver_->apply(u, rhs, r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "exception thrown " << e << std::endl;
  }
}

}} // namespace Opm, Elasticity

#endif
