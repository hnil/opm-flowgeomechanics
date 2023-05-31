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
#include <vector>
#include <algorithm>
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

    IMPL_FUNC(void, calculateStress())
{
    //assumes the grid structure is made
    const int num_neumann_faces = 0;
    num_cells_ = grid_.leafGridView().size(0); // entities of codim 0
    // assemble the mechanical system
    vem::StabilityChoice stability_choice = vem::D_RECIPE;
    //const int numdof =
    // const int tot_num_faces = accumulate(num_cell_faces_, num_cell_faces_ + num_cells_, 0);
    // const int tot_num_fcorners = accumulate(num_face_corners_, &num_face_corners_[0] + tot_num_faces, 0);
    // const int tot_num_nodes = *max_element(face_corners_, face_corners + tot_num_fcorners) + 1;
    stress_.resize(num_cells_);
    std::vector<std::array<double,6>> stress;
    stress.resize(num_cells_);

    std::vector<double> dispall;
    dispall.resize(3*grid_.leafGridView().size(3));
    {
        Vector dispalldune;
        dispalldune.resize(3*grid_.leafGridView().size(3));
        this->expandSolution(dispalldune,this->u);
        for(size_t i=0; i < dispall.size(); ++i){
            dispall[i] = dispalldune[i];//fieldvector<double,1> can be converted to double
        }
    }
    
    vem::compute_stress_3D(&coords_[0], num_cells_, &num_cell_faces_[0], &num_face_corners_[0],
                           &face_corners_[0], &ymodule_[0], &pratio_[0],
                           dispall,
                           stress,
                           stability_choice);
    // copy to dune definitions
    stress_.resize(num_cells_);
    for(size_t i=0; i < num_cells_; ++i){
        for(size_t k=0; k < 6; ++k){
            stress_[i][k] = stress[i][k];
        }
    }
        
}
    

    IMPL_FUNC(void, assemble(const Vector& pressure, bool do_matrix, bool do_vector))
{
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
        const int numdof =
            vem::assemble_mech_system_3D(&coords_[0], num_cells_, &num_cell_faces_[0], &num_face_corners_[0],
                                         &face_corners_[0], &ymodule_[0], &pratio_[0], &body_force_[0],
                                         num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                                         num_neumann_faces, nullptr, nullptr,
                                         A_entries, rhs_force_, stability_choice);
        
    
    
        this->makeDuneMatrix(A_entries);

        
        // make indexing for div operator i.e. all nodes to dofs
        std::vector<int> dof_idx(grid_.leafGridView().size(3)*3);
        std::iota(dof_idx.begin(), dof_idx.end(),0);
        
        std::set_difference(dof_idx.begin(), dof_idx.end(), fixed_dof_ixs.begin(), fixed_dof_ixs.end(),std::back_inserter(idx_free_));
        
    }
    if(do_vector){
        //NB rhs_force_ is calculated by matrix call
        vector<double> rhs_pressure;
        vector<double> std_pressure(pressure.size(), 0);
        for(size_t i = 0; i < pressure.size(); ++i){
            std_pressure[i] = pressure[i][0];
        }
        vem::potential_gradient_force_3D(&coords_[0], num_cells_, &num_cell_faces_[0], &num_face_corners_[0],
                                         &face_corners_[0], &std_pressure[0], rhs_pressure);

        //Sign is added here  i.e \div \sigma = 
        vector<double> rhs(rhs_force_);
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
