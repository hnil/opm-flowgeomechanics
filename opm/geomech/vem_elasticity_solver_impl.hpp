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



    IMPL_FUNC(void, assemble(const Vector& pressure, bool do_matrix, bool do_vector))
{
    using namespace std;
    Vector& b = A.getLoadVector();
    b = 0;
    A.getLoadVector() = 0;
    if (do_matrix)
        A.getOperator() = 0;

    std::vector<double> coords;
    std::vector<int> num_cell_faces, num_face_corners,face_corners;
    vem::getGridVectors(grid_,coords,
                        num_cell_faces,
                        num_face_corners,
                        face_corners);
    const int num_fixed_dofs = get<0>(dirichlet_);
    const vector<int>& fixed_dof_ixs = std::get<1>(dirichlet_);
    const vector<double>& fixed_dof_values = std::get<2>(dirichlet_);

  // neumann boundary conditions 
  const int num_neumann_faces = 0;
  const int num_cells = grid_.leafGridView().size(0); // entities of codim 0
  // assemble the mechanical system
  vector<tuple<int, int, double>> A_entries;
  vector<double> rhs;
  body_force_.resize(3*num_cells, 0.0);
  for(int i=0; i<num_cells; ++i){
      body_force_[3*i+2] = 2000*9.8;
  }
  const int numdof =
    vem::assemble_mech_system_3D(&coords[0], num_cells, &num_cell_faces[0], &num_face_corners[0],
                                 &face_corners[0], &ymodule_[0], &pratio_[0], &body_force_[0],
                                 num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                                 num_neumann_faces, nullptr, nullptr,
                                 A_entries, rhs);
  // adding pressure force
  vector<double> rhs_pressure;
  vector<double> std_pressure(pressure.size(), 0);
  for(size_t i = 0; i < pressure.size(); ++i){
      std_pressure[i] = pressure[i][0];
  }
  vem::potential_gradient_force_3D(&coords[0], num_cells, &num_cell_faces[0], &num_face_corners[0],
                                   &face_corners[0], &std_pressure[0],rhs_pressure);
  
 //  std::transform(rhs.begin(),rhs.end(),rhs_pressure.begin(),rhs.begin(), [](double a, double b) {return a+b;});
  // hopefully consisten with matrix
  // should be moved to initialization
  std::vector< std::set<int> > rows;
  int nrows=0;
  int ncols=0;
  for(auto matel: A_entries){
      int i = get<0>(matel);
      int j = get<1>(matel);
      nrows = std::max(nrows,i);
      ncols = std::max(ncols,j);
  }
  nrows = nrows+1;
  ncols = ncols+1;
  assert(nrows==ncols);
  rows.resize(nrows);
  for(auto matel: A_entries){
      int i = get<0>(matel);
      int j = get<1>(matel);
      rows[i].insert(j);
  }
  auto& MAT =this->A.getOperator();      
  MatrixOps::fromAdjacency(MAT, rows, nrows, ncols);
  MAT = 0;
  for(auto matel:A_entries){
      int i = get<0>(matel);
      int j = get<1>(matel);
      double val = get<2>(matel);
      A.addMatElement(i,j, val);// += val;
  }
  b.resize(rhs.size());
  // end initialization
  if(do_vector){
      for(int i=0; i < numdof; ++i){
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
