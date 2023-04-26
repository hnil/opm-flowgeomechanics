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
    using namespace Dune;
    const auto& gv = grid_.leafGridView();
  const int comp = 3+(dim-2)*3;
  static const int bfunc = 4+(dim-2)*4;
  //int loadcase = -1;
  //Dune::FieldVector<ctype,comp> eps0 = {1, 1, 1, 0, 0, 0};
  //eps0 = 0;
  Vector& b = A.getLoadVector();
  b = 0;
  A.getLoadVector() = 0;
  if (do_matrix)
    A.getOperator() = 0;


  // start VEM assembly
   // make global point coordinate vector
  vector<double> coords;
  for (const auto& v : vertices(gv)) {
    const auto c = v.geometry().corner(0);
    coords.insert(coords.end(), c.begin(), c.end());
  }

  const int num_cells = gv.size(0); // entities of codim 0
  const auto& ixset = gv.indexSet();
  // count cell faces
  vector<int> num_cell_faces;
  for (const auto& c : elements(gv)) 
    num_cell_faces.push_back(Dune::subEntities(c, Dune::Codim<1>{}).size());
  const int tot_num_cfaces = accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

  // count face corners
  vector<int> num_face_corners;
  for (const auto& c : elements(gv))
    for (const auto& f : Dune::subEntities(c, Dune::Codim<1>{}))
      num_face_corners.push_back(f.geometry().corners());

  // establish all face corners
  vector<int> face_corners;
  for (const auto& c : elements(gv))  
    for (const auto& f: Dune::subEntities(c, Dune::Codim<1>{})) 
      for (int i = 0, f_ix = 0;
           f_ix != tot_num_cfaces && i != num_face_corners[f_ix];
           ++i, f_ix += (i == num_face_corners[f_ix]))
        face_corners.push_back(ixset.subIndex(f, i, 3));

  // correct order of nodes in each face, to ensure they are mentioned in
  // clockwise or counterclockwise order. @@ This is a hack that might only work
  // for 4-faces!
  // for (int i = 0, j = 0; i != (int)num_face_corners.size(); j += num_face_corners[i++])
  //   swap(face_corners[j], face_corners[j+1]);
  
  // body force

  // dirichlet boundary conditions
  
  const int num_fixed_dofs = get<0>(dirichlet_);
  const vector<int>& fixed_dof_ixs = get<1>(dirichlet_);
  const vector<double>& fixed_dof_values = get<2>(dirichlet_);

  // neumann boundary conditions 
  const int num_neumann_faces = 0;
  
  // assemble the mechanical system
  vector<tuple<int, int, double>> A_entries;
  vector<double> rhs;
  body_force_.resize(3*num_cells, 0.0);
  const int numdof =
    vem::assemble_mech_system_3D(&coords[0], num_cells, &num_cell_faces[0], &num_face_corners[0],
                                 &face_corners[0], &ymodule_[0], &pratio_[0], &body_force_[0],
                                 num_fixed_dofs, &fixed_dof_ixs[0], &fixed_dof_values[0],
                                 num_neumann_faces, nullptr, nullptr,
                                 A_entries, rhs);
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
