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
#ifndef OPM_ELASTICITY_SOLVER_IMPL_HPP
#define OPM_ELASTICITY_SOLVER_IMPL_HPP

#include <iostream>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

namespace Opm {
namespace Elasticity {

#undef IMPL_FUNC
#define IMPL_FUNC(A,B) template<class GridType> \
                         A ElasticitySolver<GridType>::B

IMPL_FUNC(std::vector<BoundaryGrid::Vertex>, 
          extractFace(Direction dir, ctype coord))
{
  std::vector<BoundaryGrid::Vertex> result;
  const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView>  mapper(gv.leafGridView(), Dune::mcmgVertexLayout());
  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafGridView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      BoundaryGrid::Vertex v;
      v.i = mapper.index(*it);
      BoundaryGrid::extract(v.c,it->geometry().corner(0),log2(float(dir)));
      result.push_back(v);
    }
  }

  return result;
}


IMPL_FUNC(BoundaryGrid, extractMasterFace(Direction dir,
                                          ctype coord,
                                          SIDE side,
                                          bool dc))
{
  static const int V1[3][4] = {{0,2,4,6},
                               {0,1,4,5},
                               {0,1,2,3}};
  static const int V2[3][4] = {{1,3,5,7},
                               {2,3,6,7},
                               {4,5,6,7}};
  const LeafIndexSet& set = gv.leafGridView().indexSet();

  int c = 0;
  int i = log2(float(dir));
  BoundaryGrid result;
  // we first group nodes into this map through the coordinate of lower left 
  // vertex. we then split this up into pillars for easy processing later
  std::map<double, std::vector<BoundaryGrid::Quad> > nodeMap;
  for (LeafIterator cell  = gv.leafGridView().template begin<0>(); 
                    cell != gv.leafGridView().template end<0>(); ++cell, ++c) {
    std::vector<BoundaryGrid::Vertex> verts;
    int idx=0; 
    if (side == LEFT)
     idx = set.subIndex(*cell,V1[i][0],dim);
    else if (side == RIGHT)
     idx = set.subIndex(*cell,V2[i][0],dim);
    Dune::FieldVector<double, 3> pos = gv.vertexPosition(idx);
    if (isOnPlane(dir,pos,coord)) {
      for (int j=0;j<4;++j) {
        if (side == LEFT)
          idx = set.subIndex(*cell,V1[i][j],dim);
        if (side == RIGHT)
          idx = set.subIndex(*cell,V2[i][j],dim);
        pos = gv.vertexPosition(idx);
        if (!isOnPlane(dir,pos,coord))
          continue;
        BoundaryGrid::Vertex v;
        BoundaryGrid::extract(v,pos,i);
        v.i = idx;
        verts.push_back(v);
      }
    }
    if (verts.size() == 4) {
      BoundaryGrid::Quad q;
      q.v[0] = minXminY(verts);
      q.v[1] = maxXminY(verts);
      if (dc) {
        q.v[2] = minXmaxY(verts);
        q.v[3] = maxXmaxY(verts);
      } else {
        q.v[2] = maxXmaxY(verts);
        q.v[3] = minXmaxY(verts);
      }
      std::map<double, std::vector<BoundaryGrid::Quad> >::iterator it;
      for (it  = nodeMap.begin(); it != nodeMap.end(); ++it) {
        if (fabs(it->first-q.v[0].c[0]) < 1.e-7) {
          it->second.push_back(q);
          break;
        }
      }
      if (it == nodeMap.end())
        nodeMap[q.v[0].c[0]].push_back(q);

      result.add(q);
    }
  }

  int p=0;
  std::map<double, std::vector<BoundaryGrid::Quad> >::const_iterator it;
  for (it = nodeMap.begin(); it != nodeMap.end(); ++it, ++p) {
    for (size_t ii=0;ii<it->second.size();++ii)
      result.addToColumn(p,it->second[ii]);
  }

  return result;
}

// IMPL_FUNC(void, determineSideFaces(const double* min, const double* max))
// {
//   master.push_back(extractMasterFace(X,min[0]));
//   master.push_back(extractMasterFace(Y,min[1]));
//   master.push_back(extractMasterFace(Z,min[2]));

//   slave.push_back(extractFace(X,max[0]));
//   slave.push_back(extractFace(Y,max[1]));
//   slave.push_back(extractFace(Z,max[2]));
// }

IMPL_FUNC(void, findBoundaries(double* min, double* max))
{
  max[0] = max[1] = max[2] = -1e5;
  min[0] = min[1] = min[2] = 1e5;
  const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafGridView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min[i] = std::min(min[i],it->geometry().corner(0)[i]);
      max[i] = std::max(max[i],it->geometry().corner(0)[i]);
    }
  }
}

IMPL_FUNC(void, fixNodes(const std::vector<size_t>& fixed_nodes))
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid 
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(), Dune::mcmgVertexLayout());

  NodeValue zerovec;
  zerovec = 0.0;
  // iterate over vertices
  for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
      int indexi = mapper.index(*it);
      assert(indexi == gv.leafGridView().indexSet().index(it));
      bool exist = std::find(fixed_nodes.begin(), fixed_nodes.end(), indexi)
          !=
          fixed_nodes.end();
      if(exist){
          A.updateFixedNode(indexi,std::make_pair(XYZ,zerovec));
      }
  }
}

    
IMPL_FUNC(void, fixPoint(Direction dir,
                         GlobalCoordinate coord,
                         const NodeValue& value))
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid 
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(), Dune::mcmgVertexLayout());

  // iterate over vertices
  for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
    if (isOnPoint(it->geometry().corner(0),coord)) {
      int indexi = mapper.index(*it);
      A.updateFixedNode(indexi,std::make_pair(dir,value));
    }
  }
}

IMPL_FUNC(template<int comp> void,
          averageStress(Dune::BlockVector<Dune::FieldVector<ctype,comp>>& sigmacells,
                        const Vector& uarg))
{
  
  static const int bfunc = 4+(dim-2)*4;

  const LeafIterator itend = gv.leafGridView().template end<0>();

  Dune::FieldMatrix<ctype,comp,comp> C;
  Dune::FieldVector<ctype,comp> eps0;
  eps0 = 0;
  //eps0[loadcase] = 1; // NB do not understand
  int m=0;
  
  for (LeafIterator it = gv.leafGridView().template begin<0>(); it != itend; ++it) {
    materials[m]->getConstitutiveMatrix(C);
    // determine geometry type of the current element and get the matching reference element
    Dune::GeometryType gt = it->type();

    Dune::FieldVector<ctype,bfunc*dim> v;
    A.extractValues(v,uarg,it);
    Dune::FieldVector<ctype,comp> sigma;
    sigma = 0;
    double volume=0;
    // get a quadrature rule of order two for the given geometry type
    const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
        r != rule.end() ; ++r) {
      // compute the jacobian inverse transposed to transform the gradients
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
        it->geometry().jacobianInverseTransposed(r->position());

      ctype detJ = it->geometry().integrationElement(r->position());

      volume += detJ*r->weight();

      Dune::FieldMatrix<ctype,comp,dim*bfunc> lB;
      E.getBmatrix(lB,r->position(),jacInvTra);

      Dune::FieldVector<ctype,comp> s;
      E.getStressVector(s,v,eps0,lB,C);
      s *= detJ*r->weight();
      sigma += s;
    }
    sigma /= volume;
    if (Escale > 0){
        sigma /= Escale/Emin;
    }
    sigmacells[m] = sigma;
    m++;
  }
  
  
}


// IMPL_FUNC(bool, isOnPlane(Direction plane,
//                           GlobalCoordinate coord,
//                           ctype value))
// {
//   if (plane < X || plane > Z)
//     return false;
//   int p = log2(float(plane));
//   ctype delta = fabs(value-coord[p]);
//   return delta < tol;
// }

// IMPL_FUNC(void, fixLine(Direction dir,
//                         ctype x, ctype y,
//                         const NodeValue& value))
// {
//   typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
//   const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

//   // make a mapper for codim 0 entities in the leaf grid
//   using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
//   Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv, Dune::mcmgVertexLayout());

//   // iterate over vertices
//   for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
//     if (isOnLine(dir,it->geometry().corner(0),x,y)) {
//       int indexi = mapper.index(*it);
//       A.updateFixedNode(indexi,std::make_pair(XYZ,value));
//     }
//   }
// }

// IMPL_FUNC(bool, isOnLine(Direction dir,
//                          GlobalCoordinate coord,
//                          ctype x, ctype y))
// {
//   if (dir < X || dir > Z)
//     return false;
//   int ix = int(log2(float(dir))+1) % 3;
//   int iy = int(log2(float(dir))+2) % 3;
//   ctype delta = x-coord[ix];
//   if (delta > tol || delta < -tol)
//     return false;
//   delta = y-coord[iy];
//   if (delta > tol || delta < -tol)
//     return false;

//   return true;
// }

// IMPL_FUNC(bool, isOnPoint(GlobalCoordinate coord,
//                           GlobalCoordinate point))
// {
//   GlobalCoordinate delta = point-coord;
//   return delta.one_norm() < tol;
// }

    IMPL_FUNC(void, assemble(const Vector& pressure, bool matrix, bool vector))
{
  const int comp = 3+(dim-2)*3;
  static const int bfunc = 4+(dim-2)*4;
  //int loadcase = -1;
  Dune::FieldVector<ctype,comp> eps0 = {1, 1, 1, 0, 0, 0};
  //eps0 = 0;
  Vector& b = A.getLoadVector();
  b = 0;
  A.getLoadVector() = 0;
  if (matrix)
    A.getOperator() = 0;

  
  for (int i=0;i<2;++i) {
    if (color[1].size() && matrix)
      std::cout << "\tprocessing " << (i==0?"red ":"black ") << "elements" << std::endl;
#pragma omp parallel for schedule(static)
    for (size_t j=0;j<color[i].size();++j) {
      Dune::FieldMatrix<ctype,comp,comp> C;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> K;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc>* KP=0;
      Dune::FieldVector<ctype,dim*bfunc> ES;
      Dune::FieldVector<ctype,dim*bfunc>* EP=0;
      if (matrix)
        KP = &K;
      if (vector)
        EP = &ES;

      for (size_t k=0;k<color[i][j].size();++k) {
        LeafIterator it = gv.leafGridView().template begin<0>();
        for (int l=0;l<color[i][j][k];++l)
          ++it;
        size_t cell_num = color[i][j][k];
        bool valid = materials[color[i][j][k]]->getConstitutiveMatrix(C);
        if(!valid){
            OPM_THROW(std::runtime_error,"Not valid C matrix"); 
        }
        // determine geometry type of the current element and get the matching reference element
        Dune::GeometryType gt = it->type();

        Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> Aq;
        K = 0;
        ES = 0;

        // get a quadrature rule of order two for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
            r != rule.end() ; ++r) {
          // compute the jacobian inverse transposed to transform the gradients
          Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
            it->geometry().jacobianInverseTransposed(r->position());

          ctype detJ = it->geometry().integrationElement(r->position());
          if (detJ <= 1.e-5 && verbose) {
            std::cout << "cell " << color[i][j][k] << " is (close to) degenerated, detJ " << detJ << std::endl;
            double zdiff=0.0;
            for (int ii=0;ii<4;++ii)
              zdiff = std::max(zdiff, it->geometry().corner(ii+4)[2]-it->geometry().corner(ii)[2]);
            std::cout << " - Consider setting ctol larger than " << zdiff << std::endl;
          }

          Dune::FieldMatrix<ctype,comp,dim*bfunc> lB;
          E.getBmatrix(lB,r->position(),jacInvTra);

          if (matrix) {
            E.getStiffnessMatrix(Aq,lB,C,detJ*r->weight());
            K += Aq;
          }

          // load vector
          if (EP) {
            Dune::FieldVector<ctype,dim*bfunc> temp;
            //temp = Dune::FMatrixHelp::multTransposed(lB,Dune::FMatrixHelp::mult(C,eps0));
            auto Ipressure = Dune::FMatrixHelp::mult(C,eps0);
            Ipressure = eps0*pressure[cell_num][0];
            temp = Dune::FMatrixHelp::multTransposed(lB,Ipressure);
            temp *= -detJ*r->weight();
            ES += temp;
          }
        }
        A.addElement(KP,EP,it,&b); // NULL is no static forse based on the itegration point??
      }
    }
  }
}


// IMPL_FUNC(void, fixCorners(const double* min, const double* max))
// {
//   ctype c[8][3] = {{min[0],min[1],min[2]},
//                    {max[0],min[1],min[2]},
//                    {min[0],max[1],min[2]},
//                    {max[0],max[1],min[2]},
//                    {min[0],min[1],max[2]},
//                    {max[0],min[1],max[2]},
//                    {min[0],max[1],max[2]},
//                    {max[0],max[1],max[2]}};
//   for (int i=0;i<8;++i) {
//     GlobalCoordinate coord;
//     coord[0] = c[i][0]; coord[1] = c[i][1]; coord[2] = c[i][2];
//     fixPoint(XYZ,coord);
//   }
// }





IMPL_FUNC(void, solve())
{
  try {
    Dune::InverseOperatorResult r;
    Vector& rhs = A.getLoadVector();
    u.resize(rhs.size());
    u = 0;
    tsolver_->apply(u, rhs, r);
    std::cout << "\tsolution norm: " << u.two_norm() << std::endl;
  } catch (Dune::ISTLError& e) {
    std::cerr << "exception thrown " << e << std::endl;
  }
}

}} // namespace Opm, Elasticity

#endif
