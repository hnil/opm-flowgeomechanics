#include <opm/geomech/GridStretcher.hpp>
#include <opm/geomech/param_interior.hpp>
#include <array>
#include <limits>
#include <assert.h>
#include <dune/common/fvector.hh> // FieldVector
#include <dune/grid/common/mcmgmapper.hh> // for element mapper
#include <opm/geomech/convex_boundary.hpp>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace std;
using namespace Dune;

// ============================================================================
namespace // anonymous
// ============================================================================
{

using GridView = typename Opm::GridStretcher::Grid::LeafGridView;
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
using CoordType = Opm::GridStretcher::CoordType;


// ----------------------------------------------------------------------------  
double adjust_disp_to_convex(const double* const startpt,
                             const double* const endpt,
                             const double* const pt,
                             const double* const dir)
// ----------------------------------------------------------------------------  
{
  // solve for alpha in the equation:
  // alpha * dir + t * (endpt-startpt) = enpt - pt
  // This will ensure that pt + dir * alpha lies on the edge between startpt and
  // endpt.
  const double vx = endpt[0] - startpt[0];
  const double vy = endpt[1] - startpt[1];
  const double dx = endpt[0] - pt[0];
  const double dy = endpt[1] - pt[1];

  return (vy * dx - vx * dy) / (dir[0] * vy - dir[1] * vx);
}

// ----------------------------------------------------------------------------
template <typename T>
vector<T> extract_elements(const vector<T>& vec, const vector<size_t>& ixs)
// ----------------------------------------------------------------------------
{
  vector<T> result;
  for (const auto ix : ixs)
    result.push_back(vec[ix]);
  return result;
}

// ----------------------------------------------------------------------------    
double dist3D(const Dune::FieldVector<double, 3>& p1,
              const Dune::FieldVector<double, 3>& p2)
// ----------------------------------------------------------------------------    
{
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  double dz = p1[2] - p2[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

  
// ----------------------------------------------------------------------------  
vector<double> pick_nodes(const vector<size_t>& ixs, const vector<double>& coords2D)
// ----------------------------------------------------------------------------  
{
  vector<double> result;
  for (auto i : ixs) {
    result.push_back(coords2D[2*i]);
    result.push_back(coords2D[2*i+1]);
  }
  return result;
}
  
// ----------------------------------------------------------------------------
vector<double> pickCoords2D(const vector<CoordType>& coords3D)
// ----------------------------------------------------------------------------
{
  // pick the two coordinate axes with the largest spread
  const size_t N = coords3D.size(); // total number of nodes (interior and boundary)

  array<double, 3> low {coords3D[0][0], coords3D[0][1], coords3D[0][2]}, high {low};
  for (size_t i = 0; i != N; ++i) {
    for (int d = 0; d != 3; ++d) {
      low[d] = std::min(low[d], coords3D[i][d]);
      high[d] = std::max(high[d], coords3D[i][d]);
    }
  }
  const array<double, 3> span {high[0] - low[0], high[1] - low[1], high[2] - low[2] };
  const size_t min_ix = distance(span.begin(), min_element(span.begin(), span.end()));
  const size_t ix1 = (min_ix + 1) % 3;
  const size_t ix2 = (min_ix + 2) % 3;

  vector<double> coords2D;
  for (size_t i = 0; i != N; ++i) {
    coords2D.push_back(coords3D[i][ix1]);
    coords2D.push_back(coords3D[i][ix2]);
  }
  return coords2D;
}

// ----------------------------------------------------------------------------
size_t find_coord_in(const CoordType& c, const vector<CoordType>& cvec)
// ----------------------------------------------------------------------------
{
  const double TOL = 1e-3;
  
  return
    find_if(cvec.begin(),
            cvec.end(),
            [&] (const CoordType& el) {return dist3D(el, c) < TOL;})
    - cvec.begin();

}

// ----------------------------------------------------------------------------
CoordType normalize(const CoordType& vec) { return vec / vec.two_norm(); }
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
  CoordType cross(const CoordType& v1, const CoordType& v2) 
// ----------------------------------------------------------------------------    
{
  return CoordType {v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0]};
}

    
// // ----------------------------------------------------------------------------
// const vector<CoordType> compute_bnode_displacements(const vector<double>& amounts,
//                                                     const vector<size_t>& bcindices,
//                                                     const Opm::GridStretcher::CellBnodeMap& c2bix,
//                                                     const Opm::GridStretcher::Grid& grid,
//                                                     const vector<CoordType>& nodecoords,
//                                                     const BoundaryNormals& bnormals)
// // ----------------------------------------------------------------------------
// {
//   const size_t N = size(bcindices);
//   // map amounts to displacement vectors for each boundary point
//   assert(amounts.size() == N); // one scalar per boundary cell
//   vector<vector<tuple<double, CoordType>>> disp(N); // node displacements
  
//   auto aiter = amounts.begin();
//   for (size_t bc = 0; bc != amounts.size(); ++bc) {
//     // compute directional vector
//     const auto entry = c2bix.find(bcindices[bc])->second;
//     const auto elem = grid.entity(get<0>(entry));
//     const auto ccenter = elem.geometry().center();
//     const auto ecenter = (nodecoords[get<1>(entry)] + nodecoords[get<2>(entry)]) / 2;
//     const auto nvec = normalize(ecenter - ccenter);
//     disp[get<1>(entry)].push_back({*aiter, nvec});
//     disp[get<2>(entry)].push_back({*aiter++, nvec});
//   }

//   // check that every boundary node has got two displacements associated with it
    
//   // combine displacement vectors
//   vector<CoordType> result(N);
//   for (size_t i = 0; i != N; ++i) {
//     const auto& entry = disp[i];
//     const double a1 = get<0>(entry[0]);
//     const double a2 = get<0>(entry[1]);
//     const auto& v1 = get<1>(entry[0]);
//     const auto& v2 = get<1>(entry[1]);
    
//     //const auto dir = normalize(v1 + v2);
//     const auto& dir = bnnormals[i];

//     assert(dir.dot(v1) >= 0);
//     assert(dir.dot(v2) >= 0);
//     result[i] = dir * max(dir.dot(v1) * a1, dir.dot(v2) * a2);

//   }
//   return result;
// }

    
}; // end anonymous namespace

// ============================================================================
namespace Opm
// ============================================================================
//#include <dune/geometry/referenceelements.hh>
{    
vector<size_t> GridStretcher::boundary_node_indices(const Grid& grid){
  //typedef ReferenceElement< typename Grid::ctype, dimGrid > RefElement;
  //typedef ReferenceElements< typename Grid::ctype, dimGrid > RefElements;
  vector<size_t> bix;
  auto gv = grid.leafGridView();
  for(const auto& el : Dune::elements(gv)){
    if(!el.hasBoundaryIntersections()) 
       continue;
    const auto& refEl = Dune::referenceElement(el);
    //const auto dimW = 2;
    for(auto& is : Dune::intersections(gv,el)){
      if(is.boundary()){
        auto inside = is.indexInInside();
        auto faceSize = refEl.size(inside,/*face*/ 1,/*node*/ 2);
        for(int i = 0; i < faceSize; ++i){
          // this is a 2 grid where faces=cells=codim 1 nodes is of codim 2
          auto corner = refEl.subEntity(/*facenr*/inside,/*face*/1, /*nodenum*/i, /*node*/2);
          auto cornerIndex = gv.indexSet().subIndex(el,corner,2);
          bix.push_back(cornerIndex);
        }
      }
    }
  }
  // make unique
  std::sort(bix.begin(), bix.end());
  bix.erase( unique( bix.begin(), bix.end() ), bix.end() );

  // ensure correct order (we use convex hull algorithm for that)
  const vector<CoordType> ncoords = node_coordinates(grid);
  std::array<int, 2> plane = best2Dplane(ncoords);
  vector<CoordType> bcoords;
  for(auto i : bix)
    bcoords.push_back(ncoords[i]);

  vector<double> bpts2D = projectPointsTo2DPlane(bcoords, plane);
  vector<size_t> bix_ordered = convex_hull(bpts2D);
  assert(bix_ordered.size() == bix.size()); // should be the case if boundary is convex
  for (auto& i : bix_ordered)
    i = bix[i];
  
  return bix_ordered;
}
  
  // ----------------------------------------------------------------------------
vector<size_t> GridStretcher::boundary_node_indices_old(const Grid& grid)
// ----------------------------------------------------------------------------  
{
  const auto view = grid.leafGridView();


  // register all node coordinates
  vector<CoordType> vcoords(node_coordinates(grid));

  // determine boundary edges
  vector<int> count(view.size(1), 0); // number of edges (codim 1)

  for (const auto& elem : elements(view))
    for (size_t i = 0; i != elem.subEntities(1); ++i)
      count[view.indexSet().index(elem.subEntity<1>(i))] += 1;

  // @@ This is suboptimal and brittle - there must be a better way to identify
  // boundary nodes than by geometric comparison!

  set<size_t> bix;
  for(const auto& ep : Dune::edges(view)){
  //for (auto ep = view.begin<1>(); ep != view.end<1>(); ++ep, ++ix)
    int edge_index = view.indexSet().index(ep);
    if (count[edge_index] < 2){
      for (int i = 0; i != 2; ++i) {
        const auto pt = ep.geometry().corner(i);
        const size_t pt_ix = find_coord_in(pt, vcoords);

        assert(pt_ix < vcoords.size());
        bix.insert(pt_ix);
      }
    }
  }
  return vector<size_t>(bix.begin(), bix.end());

}

// ----------------------------------------------------------------------------  
vector<size_t> GridStretcher::complement_of(const vector<size_t>& vec,
                                            const size_t N)
// ----------------------------------------------------------------------------  
{
  // indicator vector to identify the indices mentioned in 'vec'
  vector<int> tmp(N, 0);
  for (const auto v : vec)
    tmp[v] = 1;

  // indices not mentioned in 'vec' are stored in 'result'
  vector<size_t> result;
  for (size_t ix = 0; ix != N; ++ix)
    if (tmp[ix] == 0)
      result.push_back(ix);

  return result;
}

// ----------------------------------------------------------------------------
vector<CoordType> GridStretcher::node_coordinates(const Grid& g)
// ----------------------------------------------------------------------------
{
  vector<CoordType> result;
  for (const auto& v : vertices(g.leafGridView()))
    result.push_back(v.geometry().corner(0));
  return result;
}

// ----------------------------------------------------------------------------
BoundaryNormals GridStretcher::boundary_normals(const Grid& grid,
                                                const CellBnodeMap& c2bix,
                                                const std::vector<size_t>& bcindices,
                                                const std::vector<size_t>& bnindices,
                                                const vector<CoordType>& nodecoords)
// ----------------------------------------------------------------------------
{ 
  const size_t N = size(bcindices);
  // map amounts to displacement vectors for each boundary point

  BoundaryNormals result {vector<CoordType> {N, {0, 0, 0}},    // nodenormals
                          vector<CoordType> {N, {0, 0, 0}},    // cell "normals"
                          vector<CoordType> {N, {0, 0, 0}}};   // edge normals
  int pos = 0;
  for (size_t bix : bcindices) {
    // compute directional vector from cell centroid to edge centroid
    const auto entry = c2bix.find(bix)->second;
    const auto elem = grid.entity(get<0>(entry));
    const auto ccenter = elem.geometry().center();
    const auto ecenter = (nodecoords[bnindices[get<1>(entry)]] +
                          nodecoords[bnindices[get<2>(entry)]]) / 2;
    const auto dvec = normalize(ecenter - ccenter); // cell centroid to edge centroid
    result.bcell_normals[pos] = dvec;

    // compute outward edge normal
    const auto tangent = nodecoords[bnindices[get<1>(entry)]] -
                         nodecoords[bnindices[get<2>(entry)]];
    const auto enormal = normalize(cross(cross(tangent, dvec), tangent));
    result.bedge_normals[pos] = enormal;

    // accumulate result for average node normal
    result.bnode_normals[get<1>(entry)] += enormal;
    result.bnode_normals[get<2>(entry)] += enormal;

    pos++;
  }

  // each entry is the sum of two normals, so we need to normalize again
  transform(result.bnode_normals.begin(), result.bnode_normals.end(),
            result.bnode_normals.begin(),
            [](const CoordType& c) {return normalize(c);});
  
  return result;
}

// ----------------------------------------------------------------------------  
std::vector<double>
GridStretcher::bcentroid_param_mat(const Grid& grid,
                                   const std::vector<size_t>& bnindices,
                                   const std::vector<size_t>& iindices,
                                   const std::vector<double>& iparam,
                                   const CellBnodeMap& c2bix)
// ----------------------------------------------------------------------------    
{
  const auto view = grid.leafGridView();

  // helper function to identify the internal node and two boundary nodes of a
  // given boundary triangle
  auto corner_nodes = [&] (const CellBnodeMap::value_type& pair) {
    const size_t bn1 = bnindices[get<1>(pair.second)]; // boundary node 1
    const size_t bn2 = bnindices[get<2>(pair.second)]; // boundary node 2
    const auto elem = grid.entity(get<0>(pair.second));
    
    assert(elem.subEntities(2) == 3); // should have three corners
    size_t in(0);
    for (int i = 0; i != 3; ++i) {
      in = view.indexSet().index(elem.subEntity<2>(i));
      if (in != bn1 && in != bn2) break;
    }
    return array<size_t, 3> {in, bn1, bn2};
  };

  // computing the bcentroid parameter matrix
  vector<double> result;
  const size_t Nb = bnindices.size();
  for (const auto& pair : c2bix) {

    const auto cnodes = corner_nodes(pair);
    const size_t iix =
      find(iindices.begin(), iindices.end(), cnodes[0]) - iindices.begin();

    // write in one-third of the parameterization of the internal point
    transform(&iparam[iix * Nb], &iparam[(iix+1) * Nb], back_inserter(result),
              [] (const double p) {return p/3;});

    // add one-third of each of the two boundary points
    for (int i = 1; i != 3; ++i)
      result[result.size() - Nb + cnodes[i]] += 1.0/3.0;
    
  }
  return result;
}

// ----------------------------------------------------------------------------
double GridStretcher::objective(const std::vector<double>& bndisp,
                                const std::vector<double>& dtarget,
                                std::vector<double>& grad,
                                bool fixed_cell_centroids)
// ----------------------------------------------------------------------------
{
  const vector<CoordType>& ncoords = nodecoords();
  const vector<CoordType>& normals = bnodenormals();
  const size_t Nb = bnindices_.size(); // number of boundary nodes (and cells)
  assert(bndisp.size() == Nb);
  
  // compute boundary node positions
  vector<CoordType> bnodes, bnodes0;
  auto disp_iter = bndisp.begin();
  for (auto bix : bnindices_) {
    bnodes0.push_back(ncoords[bix]);
    bnodes.push_back(ncoords[bix] + normals[bix] * *disp_iter++);
  }
  const vector<CoordType>& bnodes_for_centroids = fixed_cell_centroids ? bnodes0 : bnodes;
  // compute cell and face centroid positions, and distance vector
  vector<CoordType> cell_centroids;
  vector<CoordType> face_centroids;
  vector<CoordType> distance;
  auto cpar_iter = bcentroid_param_.begin();
  for (const auto& pair : c2bix_) {
    const size_t bn1 = get<1>(pair.second); // boundary node 1
    const size_t bn2 = get<2>(pair.second); // boundary node 2

    // compute cell centroid as a linear combination of all boundary nodes
    CoordType cc {0, 0, 0};
    for (size_t i = 0; i != Nb; ++i)
      cc += bnodes_for_centroids[i] * *cpar_iter++;
    cell_centroids.push_back(cc);

    // compute face centroid
    face_centroids.push_back(0.5 * (bnodes[bn1] + bnodes[bn2]));

    // compute distance vector between face and cell centroid
    distance.push_back(face_centroids.back() - cell_centroids.back());
  }

  // compute objective value
  double objval = 0;
  for (size_t i = 0; i != distance.size(); ++i)
    // o = o + (d^2 - dtarget^2) ^ 2
    objval += pow(distance[i].two_norm() - dtarget[i], 2);

  objval /= 2; 

  // compute gradient
  grad.resize(Nb);
  fill(grad.begin(), grad.end(), 0);
  auto c2bix_iter = c2bix_.begin();
  for (size_t i = 0; i != Nb; ++i, ++c2bix_iter) { // loop over cells
    const double dfac = (distance[i].two_norm() - dtarget[i]) / distance[i].two_norm();
    for (size_t j = 0; j != Nb; ++j) {// loop over boundary nodes
      // efac is zero unless this boundary node is part of the current cell
      const double efac = (j == get<1>(c2bix_iter->second) ||
                           j == get<2>(c2bix_iter->second)) ? 0.5 : 0.0;
      // derivative of centroid position with respect to boundary node 'j'
      const double cpar = fixed_cell_centroids ? 0 : bcentroid_param_[i*Nb + j]; //d(d_i)/d(b_j);
      const double m = (efac - cpar);
      for (size_t d = 0; d != 3; ++d) // loop over dimensions
        grad[j] += dfac * distance[i][d] * m * normals[j][d];
    }
  }
  
  return objval;
}

  
// ----------------------------------------------------------------------------  
vector<double> GridStretcher::interior_parametrization(const vector<size_t>& bindices,
                                                              const vector<size_t>& iindices,
                                                              const vector<CoordType>& coords3D)
// ----------------------------------------------------------------------------  
{
  // pick two out of the three coordinate axes to use for computing the 2D parametrization
  const vector<double> coords2D(pickCoords2D(coords3D));

  // sort node coordinates into boundary and interior
  const vector<double> bcoords(pick_nodes(bindices, coords2D));
  const vector<double> icoords(pick_nodes(iindices, coords2D));

  // call parametrization function
  vector<double> result;
  parametrize_interior_2D(bcoords, icoords, result);

  return result;
}

// ----------------------------------------------------------------------------  
GridStretcher::CellBnodeMap
GridStretcher::compute_cell_2_bindices_mapping(const Grid& grid,
                                               const vector<size_t>& bindices)
// ----------------------------------------------------------------------------
{
  CellBnodeMap result;
  const auto gview = grid.leafGridView();

  // coordinates of boundary indices @@ We really ought to be able to do this
  // without resorting to geometric calculations!
  // register all node coordinates
  const vector<CoordType> ncoords(node_coordinates(grid));
  const vector<CoordType> bcoords(extract_elements(ncoords, bindices));
  
  const ElementMapper mapper(gview, Dune::mcmgElementLayout());
  for (const auto& elem: elements(gview))
    for (const auto& is : intersections(gview, elem)) 
      if (is.boundary()) 
        // this is a boundary cell
        result[mapper.index(elem)] =
          {elem.seed(),
           find_coord_in(is.geometry().corner(0), bcoords),
           find_coord_in(is.geometry().corner(1), bcoords)};
  return result;
}

// ----------------------------------------------------------------------------
vector<double>
GridStretcher::computeBoundaryNodeDisplacements(const vector<double>& amounts,
                                                const vector<CoordType> bnodenormals) const
// ----------------------------------------------------------------------------
{
  const size_t N = bcindices_.size();

  const vector<CoordType>& bnnorm = bnodenormals.empty() ?
                                    boundary_normals_.bnode_normals :
                                    bnodenormals;
  assert(bnnorm.size() == N);
  
  vector<double> node_amounts(N, -1 * std::numeric_limits<double>::infinity());

  size_t cur_bcell_ix = 0;
  for (const size_t bix : bcindices_) {
    const auto entry = c2bix_.find(bix)->second;
    const size_t bnodeix1 = get<1>(entry);
    const size_t bnodeix2 = get<2>(entry);
    
    const auto& enorm = boundary_normals_.bedge_normals[cur_bcell_ix];
    const auto& cnorm = boundary_normals_.bcell_normals[cur_bcell_ix];
    const auto& bn1 = bnnorm[bnodeix1];
    const auto& bn2 = bnnorm[bnodeix2];
    
    const double l1 = amounts[cur_bcell_ix] * enorm.dot(cnorm) / enorm.dot(bn1);
    const double l2 = amounts[cur_bcell_ix] * enorm.dot(cnorm) / enorm.dot(bn2);
    
    node_amounts[bnodeix1] = max(node_amounts[bnodeix1], l1); // @@ correct for negative values?
    node_amounts[bnodeix2] = max(node_amounts[bnodeix2], l2);
    cur_bcell_ix++;
  }

  return node_amounts;
}
  
// ----------------------------------------------------------------------------  
void GridStretcher::expandBoundaryCells(const vector<double>& amounts)
                                        
// ----------------------------------------------------------------------------
{
  vector<double> distance =
    computeBoundaryNodeDisplacements(amounts,
                                    boundary_normals_.bnode_normals);

  vector<CoordType> displacements;
  for (size_t i = 0; i != amounts.size(); ++i)
    displacements.push_back(distance[i] * boundary_normals_.bnode_normals[i]);

  applyBoundaryNodeDisplacements(displacements);

}

// ----------------------------------------------------------------------------
void GridStretcher::applyBoundaryNodeDisplacements(const vector<CoordType>& disp)
// ----------------------------------------------------------------------------  
{
  assert(disp.size() == bnindices_.size());

  vector<CoordType> ncoords(node_coordinates(grid_)); 
  
  // compute new boundary node coordinates
  for (size_t i = 0; i != bnindices_.size(); ++i) 
    ncoords[bnindices_[i]] += disp[i];

  // reset all internal nodes to zero
  for (size_t iix : iindices_)
    ncoords[iix] = 0; 
  
  // compute new coordinates for internal nodes
  auto ipiter = iparam_.begin();
  for (size_t iix : iindices_)
    for (size_t bix = 0; bix != bnindices_.size(); ++bix)
      ncoords[iix] +=
        ncoords[bnindices_[bix]] * (*ipiter++);

  // write all new node positions to mesh
  size_t vcount = 0;
  for (auto& vertex : vertices(grid_.leafGridView()))
    grid_.setPosition(vertex, ncoords[vcount++]);

  // recompute stored geometric information
  nodecoords_ = node_coordinates(grid_);
  boundary_normals_ = boundary_normals(grid_, c2bix_, bcindices_, bnindices_, nodecoords_);
}


// ----------------------------------------------------------------------------
std::vector<double> GridStretcher::centroidEdgeDist() const
// ----------------------------------------------------------------------------  
{
  const size_t N = size(bcindices_);
  vector<double> result(N);
  
  for (size_t bc = 0; bc != N; ++bc) { 
    const auto entry = c2bix_.find(bcindices_[bc])->second;
    const auto elem = grid_.entity(get<0>(entry));
    const auto ccenter = elem.geometry().center();
    const auto ecenter = (nodecoords_[bnindices_[get<1>(entry)]] +
                          nodecoords_[bnindices_[get<2>(entry)]]) / 2;
    result[bc] = (ccenter - ecenter).two_norm();
  }
  return result;
}

// ----------------------------------------------------------------------------
double GridStretcher::maxBoxLength() const
// ----------------------------------------------------------------------------  
{
  assert(nodecoords_.size() > 0);
  std::array<double, 3> low {nodecoords_[0][0], nodecoords_[0][1], nodecoords_[0][2]};
  std::array<double, 3> high {low};

  for(const auto n : nodecoords_)
    for (int d = 0; d != 3; ++d) {
      low[d] = std::min(low[d], n[d]);
      high[d] = std::max(high[d], n[d]);
    }

  for (int d = 0; d != 3; ++d)
    high[d] = high[d] - low[d];

  return *std::max_element(high.begin(), high.end());
}

// ----------------------------------------------------------------------------
void GridStretcher::adjustToConvex(std::vector<double>& disp,
                                   std::vector<double>& total_disp,
                                   const std::vector<CoordType>& dirs) const
// ----------------------------------------------------------------------------
{
  const size_t N = bnindices_.size();
  assert(dirs.size() == N);
  assert(disp.size() == N);
  assert(total_disp.size() == N);

  // compute new boundary node coordinates, before convexity is enforced
  vector<CoordType> old_bcoords(dirs.size()), new_bcoords(dirs.size());
  for (size_t i = 0; i != N; ++i) {
    old_bcoords[i] = nodecoords_[bnindices_[i]];
    new_bcoords[i] = old_bcoords[i] + dirs[i] * disp[i];
  }

  // reduce geometry to 2D, to allow convex hull algorithm to work
  const array<int, 2> plane = best2Dplane(new_bcoords);
  const vector<double> bpts2D = projectPointsTo2DPlane(new_bcoords, plane);
  const vector<double> bpts2D_old = projectPointsTo2DPlane(old_bcoords, plane);
  const vector<double> dirs2D = projectPointsTo2DPlane(dirs, plane);

  // compute convex hull of boundary points
  const vector<size_t> cvhull_pts = convex_hull(bpts2D);
  assert(cvhull_pts.size() > 2);

  // adjust positions of points that are inside the convex hull so that they
  // reach the convex hull boundary, ensuring that the fracture outline will be
  // convex
  for (size_t i = 0; i != cvhull_pts.size(); ++i) {
    const size_t cstart = cvhull_pts[i];
    const size_t cend = cvhull_pts[(i+1) % cvhull_pts.size()];

    // adjust boundary points in the interval between bp_start_ix and bp_end_ix
    for (size_t j = (cstart + 1)%N; j != cend; j = (j + 1) % N) {
      const double alpha = adjust_disp_to_convex(&bpts2D[2*cstart],
                                                 &bpts2D[2*cend],
                                                 &bpts2D_old[2*j],
                                                 &dirs2D[2*j]);
      const double diff = alpha - disp[j];
      disp[j] = alpha;
      total_disp[j] += diff;
    }
  }
}


// ----------------------------------------------------------------------------
void GridStretcher::dumpToVTK(const char* filename, const std::vector<std::vector<double>> data) const
// ----------------------------------------------------------------------------  
{
  //std::cout << "Hello" << std::endl;
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_.leafGridView(),
                                                          Dune::VTK::nonconforming);

  for (int i = 0; i != data.size(); ++i) 
    vtkwriter->addCellData(data[i], "data" + std::to_string(i));
      
  vtkwriter->write(filename);
}

  
}; // end namespace Opm

