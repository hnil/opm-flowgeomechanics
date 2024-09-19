#include <opm/geomech/GridStretcher.hpp>
#include <opm/geomech/param_interior.hpp>
#include <array>
#include <assert.h>
#include <dune/common/fvector.hh> // FieldVector
#include <dune/grid/common/mcmgmapper.hh> // for element mapper

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
double dist2D(const Dune::FieldVector<double, 3>& p1,
              const Dune::FieldVector<double, 3>& p2)
// ----------------------------------------------------------------------------    
{
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  return sqrt(dx*dx + dy*dy);
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
            [&] (const CoordType& el) {return dist2D(el, c) < TOL;})
    - cvec.begin();

}

// ----------------------------------------------------------------------------
CoordType normalize(const CoordType& vec) { return vec / vec.two_norm(); }
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
const vector<CoordType> compute_bnode_displacements(const vector<double>& amounts,
                                                    const vector<size_t>& bcindices,
                                                    const Opm::GridStretcher::CellBnodeMap& c2bix,
                                                    const Opm::GridStretcher::Grid& grid,
                                                    const vector<CoordType>& nodecoords)
// ----------------------------------------------------------------------------
{
  const size_t N = size(bcindices);
  // map amounts to displacement vectors for each boundary point
  assert(amounts.size() == N); // one scalar per boundary cell
  vector<vector<tuple<double, CoordType>>> disp(N); // node displacements
  
  auto aiter = amounts.begin();
  for (size_t bc = 0; bc != amounts.size(); ++bc) {
    // compute directional vector
    const auto entry = c2bix.find(bcindices[bc])->second;
    const auto elem = grid.entity(get<0>(entry));
    const auto ccenter = elem.geometry().center();
    const auto ecenter = (nodecoords[get<1>(entry)] + nodecoords[get<2>(entry)]) / 2;
    const auto nvec = normalize(ecenter - ccenter);
    disp[get<1>(entry)].push_back({*aiter, nvec});
    disp[get<2>(entry)].push_back({*aiter++, nvec});
  }

  // check that every boundary node has got two displacements associated with it
    
  // combine displacement vectors
  vector<CoordType> result(N);
  for (size_t i = 0; i != N; ++i) {
    const auto& entry = disp[i];
    const double a1 = get<0>(entry[0]);
    const double a2 = get<0>(entry[1]);
    const auto& v1 = get<1>(entry[0]);
    const auto& v2 = get<1>(entry[1]);
    
    const auto dir = normalize(v1 + v2);

    //result[i] = dir * max(dir.dot(v1) * a1, dir.dot(v2) * a2);
    result[i] = dir *
      ((fabs(dir.dot(v1) * a1) > fabs(dir.dot(v2) * a2)) ? dir.dot(v1) * a1 :
                                                           dir.dot(v2) * a2);
  }
  return result;
}
  
}; // end anonymous namespace

// ============================================================================
namespace Opm
// ============================================================================
{    

  // ----------------------------------------------------------------------------
vector<size_t> GridStretcher::boundary_node_indices(const Grid& grid)
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

  int ix = 0;
  set<size_t> bix;
  for (auto ep = view.begin<1>(); ep != view.end<1>(); ++ep, ++ix)
    if (count[ix] < 2)
      for (int i = 0; i != 2; ++i) {
        const auto pt = ep->geometry().corner(i);
        const size_t pt_ix = find_coord_in(pt, vcoords);

        if (pt_ix < vcoords.size())
          bix.insert(pt_ix);
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
vector<CoordType> GridStretcher::boundary_normals(const Grid& grid,
                                                  const CellBnodeMap& c2bix,
                                                  const std::vector<size_t>& bcindices,
                                                  const vector<CoordType>& nodecoords)
// ----------------------------------------------------------------------------
{ 
  const size_t N = size(bcindices);
  // map amounts to displacement vectors for each boundary point

  vector<CoordType> result {N, {0, 0, 0}};
  
  for (size_t bix : bcindices) {
    // compute directional vector
    const auto entry = c2bix.find(bix)->second;
    const auto elem = grid.entity(get<0>(entry));
    const auto ccenter = elem.geometry().center();
    const auto ecenter = (nodecoords[get<1>(entry)] + nodecoords[get<2>(entry)]) / 2;
    const auto nvec = normalize(ecenter - ccenter);
      
    result[get<1>(entry)] += nvec;
    result[get<2>(entry)] += nvec; 
  }

  // each entry is the sum of two normals, so we need to normalize again
  transform(result.begin(), result.end(), result.begin(),
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
    const size_t bn1 = get<1>(pair.second); // boundary node 1
    const size_t bn2 = get<2>(pair.second); // boundary node 2
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
void GridStretcher::expandBoundaryCells(const vector<double>& amounts)
// ----------------------------------------------------------------------------
{
  applyBoundaryNodeDisplacements(compute_bnode_displacements(amounts,
                                                             bcindices_,
                                                             c2bix_,
                                                             grid_,
                                                             nodecoords_));
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
  bnode_normals_ = boundary_normals(grid_, c2bix_, bcindices_, nodecoords_);
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
    const auto ecenter = (nodecoords_[get<1>(entry)] + nodecoords_[get<2>(entry)]) / 2;
    result[bc] = (ccenter - ecenter).two_norm();
  }
  return result;
}


}; // end namespace Opm

