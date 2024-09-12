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
vector<CoordType> node_coordinates(const Opm::GridStretcher::Grid& g)
// ----------------------------------------------------------------------------
{
  vector<CoordType> result;
  for (auto v : vertices(g.leafGridView()))
    result.push_back(v.geometry().corner(0));
  return result;
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
vector<double> pickCoords2D(const vector<double>& coords3D)
// ----------------------------------------------------------------------------
{
  // pick the two coordinate axes with the largest spread
  const size_t N = coords3D.size() / 3; // total number of nodes (interior and boundary)

  array<double, 3> low {coords3D[0], coords3D[1], coords3D[2]}, high {low};
  for (size_t i = 0; i != N; ++i) {
    for (int d = 0; d != 3; ++d) {
      low[d] = std::min(low[d], coords3D[3*i + d]);
      high[d] = std::max(high[d], coords3D[3*i + d]);
    }
  }
  const array<double, 3> span {high[0] - low[0], high[1] - low[1], high[2] - low[2] };
  const size_t min_ix = distance(span.begin(), min_element(span.begin(), span.end()));
  const size_t ix1 = (min_ix + 1) % 3;
  const size_t ix2 = (min_ix + 2) % 3;

  vector<double> coords2D;
  for (size_t i = 0; i != N; ++i) {
    coords2D.push_back(coords3D[i*3 + ix1]);
    coords2D.push_back(coords3D[i*3 + ix2]);
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
                                                    const Opm::GridStretcher::Grid& grid)
// ----------------------------------------------------------------------------
{
  const size_t N = size(bcindices);
  // map amounts to displacement vectors for each boundary point
  assert(amounts.size() == N); // one scalar per boundary cell
  vector<vector<tuple<double, CoordType>>> disp(N); // node displacements
  
  const auto nodecoords = node_coordinates(grid); // @@
  for (size_t bc = 0; bc != amounts.size(); ++bc) {
    // compute directional vector
    const auto entry = c2bix.find(bcindices[bc])->second;
    const auto elem = grid.entity(get<0>(entry));
    const auto ccenter = elem.geometry().center();
    const auto ecenter = (nodecoords[get<1>(entry)] + nodecoords[get<2>(entry)]) / 2;
    const auto nvec = normalize(ecenter - ccenter);
    disp[get<1>(entry)].push_back({amounts[bc], nvec});
    disp[get<2>(entry)].push_back({amounts[bc], nvec});
  }

  // check that every boundary node has got two displacements associated with it
  assert(std::all_of(disp.begin(), disp.end(),
                     [](const vector<tuple<double, CoordType>>& v)
                     { return v.size() == 2;}));
    
  // combine displacement vectors
  vector<CoordType> result(N);
  return result;
  for (size_t i = 0; i != N; ++i) {
    const auto& entry = disp[i];
    const double a1 = get<0>(entry[0]);
    const double a2 = get<0>(entry[1]);
    const auto& v1 = get<1>(entry[0]);
    const auto& v2 = get<1>(entry[1]);
    
    const auto dir = normalize(v1 + v2);

    result[i] = dir * max(dir.dot(v1) * a1, dir.dot(v2) * a2);
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
vector<double> GridStretcher::nodecoords(const Grid& grid)
// ----------------------------------------------------------------------------  
{
  // store coordinates as [x1, y1, z1, ...., xn, yn, zn, .... xlast, ylast, zlast]
  vector<double> coords;
  for (const auto& vertex : vertices(grid.leafGridView()))
    for (int dim = 0; dim != 3; ++dim)
      coords.push_back(vertex.geometry().corner(0)[dim]);

  return coords;
}
  

// ----------------------------------------------------------------------------  
vector<double> GridStretcher::interior_parametrization(const vector<size_t>& bindices,
                                                              const vector<size_t>& iindices,
                                                              const vector<double>& coords3D)
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
                                                             grid_));
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
  for (size_t iix : iindices_)
    for (size_t bix = 0; bix != bnindices_.size(); ++bix)
      ncoords[iix] +=
        ncoords[bnindices_[bix]] * iparam_[iix * bnindices_.size() + bnindices_[bix]];

  // write all new node positions to mesh
  size_t vcount = 0;
  for (auto& vertex : vertices(grid_.leafGridView()))
    grid_.setPosition(vertex, ncoords[vcount++]);
}

  
}; // end namespace Opm
