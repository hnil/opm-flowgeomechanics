#include <opm/geomech/GridStretcher.hpp>
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
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

  
// ----------------------------------------------------------------------------  
vector<double> pick_nodes(const vector<size_t>& ixs, vector<double>& coords2D)
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

  array<size_t, 3> low {coords[0], coords[1], coords[2]}, high {low};
  for (size_t i = 0; i != N ++i) {
    for (int d = 0; d != 3; ++d) {
      low[d] = std::min(low[d], coords[3*i + d]);
      high[d] = std::max(high[d], coords[3*i + d]);
    }
  }
  const array<size_t, 3> span {high[0] - low[0], high[1] - low[1], high[2] - low[2] };
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
  
}; // end anonymous namespace

// ============================================================================
namespace Opm
// ============================================================================
{    

  // ----------------------------------------------------------------------------
static vector<size_t> GridStretcher::boundary_node_indices(const Grid& grid)
// ----------------------------------------------------------------------------  
{
  const auto view = g.leafGridView();

  //using FVec = Dune::FieldVector<double, 3>;
  using CoordType = Dune::FieldVector<double, 3>;

  // register all node coordinates
  vector<CoordType> vcoords;
  for (auto v : vertices(view)) 
    vcoords.push_back(v.geometry().corner(0));

  // determine boundary edges
  vector<int> count(view.size(1), 0); // number of edges (codim 1)

  for (const auto& elem : elements(view))
    for (int i = 0; i != elem.subEntities(1); ++i)
      count[view.indexSet().index(elem.subEntity<1>(i))] += 1;

  // @@ This is suboptimal and brittle - there must be a better way to identify
  // boundary nodes than by geometric comparison!

  int ix = 0;
  set<int> bix;
  const double TOL = 1e-3;
  for (auto ep = view.begin<1>(); ep != view.end<1>(); ++ep, ++ix)
    if (count[ix] < 2)
      for (int i = 0; i != 2; ++i) {
        auto pt = ep->geometry().corner(i);
        bool found = false;
        for (int j = 0; j != vcoords.size(); ++j && !found)
          if (dist2D(pt, vcoords[j]) < TOL) {
            found = true;
            bix.insert(j);
          }
      }

  return vector(bix.begin(), bix.end());

}

// ----------------------------------------------------------------------------  
static vector<size_t> GridStretcher::complement_of(const vector<size_t>& vec,
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
static vector<double> GridStretcher::nodecoords(const Grid& grid)
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
static vector<double> GridStretcher::interior_parametrization(const vector<size_t>& bindices,
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
static CellBnodeMap GridStretcher::compute_cell_2_bindices_mapping(
                                                       const Grid& grid,
                                                       const vector<size_t>& bindices)
// ----------------------------------------------------------------------------
{
  CellBnodeMap result;
  const auto gview = grid.leafGridView();

  const ElementMapper mapper(gview, Dune::mcmgElementLayout());
  for (const auto& elem: elements(gview)){
    for (const auto& is : intersections(gview, elem)) {
      if (is.boundary()) {
        // this is boundary cell
        const size_t cell_ix = mapper.index(elem);
        const size_t n1_ix = mapper.index(is.geometry().corner(0));
        const size_t n2_ix = mapper.index(is.geometry().corner(1));

        // @@ currently, only one boundary edge per cell is supported
        assert(result.find(cell_ix) == result.end()); 

        // store the association of this cell index with these two corner indices
        result[cell_ix] = make_pair(n1_ix, n2_ix);
      }
    }
  }
  return result;
}

// ----------------------------------------------------------------------------  
void GridStretcher::expandBoundaryCells(const vector<double>& amounts)
// ----------------------------------------------------------------------------
{

}
  
}; // end namespace Opm
