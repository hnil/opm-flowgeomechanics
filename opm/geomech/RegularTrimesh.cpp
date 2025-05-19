#include "RegularTrimesh.hpp"
#include "GeometryHelpers.hpp"
#include <algorithm>
#include <fstream>
#include <opm/geomech/RegularTrimesh.hpp>
#include <set>
#include <tuple>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Opm;
using namespace std;

namespace
{

array<EdgeRef, 3>
cell2edges(const CellRef& cell)
{
    return (cell[2] == 0) ? array<EdgeRef, 3> {EdgeRef {cell[0], cell[1], 0},
                                               EdgeRef {cell[0], cell[1], 1},
                                               EdgeRef {cell[0], cell[1], 2}}
                          : array<EdgeRef, 3> {EdgeRef {cell[0], cell[1] + 1, 0},
                                               EdgeRef {cell[0] + 1, cell[1], 1},
                                               EdgeRef {cell[0], cell[1], 2}};
}

array<CellRef, 2>
edge2cells(const EdgeRef& edge)
{
    return (edge[2] == 0) ? array<CellRef, 2> {CellRef {edge[0], edge[1], 0}, CellRef {edge[0], edge[1] - 1, 1}}
        : (edge[2] == 1)  ? array<CellRef, 2> {CellRef {edge[0], edge[1], 0}, CellRef {edge[0] - 1, edge[1], 1}}
                          : array<CellRef, 2> {CellRef {edge[0], edge[1], 1}, CellRef {edge[0], edge[1], 0}};
}

array<CellRef, 3>
cellNeighbors(const CellRef& cell)
{
    return (cell[2] == 0) ? array<CellRef, 3> {CellRef {cell[0], cell[1], 1},
                                               CellRef {cell[0] - 1, cell[1], 1},
                                               CellRef {cell[0], cell[1] - 1, 1}}
                          : array<CellRef, 3> {CellRef {cell[0], cell[1], 0},
                                               CellRef {cell[0] + 1, cell[1], 0},
                                               CellRef {cell[0], cell[1] + 1, 0}};
};

// += operator for Opm::Coord3D
Opm::Coord3D&
operator+=(Opm::Coord3D& lhs, const Opm::Coord3D& rhs)
{
    for (int i = 0; i != 3; ++i)
        lhs[i] += rhs[i];
    return lhs;
}

// *= operator for Opm::Coord3D
Opm::Coord3D&
operator*=(Opm::Coord3D& lhs, const double rhs)
{
    for (int i = 0; i != 3; ++i)
        lhs[i] *= rhs;
    return lhs;
}

// /= operator for Opm::Coord3D
Opm::Coord3D&
operator/=(Opm::Coord3D& lhs, const double rhs)
{
    for (int i = 0; i != 3; ++i)
        lhs[i] /= rhs;
    return lhs;
}

// + operator for Opm::Coord3D
Opm::Coord3D
operator+(const Opm::Coord3D& lhs, const Opm::Coord3D& rhs)
{
    Opm::Coord3D result(lhs);
    result += rhs;
    return result;
}

// / operator for Opm::Coord3D
Opm::Coord3D
operator/(const Opm::Coord3D& lhs, const double rhs)
{
    Opm::Coord3D result(lhs);
    result /= rhs;
    return result;
}

// * operator for Opm::Coord3D
Opm::Coord3D
operator*(const Opm::Coord3D& lhs, const double rhs)
{
    Opm::Coord3D result(lhs);
    result *= rhs;
    return result;
}

// * operator for Opm::Coord3D with double on lhs
Opm::Coord3D
operator*(const double lhs, const Opm::Coord3D& rhs)
{
    return rhs * lhs;
};

// == operator for Opm::NodeRef
bool
operator==(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs)
{
    return lhs[0] == rhs[0] && lhs[1] == rhs[1];
}

// < operator for Opm::NodeRef
bool
operator<(const Opm::NodeRef& lhs, const Opm::NodeRef& rhs)
{
    return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

// == operator for Opm::EdgeRef
bool
operator==(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs)
{
    return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
}

// < operator for Opm::EdgeRef
bool
operator<(const Opm::EdgeRef& lhs, const Opm::EdgeRef& rhs)
{
    return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1])
        || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
}

// ----------------------------------------------------------------------------
array<EdgeRef, 3>
half_edges(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    return (cell[2] == 0) ? array<EdgeRef, 3> {EdgeRef {cell[0], cell[1], 0},
                                               EdgeRef {cell[0], cell[1], 1},
                                               EdgeRef {cell[0], cell[1], 2}}
                          : array<EdgeRef, 3> {EdgeRef {cell[0], cell[1] + 1, 0},
                                               EdgeRef {cell[0] + 1, cell[1], 1},
                                               EdgeRef {cell[0], cell[1], 2}};
}

// ----------------------------------------------------------------------------
array<CellRef, 3>
neigh_cells(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    return (cell[2] == 0) ? array<CellRef, 3> {CellRef {cell[0], cell[1] - 1, 1},
                                               CellRef {cell[0] - 1, cell[1], 1},
                                               CellRef {cell[0], cell[1], 1}}
                          :

                          array<CellRef, 3> {CellRef {cell[0], cell[1] + 1, 0},
                                             CellRef {cell[0] + 1, cell[1], 0},
                                             CellRef {cell[0], cell[1], 0}};
}


// ----------------------------------------------------------------------------
array<bool, 3>
identify_boundary(const CellRef& cell, const RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    // return 'true' for the triangle edges that lie on the boundary
    const auto ncells = neigh_cells(cell);
    return {!mesh.isActive(ncells[0]), !mesh.isActive(ncells[1]), !mesh.isActive(ncells[2])};
}

// ----------------------------------------------------------------------------
RegularTrimesh
remove_mesh_corners(const RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    RegularTrimesh result(mesh);

    bool redo = true;

    while (redo) {
        redo = false;
        for (const auto& cell : result.boundaryCells()) {
            const auto bnd = identify_boundary(cell, result);
            if (count(bnd.begin(), bnd.end(), true) > 1) {
                result.setInactive(cell);
                redo = true;
            }
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
vector<array<NodeRef, 3>>
single_split(const CellRef& cell, const array<bool, 3>& bnd)
// ----------------------------------------------------------------------------
{
    assert(bnd[0] + bnd[1] + bnd[2] == 1);

    const int i = cell[0];
    const int j = cell[1];

    const NodeRef A = (cell[2] == 0) ? NodeRef {2 * i, 2 * j} : NodeRef {2 * (i + 1), 2 * (j + 1)};
    const NodeRef B = (cell[2] == 0) ? NodeRef {2 * (i + 1), 2 * j} : NodeRef {2 * i, 2 * (j + 1)};
    const NodeRef C = (cell[2] == 0) ? NodeRef {2 * i, 2 * (j + 1)} : NodeRef {2 * (i + 1), 2 * j};

    const array<NodeRef, 3> tri = bnd[0] ? array<NodeRef, 3> {A, B, C}
        : bnd[1]                         ? array<NodeRef, 3> {C, A, B}
                                         : array<NodeRef, 3> {B, C, A};
    const NodeRef X = {(tri[0][0] + tri[1][0]) / 2, (tri[0][1] + tri[1][1]) / 2};

    return vector<array<NodeRef, 3>> {array<NodeRef, 3> {tri[0], X, tri[2]}, array<NodeRef, 3> {X, tri[1], tri[2]}};
}



// ----------------------------------------------------------------------------
vector<array<NodeRef, 3>>
tesselate_coarsecell(const CellRef& cell, const RegularTrimesh& mesh, const bool skip = false)
// ----------------------------------------------------------------------------
{
    if (skip) {
        const auto& nodes = mesh.cellNodes(cell);
        return vector<array<NodeRef, 3>> {array<NodeRef, 3> {RegularTrimesh::coarse_to_fine(nodes[0], 1),
                                                             RegularTrimesh::coarse_to_fine(nodes[1], 1),
                                                             RegularTrimesh::coarse_to_fine(nodes[2], 1)}};
    }

    // tessellate this cell as if its boundary intersects with a refined triangulation
    const array<bool, 3> bnd = identify_boundary(cell, mesh);

    const int sum_bnd = bnd[0] + bnd[1] + bnd[2];
    assert(sum_bnd <= 1); // before calling this function , the grid should have been
                          // 'rounded' so that there are no boundary cells with more than
                          // one boundary edge.

    if (sum_bnd == 1)
        // triangle split in the middle
        return single_split(cell, bnd);

    // if function is called on interior cell, just return the triangle associated
    // with that cell
    const auto& nodes = mesh.cellNodes(cell);
    return vector<array<NodeRef, 3>> {array<NodeRef, 3> {RegularTrimesh::coarse_to_fine(nodes[0]),
                                                         RegularTrimesh::coarse_to_fine(nodes[1]),
                                                         RegularTrimesh::coarse_to_fine(nodes[2])}};
}

// ----------------------------------------------------------------------------
bool
is_boundary_cell(const CellRef& cell, const RegularTrimesh& mesh)
// ----------------------------------------------------------------------------
{
    const auto bnd = identify_boundary(cell, mesh);
    return count(bnd.begin(), bnd.end(), true) > 0;
}

}; // namespace

namespace Opm
{
// ----------------------------------------------------------------------------
RegularTrimesh::RegularTrimesh(const int layers,
                               const array<double, 3>& origin,
                               const array<double, 3>& axis1,
                               const array<double, 3>& axis2,
                               const array<double, 2>& edgelen)
    // ----------------------------------------------------------------------------
    : origin_(origin)
    , axis1_(RegularTrimesh::normalize(axis1))
    , axis2_(RegularTrimesh::normalize(axis2))
    , edgelen_(edgelen)
{
    cellinfo_[{0, 0, 0}] = CellAttributes(); // set a single seed cell
    for (int i = 0; i != layers; ++i)
        expandGrid();
}

// ----------------------------------------------------------------------------
RegularTrimesh::RegularTrimesh(const double radius,
                               const array<double, 3>& origin,
                               const array<double, 3>& axis1,
                               const array<double, 3>& axis2,
                               const array<double, 2>& edgelen)
    // ----------------------------------------------------------------------------
    : origin_(origin)
    , axis1_(RegularTrimesh::normalize(axis1))
    , axis2_(RegularTrimesh::normalize(axis2))
    , edgelen_(edgelen)
{
    const double R2 = radius * radius;
    const double denom2 = (5.0 / 4.0) - (axis1_[0] * axis2_[0] + axis1_[1] * axis2_[1] + axis1_[2] * axis2_[2]);

    // const double denom2 = 2 * (1 - axis1_[0] * axis2_[0] -
    //                                axis1_[1] * axis2_[1] -
    //                                axis1_[2] * axis2_[2]);

    const int c = ceil(max(sqrt(R2 / denom2), radius));

    for (int i = (-c - 1); i <= c; ++i)
        for (int j = (-c - 1); j <= c; ++j)
            for (int k = 0; k != 2; ++k)
                if (dist(origin_, cellCentroid({i, j, k})) < radius)
                    setActive(CellRef {i, j, k});

    auto tmp = remove_mesh_corners(*this);
    swap(tmp);
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::swap(RegularTrimesh& other)
// ----------------------------------------------------------------------------
{
    std::swap(cellinfo_, other.cellinfo_);
    std::swap(origin_, other.origin_);
    std::swap(axis1_, other.axis1_);
    std::swap(axis2_, other.axis2_);
    std::swap(edgelen_, other.edgelen_);
}

// ----------------------------------------------------------------------------
vector<CellRef>
RegularTrimesh::cellIndices() const
// ----------------------------------------------------------------------------
{
    vector<CellRef> indices;
    for (const auto& cell : cellinfo_)
        indices.push_back(cell.first);

    sort(indices.begin(), indices.end());
    return indices;
}

// ----------------------------------------------------------------------------
vector<EdgeRef>
RegularTrimesh::edgeIndices() const
// ----------------------------------------------------------------------------
{
    vector<EdgeRef> all_edges = all_half_edges_();
    // keep only unique elements
    const auto last = unique(all_edges.begin(), all_edges.end());

    // remove duplicates from all_edges
    all_edges.erase(last, all_edges.end());

    return all_edges;
}

// ----------------------------------------------------------------------------
vector<NodeRef>
RegularTrimesh::nodeIndices() const
// ----------------------------------------------------------------------------
{
    vector<NodeRef> indices;

    for (const auto& entry : cellinfo_) {
        const auto& cell = entry.first;
        for (int i = 0; i != 3; ++i)
            indices.push_back(i == 0 ? (cell[2] == 0 ? NodeRef {cell[0], cell[1]} : NodeRef {cell[0] + 1, cell[1] + 1})
                                  : i == 1 ? NodeRef {cell[0] + 1, cell[1]}
                                           : NodeRef {cell[0], cell[1] + 1});
    }

    sort(indices.begin(), indices.end());
    const auto last = unique(indices.begin(), indices.end());

    return vector<NodeRef>(indices.begin(), last);
}

// ----------------------------------------------------------------------------
vector<EdgeRef>
RegularTrimesh::boundaryEdges() const
// ----------------------------------------------------------------------------
{
    // make a vector of all edges.  Count internal edges twice
    const vector<EdgeRef> all_edges = all_half_edges_();

    // boundary edges are those that are not duplicated
    vector<EdgeRef> result;
    for (auto it = all_edges.begin(); it != all_edges.end(); ++it)
        if (*it != *(it + 1))
            result.push_back(*it);
        else
            ++it; // skip the next one, which we already know is duplicated

    // sort entries
    sort(result.begin(), result.end());

    return result;
}

// ----------------------------------------------------------------------------
vector<CellRef>
RegularTrimesh::boundaryCells() const
// ----------------------------------------------------------------------------
{
    vector<CellRef> result;
    for (const auto& edge : boundaryEdges())
        for (const auto& cell : edge2cells(edge))
            if (isActive(cell))
                result.push_back(cell);

    // remove any duplicates
    sort(result.begin(), result.end());
    const auto last = unique(result.begin(), result.end());

    // shrink result to remove duplicates
    result.erase(last, result.end());

    // sort entries
    sort(result.begin(), result.end());

    return result;
}

// ----------------------------------------------------------------------------
vector<CellRef>
RegularTrimesh::interiorCells() const
// ----------------------------------------------------------------------------
{
    const vector<CellRef> bcells = boundaryCells();
    const vector<CellRef> allcells = cellIndices();

    vector<CellRef> result;
    set_difference(allcells.begin(), allcells.end(), bcells.begin(), bcells.end(), back_inserter(result));
    return result;
}

// ----------------------------------------------------------------------------
vector<pair<array<unsigned int, 3>, array<CellRef, 2>>>
RegularTrimesh::boundary_smoothing_triangles_() const
// ----------------------------------------------------------------------------
{
    // identify all 'internal' edges within boundary cells
    const auto bcells = boundaryCells();
    const auto bedges = boundaryEdges(); 
    set<EdgeRef> identified_edges(bedges.begin(), bedges.end()); // better for search?
    array<vector<EdgeRef>, 3> internal_edges;
    for (const auto& cell : bcells) {
        const auto cell_edges = cell2edges(cell);
        for (const auto& e : cell_edges)
            if (identified_edges.find(e) == identified_edges.end()) 
                internal_edges[e[2]].push_back(e);
    }
    
    // identify internal edges that 'line up' along one of the three cardinal grid directions
    const array<int, 3> ioffsets {-1, 2, 1}, joffsets {2, -1, 1};
    array<vector<EdgeRef>, 3> candidate_sites; // candidates for where to place a smoothing triangle

    for (int i = 0; i != 3; ++i)
        for (int k = 0; k != internal_edges[i].size(); ++k)
            if (find(internal_edges[i].begin(), internal_edges[i].end(),
                     EdgeRef { internal_edges[i][k][0] + ioffsets[i],
                               internal_edges[i][k][1] + joffsets[i], i}) != internal_edges[i].end()) 
                candidate_sites[i].push_back(internal_edges[i][k]);
    
    // determine smoothing triangles
    vector<NodeRef> corners;
    vector<CellRef> neigh_cells;
    const array<int, 3> norient {1, 1, 0}; // orientation of "opposing" boundary triangle
    const array<array<int, 2>, 3> ocell {{ {-1, 1}, {1 , -1}, {1, 1 } }}; // location of 'opposing' cell    
    const array<array<int, 3>, 3> n1 {{ {-1, 0}, {0, -1}, {0, 1} }};
    const array<array<int, 3>, 3> n2 {{ {0, 0}, {0, 0}, {1, 0} }};
    const array<array<array<int, 2>, 3>, 3>n1corner {{ {{ {0, 0}, {0, 1}, {-1, 2} }},
                                                       {{ {0, 0}, {2,-1}, {1, 0}  }},
                                                       {{ {0, 1}, {1, 1}, {1, 2}  }} }};
    const array<array<array<int, 2>, 3>, 3>n2corner {{ {{ {1, 0}, {0, 2}, {0, 1} }},
                                                       {{ {0, 1}, {1, 0}, {2, 0} }},
                                                       {{ {1, 0}, {2, 1}, {1, 1} }} }};
    for (int dir = 0; dir != 3; ++dir) 
        for (const auto& edge : candidate_sites[dir]) {
            if (! isActive( {edge[0] + n1[dir][0], edge[1] + n1[dir][1], norient[dir]} )) {
                corners.insert(corners.end(), { {edge[0] + n1corner[dir][0][0], edge[1] + n1corner[dir][0][1]},
                                                {edge[0] + n1corner[dir][1][0], edge[1] + n1corner[dir][1][1]},
                                                {edge[0] + n1corner[dir][2][0], edge[1] + n1corner[dir][2][1]} });
                neigh_cells.insert(neigh_cells.end(), { {edge[0], edge[1], (norient[dir]+1)%2},
                                                        {edge[0] + ocell[dir][0], edge[1] + ocell[dir][1], norient[dir]}});
            }
            if (! isActive( {edge[0] + n2[dir][0], edge[1] + n2[dir][1], norient[dir]} )) {
                corners.insert(corners.end(), { {edge[0] + n2corner[dir][0][0], edge[1] + n2corner[dir][0][1]},
                                                {edge[0] + n2corner[dir][1][0], edge[1] + n2corner[dir][1][1]},
                                                {edge[0] + n2corner[dir][2][0], edge[1] + n2corner[dir][2][1]} });
                neigh_cells.insert(neigh_cells.end(), { {edge[0], edge[1], (norient[dir]+1)%2},
                                                        {edge[0] + ocell[dir][0], edge[1] + ocell[dir][1], norient[dir]}});
            }
        }
    // prepare result by changing NodeRefs to indices
    const vector<unsigned int> node_ixs = noderefs_to_indices_(corners);
    assert(node_ixs.size() == 3 * (neigh_cells.size() / 2)); // each traingle has 3 nodes and 2 cell neighbors
    const int N = node_ixs.size() / 3;
    vector<pair<array<unsigned int, 3>, array<CellRef, 2>>> result;
    for (int i = 0; i != N; ++i) 
        result.push_back({{node_ixs[3 * i],        // triangle corner 1
                           node_ixs[3 * i + 1],    // triangle corner 2
                           node_ixs[3 * i + 2]},   // triangle corner 3
                          {neigh_cells[2 * i],     // cell neighbor 1
                           neigh_cells[2 * i + 1]}}); // cell neighbor 2
    return result;
}

// ----------------------------------------------------------------------------
vector<EdgeRef>
RegularTrimesh::all_half_edges_() const
// ----------------------------------------------------------------------------
{
    // make a vector of all edges.  Count internal edges twice.  Sort the result
    vector<EdgeRef> all_edges;
    for (const auto& entry : cellinfo_)
        for (const auto& edge : half_edges(entry.first))
            all_edges.push_back(edge);

    sort(all_edges.begin(), all_edges.end());
    return all_edges;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::nodeCoord(const NodeRef& node) const
// ----------------------------------------------------------------------------
{
    Coord3D result(origin_);
    result += node[0] * edgelen_[0] * axis1_;
    result += node[1] * edgelen_[1] * axis2_;
    return result;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::cellCentroid(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    Coord3D result = cell[2] == 0 ? nodeCoord({cell[0], cell[1]}) : nodeCoord({cell[0] + 1, cell[1] + 1});
    result += nodeCoord({cell[0] + 1, cell[1]});
    result += nodeCoord({cell[0], cell[1] + 1});

    result /= 3;
    return result;
}

// ----------------------------------------------------------------------------
Coord3D
RegularTrimesh::edgeCentroid(const EdgeRef& edge) const
// ----------------------------------------------------------------------------
{
    return (edge[2] == 0) ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0] + 1, edge[1]})) / 2
        : (edge[2] == 1)  ? (nodeCoord({edge[0], edge[1]}) + nodeCoord({edge[0], edge[1] + 1})) / 2
                          : (nodeCoord({edge[0], edge[1] + 1}) + nodeCoord({edge[0] + 1, edge[1]})) / 2;
}

// ----------------------------------------------------------------------------
vector<Coord3D>
RegularTrimesh::cellCentroids() const
// ----------------------------------------------------------------------------
{
    vector<Coord3D> result;
    for (const auto& entry : cellinfo_)
        result.push_back(cellCentroid(entry.first));
    return result;
}

// ----------------------------------------------------------------------------
vector<Coord3D>
RegularTrimesh::edgeCentroids() const
// ----------------------------------------------------------------------------
{
    vector<Coord3D> result;
    for (const auto& edge : edgeIndices())
        result.push_back(edgeCentroid(edge));
    return result;
}

// ----------------------------------------------------------------------------
vector<Coord3D>
RegularTrimesh::nodeCoords() const
// ----------------------------------------------------------------------------
{
    vector<Coord3D> result;
    for (const auto& node : nodeIndices())
        result.push_back(nodeCoord(node));
    return result;
}

//------------------------------------------------------------------------------
pair<vector<array<unsigned int, 3>>, vector<int>>
RegularTrimesh::getTriangles() const
//------------------------------------------------------------------------------
{
    vector<array<unsigned int, 3>> result;
    for (const auto& cell : cellNodesLinear())
        result.push_back({(unsigned int)cell[0], (unsigned int)cell[1], (unsigned int)cell[2]});

    vector<int> trivialmap(result.size());
    for (size_t i = 0; i != result.size(); ++i)
        trivialmap[i] = i;
    return make_pair(result, trivialmap);
}

//------------------------------------------------------------------------------
pair<vector<array<unsigned int, 3>>, vector<int>>
RegularTrimesh::getMultiresTriangles(const vector<CellRef>& fixed_cells, const int max_levels) const
//------------------------------------------------------------------------------
{
    vector<NodeRef> triangles; // 3 consecutive entries make up a triangle
    vector<int> cellmap;
    RegularTrimesh mesh(*this);
    const int CELLNUM_THRESHOLD = 4;
    int level = 0;

    while (mesh.numCells() > CELLNUM_THRESHOLD && level < max_levels) {

        // create new coarsened level
        auto coarsened = coarsen_mesh(mesh, level == 0 ? fixed_cells : vector<CellRef>());
        auto& new_mesh = coarsened.first;
        auto& uncoarsened_cells = coarsened.second;

        // add uncoarsened cells into result
        for (const auto& cell : uncoarsened_cells)
            for (const auto& tri : tesselate_coarsecell(cell, mesh, level == 0))
                for (int i = 0; i != 3; ++i)
                    triangles.push_back(coarse_to_fine(tri[i], level - 1));

        // add corresponding entries to cellmap
        if (level == 0)
            for (const auto& cell : uncoarsened_cells)
                cellmap.push_back(linearCellIndex(cell));
        else
            cellmap.insert(cellmap.end(), uncoarsened_cells.size(), -1);

        // increment level and swap meshes
        level++;
        mesh.swap(new_mesh);
    }

    // add remaining cells into result
    const int num_tri_before = triangles.size() / 3;
    for (const auto& cell : mesh.cellIndices())
        for (const auto& tri : tesselate_coarsecell(cell, mesh, level == 0))
            for (int i = 0; i != 3; ++i)
                triangles.push_back(coarse_to_fine(tri[i], level - 1));

    const int new_tri = (triangles.size() / 3) - num_tri_before;
    cellmap.insert(cellmap.end(), new_tri, -1);

    // return result after having translated NodeRefs to indices
    const auto triangles_ixs = noderefs_to_indices_(triangles);
    vector<array<unsigned int, 3>> result_triangles;
    for (size_t i = 0; i != triangles_ixs.size(); i += 3)
        result_triangles.push_back(
            {(unsigned int)triangles_ixs[i], (unsigned int)triangles_ixs[i + 1], (unsigned int)triangles_ixs[i + 2]});

    return make_pair(result_triangles, cellmap);
}

//------------------------------------------------------------------------------
pair<RegularTrimesh, vector<CellRef>>
RegularTrimesh::coarsen_mesh(const RegularTrimesh& mesh, const vector<CellRef>& fixed_cells)
//------------------------------------------------------------------------------
{
    RegularTrimesh tmp_mesh(mesh);
    for (const auto& cell : fixed_cells)
        tmp_mesh.setInactive(cell);

    tmp_mesh.contractGrid(); // exclude boundary cells from coarsening

    const RegularTrimesh coarsened = remove_mesh_corners(tmp_mesh.coarsen(true));

    // identify which fine-scale cells have been covered by the coarse grid
    vector<CellRef> covered_finecells; // will be gradually filled below
    for (const auto& cell : coarsened.cellIndices())
        for (const auto& cref : coarse_to_fine(cell))
            covered_finecells.push_back(cref);

    vector<CellRef> all_finecells = mesh.cellIndices(); // these are sorted
    vector<CellRef> uncovered_finecells;
    sort(covered_finecells.begin(), covered_finecells.end()); // required by set_difference
    set_difference(all_finecells.begin(),
                        all_finecells.end(),
                        covered_finecells.begin(),
                        covered_finecells.end(),
                        back_inserter(uncovered_finecells));

    return make_pair(coarsened, uncovered_finecells);
}

//------------------------------------------------------------------------------
vector<unsigned int>
RegularTrimesh::noderefs_to_indices_(const vector<NodeRef>& noderefs) const
//------------------------------------------------------------------------------
{
    map<NodeRef, unsigned int> nodemap;
    vector<NodeRef> nix = nodeIndices();
    for (size_t i = 0; i != nix.size(); ++i)
        nodemap[nix[i]] = i;

    vector<unsigned int> result;
    for (const auto& node : noderefs)
        result.push_back(nodemap[node]);

    return result;
}

//------------------------------------------------------------------------------
vector<unsigned int>
RegularTrimesh::cellrefs_to_indices_(const vector<CellRef>& cellrefs) const
//------------------------------------------------------------------------------
{
    map<CellRef, unsigned int> cellmap;
    vector<CellRef> cix = cellIndices();
    for (size_t i = 0; i != cix.size(); ++i)
        cellmap[cix[i]] = i;

    vector<unsigned int> result;
    for (const auto& cell : cellrefs)
        result.push_back(cellmap[cell]);

    return result;
}

//------------------------------------------------------------------------------
std::tuple<std::unique_ptr<Grid>, std::vector<int>, std::map<int, int>>    
RegularTrimesh::createDuneGrid(const int coarsen_levels,
                               const vector<CellRef>& fixed_cells,
                               const bool add_smoothing_triangles) const
//------------------------------------------------------------------------------
{
    Dune::GridFactory<Grid> factory;

    // define points
    for (const auto& node : nodeCoords())
        factory.insertVertex(Dune::FieldVector<double, 3> {node[0], node[1], node[2]});

    // define triangles
    const auto tmp = coarsen_levels > 0 ? getMultiresTriangles(fixed_cells, coarsen_levels) : getTriangles();
    const auto& fsmap = tmp.second;
    
    for (const auto& tri : tmp.first)
        factory.insertElement(Dune::GeometryTypes::simplex(2),
                              vector<unsigned int> {tri[0], tri[1], tri[2]});

    // create map from trimesh boundary cells to Dune grid cells
    std::map<int, int> boundary_map;
    std::vector<size_t> fsmap_inv(numCells(), -1);
    for (int i = 0; i != fsmap.size(); ++i)
        if (fsmap[i] != -1)
            fsmap_inv[fsmap[i]] = i; // mapping from Trimesh to Dune grid cells
    const auto bcells = cellrefs_to_indices_(boundaryCells());
    for (const auto& cell : bcells)
        boundary_map[cell] = fsmap_inv[cell]; // mapping from Dune grid cells to Trimesh cells
        
    
    if (add_smoothing_triangles) {
        const auto smoothing_triangles = boundary_smoothing_triangles_();
        int smoothing_cell_index = fsmap.size();
        for (const auto& tri : smoothing_triangles) {
            factory.insertElement(Dune::GeometryTypes::simplex(2),
                                  vector<unsigned int> {tri.first[0], tri.first[1], tri.first[2]});
            // make the associated boundary cells in the TriMesh map to the
            // smoothing triangle in the Dune grid
            boundary_map[linearCellIndex(tri.second[0])] = smoothing_cell_index;
            boundary_map[linearCellIndex(tri.second[1])] = smoothing_cell_index;
            ++smoothing_cell_index;
        }
    }
    
    return make_tuple(factory.createGrid(), fsmap, boundary_map);
}

// ----------------------------------------------------------------------------
pair<NodeRef, NodeRef>
RegularTrimesh::edgeNodes(const EdgeRef& e) const
// ----------------------------------------------------------------------------
{
    return e[2] == 0 ? make_pair(NodeRef {e[0], e[1]}, NodeRef {e[0] + 1, e[1]})
        : e[2] == 1  ? make_pair(NodeRef {e[0], e[1]}, NodeRef {e[0], e[1] + 1})
                     : make_pair(NodeRef {e[0], e[1] + 1}, NodeRef {e[0] + 1, e[1]});
}

// ----------------------------------------------------------------------------
vector<pair<size_t, size_t>>
RegularTrimesh::edgeNodeIndices(bool only_boundary) const
// ----------------------------------------------------------------------------
{
    const auto nodeindices = nodeIndices();
    const auto edgeindices = only_boundary ? boundaryEdges() : edgeIndices();

    // make mapping from node to index
    map<NodeRef, size_t> node2index;
    for (size_t i = 0; i != nodeindices.size(); ++i)
        node2index[nodeindices[i]] = i;

    // map back to indices, for all edges in the mesh
    vector<pair<size_t, size_t>> result;
    for (const auto& edge : edgeindices) {
        const auto nodes = edgeNodes(edge);
        result.push_back({node2index[nodes.first], node2index[nodes.second]});
    }

    return result;
}

// ----------------------------------------------------------------------------
pair<Coord3D, Coord3D>
RegularTrimesh::edgeNodeCoords(const EdgeRef& edge) const
// ----------------------------------------------------------------------------
{
    return make_pair(nodeCoord(edgeNodes(edge).first), nodeCoord(edgeNodes(edge).second));
}

// ----------------------------------------------------------------------------
vector<pair<Coord3D, Coord3D>>
RegularTrimesh::edgeNodeCoords() const
// ----------------------------------------------------------------------------
{
    const auto edgeindices = edgeIndices();
    vector<pair<Coord3D, Coord3D>> result;
    for (const auto& edge : edgeindices)
        result.push_back(edgeNodeCoords(edge));
    return result;
}

// ----------------------------------------------------------------------------
array<NodeRef, 3>
RegularTrimesh::cellNodes(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    return cell[2] == 0
        ? array<NodeRef, 3> {NodeRef {cell[0], cell[1]}, NodeRef {cell[0] + 1, cell[1]}, NodeRef {cell[0], cell[1] + 1}}
        : array<NodeRef, 3> {
            NodeRef {cell[0] + 1, cell[1]}, NodeRef {cell[0], cell[1] + 1}, NodeRef {cell[0] + 1, cell[1] + 1}};
}

// ----------------------------------------------------------------------------
vector<array<size_t, 3>>
RegularTrimesh::cellNodesLinear() const
// ----------------------------------------------------------------------------
{
    vector<array<size_t, 3>> result;
    const auto nodeindices = nodeIndices();
    const auto cellindices = cellIndices();

    const auto findnode = [&nodeindices](const NodeRef& node) {
        return size_t(find(nodeindices.begin(), nodeindices.end(), node) - nodeindices.begin());
    };

    for (const auto& cell : cellindices) {
        const auto noderefs = cellNodes(cell);
        result.push_back({findnode(noderefs[0]), findnode(noderefs[1]), findnode(noderefs[2])});
    }
    return result;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::writeMatlabTriangulation(ostream& out) const
// ----------------------------------------------------------------------------
{
    const auto nodeindices = nodeIndices();
    const auto nodecoords = nodeCoords();

    // define the vertices of the triangulation
    out << "vertices = [";
    for (const auto& node : nodecoords) {
        out << node[0] << " " << node[1] << " " << node[2] << "; ";
    }
    out << "];\n";

    // define the triangles of the triangulation
    const auto cellnodes = cellNodesLinear();
    out << "triangles = [";
    for (const auto& cell : cellnodes) {
        out << cell[0] + 1 << " " << cell[1] + 1 << " " << cell[2] + 1 << "; ";
    }
    out << "];\n";

    // create a triangulation object
    out << "tri = triangulation(triangles, vertices);\n";

    // plot the triangulation, using triplot
    out << "figure;\n";
    out << "triplot(tri);\n";

    // set axis equal
    out << "axis equal;\n";
}

// ----------------------------------------------------------------------------
vector<CellRef>
RegularTrimesh::inner_ring_cells() // static function
// ----------------------------------------------------------------------------
{
    // utility function to quickly construct a vector with reference to all
    // cells having the origin as a corner
    return vector<CellRef> {{0, 0, 0}, {-1, 0, 1}, {-1, 0, 0}, {-1, -1, 1}, {0, -1, 0}, {0, -1, 1}};
}

// ----------------------------------------------------------------------------
void
writeMeshToVTK(const RegularTrimesh& mesh,
               const char* const filename,
               const int coarsen_levels,
               const vector<CellRef>& fixed_cells,
               const bool add_smoothing_triangles)
// ----------------------------------------------------------------------------
{
    const auto grid = mesh.createDuneGrid(coarsen_levels, fixed_cells, add_smoothing_triangles);

    // write grid to file
    auto vtkwriter
        = make_unique<Dune::VTKWriter<Grid::LeafGridView>>(get<0>(grid)->leafGridView(), Dune::VTK::nonconforming);
    // write flag to file
    if (coarsen_levels == 0) {
        vector<int> flags = mesh.getCellFlags();
        if (add_smoothing_triangles)
            flags.insert(flags.end(), vector<int>(grid->
        vtkwriter->addCellData(flags, "flag");

            
        vtkwriter->write(filename);
    } else {
        vtkwriter->write(filename);
    }
}

// ----------------------------------------------------------------------------
void
writeMeshToVTKDebug(const RegularTrimesh& mesh, const char* const filename, const int coarsen_levels)
// ----------------------------------------------------------------------------
{
    writeMeshToVTK(mesh, filename, coarsen_levels);
}

// ----------------------------------------------------------------------------
void
writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename)
// ----------------------------------------------------------------------------
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    const vector<Coord3D> points = mesh.nodeCoords();
    const vector<pair<size_t, size_t>> edges = mesh.edgeNodeIndices(true);

    // header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << edges.size() << "\">\n";
    // points
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& point : points)
        file << "          " << point[0] << " " << point[1] << " " << point[2] << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // cells (edges)
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& edge : edges)
        file << "          " << edge.first << " " << edge.second << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= edges.size(); ++i)
        file << "          " << i * 2 << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < edges.size(); ++i)
        file << "          3\n";
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Footer
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}

// ----------------------------------------------------------------------------
size_t
RegularTrimesh::numActive() const
// ----------------------------------------------------------------------------
{
    return cellinfo_.size();
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::isActive(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    return cellinfo_.find(cell) != cellinfo_.end();
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::setInactive(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    if (!isActive(cell))
        return false;
    // remove this cell from cellinfo
    cellinfo_.erase(cell);
    return true;
}

// ----------------------------------------------------------------------------
bool
RegularTrimesh::setActive(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    if (isActive(cell))
        return false;
    cellinfo_.insert({cell, CellAttributes()});

    assert(cellinfo_[cell].flag == 0);
    return true;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    int result = 0;
    const auto edges = cell2edges(cell);
    for (const auto& e : edges)
        for (const auto& c : edge2cells(e))
            result += setActive(c);

    return result;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::contractGrid()
// ----------------------------------------------------------------------------
{
    // identify boundary cells and remove them
    const vector<CellRef> bcells = boundaryCells();

    // identify all cells that are not boundary cells, using set_difference
    for (const auto& bc : bcells)
        cellinfo_.erase(bc);

    return (int)bcells.size();
}


// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid(const vector<CellRef>& cells)
// ----------------------------------------------------------------------------
{
    int result = 0;
    for (const auto& c : cells)
        result += expandGrid(c);
    return result;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::expandGrid()
// ----------------------------------------------------------------------------
{
    // regular expansion in all directions
    return expandGrid(boundaryCells());
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::removeSawtooths()
// ----------------------------------------------------------------------------
{
    const auto bedges = boundaryEdges();
    vector<CellRef> candidates;
    for (const auto& edge : bedges)
        for (const auto& cell : edge2cells(edge))
            if (!isActive(cell))
                candidates.push_back(cell);

    sort(candidates.begin(), candidates.end());

    // inactive cells adjacent to more than one boundary edge should be activated
    for (auto it = candidates.begin(); it != candidates.end(); ++it) {
        const auto range = equal_range(it, candidates.end(), *it);
        if (range.second - range.first > 1) {
            setActive(*it);
            it = range.second - 1;
        }
    }
}

// ----------------------------------------------------------------------------
size_t
RegularTrimesh::linearCellIndex(const CellRef& cell) const
{
    // cellinfo is a map from CellRef to CellAttributes.
    assert(cellinfo_.find(cell) != cellinfo_.end());
    return distance(cellinfo_.begin(), cellinfo_.find(cell));
}

// ----------------------------------------------------------------------------
CellRef
RegularTrimesh::cellIndex(const size_t index) const
{
    auto it = cellinfo_.begin();
    advance(it, index);
    return it->first;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setCellFlag(const CellRef& cell, const int value)
// ----------------------------------------------------------------------------
{
    cellinfo_[cell].flag = value;
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setCellFlags(const vector<CellRef>& cells, const int value)
// ----------------------------------------------------------------------------
{
    for (const auto& cell : cells)
        setCellFlag(cell, value);
}

// ----------------------------------------------------------------------------
void
RegularTrimesh::setAllFlags(const int value)
// ----------------------------------------------------------------------------
{
    for (auto& it : cellinfo_)
        it.second.flag = value;
}

// ----------------------------------------------------------------------------
int
RegularTrimesh::getCellFlag(const CellRef& cell) const
// ----------------------------------------------------------------------------
{
    const auto it = cellinfo_.find(cell);
    return it->second.flag;
}

// ----------------------------------------------------------------------------
vector<int>
RegularTrimesh::getCellFlags() const
// ----------------------------------------------------------------------------
{
    vector<int> result;
    for (const auto& el : cellinfo_)
        result.push_back(el.second.flag);
    return result;
}

// ----------------------------------------------------------------------------
RegularTrimesh
RegularTrimesh::refine() const
// ----------------------------------------------------------------------------
{
    map<CellRef, CellAttributes> new_cells;
    for (const auto& e : cellinfo_)
        for (const auto& c : coarse_to_fine(e.first))
            new_cells[c] = e.second;

    return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0] / 2, edgelen_[1] / 2}};
}

// ----------------------------------------------------------------------------
RegularTrimesh
RegularTrimesh::coarsen(bool strict) const
// ----------------------------------------------------------------------------
{
    map<CellRef, CellAttributes> new_cells;
    for (const auto& e : cellinfo_) {
        const CellRef& c = e.first;
        if (!strict || !is_boundary_cell(c, *this))
            if ((abs(c[0]) % 2 == 1 && abs(c[1]) % 2 == 1 && c[2] == 0)
                || (abs(c[0]) % 2 == 0 && abs(c[1]) % 2 == 0 && c[2] == 1))
                new_cells[fine_to_coarse(c)] = e.second;
    }

    return RegularTrimesh {new_cells, origin_, axis1_, axis2_, {edgelen_[0] * 2, edgelen_[1] * 2}};
}


// ----------------------------------------------------------------------------
array<CellRef, 4>
RegularTrimesh::coarse_to_fine(const CellRef& c)
// ----------------------------------------------------------------------------
{
    return (c[2] == 0) ? array<CellRef, 4> {CellRef {2 * c[0], 2 * c[1], 0},
                                            CellRef {2 * c[0] + 1, 2 * c[1], 0},
                                            CellRef {2 * c[0], 2 * c[1] + 1, 0},
                                            CellRef {2 * c[0], 2 * c[1], 1}}
                       : array<CellRef, 4> {CellRef {2 * c[0] + 1, 2 * c[1] + 1, 0},
                                            CellRef {2 * c[0] + 1, 2 * c[1], 1},
                                            CellRef {2 * c[0], 2 * c[1] + 1, 1},
                                            CellRef {2 * c[0] + 1, 2 * c[1] + 1, 1}};
}

// ----------------------------------------------------------------------------
CellRef
RegularTrimesh::fine_to_coarse(const CellRef& cell)
// ----------------------------------------------------------------------------
{
    const int i = cell[0] >= 0 ? cell[0] : cell[0] - 1;
    const int j = cell[1] >= 0 ? cell[1] : cell[1] - 1;
    return {i / 2, j / 2, cell[2] == 0 ? 1 : 0};
}

// ----------------------------------------------------------------------------
NodeRef
RegularTrimesh::coarse_to_fine(const NodeRef& node, const int level)
// ----------------------------------------------------------------------------
{
    if (level < 0)
        return {node[0] / (1 << abs(level)), node[1] / (1 << abs(level))};
    else
        return {node[0] * (1 << level), node[1] * (1 << level)};
    // return {2*node[0], 2*node[1]};
}

// ----------------------------------------------------------------------------
vector<CellRef>
RegularTrimesh::interior_coarsegrid_() const
// ----------------------------------------------------------------------------
{
    map<CellRef, int> cell_count;
    for (const auto& e : cellinfo_) {
        const CellRef coarse_cell = fine_to_coarse(e.first);
        if (cell_count.find(coarse_cell) == cell_count.end())
            cell_count[coarse_cell] = 1;
        else
            cell_count[coarse_cell]++;
    }
    vector<CellRef> result;
    for (const auto& e : cell_count)
        if (e.second == 4) // this cell is fully covered by all its fine cells
            result.push_back(e.first);

    return result;
}

// ----------------------------------------------------------------------------
RegularTrimesh
expand_to_criterion(const RegularTrimesh& mesh,
                    function<vector<double>(const RegularTrimesh&)> score_function,
                    double threshold)
// ----------------------------------------------------------------------------
{
    // @@ for the moment, an extremely simple expansion algorithm
    RegularTrimesh working_mesh = mesh; // make a working copy of the mesh;

    while (true) { // keep looping as long as grid still needs expansion
        const vector<double> bnd_scores = score_function(working_mesh);
        const vector<CellRef> bnd_cells = mesh.boundaryCells();

        assert(bnd_scores.size() == bnd_cells.size());

        vector<CellRef> expand_cells;
        for (size_t i = 0; i != bnd_scores.size(); ++i)
            if (bnd_scores[i] > threshold)
                expand_cells.push_back(bnd_cells[i]);

        if (expand_cells.size() == 0)
            break;

        // @@ debugging info
        // working_mesh.setAllFlags(0);
        // working_mesh.setCellFlags(expand_cells, 1);
        // writeMeshToVTK(working_mesh, "before", false);

        working_mesh.expandGrid(expand_cells);
        working_mesh.removeSawtooths();

        // writeMeshToVTK(working_mesh, "after", false);
    }
    return working_mesh;
}

} // namespace Opm
