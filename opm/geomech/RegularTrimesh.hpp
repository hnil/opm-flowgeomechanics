#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Opm
{

typedef std::array<int, 3> CellRef; // (i, j, [0 | 1])
typedef std::array<int, 3> EdgeRef; // (i, j, [0 | 1 | 2])
typedef std::array<int, 2> NodeRef; // (i, j)
typedef std::array<double, 3> Coord3D;
using Grid = Dune::FoamGrid<2, 3>;
// typedef Dune::FoamGrid<2, 3> Grid;
// typedef int Grid; // @@ dummy


struct CellAttributes {
    CellAttributes()
        : flag(0) {};
    int flag;
}; // fill this with contents as needed



class RegularTrimesh
{
public:
    template <typename Iterator>
    RegularTrimesh(const Iterator cells_begin,
                   const Iterator cells_end,
                   const std::array<double, 3>& origin = {0, 0, 0},
                   const std::array<double, 3>& axis1 = {1, 0, 0},
                   const std::array<double, 3>& axis2 = {0.5, std::sqrt(3) / 2, 0},
                   const std::array<double, 2>& edgelen = {1, 1})
        : origin_(origin)
        , axis1_(RegularTrimesh::normalize(axis1))
        , axis2_(RegularTrimesh::normalize(axis2))
        , edgelen_(edgelen)
    {
        for (auto it = cells_begin; it != cells_end; ++it)
            cellinfo_[*it] = CellAttributes();
    }

    RegularTrimesh(const std::map<CellRef, CellAttributes>& cells,
                   const std::array<double, 3>& origin = {0, 0, 0},
                   const std::array<double, 3>& axis1 = {1, 0, 0},
                   const std::array<double, 3>& axis2 = {0.5, std::sqrt(3) / 2, 0},
                   const std::array<double, 2>& edgelen = {1, 1})
        : cellinfo_(cells)
        , origin_(origin)
        , axis1_(RegularTrimesh::normalize(axis1))
        , axis2_(RegularTrimesh::normalize(axis2))
        , edgelen_(edgelen)
    {
    }

    RegularTrimesh(const int layers = 1,
                   const std::array<double, 3>& origin = {0, 0, 0},
                   const std::array<double, 3>& axis1 = {1, 0, 0},
                   const std::array<double, 3>& axis2 = {0.5, std::sqrt(3) / 2, 0},
                   const std::array<double, 2>& edgelen = {1, 1});

    RegularTrimesh(const double radius,
                   const std::array<double, 3>& origin = {0, 0, 0},
                   const std::array<double, 3>& axis1 = {1, 0, 0},
                   const std::array<double, 3>& axis2 = {0.5, std::sqrt(3) / 2, 0},
                   const std::array<double, 2>& edgelen = {1, 1});
    
    // --------------------- Functions for inspecting the grid ---------------------
    size_t numCells() const
    {
        return cellinfo_.size();
    }
    std::vector<CellRef> cellIndices() const; // result is sorted
    std::vector<EdgeRef> edgeIndices() const; // result is sorted
    std::vector<NodeRef> nodeIndices() const; // result is sorted

    std::vector<EdgeRef> boundaryEdges() const; // result is sorted
    std::vector<CellRef> boundaryCells() const; // result is sorted
    std::vector<CellRef> interiorCells() const; // result is sorted

    Coord3D cellCentroid(const CellRef& cell) const;
    Coord3D edgeCentroid(const EdgeRef& edge) const;
    Coord3D nodeCoord(const NodeRef& node) const;

    std::vector<Coord3D> cellCentroids() const;
    std::vector<Coord3D> edgeCentroids() const;
    std::vector<Coord3D> nodeCoords() const;

    static std::array<NodeRef, 3> cellNodes(const CellRef& cell);
    std::vector<std::array<size_t, 3>> cellNodesLinear() const;

    std::pair<NodeRef, NodeRef> edgeNodes(const EdgeRef& edge) const;
    std::vector<std::pair<size_t, size_t>> edgeNodeIndices(bool only_boundary = false) const;
    std::pair<Coord3D, Coord3D> edgeNodeCoords(const EdgeRef& edge) const;
    std::vector<std::pair<Coord3D, Coord3D>> edgeNodeCoords() const;

    bool isActive(const CellRef& cell) const;
    size_t numActive() const;
    size_t linearCellIndex(const CellRef& cell) const;
    CellRef cellIndex(const size_t idx) const;
    int getCellFlag(const CellRef& cell) const;
    std::vector<int> getCellFlags() const;
    std::vector<CellRef> activeNeighborCells(const std::vector<CellRef>& cells) const;

    // --------------------- Functions for modifying the grid ---------------------
    void setAllFlags(const int value);
    void setCellFlag(const CellRef& cell, const int value);
    void setCellFlags(const std::vector<CellRef>& cells, const int value);
    bool setActive(const CellRef& cell);
    bool setInactive(const CellRef& cell);
    int expandGrid(const CellRef& cell);
    int expandGrid(const std::vector<CellRef>& cells);
    int expandGrid(); // uniform expansion all directions
    int contractGrid(); // uniform contraction in all directions
    void removeSawtooths();
    void swap(RegularTrimesh& other);
    // ---------------------- Functions for outputting other grid types -------------

    // first return value is the Dune grid, second is a map from Dune cells to
    // Trimesh cells, third is a map from TriMesh boundary cells to the relevant
    // boundary cell in the Dune grid (not always trivial, due to smoothing
    // cells)
    std::tuple<std::unique_ptr<Grid>, std::vector<std::vector<CellRef>>, std::map<CellRef, int>>
    createDuneGrid(const int coarsen_levels = 0,
                   const std::vector<CellRef>& fixed_cells = std::vector<CellRef>(),
                   const bool add_smoothing_triangles = true,
                   const int cellnum_threshold = 10) const;
    void writeMatlabTriangulation(std::ostream& out) const;

    // ------------- Functions for creating new RegularTrimesh objects -------------
    RegularTrimesh refine() const;
    RegularTrimesh coarsen(bool strict = false) const;

    // -------------- Functions for getting the triangles of the mesh --------------

    // Returns a pair of vectors, first is the triangles, second is the set of
    // triangles in the Trimesh covered by each triangle (for
    // `getTriangles`, this is a trivial one-to-one mapping, but not for its
    // sister function `getMultiresTriangles`).
    std::pair<std::vector<std::array<NodeRef, 3>>,std::vector<std::vector<CellRef>>>
    getTriangles() const;

    // Returns a pair of vectors, first is the triangles, second is the set of
    // triangles in the Trimesh covered by each triangle (note that this is a
    // many-to-many relationship, since when a coarse triangle is split in two
    // (for conformity reasons), each of the two resulting triangles will be
    // considered "covering" all the Trimesh triangles covered by the coarse
    // triangle before splitting.
    std::pair<std::vector<std::array<NodeRef, 3>>, std::vector<std::vector<CellRef>>>
    getMultiresTriangles(const std::vector<CellRef>& fixed_cells = std::vector<CellRef>(),
                        const int max_levels = 5, const int cellnum_threshold = 10) const;
//                        const int max_levels, const int cellnum_threshold ) const;
//                         

    // static functions
    static std::vector<CellRef> coarse_to_fine(const CellRef& cell, const int levels=1);
    static NodeRef coarse_to_fine(const NodeRef& node, const int levels=1);
    static CellRef fine_to_coarse(const CellRef& cell, const int levels=1);
    static std::vector<CellRef> inner_ring_cells();
    static std::pair<RegularTrimesh, std::vector<CellRef>>
    coarsen_mesh(const RegularTrimesh& mesh, const std::vector<CellRef>& fixed_cells);
                                                                             

private:
    // helper functions
    std::vector<EdgeRef> all_half_edges_() const; // internal edges are duplicated
    std::vector<CellRef> interior_coarsegrid_() const; // all coarse cells fully covered by fine ones
    std::vector<unsigned int> noderefs_to_indices_(const std::vector<NodeRef>& noderefs) const;
    std::vector<unsigned int> cellrefs_to_indices_(const std::vector<CellRef>& cellrefs) const;
    std::vector<std::array<unsigned int, 3>> ref2index_(const std::vector<std::array<NodeRef, 3>>& vec) const;
    
    std::vector<std::pair<std::array<unsigned int, 3>,
                          std::array<CellRef, 2>>> boundary_smoothing_triangles_() const;

    static double norm(const Coord3D& v)
    {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    static double dist(const Coord3D& v1, const Coord3D& v2) {
        return std::sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) +
                         (v1[1] - v2[1]) * (v1[1] - v2[1]) +
                         (v1[2] - v2[2]) * (v1[2] - v2[2]));
    }
        
    static Coord3D normalize(const Coord3D& v)
    {
        const double n = norm(v); //std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        return {v[0] / n, v[1] / n, v[2] / n};
    }

    static void rotate60(NodeRef& node) // rotate node 60 degrees clockwise
    {
        NodeRef tmp {node[0] + node[1], -node[0]};
        std::swap(node, tmp);
    }

    static void rotate60(CellRef& cell) // rotate cell 60 degrees clockwise
    {
        CellRef tmp {cell[0] + cell[1] + cell[2], -cell[0]  - 1, (cell[2] + 1) % 2} ;
        std::swap(cell, tmp);
    }

    static NodeRef edge2node(const EdgeRef& edge) // convert edge to corner node
    {
        return { edge[0], edge[1] }; 
    }

    static CellRef node2cell(const NodeRef& node, bool second) // make cell from corner node and orientation
    {
        return { node[0], node[1], int(second) };
    }
    
    // data members
    std::map<CellRef, CellAttributes> cellinfo_;
    Coord3D origin_;
    Coord3D axis1_;
    Coord3D axis2_;
    std::array<double, 2> edgelen_;
};

void writeMeshToVTK(const RegularTrimesh& mesh, const char* const filename, const int coarsen_levels = 0,
                    const std::vector<CellRef>& fixed_cells = std::vector<CellRef>(), const bool add_smoothing_triangles = true);
// @@ just because debugger doesnt handle the default arg.    
void writeMeshToVTKDebug(const RegularTrimesh& mesh, const char* const filename, const int coarsen_levels = 0, const bool add_smoothing_triangles = true);

void writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename);


std::tuple<RegularTrimesh, int>
expand_to_criterion(const RegularTrimesh& mesh,
        std::function<std::vector<double>(const RegularTrimesh&,
                                          const int level)> score_function,
        double threshold, const std::vector<CellRef>& fixed_cells,
        const int target_cellcount, // target number of cells in final mesh
        const int cellcount_threshold // target number of cells in initial mesh to expand (start level will be determined by this)
        );
} // namespace Opm
