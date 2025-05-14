#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

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
    std::pair<std::unique_ptr<Grid>, std::vector<int>>
    createDuneGrid(const int coarsen_levels = 0,
                   const std::vector<CellRef>& fixed_cells = std::vector<CellRef>()) const;
    void writeMatlabTriangulation(std::ostream& out) const;

    // ------------- Functions for creating new RegularTrimesh objects -------------
    RegularTrimesh refine() const;
    RegularTrimesh coarsen(bool strict = false) const;

    // -------------- Functions for getting the triangles of the mesh --------------
    std::pair<std::vector<std::array<unsigned int, 3>>, std::vector<int>> getTriangles() const;
    std::pair<std::vector<std::array<unsigned int, 3>>, std::vector<int>>
    getMultiresTriangles(const std::vector<CellRef>& fixed_cells = std::vector<CellRef>(),
                         const int max_levels = 5) const;

    // static functions
    static std::array<CellRef, 4> coarse_to_fine(const CellRef& cell);
    static NodeRef coarse_to_fine(const NodeRef& node, const int level=1);
    static CellRef fine_to_coarse(const CellRef& cell);
    static std::vector<CellRef> inner_ring_cells();
    static std::pair<RegularTrimesh, std::vector<CellRef>>
    coarsen_mesh(const RegularTrimesh& mesh, const std::vector<CellRef>& fixed_cells);
                                                                             

private:
    // helper functions
    std::vector<EdgeRef> all_half_edges_() const; // internal edges are duplicated
    std::vector<CellRef> interior_coarsegrid_() const; // all coarse cells fully covered by fine ones
    std::vector<unsigned int> noderefs_to_indices(const std::vector<NodeRef>& noderefs) const;
    std::vector<unsigned int> cellrefs_to_indices(const std::vector<CellRef>& cellrefs) const;
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

    // data members
    std::map<CellRef, CellAttributes> cellinfo_;
    Coord3D origin_;
    Coord3D axis1_;
    Coord3D axis2_;
    std::array<double, 2> edgelen_;
};

void writeMeshToVTK(const RegularTrimesh& mesh, const char* const filename, const int coarsen_levels = 0,
                    const std::vector<CellRef>& fixed_cells = std::vector<CellRef>());
// @@ just because debugger doesnt handle the default arg.    
void writeMeshToVTKDebug(const RegularTrimesh& mesh, const char* const filename, const int coarsen_levels = 0);

void writeMeshBoundaryToVTK(const RegularTrimesh& mesh, const char* const filename);


RegularTrimesh expand_to_criterion(const RegularTrimesh& mesh,
                                   std::function<std::vector<double>(const RegularTrimesh&)> score_function,
                                   double threshold);
} // namespace Opm
