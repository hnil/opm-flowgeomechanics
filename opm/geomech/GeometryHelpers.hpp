#pragma once
#include <array>
#include <vector>
#include <unordered_map>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <opm/grid/CpGrid.hpp>

// for resinsight search
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>

#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractor.hpp>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>
#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>
#include <external/resinsight/CommonCode/cvfStructGrid.h>
#include <external/resinsight/LibGeometry/cvfBoundingBox.h>

namespace Opm{
    // copy from base vangaurd
    class GeometryHelper{
        //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        // static const int dimension = Grid::dimension;
        // static const int dimensionworld = Grid::dimensionworld;
        // using Element = typename GridView::template Codim<0>::Entity;
        // using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
        //using type = Dune::MultipleCodimMultipleGeomTypeMapper<GetPropType<TypeTag, Properties::GridView>>;
    public:
        using Point3D = Dune::FieldVector<double,3>;
        template<class Grid>
        GeometryHelper(const Grid& grid){
            //cartDim_ = cartdim;
            this->init(grid);
        }
        // static constexpr int Dimension = 3;
        // using Point3D = Dune::FieldVector<double,3>;
        // unsigned cartesianIndex(const std::array<int, Dimension>& coords) const
        // {
        //     unsigned cartIndex = coords[0];
        //     int factor = cartDim_[0];
        //     for (unsigned i = 1; i < cartDim_.size(); ++i) {
        //         cartIndex += coords[i]*factor;
        //         factor *= cartDim_[i];
        //     }

        //     return cartIndex;
        // }
        int compressedIndex(int cartesianCellIdx) const
        {
            //using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
            auto index_pair = cartesianToCompressed_.find(cartesianCellIdx);
            if (index_pair!=cartesianToCompressed_.end())
                return index_pair->second;
            else
                return -1;
        }

        template<class Grid>
        void init(const Grid& grid)
        {
            using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
            CartesianIndexMapper cartmapper(grid);
            std::size_t num_cells = grid.leafGridView().size(0);
            is_interior_.resize(num_cells);
            centroids_.resize(num_cells);
            using GridView = typename Grid::LeafGridView;
            using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
            ElementMapper elemMapper(grid.leafGridView(), Dune::mcmgElementLayout());
            for (const auto& element : elements(grid.leafGridView()))
            {
                const auto elemIdx = elemMapper.index(element);
                unsigned cartesianCellIdx = cartmapper.cartesianIndex(elemIdx);
                auto geom = element.geometry();
                cartesianToCompressed_[cartesianCellIdx] = elemIdx;
                if (element.partitionType() == Dune::InteriorEntity)
                {
                    is_interior_[elemIdx] = 1;
                }
                else
                {
                    is_interior_[elemIdx] = 0;
                }
                centroids_[elemIdx] = geom.center();
            }

        }
        Point3D centroid(int cell_index){return centroids_[cell_index];}
    private:
        std::unordered_map<int,int> cartesianToCompressed_;
        std::vector<Point3D> centroids_;
        std::vector<bool> is_interior_;
        //std::array<int,3> cartDim_;
    };

    template<class Grid>
    int findCell(const Grid& grid, const Dune::FieldVector<double,3> point){
        // bad brut forse search
        using GridView = typename Grid::LeafGridView;
        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
        ElementMapper elemMapper(grid.leafGridView(), Dune::mcmgElementLayout());
        for (const auto& element : elements(grid.leafGridView()))
        {
            const auto elemIdx = elemMapper.index(element);
            auto geom = element.geometry();
            auto local = geom.local(point);
            auto ref_el = Dune::referenceElement(geom);
            if(ref_el.checkInside(local)){
                return elemIdx;
            }

        }
        return -1;
    }
}
namespace external{
    std::vector<size_t>
    findCloseCellIndices(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
                         cvf::BoundingBox& bb);

    // copy from RigEclipseWellLogExtractor
    void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Opm::EclipseGrid& m_grid);
    void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Dune::CpGrid& grid);

}
