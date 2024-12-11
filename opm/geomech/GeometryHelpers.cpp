#include "config.h"
#include "GeometryHelpers.hpp"
namespace external{
    std::vector<size_t>
    findCloseCellIndices(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
                         cvf::BoundingBox& bb){
        //cvf::BoundingBox bb;
        //bb.add( p1 );
        //bb.add( p2 );
        std::vector<size_t> closeCells;
        //this->findIntersectingCells( bb, &closeCells );
        m_cellSearchTree->findIntersections( bb, &closeCells);
        return closeCells;
    }
  void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, 
                            const Dune::CpGrid& grid){
        using GridView = Dune::CpGrid::LeafGridView;
        const auto& gv = grid.leafGridView();
        //auto nx = m_grid.getNX();
        //auto ny = m_grid.getNY();
        //auto nz = m_grid.getNZ();
        //auto gv = grid.leafGridView();
        size_t cellCount = gv.size(0);
        std::vector<size_t>           cellIndicesForBoundingBoxes;
        std::vector<cvf::BoundingBox> cellBoundingBoxes;

        std::array<double, 3> cornerPointArray;
        cvf::Vec3d cornerPoint;
        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
        ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interally
	    for (const auto& element : Dune::elements(gv))
        {
            int index = mapper.index(element);
            auto geom = element.geometry();
            assert(geom.corners() == 8);
            cvf::BoundingBox cellBB;
            cvf::Vec3d cornerPoint;
            //NB order should not matter when adding to bounding box: dune ordring and resinsight ordering is different
            // dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
            for (std::size_t l = 0; l < 8; l++) {
                auto cornerPointArray = geom.corner(l);
                cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);     
                cellBB.add(cornerPoint);
            }
            cellIndicesForBoundingBoxes.emplace_back(index);
            cellBoundingBoxes.emplace_back( cellBB);
        }
        m_cellSearchTree = new cvf::BoundingBoxTree;
        m_cellSearchTree->buildTreeFromBoundingBoxes( cellBoundingBoxes, &cellIndicesForBoundingBoxes );
    }

      void buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Opm::EclipseGrid& m_grid){

        auto nx = m_grid.getNX();
        auto ny = m_grid.getNY();
        auto nz = m_grid.getNZ();

        size_t cellCount = nx * ny * nz;
        std::vector<size_t>           cellIndicesForBoundingBoxes;
        std::vector<cvf::BoundingBox> cellBoundingBoxes;

// #pragma omp parallel
//         {
            size_t threadCellCount = cellCount;

            std::vector<size_t>           threadIndicesForBoundingBoxes;
            std::vector<cvf::BoundingBox> threadBoundingBoxes;

            threadIndicesForBoundingBoxes.reserve( threadCellCount );
            threadBoundingBoxes.reserve( threadCellCount );

            std::array<double, 3> cornerPointArray;
            cvf::Vec3d cornerPoint;

// #pragma omp for
            for (int cIdx = 0; cIdx < (int)cellCount; ++cIdx) {
                const auto[i,j,k] = m_grid.getIJK(cIdx);
                cvf::BoundingBox cellBB;
                for (std::size_t l = 0; l < 8; l++) {
                     cornerPointArray = m_grid.getCornerPos(i,j,k,l);
                     cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
                     cellBB.add(cornerPoint);
                }

                if (cellBB.isValid()) {
                    threadIndicesForBoundingBoxes.emplace_back(cIdx);
                    threadBoundingBoxes.emplace_back(cellBB);
                }
            }

            threadIndicesForBoundingBoxes.shrink_to_fit();
            threadBoundingBoxes.shrink_to_fit();
// #pragma omp critical
//          {
             cellIndicesForBoundingBoxes.insert( cellIndicesForBoundingBoxes.end(),
                                                 threadIndicesForBoundingBoxes.begin(),
                                                 threadIndicesForBoundingBoxes.end() );

             cellBoundingBoxes.insert( cellBoundingBoxes.end(), threadBoundingBoxes.begin(), threadBoundingBoxes.end() );
        //  }
    // } #pragma omp parallel
         m_cellSearchTree = new cvf::BoundingBoxTree;
         m_cellSearchTree->buildTreeFromBoundingBoxes( cellBoundingBoxes, &cellIndicesForBoundingBoxes );
    }
}
namespace Opm{


}
