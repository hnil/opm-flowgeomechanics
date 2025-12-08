#include "GeometryHelpers.hpp"
#include <opm/geomech/CGAL_helper.hpp>
#include "config.h"
namespace external
{
std::vector<size_t>
findCloseCellIndices(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, cvf::BoundingBox& bb)
{
    // cvf::BoundingBox bb;
    // bb.add( p1 );
    // bb.add( p2 );
    std::vector<size_t> closeCells;
    // this->findIntersectingCells( bb, &closeCells );
    m_cellSearchTree->findIntersections(bb, &closeCells);
    return closeCells;
}
void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seeds, const Dune::CpGrid& grid)
{
    using GridView = Dune::CpGrid::LeafGridView;
    const auto& gv = grid.leafGridView();
    // auto nx = m_grid.getNX();
    // auto ny = m_grid.getNY();
    // auto nz = m_grid.getNZ();
    // auto gv = grid.leafGridView();
    size_t cellCount = gv.size(0);
    std::vector<size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    //    std::array<double, 3> cornerPointArray;
    //cvf::Vec3d cornerPoint;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interally
    entity_seeds.resize(cellCount);
    for (const auto& element : Dune::elements(gv)) {
        int index = mapper.index(element);
        entity_seeds[index] = element.seed();
        auto geom = element.geometry();
        assert(geom.corners() == 8);
        cvf::BoundingBox cellBB;
        cvf::Vec3d cornerPoint;
        // NB order should not matter when adding to bounding box: dune ordring and resinsight ordering
        // is different
        //  dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
        for (std::size_t l = 0; l < 8; l++) {
            auto cornerPointArray = geom.corner(l);
            cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
            cellBB.add(cornerPoint);
        }
        cellIndicesForBoundingBoxes.emplace_back(index);
        cellBoundingBoxes.emplace_back(cellBB);
    }
    m_cellSearchTree = new cvf::BoundingBoxTree;
    m_cellSearchTree->buildTreeFromBoundingBoxes(cellBoundingBoxes, &cellIndicesForBoundingBoxes);
    // check bounding box tree
    for (const auto& element : Dune::elements(gv)) {
         int index = mapper.index(element);
        auto geom = element.geometry();
        //external::cvf::BoundingBox bb;
        using Vec3d = external::cvf::Vec3d;
        auto vertex = geom.center();
        Vec3d point(vertex[0], vertex[1], vertex[2]);
        //bb.add(point);
        if (element.partitionType() == Dune::InteriorEntity){
            int cell = external::cellOfPoint(m_cellSearchTree, grid, entity_seeds, point);
            assert(cell == index);  
        }else{
            // may get the neares cell in partition
            //assert(cell == -1);  
        }
    }
}

int cellOfPoint(const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
                const Dune::CpGrid& grid,
                const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seed,
                const cvf::Vec3d& point){
        external::cvf::BoundingBox bb;
        bb.add(point);
        std::vector<size_t> cells = external::findCloseCellIndices(m_cellSearchTree, bb);
        //int count =0;
        int outcell = -1;
        double min_distance =1e99;
        Dune::FieldVector<double,3> local;
        for(const auto& cell: cells){
            auto element = grid.entity(entity_seed[cell]);
            if (element.partitionType() != Dune::InteriorEntity){
                // outside of parition is considered outside model
                continue;
            }
            auto geom = element.geometry();
            auto center = geom.center();
            Dune::FieldVector<double,3> vec;
            vec[0] = point[0];
            vec[1] = point[1];
            vec[2] = point[2];
            auto d_vec = vec-center;
            if(d_vec.two_norm() < min_distance){
                min_distance = d_vec.two_norm();
                outcell = cell;
            }
            //NB should probably point is inside model
            // // consistent inside is difficult due to round off errors
            // bool check_inside = false; 
            // if (check_inside)
            // {
            //     local = geom.local(vec);
            //     bool inside = true;
            //     for (int i = 0; i < 3; ++i) {
            //         if (local[i] < 0) {
            //             inside = false;
            //         }
            //         if (local[i] > 1) {
            //             inside = false;
            //         }
            //     }
            //     auto refEl = Dune::ReferenceElements<double, 3>::cube();
            //     assert(inside == refEl.checkInside(local));
            //     if (inside) {
            //         count += 1;
            //         outcell = cell;
            //     }
            // }
     }
     if(outcell>=0){//} && count ==1){
        return outcell;
     }else{
        // std::cout << "outcell" << outcell << std::endl;
        // std::cout << "count" << count << std::endl;
        // std::cout << "local" << local[0] << " " << local[1] << " " << local[2] << " " << std::endl;
        //assert(false);
        //assert(cells.size()==0);
        //std::stringstream os;
        std::cout << "Point outside bounding box or overlap cells?" << std::endl;
        //OpmLog::info(os.str())
        return -1;
     }
}

void
cellsOfTri(std::vector<int>& cells,
           std::vector<double>& areas,
           std::vector<Dune::FieldVector<double,3>>& centroids,
           const cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree,
           const Dune::CpGrid& grid,
           const std::vector<Dune::CpGrid::Codim<0>::Entity::EntitySeed>& entity_seed,
           const std::vector<std::array<double, 3>>& tri_corners)
{
    external::cvf::BoundingBox bb;
    for (const auto& tpoint : tri_corners) {
        cvf::Vec3d point(tpoint[0], tpoint[1], tpoint[2]);
        bb.add(point);
    }
    cells.resize(0);
    areas.resize(0);
    std::vector<size_t> cells_tmp = external::findCloseCellIndices(m_cellSearchTree, bb);
    // int count =0;
    int outcell = -1;
    double min_distance = 1e99;
    Dune::FieldVector<double, 3> local;
    for (const auto& cell : cells_tmp) {
        auto element = grid.entity(entity_seed[cell]);
        if (element.partitionType() != Dune::InteriorEntity) {
            // outside of parition is considered outside model
            continue;
        }
        std::vector<std::array<double, 3>> hex_corners;
        {
            auto geom_hex = entity_seed[cell].geometry();
            for (size_t i = 0; i < size_t(geom_hex.corners()); ++i) {
                auto corner = geom_hex.corner(i);
                std::array<double, 3> tmp({corner[0], corner[1], corner[2]});
                hex_corners.push_back(tmp);
            }
        }
        double area = area_of_intersection(hex_corners, tri_corners);
        if (area > 1e-12) {
            cells.push_back(cell);
            areas.push_back(area);
            auto centroid = element.geometry().center();
            centroids.push_back(centroid);
        }
    }
}





void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const Opm::EclipseGrid& m_grid)
{

    auto nx = m_grid.getNX();
    auto ny = m_grid.getNY();
    auto nz = m_grid.getNZ();

    size_t cellCount = nx * ny * nz;
    std::vector<size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    // #pragma omp parallel
    //         {
    size_t threadCellCount = cellCount;

    std::vector<size_t> threadIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> threadBoundingBoxes;

    threadIndicesForBoundingBoxes.reserve(threadCellCount);
    threadBoundingBoxes.reserve(threadCellCount);

    std::array<double, 3> cornerPointArray;
    cvf::Vec3d cornerPoint;

    // #pragma omp for
    for (int cIdx = 0; cIdx < (int)cellCount; ++cIdx) {
        const auto [i, j, k] = m_grid.getIJK(cIdx);
        cvf::BoundingBox cellBB;
        for (std::size_t l = 0; l < 8; l++) {
            cornerPointArray = m_grid.getCornerPos(i, j, k, l);
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
    cellIndicesForBoundingBoxes.insert(cellIndicesForBoundingBoxes.end(),
                                       threadIndicesForBoundingBoxes.begin(),
                                       threadIndicesForBoundingBoxes.end());

    cellBoundingBoxes.insert(
        cellBoundingBoxes.end(), threadBoundingBoxes.begin(), threadBoundingBoxes.end());
    //  }
    // } #pragma omp parallel
    m_cellSearchTree = new cvf::BoundingBoxTree;
    m_cellSearchTree->buildTreeFromBoundingBoxes(cellBoundingBoxes, &cellIndicesForBoundingBoxes);
}
} // namespace external
namespace Opm
{


}
