#ifndef VEM_UTILS_HPP
#define VEM_UTILS_HPP

#include <array>
#include <tuple>
#include <vector>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/CpGrid.hpp>
#ifdef HAVE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <opm/geomech/boundaryutils.hh>
namespace vem
{
    using PolyGrid = Dune::PolyhedralGrid<3, 3>;
    #ifdef HAVE_ALUGRID
    using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming >;
    #endif
    void getGridVectors(const PolyGrid& grid, std::vector<double>& coords,
                        std::vector<int>& num_cell_faces,
                        std::vector<int>& num_face_corners,
                        std::vector<int>& face_corners);

    void getGridVectorsDune(const PolyGrid& grid, std::vector<double>& coords,
                            std::vector<int>& num_cell_faces,
                            std::vector<int>& num_face_corners,
                            std::vector<int>& face_corners);

    template<class GridType>
    void getGridVectors(const GridType& grid, std::vector<double>& coords,
                std::vector<int>& num_cell_faces,
                std::vector<int>& num_face_corners,
                std::vector<int>& face_corners)
{
   static constexpr int dim = GridType::dimension;
    using namespace std;
    using namespace Dune;
    const auto& gv = grid.leafGridView();

  for (const auto& v : vertices(gv)) {
    const auto c = v.geometry().corner(0);
    coords.insert(coords.end(), c.begin(), c.end());
  }

  for (const auto& cell : elements(gv)){
      num_cell_faces.push_back(6);
      for(int i = 0; i < 6; ++i){
          auto faceDir = Opm::Elasticity::faceToFaceDir(i);
          std::array<int,4> nodes = Opm::Elasticity::faceDirToNodes(faceDir);
          num_face_corners.push_back(nodes.size());
          for(auto nind: nodes){
              auto global_ind = gv.indexSet().subIndex(cell,nind,dim);
              face_corners.push_back(global_ind);
          }
      }
  }
}


void getGridVectors(const Dune::CpGrid& grid, std::vector<double>& coords,
                       std::vector<int>& num_cell_faces,
                       std::vector<int>& num_face_corners,
                       std::vector<int>& face_corners);

}
#endif
