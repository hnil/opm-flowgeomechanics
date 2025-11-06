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
#include <opm/geomech/vem/vem.hpp>
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
  int assemble_mech_system_3D_dune(const Dune::CpGrid& grid, const double* const points,
                        const int num_cells,
                        const int* const num_cell_faces, // cell faces per cell
                        const int* const num_face_corners, // corners per face
                        const int* const face_corners,
                        const double* const young,
                        const double* const poisson,
                        const double* const body_force, // 3 * number of cells
                        const int num_fixed_dofs, // dirichlet
                        const int* const fixed_dof_ixs, // indices must be sorted
                        const double* const fixed_dof_values,
                        const int num_neumann_faces,
                        const int* const neumann_faces,
                        const double* const neumann_forces, // 3 * number of neumann faces
                        std::vector<std::tuple<int, int, double>>& A_entries,
                        std::vector<double>& b,
                             const vem::StabilityChoice stability_choice,
                             bool reduce_boundary);

void
potential_gradient_force_3D_dune(const Dune::CpGrid& grid, const double* const points,
                            const int num_cells,
                            const int* const num_cell_faces, // cell faces per cell
                            const int* const num_face_corners, // corners per cellface
                            const int* const face_corners,
                            const double* const field,
                                 std::vector<double>& fgrad,
                                 std::vector<std::tuple<int, int, double>>& div,
                                 bool get_matrix);

Dune::BlockVector<Dune::FieldVector<double,1>> smoothCellVector(const Dune::CpGrid& grid,const Dune::BlockVector<Dune::FieldVector<double,1>>& cell_vector);                                

                 
  
}
#endif
