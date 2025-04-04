#include "vemutils.hpp"
namespace vem{
  using PolyGrid = Dune::PolyhedralGrid<3, 3>;
#ifdef HAVE_ALUGRID
  using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming >;
#endif
  void getGridVectors(const PolyGrid& grid, std::vector<double>& coords,
                      std::vector<int>& num_cell_faces,
                      std::vector<int>& num_face_corners,
                      std::vector<int>& face_corners)
  {
    const UnstructuredGrid& ungrid = grid;
    static constexpr int dim = PolyGrid::dimension;
    coords.resize(ungrid.number_of_nodes*dim);
    for(int i=0; i< ungrid.number_of_nodes; ++i){
      for(int j=0; j < dim; ++j){
        coords[3*i+j] = ungrid.node_coordinates[3*i+j];
      }
    }
    num_cell_faces.resize(ungrid.number_of_cells);
    for(int i=0; i< ungrid.number_of_cells; ++i){
      num_cell_faces[i]=ungrid.cell_facepos[i+1]-ungrid.cell_facepos[i];
    }
    num_face_corners.resize(ungrid.cell_facepos[ungrid.number_of_cells]);
    int tot_num_face_corners = 0;
    //for(int i=0; i< num_cellfaces; ++i){
    for(int cell=0; cell < ungrid.number_of_cells; ++cell){
      for(grid_size_t hface = ungrid.cell_facepos[cell]; hface < ungrid.cell_facepos[cell+1]; hface ++){
        int face = ungrid.cell_faces[hface];//ungrid.cell_facepos[i]];
        int num_local_corners = ungrid.face_nodepos[face+1] - ungrid.face_nodepos[face];
        num_face_corners[hface] = num_local_corners;
        // NB maybe we should order with outwards normal
        if(cell == ungrid.face_cells[2*face]){
          for(grid_size_t j  = ungrid.face_nodepos[face]; j <  ungrid.face_nodepos[face+1]; ++j){
            face_corners.push_back(ungrid.face_nodes[j]);
          }
        }else{
          // flip orientation for hface
          for(grid_size_t j = ungrid.face_nodepos[face+1]-1; j >=  ungrid.face_nodepos[face]; --j){
            face_corners.push_back(ungrid.face_nodes[j]);
          }
        }
        tot_num_face_corners += num_local_corners;
      }
    }
    assert(face_corners.size() == static_cast<std::size_t>(tot_num_face_corners));
  }
  
  void getGridVectorsDune(const PolyGrid& grid, std::vector<double>& coords,
                          std::vector<int>& num_cell_faces,
                          std::vector<int>& num_face_corners,
                          std::vector<int>& face_corners)
  {
    using namespace std;
    using namespace Dune;
    const auto& gv = grid.leafGridView();
    //int loadcase = -1;
    //Dune::FieldVector<ctype,comp> eps0 = {1, 1, 1, 0, 0, 0};
    //eps0 = 0;

    // start VEM assembly
    // make global point coordinate vector
    //vector<double> coords;
    for (const auto& v : vertices(gv)) {
      const auto c = v.geometry().corner(0);
      coords.insert(coords.end(), c.begin(), c.end());
    }

    const auto& ixset = gv.indexSet();
    // count cell faces
    //vector<int> num_cell_faces;
    for (const auto& c : elements(gv))
      num_cell_faces.push_back(Dune::subEntities(c, Dune::Codim<1>{}).size());
    const int tot_num_cfaces = accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

    // count face corners
    //vector<int> num_face_corners;
    for (const auto& c : elements(gv))
      for (const auto& f : Dune::subEntities(c, Dune::Codim<1>{}))
        num_face_corners.push_back(f.geometry().corners());

    // establish all face corners
    //vector<int> face_corners;
    for (const auto& c : elements(gv))
      for (const auto& f: Dune::subEntities(c, Dune::Codim<1>{}))
        for (int i = 0, f_ix = 0;
             f_ix != tot_num_cfaces && i != num_face_corners[f_ix];
             ++i, f_ix += (i == num_face_corners[f_ix]))
          face_corners.push_back(ixset.subIndex(f, i, 3));

    // correct order of nodes in each face, to ensure they are mentioned in
    // clockwise or counterclockwise order. @@ This is a hack that might only work
    // for 4-faces!
    // for (int i = 0, j = 0; i != (int)num_face_corners.size(); j += num_face_corners[i++])
    //   swap(face_corners[j], face_corners[j+1]);

    // body force

    // dirichlet boundary conditions

  }



  void getGridVectors(const Dune::CpGrid& grid, std::vector<double>& coords,
                      std::vector<int>& num_cell_faces,
                      std::vector<int>& num_face_corners,
                      std::vector<int>& face_corners)
  {
    //assert(false);
    // std::vector<double> coords_tmp;
    // std::vector<int> num_cell_faces_tmp;
    // std::vector<int> num_face_corners_tmp;
    // std::vector<int> face_corners_tmp;
    // getGridVectors<Dune::CpGrid>(grid, coords,
    //                               num_cell_faces,
    //                               num_face_corners,
    //                               face_corners);
    // return;
    // coords = coords_tmp;
    // num_cell_faces = num_cell_faces_tmp;
    // num_face_corners =  num_face_corners_tmp;
    // face_corners = face_corners_tmp;
    //return;
    //using GridType = Dune::CpGrid;
    //static constexpr int dim = GridType::dimension;
    using namespace std;
    using namespace Dune;
    const auto& gv = grid.leafGridView();

    for (const auto& v : vertices(gv)) {
      const auto c = v.geometry().corner(0);
      coords.insert(coords.end(), c.begin(), c.end());
    }
    //const auto& faces = grid.getFaceToPoint();
    //const auto& facePos = grid.getFacePos();
    for (const auto& cell : elements(gv)){
      int cellIdx = gv.indexSet().index(cell);
      int nf = grid.numCellFaces(cellIdx);
      //assert(nf==6);
      num_cell_faces.push_back(nf);
      for(int f=0; f < nf; ++f){
        auto face = grid.cellFace(cellIdx,f);
        auto faceSize = grid.numFaceVertices(face);
        num_face_corners.push_back(faceSize);
        //assert(faceSize == 4);
        auto out_cell = grid.faceCell(face,1);
        auto in_cell = grid.faceCell(face,0);
        if(out_cell!=cellIdx){
          assert(in_cell==cellIdx);
          for(int v = 0; v < faceSize; ++v){
            int fv = grid.faceVertex(face,v);
            //const auto& point = cpgrid::Entity<3>(*currentData().back(), fv, true);
            //const auto& localId = currentData().back()->localIdSet().id(point);
          //const auto& globalId = currentData().back()->globalIdSet().id(point);
            face_corners.push_back(fv);
          }
        }else{
          //assert(in_cell==cellIdx);
          for(int v = faceSize-1; v >  -1; --v){
            int fv = grid.faceVertex(face,v);
            //const auto& point = cpgrid::Entity<3>(*currentData().back(), fv, true);
            //const auto& localId = currentData().back()->localIdSet().id(point);
          //const auto& globalId = currentData().back()->globalIdSet().id(point);
            face_corners.push_back(fv);
          }
          
        }
      }

    }
    //coords = coords_tmp;
    //num_cell_faces = num_cell_faces_tmp;
    //num_face_corners =  num_face_corners_tmp;
    //face_corners = face_corners_tmp;
    // for(size_t i=0; i < face_corners.size(); ++i){
    //   std::cout << face_corners[i] << " " << face_corners_tmp[i] << std::endl; 
    // }
    //    assert(false);
   
  }

}
