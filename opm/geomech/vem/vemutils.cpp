#include "config.h"
#include "vemutils.hpp"
#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/materials.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/grid/common/mcmgmapper.hh>
namespace vem
{
using PolyGrid = Dune::PolyhedralGrid<3, 3>;
#ifdef HAVE_ALUGRID
using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
#endif
void
getGridVectors(const PolyGrid& grid,
               std::vector<double>& coords,
               std::vector<int>& num_cell_faces,
               std::vector<int>& num_face_corners,
               std::vector<int>& face_corners)
{
    const UnstructuredGrid& ungrid = grid;
    static constexpr int dim = PolyGrid::dimension;
    coords.resize(ungrid.number_of_nodes * dim);
    for (int i = 0; i < ungrid.number_of_nodes; ++i) {
        for (int j = 0; j < dim; ++j) {
            coords[3 * i + j] = ungrid.node_coordinates[3 * i + j];
        }
    }
    num_cell_faces.resize(ungrid.number_of_cells);
    for (int i = 0; i < ungrid.number_of_cells; ++i) {
        num_cell_faces[i] = ungrid.cell_facepos[i + 1] - ungrid.cell_facepos[i];
    }
    num_face_corners.resize(ungrid.cell_facepos[ungrid.number_of_cells]);
    int tot_num_face_corners = 0;
    // for(int i=0; i< num_cellfaces; ++i){
    for (int cell = 0; cell < ungrid.number_of_cells; ++cell) {
        for (grid_size_t hface = ungrid.cell_facepos[cell]; hface < ungrid.cell_facepos[cell + 1];
             hface++) {
            int face = ungrid.cell_faces[hface]; // ungrid.cell_facepos[i]];
            int num_local_corners = ungrid.face_nodepos[face + 1] - ungrid.face_nodepos[face];
            num_face_corners[hface] = num_local_corners;
            // NB maybe we should order with outwards normal
            if (cell == ungrid.face_cells[2 * face]) {
                for (grid_size_t j = ungrid.face_nodepos[face]; j < ungrid.face_nodepos[face + 1]; ++j) {
                    face_corners.push_back(ungrid.face_nodes[j]);
                }
            } else {
                // flip orientation for hface
                for (grid_size_t j = ungrid.face_nodepos[face + 1] - 1; j >= ungrid.face_nodepos[face];
                     --j) {
                    face_corners.push_back(ungrid.face_nodes[j]);
                }
            }
            tot_num_face_corners += num_local_corners;
        }
    }
    assert(face_corners.size() == static_cast<std::size_t>(tot_num_face_corners));
}

void
getGridVectorsDune(const PolyGrid& grid,
                   std::vector<double>& coords,
                   std::vector<int>& num_cell_faces,
                   std::vector<int>& num_face_corners,
                   std::vector<int>& face_corners)
{
    using namespace std;
    using namespace Dune;
    const auto& gv = grid.leafGridView();
    // int loadcase = -1;
    // Dune::FieldVector<ctype,comp> eps0 = {1, 1, 1, 0, 0, 0};
    // eps0 = 0;

    // start VEM assembly
    // make global point coordinate vector
    // vector<double> coords;
    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().corner(0);
        coords.insert(coords.end(), c.begin(), c.end());
    }

    const auto& ixset = gv.indexSet();
    // count cell faces
    // vector<int> num_cell_faces;
    for (const auto& c : elements(gv))
        num_cell_faces.push_back(Dune::subEntities(c, Dune::Codim<1> {}).size());
    const int tot_num_cfaces = accumulate(num_cell_faces.begin(), num_cell_faces.end(), 0);

    // count face corners
    // vector<int> num_face_corners;
    for (const auto& c : elements(gv))
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {}))
            num_face_corners.push_back(f.geometry().corners());

    // establish all face corners
    // vector<int> face_corners;
    for (const auto& c : elements(gv))
        for (const auto& f : Dune::subEntities(c, Dune::Codim<1> {}))
            for (int i = 0, f_ix = 0; f_ix != tot_num_cfaces && i != num_face_corners[f_ix];
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



void
getGridVectors(const Dune::CpGrid& grid,
               std::vector<double>& coords,
               std::vector<int>& num_cell_faces,
               std::vector<int>& num_face_corners,
               std::vector<int>& face_corners)
{
    // assert(false);
    //  std::vector<double> coords_tmp;
    //  std::vector<int> num_cell_faces_tmp;
    //  std::vector<int> num_face_corners_tmp;
    //  std::vector<int> face_corners_tmp;
    //  getGridVectors<Dune::CpGrid>(grid, coords,
    //                                num_cell_faces,
    //                                num_face_corners,
    //                                face_corners);
    //  return;
    //  coords = coords_tmp;
    //  num_cell_faces = num_cell_faces_tmp;
    //  num_face_corners =  num_face_corners_tmp;
    //  face_corners = face_corners_tmp;
    // return;
    // using GridType = Dune::CpGrid;
    // static constexpr int dim = GridType::dimension;
    using namespace std;
    using namespace Dune;
    const auto& gv = grid.leafGridView();
    int nc = grid.numCells();
    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().corner(0);
        coords.insert(coords.end(), c.begin(), c.end());
    }
    // const auto& faces = grid.getFaceToPoint();
    // const auto& facePos = grid.getFacePos();
    for (const auto& cell : elements(gv)) {
        int cellIdx = gv.indexSet().index(cell);
        int nf = grid.numCellFaces(cellIdx);
        // assert(nf==6);
        num_cell_faces.push_back(nf);
        for (int f = 0; f < nf; ++f) {
            auto face = grid.cellFace(cellIdx, f);
            auto faceSize = grid.numFaceVertices(face);
            num_face_corners.push_back(faceSize);
            // assert(faceSize == 4);
            // auto numface_cells = grid.numCellFaces();
            // assert(numface_cells == 2);
            auto out_cell = grid.faceCell(face, 1);
            auto in_cell = grid.faceCell(face, 0);
            assert((out_cell < nc) && (in_cell < nc));
            assert(out_cell != in_cell);
            if (out_cell != cellIdx) {
                assert(in_cell == cellIdx);
                for (int v = 0; v < faceSize; ++v) {
                    int fv = grid.faceVertex(face, v);
                    // const auto& point = cpgrid::Entity<3>(*currentData().back(), fv, true);
                    // const auto& localId = currentData().back()->localIdSet().id(point);
                    // const auto& globalId = currentData().back()->globalIdSet().id(point);
                    face_corners.push_back(fv);
                }
            } else {
                // assert(in_cell==cellIdx);
                for (int v = faceSize - 1; v > -1; --v) {
                    int fv = grid.faceVertex(face, v);
                    // const auto& point = cpgrid::Entity<3>(*currentData().back(), fv, true);
                    // const auto& localId = currentData().back()->localIdSet().id(point);
                    // const auto& globalId = currentData().back()->globalIdSet().id(point);
                    face_corners.push_back(fv);
                }
            }
        }
    }



    // coords = coords_tmp;
    // num_cell_faces = num_cell_faces_tmp;
    // num_face_corners =  num_face_corners_tmp;
    // face_corners = face_corners_tmp;
    //  for(size_t i=0; i < face_corners.size(); ++i){
    //    std::cout << face_corners[i] << " " << face_corners_tmp[i] << std::endl;
    //  }
    //     assert(false);
}

bool cellFemCell(const Dune::CpGrid& grid, int cellIdx,
        Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapper){
        bool is_ok = true;
        int nf = grid.numCellFaces(cellIdx);
        // if(nf != 6){
        //     //not needed
        //     is_ok = false;
        //     return false;
        // }
        //auto cartMapper = Dune::CartesianIndexMapper<Dune::CpGrid>(grid);
        // using ElementMapper =   
        //   Dune::MultipleCodimMultipleGeomTypeMapper<typename Dune::CpGrid::LeafGridView>;
        
        for (int f = 0; f < nf; ++f) {
            auto face = grid.cellFace(cellIdx, f);
            auto faceSize = grid.numFaceVertices(face);
            // if(faceSize !=4){
            //     is_ok = false;
            //     return false;
            // }
            auto out_cell = grid.faceCell(face, 1);
            auto in_cell = grid.faceCell(face, 0);

            if(in_cell<0 || out_cell <0){
              // assume boundary is ok
                //is_ok = false;
                //return false;
                continue;

            }
            std::array<int,3> cartCoord1;
            cartMapper.cartesianCoordinate(in_cell, cartCoord1);
            std::array<int,3> cartCoord2;
            cartMapper.cartesianCoordinate(out_cell, cartCoord2);
            if(cartCoord1[0] ==cartCoord2[0] && cartCoord1[1] ==cartCoord2[1]){
                // same column
                continue;
            }
            if(std::abs(cartCoord1[2]-cartCoord2[2]) != 0){
                // not normal
                is_ok = false;
                return false;
            }
        }
        assert(is_ok==true);
        return true;    
    }

    int
assemble_mech_system_3D_dune(const Dune::CpGrid& grid,const double* const points,
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
                        bool reduce_boundary)
// ----------------------------------------------------------------------------
{
    // preliminary computations
  const int tot_num_cell_faces = std::accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
    const int tot_num_face_corners
      = std::accumulate(num_face_corners, num_face_corners + tot_num_cell_faces, 0);
    const int num_points = *std::max_element(face_corners, face_corners + tot_num_face_corners) + 1;

    // assemble full system matrix
    A_entries.clear();
    b.resize(num_points * 3);
    std::fill(b.begin(), b.end(), 0);

    // loop over cells and assemble system matrix
    std::vector<int> loc_indexing;
    std::vector<double> loc; // use as local 'scratch' vector
    std::array<double, 3> centroid; // will contain the centroid for the currently treated cell
    int cf_ix = 0; // index to first cell face for the current cell
    int fcorners_start = 0; // index to first cell corner for current cell

    // cout << "Starting assembly" << endl;
    //for (int c = 0; c != num_cells; ++c) {
    // for FEM ---------------------------------
     static constexpr int dim = 3;//GridType::dimension;
     static constexpr int comp = 3 + (dim - 2) * 3;
     static constexpr int bfunc = 4 + (dim - 2) * 4;
     static constexpr int esize = dim*bfunc;
     //int loadcase = -1;
     //Dune::FieldVector<double,comp> eps0 = {1, 1, 1, 0, 0, 0};
  //eps0 = 0;
  //ASMHandler<GridType> A(grid_);//need to set this.
  //A.getOperator() = 0;
 
     using ctype = double;
      Dune::FieldMatrix<ctype,comp,comp> C;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> K;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> Aq;
      //Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc>* KP=0;
      //Dune::FieldVector<ctype,dim*bfunc> ES;
      //Dune::FieldVector<ctype,dim*bfunc>* EP=0;
      Opm::Elasticity::Elasticity E(grid);

    //----------------------------------------------
    // find if cells are of fem type
      const auto gv = grid.leafGridView();
      std::vector<bool> cell_fem_type(gv.size(0),false);
      bool use_fem_type = StabilityChoice::FEM == stability_choice;
      if(use_fem_type){
        //const auto& idx  = grid.getCellIndexSet();
        //const auto& rIdx = grid.getCellRemoteIndices();
        // Dune::Interface cominterface;
        // using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::OwnerSet;//Dune::EnumItem<AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>;
        // using AllSet =  Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;//Dune::AllSet<AttributeSet>;
        // OwnerSet soureFlags;
        // //AllSet destFlags;
        // AllSet destFlags;
        // cominterface.build(rIdx, soureFlags, destFlags);
        // Dune::BufferedCommunicator cell_cell_comm;
        // using Vector = std::vector<bool>;
        // cell_cell_comm.template build<Vector>(cominterface);
        auto cartMapper = Dune::CartesianIndexMapper<Dune::CpGrid>(grid);
        Dune::MultipleCodimMultipleGeomTypeMapper<typename Dune::CpGrid::LeafGridView> elementMapper(gv,Dune::mcmgElementLayout());
        for(const auto elem : elements(gv)){
          int c = elementMapper.index(elem);
          //int c = gv.indexSet().index(elem);
          bool is_femcell = cellFemCell(grid,c, cartMapper);
          cell_fem_type[c] = is_femcell;
        }
        //cell_cell_comm.template forward<Dune::CopyGatherScatter<Vector>>(cell_fem_type, cell_fem_type);
        // couldhave used this from linear solver interfae
        //cell_cell_comm.copyOwnerToAll(fem_type, fem_type);
        if(grid.comm().size()>1){
           //grid.cellScatterGatherInterface() const;
           //gridpointScatterGatherInterface()
          const auto& cellcomm = grid.cellCommunication();
          cellcomm.copyOwnerToAll(cell_fem_type,cell_fem_type);
        }
      }
      
      


    for(const auto elem : elements(gv)){
        int c = gv.indexSet().index(elem);  
        // computing local stiffness matrix, writing its entries into the global matrix
       
        auto loc_stability_choice =  stability_choice;
        bool vem_type = StabilityChoice::FEM != stability_choice;
        if(!vem_type){
          // check if we can do vem on this cell
          //auto cartMapper = Dune::CartesianIndexMapper<Dune::CpGrid>(grid);
          // maybe move out the desition for parallel case
          bool is_femcell = cell_fem_type[c];//cellFemCell(grid,c, cartMapper);
          if(!is_femcell){
            // if fem and not a "normal" cell switch to vem
            vem_type = true;
            loc_stability_choice = StabilityChoice::D_RECIPE;
          }
        }

        if(vem_type){
          assemble_stiffness_matrix_3D(points,
                                     &face_corners[fcorners_start],
                                     &num_face_corners[cf_ix],
                                     num_cell_faces[c],
                                     young[c],
                                     poisson[c],
                                     loc_stability_choice,
                                     centroid,
                                     loc_indexing,
                                     loc);
          const int ncv = int(loc_indexing.size()); // number of cell vertices
          for (int i = 0; i != 3 * ncv; ++i) {
            for (int j = 0; j != 3 * ncv; ++j) {
              // using integer division to go from 'i' and 'j' to corresponding
              // local point indices, which are subsequently converted to global indexing
              const int I = 3 * loc_indexing[i / 3] + i % 3;
              const int J = 3 * loc_indexing[j / 3] + j % 3;
              const double val = loc[i * 3 * ncv + j]; // local matrix entry
              A_entries.push_back(std::tuple<int, int, double> {I, J, val});
              // const auto it = A_entrymap.find(make_pair(I, J));
              // (it != A_entrymap.end()) ? it->second += val :
                //                            A_entrymap[make_pair(I, J)] = val;
            }
          }
        }else{
          // do fem
          int cellIdx = gv.indexSet().index(elem);
        // determine geometry type of the current element and get the matching reference element
          K=0;
          int  id = 1;
          double rho = 0.0;
          Opm::Elasticity::Isotropic  material(id, young[cellIdx], poisson[cellIdx], rho);
          bool valid = material.getConstitutiveMatrix(C);
          if(!valid){
            OPM_THROW(std::runtime_error,"Not valid C matrix"); 
        }
          Dune::GeometryType gt = elem.type();
          // get a quadrature rule of order two for the given geometry type
          const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
          for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
               r != rule.end() ; ++r) {
            // compute the jacobian inverse transposed to transform the gradients
            Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
              elem.geometry().jacobianInverseTransposed(r->position());
            bool verbose = true;
            ctype detJ = elem.geometry().integrationElement(r->position());
            assert(elem.geometry().corners() == 8);
            if (detJ <= 1.e-5 && verbose) {
              double zdiff=0.0;
              for (int ii=0;ii<4;++ii)
                zdiff = std::max(zdiff, elem.geometry().corner(ii+4)[2]-elem.geometry().corner(ii)[2]);
              std::cout << " - Consider setting ctol larger than " << zdiff << std::endl;
            }

            Dune::FieldMatrix<ctype,comp,dim*bfunc> lB;
            E.getBmatrix(lB,r->position(),jacInvTra);
            
            E.getStiffnessMatrix(Aq,lB,C,detJ*r->weight());
            K += Aq;
          }  
          //for (int li =0;li<esize/dim;++li) {//nomber for nodes of reference element
          for (int li =0;li<8;++li) {//nomber for nodes of reference element
            int node_index1 = gv.indexSet().subIndex(elem,li,dim);
            for (int ki=0;ki<dim;++ki) {
              int ldof1 = li*dim+ki;
              int gdof1 = node_index1*dim+ki;
              for (int lj=0;lj<esize/dim;++lj) {
                int node_index2 = gv.indexSet().subIndex(elem,lj,dim);
                for (int kj=0;kj<dim;++kj) {
                  int ldof2 = lj*dim+kj;
                  int gdof2 = node_index2*dim+kj;
                  A_entries.push_back(std::tuple<int, int, double> {gdof1, gdof2, K[ldof1][ldof2]});
                }
              }
            }
          }

        }
        

        
        // add contribution to right hand side from body forces.  These are
        // written directly into the global right-hand side vector
        vem::compute_bodyforce_3D(points,
                             &face_corners[fcorners_start],
                             &num_face_corners[cf_ix],
                             num_cell_faces[c],
                             &centroid[0],
                             &body_force[3 * c],
                             &b[0]);
        // increment local indexing
        fcorners_start
          += std::accumulate(&num_face_corners[cf_ix], &num_face_corners[cf_ix + num_cell_faces[c]], 0);
        cf_ix += num_cell_faces[c];
    }
    // cout << "Applying forces" << endl;

    // add contribution to right hand side from applied forces.  Values are written
    // directly into the global right-hand-side vector
    vem::compute_applied_forces_3D(
        points, num_face_corners, face_corners, num_neumann_faces, neumann_faces, neumann_forces, &b[0]);

    if (reduce_boundary) {
        // cout << "Set boundary conditions: Reducing system" << endl;
        vem::reduce_system(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    } else {
        // cout << "Set boundary conditions: Set it directly in system" << endl;
        vem::set_boundary_conditions(A_entries, b, num_fixed_dofs, fixed_dof_ixs, fixed_dof_values);
    }

    // cout << "Finished assembly.  Returning." << endl;

    return b.size();
}

void
potential_gradient_force_3D_dune(const Dune::CpGrid& grid, const double* const points,
                            const int num_cells,
                            const int* const num_cell_faces, // cell faces per cell
                            const int* const num_face_corners, // corners per cellface
                            const int* const face_corners,
                            const double* const field,
                                 std::vector<double>& fgrad,
                                 std::vector<std::tuple<int, int, double>>& div,
                            bool get_matrix)
// ----------------------------------------------------------------------------
{
    // initializing result vector
  const int tot_num_faces = std::accumulate(num_cell_faces, num_cell_faces + num_cells, 0);
  const int tot_num_fcorners = std::accumulate(num_face_corners, num_face_corners + tot_num_faces, 0);
  const int tot_num_nodes = *std::max_element(face_corners, face_corners + tot_num_fcorners) + 1;
  fgrad = std::vector<double>(3 * tot_num_nodes, 0.0);

    // The contribution from a given cell to the gradient of pressure evaluated at
    // one of its corners can be obtained from vem by considering the first three componnets
    // of `volume * P * Wc`, where `P = [p p p 0 0 0]`.  Therefore, only the three first columns
    // of Wc matter, which is `diag([2*q1, 2*q2, 2*q3])`.  Moreover, `q_i = 1/(2*volume) z_i`,
    // where z_i is the face integral of the phi_i basis function of the node in question.  As such,
    // we can compute the contribution of this cell to the gradient of pressure in the given node
    // directly as `p * z_i`, and bypass the volume computation entirely.

    // Below, we compute the contribution from all corners of all faces to the global 'fgrad' vector.
    bool vem_method = false;
    std::cout << "Using FEM assembly of source term" << std::endl;
  std::vector<int> indexing;
    int cur_fcor_start = 0;
    int cur_cface_start = 0;
    //    for (int cell = 0; cell != num_cells; ++cell) {
    
    Opm::Elasticity::Elasticity elasticity(grid);
    const auto gv = grid.leafGridView();
    for(const auto elem : elements(gv)){
        int cell = gv.indexSet().index(elem);  
  
        const auto reindex = global_to_local_indexing(&face_corners[cur_fcor_start],
                                                      &num_face_corners[cur_cface_start],
                                                      num_cell_faces[cell],
                                                      indexing);
        const int num_corners = int(indexing.size());
        const auto corners_loc = pick_points<3>(points, &indexing[0], num_corners);

        const int tot_num_cellface_corners
          = std::accumulate(&num_face_corners[cur_cface_start],
                         &num_face_corners[cur_cface_start] + num_cell_faces[cell],
                         0);

        std::vector<int> faces_loc(tot_num_cellface_corners);
        std::transform(&face_corners[cur_fcor_start],
                  &face_corners[cur_fcor_start] + tot_num_cellface_corners,
                  faces_loc.begin(),
                  [&reindex](const int I) { return reindex.find(I)->second; });

        double volume; // to be computed below
        std::vector<double> outward_normals, face_centroids; // to be computed below
        std::array<double, 3> star_point, cell_centroid; // to be computed below
        vem::compute_cell_geometry(&corners_loc[0],
                              num_corners,
                              &faces_loc[0],
                              &num_face_corners[cur_cface_start],
                              num_cell_faces[cell],
                              outward_normals,
                              face_centroids,
                              cell_centroid,
                              star_point,
                              volume);

        // Since computing q involves dividing by volume, we are really computing q * volume
        // in the below call, by passing the value '1' for volume.  In other words, we are computing
        // (1/2) z_i.
        const auto qv = compute_q_3D(&corners_loc[0],
                                     num_corners,
                                     &faces_loc[0],
                                     &num_face_corners[cur_cface_start],
                                     num_cell_faces[cell],
                                     1,
                                     outward_normals);
        
        if (vem_method) {
            // fill in entries in global fgrad vector
            for (int c = 0; c != num_corners; ++c) {
                for (int d = 0; d != 3; ++d) {
                    fgrad[3 * indexing[c] + d] += 2 * field[cell] * qv[3 * c + d];
                    if (get_matrix) {
                        int I = 3 * indexing[c] + d;
                        int J = cell;
                        double val = 2 * qv[3 * c + d];
                        div.push_back(std::tuple<int, int, double> {I, J, val});
                    }
                }
            }
        } else {
            static constexpr int dim = 3; // GridType::dimension;
            static constexpr int comp = 3 + (dim - 2) * 3;
            static constexpr int bfunc = 4 + (dim - 2) * 4;
            //static constexpr int esize = dim * bfunc;
            Dune::GeometryType gt = elem.type();
            using ctype = double;
            Dune::FieldVector<ctype, dim * bfunc> ES;
            Dune::FieldMatrix<ctype, comp, comp> dummyC;
            dummyC = 0;
            // Opm::Elasticity::Isotropic  material(id, young[cellIdx], poisson[cellIdx], rho);
            // bool valid = material.getConstitutiveMatrix(C);
            // if(!valid){
            //   OPM_THROW(std::runtime_error,"Not valid C matrix");
            // }

            const Dune::QuadratureRule<ctype, dim>& rule
                = Dune::QuadratureRules<ctype, dim>::rule(gt, 2);
            ES = 0;
            for (typename Dune::QuadratureRule<ctype, dim>::const_iterator r = rule.begin();
                 r != rule.end();
                 ++r) {
                // compute the jacobian inverse transposed to transform the gradients
                Dune::FieldMatrix<ctype, dim, dim> jacInvTra
                    = elem.geometry().jacobianInverseTransposed(r->position());

                ctype detJ = elem.geometry().integrationElement(r->position());
                Dune::FieldMatrix<ctype, comp, dim * bfunc> lB;
                elasticity.getBmatrix(lB, r->position(), jacInvTra);
                Dune::FieldVector<ctype, dim * bfunc> temp;
                Dune::FieldVector<ctype, comp> eps0
                    = {1, 1, 1, 0, 0, 0}; // since pressure has one dof we can suse scalar
                // temp = Dune::FMatrixHelp::multTransposed(lB,Dune::FMatrixHelp::mult(C,eps0));
                auto Ipressure = Dune::FMatrixHelp::mult(dummyC, eps0); // just to get correct size??
                // Dune::FieldMatrix<ctype,dim*bfunc,1> Ipressure = eps0;
                Ipressure = eps0 * 1.0; //*pressure[cell_num][0];
                auto vec = Dune::FMatrixHelp::multTransposed(lB, Ipressure);
                temp = vec;
                temp *= detJ * r->weight();
                ES += temp;
            }
            assert(elem.geometry().corners() == 8);
            for (int li = 0; li < 8; ++li) { // nomber for nodes of reference element
                int node_index1 = gv.indexSet().subIndex(elem, li, dim);
                for (int ki = 0; ki < dim; ++ki) {
                    int ldof1 = li * dim + ki;
                    int gdof1 = node_index1 * dim + ki;
                    int J = cell;
                    div.push_back(std::tuple<int, int, double> {gdof1, J, ES[ldof1]});
                }
            }
        }

        // keep track of our current position in the `face_corners` array
        cur_fcor_start += tot_num_cellface_corners;
        cur_cface_start += num_cell_faces[cell];
    }
}
Dune::BlockVector<Dune::FieldVector<double,1>> smoothCellVector(const Dune::CpGrid& grid,const Dune::BlockVector<Dune::FieldVector<double,1>>& cell_vector){
  
  const auto& gv = grid.leafGridView();
  Dune::BlockVector<Dune::FieldVector<double,1>> smoothed_cell_vector(cell_vector.size());
  for(const auto& elem : elements(gv)){
    int cellIdx = gv.indexSet().index(elem);
    double sum = 0.0;
    int count = 0;
    for (const auto& intersection : intersections(gv, elem)) {
       if (intersection.boundary()) {
         continue;
        }
        auto neigh = intersection.outside();
        int cell_outside = gv.indexSet().index(neigh);
        sum += cell_vector[cell_outside][0];
        count += 1;
    }
    // include self
    sum += cell_vector[cellIdx][0];
    count += 1;
    smoothed_cell_vector[cellIdx][0] = sum / double(count);
  }
  return smoothed_cell_vector;
}

Dune::BlockVector<Dune::FieldVector<double,6>>  computeStressFem(const Dune::CpGrid& grid,
                                                                const Dune::BlockVector<Dune::FieldVector<double,1>>& disp,
                                                                const std::vector <double>& young,
                                                                const std::vector <double>& poisson){
        std::cout << "Compute stress via FEM" <<std::endl;                                                          
        static constexpr int dim = 3; // GridType::dimension;
        static constexpr int comp = 3 + (dim - 2) * 3;
        static constexpr int bfunc = 4 + (dim - 2) * 4;
        //        static constexpr int esize = dim * bfunc;
      Opm::Elasticity::Elasticity elasticity(grid);
       const auto& gv = grid.leafGridView();
        using ctype = double;
        Dune::FieldMatrix<ctype, comp, comp> C;
        Dune::FieldVector<ctype, comp> eps0;
        eps0 = 0;
        // eps0[loadcase] = 1; // NB do not understand
        Dune::BlockVector<Dune::FieldVector<double,6>> sigmacells(grid.numCells());                                                       
        for (const auto& elem: elements(gv)){
            int cellIdx = gv.indexSet().index(elem);
            double rho = 0.0;
            Opm::Elasticity::Isotropic  material(1.0, young[cellIdx], poisson[cellIdx], rho);
             bool valid = material.getConstitutiveMatrix(C);
            // determine geometry type of the current element and get the matching reference element
            Dune::GeometryType gt = elem.type();

            Dune::FieldVector<ctype, bfunc * dim> v;
            for (int li = 0; li < 8; ++li) { // nomber for nodes of reference element
                int node_index1 = gv.indexSet().subIndex(elem, li, dim);
                for (int ki = 0; ki < dim; ++ki) {
                    int ldof1 = li * dim + ki;
                    v[ldof1] = disp[node_index1*dim+ki][0];
                }
            }
            
            Dune::FieldVector<ctype, comp> sigma;
            sigma = 0;
            double volume = 0;
            // get a quadrature rule of order two for the given geometry type
            const Dune::QuadratureRule<ctype, dim>& rule = Dune::QuadratureRules<ctype, dim>::rule(gt, 2);
            for (typename Dune::QuadratureRule<ctype, dim>::const_iterator r = rule.begin(); r != rule.end(); ++r) {
                // compute the jacobian inverse transposed to transform the gradients
                Dune::FieldMatrix<ctype, dim, dim> jacInvTra = elem.geometry().jacobianInverseTransposed(r->position());

                ctype detJ = elem.geometry().integrationElement(r->position());

                volume += detJ * r->weight();

                Dune::FieldMatrix<ctype, comp, dim * bfunc> lB;
                elasticity.getBmatrix(lB, r->position(), jacInvTra);

                Dune::FieldVector<ctype, comp> s;
                elasticity.getStressVector(s, v, eps0, lB, C);
                s *= detJ * r->weight();
                sigma += s;
            }
            sigma /= volume;
            // if (Escale > 0) {
            //     sigma /= Escale / Emin;
            // }
            // switch to voits notation
            std::swap(sigma[4],sigma[5]);
            
            sigmacells[cellIdx] = sigma;
        }
        return sigmacells;
    }
   

} // namespace vem
