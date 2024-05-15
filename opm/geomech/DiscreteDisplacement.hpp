#pragma once
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include "CutDe.hpp"
namespace ddm
{
template <class Element>
std::array<Real3, 3>
getTri(const Element& elem){
    // dim = Dune::Codim<Grid::dimension>{}
    int codim = 2;//vertices in this cases
    std::array<Real3, 3> tri;
    int i = 0;
    for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<2>{})){
        // assume order is ok and orientation is fine;
        auto geom = vertex.geometry();
        auto corner  = geom.center();
        Real3 corner_tmp = make3(corner[0],corner[1],corner[2]);
        tri[i] = corner_tmp;
        ++i;
    }
    return tri;
}
template <class Element>
Dune::FieldVector<double, 3>
TDDispFS(const Dune::FieldVector<double, 3>& obs,
         const Element& elem,
         const Dune::FieldVector<double, 3>& slip,
         double nu)
{
    std::array<Real3, 3> tri = getTri(elem);
    Real3 obs_tmp = make3(obs[0], obs[1], obs[2]);
    Real3 slip_tmp = make3(slip[0], slip[1], slip[2]);
    Real3 disp_tmp = disp_fs(obs_tmp, tri, slip_tmp, nu);
    Dune::FieldVector<double, 3> disp;
    disp[0] = disp_tmp.x;
    disp[1] = disp_tmp.y;
    disp[2] = disp_tmp.z;
    return disp;
}
template <class Element>
const Dune::FieldVector<double, 6>
TDStrainFS(const Dune::FieldVector<double, 3>& obs,
           const Element& elem,
           const Dune::FieldVector<double, 3>& slip,
           double nu)
{
    std::array<Real3, 3> tri = getTri(elem);
    Real3 obs_tmp = make3(obs[0], obs[1], obs[2]);
    Real3 slip_tmp = make3(slip[0], slip[1], slip[2]);
    Real6 strain_tmp = strain_fs(obs_tmp, tri, slip_tmp, nu);
    Dune::FieldVector<double, 6> strain;
    strain[0] = strain_tmp.x;
    strain[1] = strain_tmp.y;
    strain[2] = strain_tmp.z;
    strain[3] = strain_tmp.a;
    strain[4] = strain_tmp.b;
    strain[5] = strain_tmp.c;
    return strain;
}
void
//assembleMatrix(Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)
assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)

{
    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
    for (auto elem1 : elements(grid.leafGridView())) {
        int idx1 = mapper.index(elem1);
        for (auto elem2 : elements(grid.leafGridView())) {
            int idx2 = mapper.index(elem2);
            auto geom = elem1.geometry();
            auto center = geom.center();
            // check if this is defined in relative coordinates
            Dune::FieldVector<double,3> slip;// = make3(1.0,0.0, 0.0);
            slip[0]=1;slip[1]=1;slip[2]=0;
            // symmetric stress voit notation
            Dune::FieldVector<double,6> strain = TDStrainFS(center, elem2, slip, nu);
            /*
            Dune::FieldVector<double,6> stress = strainToStress(E,nu,strain);
            double ntraction = normalTraction(stress,normal);
            */
            double ntraction = strain[0];// strain[0][0]
            //  auto traction
            matrix[idx1][idx2] = ntraction;
        }
    }
}

} // namespace ddm
