#pragma once
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include "CutDe.hpp"
#include <opm/geomech/Math.hpp>
namespace ddm
{
template <class Element>
Dune::FieldVector<double, 3>
normalOfElement(const Element& elem){
    // dim = Dune::Codim<Grid::dimension>{}
    int codim = 2;//vertices in this cases
    std::array<Dune::FieldVector<double, 3>, 3> tri;
    int i = 0;
    for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<2>{})){
        // assume order is ok and orientation is fine;
        auto geom = vertex.geometry();
        auto corner  = geom.center();
        tri[i] = corner;
        ++i;
    }
    // NB sign should not matter as long as tensile fracture
    Dune::FieldVector<double, 3> e01 = tri[1]-tri[0];
    Dune::FieldVector<double, 3> e02 = tri[2]-tri[0];
    Dune::FieldVector<double, 3> normal = Opm::crossProduct(e01, e02);
    normal /= normal.two_norm();
    return normal;
}

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

inline double traceSymTensor(const Dune::FieldVector<double,6> symtensor){
  double trace = 0;
  for(int i=0; i < 3; ++i){
    trace += symtensor[i];
  }
  return trace;
}

// compute the stress tensor (6 components) from the strain tensor (6 components)
inline Dune::FieldVector<double,6>
strainToStress(const double E,const double nu,const Dune::FieldVector<double,6> strain){
    Dune::FieldVector<double,6> stress;
    double volume_strain = traceSymTensor(strain);
    double mu = E/(2*(1+nu));//??
    double lambda = 2 * mu * nu / (1 - 2 * nu);
    for(int i=0; i < 3; ++i){
        stress[i]   =  2*mu*strain[i] + (lambda * volume_strain);
        stress[i+3] += 2*mu*strain[i+3]; // [xy,xz,yz]
    }
    return stress;
}

inline double
tractionSymTensor(const Dune::FieldVector<double,6> symtensor, Dune::FieldVector<double,3> normal){
    const Dune::FieldVector<double,6> stress;
    assert(std::abs(normal.two_norm()-1) < 1e-13);
    double traction=0.0;
    for(int i=0; i < 3; ++i){
        traction += symtensor[i]*normal[i]*normal[i];
    }
    traction += 2*symtensor[3]*normal[0]*normal[1];// xy*nx*ny;
    traction += 2*symtensor[4]*normal[0]*normal[2];// xz*nx*nz
    traction += 2*symtensor[6]*normal[1]*normal[2];// yz*ny*nz
    return traction;
}

inline void
//assembleMatrix(Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)
assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)

{
    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
    for (auto elem1 : elements(grid.leafGridView())) {
        int idx1 = mapper.index(elem1);
        auto geom = elem1.geometry();
        auto center = geom.center();
        auto normal = normalOfElement(elem1);
        for (auto elem2 : elements(grid.leafGridView())) {
            int idx2 = mapper.index(elem2);
            // check if this is defined in relative coordinates
            Dune::FieldVector<double,3> slip;// = make3(1.0,0.0, 0.0);
            slip[0]=1;slip[1]=0;slip[2]=0;

            // symmetric stress voit notation
            Dune::FieldVector<double,6> strain = TDStrainFS(center, elem2, slip, nu);
            Dune::FieldVector<double,6> stress = strainToStress(E,nu,strain);

            double ntraction = tractionSymTensor(stress,normal);
            // matrix relate to pure traction not area weighted
            matrix[idx1][idx2] = ntraction;
        }
    }
}

} // namespace ddm


