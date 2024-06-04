#pragma once
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include "CutDe.hpp"
#include <opm/geomech/Math.hpp>
#include <dune/grid/common/entity.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
namespace ddm
{
template <class Element>
Dune::FieldVector<double, 3>
normalOfElement(const Element& elem){
    // dim = Dune::Codim<Grid::dimension>{}
    //int codim = 2;//vertices in this cases
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

double fractureK1(double dist,double width, double E, double nu);

//double fractureK2(double dist,double edist, double E, double nu);
//double fractureK3(double dist,double edist, double E, double nu);

template <class Element>
std::array<Real3, 3>
getTri(const Element& elem){
    // dim = Dune::Codim<Grid::dimension>{}
    //int codim = 2;//vertices in this cases
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

double traceSymTensor(const Dune::FieldVector<double,6> symtensor);

Dune::FieldMatrix<double,3,3> symTensor2Matrix(const Dune::FieldVector<double,6> symtensor);

// compute the stress tensor (6 components) from the strain tensor (6 components)
Dune::FieldVector<double,6>
strainToStress(const double E,const double nu,const Dune::FieldVector<double,6> strain);

double
tractionSymTensor(const Dune::FieldVector<double,6> symtensor, Dune::FieldVector<double,3> normal);

//assembleMatrix(Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)
void assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid);

} // namespace ddm


