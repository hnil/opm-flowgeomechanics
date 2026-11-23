/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <array>
#include <map>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include "CutDe.hpp"
#include <opm/geomech/FractureMesh.hpp>
#include <opm/geomech/Math.hpp>
#include <dune/grid/common/entity.hh>
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif
#include <dune/grid/common/mcmgmapper.hh>
namespace ddm
{
template <class Element>
Dune::FieldVector<double, 3>
normalOfElement(const Element& elem){
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


template <class Element>
std::array<Real3, 3>
getTri(const Element& elem){
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
Dune::FieldVector<double, 3>
TDDispFSTri(const Dune::FieldVector<double, 3>& obs,
            const std::array<Real3, 3>& tri,
            const Dune::FieldVector<double, 3>& slip,
            double nu);

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

const Dune::FieldVector<double, 6>
TDStrainFSTri(const Dune::FieldVector<double, 3>& obs,
              const std::array<Real3, 3>& tri,
              const Dune::FieldVector<double, 3>& slip,
              double nu);

template <class Element>
const Dune::FieldVector<double, 6>
TDStrainFS(const Dune::FieldVector<double, 3>& obs,
           const Element& elem,
           const Dune::FieldVector<double, 3>& slip,
           double nu)
{
    std::array<Real3, 3> tri = getTri(elem);
    return TDStrainFSTri(obs, tri,slip,nu);
}



double traceSymTensor(const Dune::FieldVector<double,6> symtensor);

Dune::FieldMatrix<double,3,3> symTensor2Matrix(const Dune::FieldVector<double,6> symtensor);

// compute the stress tensor (6 components) from the strain tensor (6 components)
Dune::FieldVector<double,6>
strainToStress(const double E,const double nu,const Dune::FieldVector<double,6> strain);

// Normal traction of a tensor ordered [xx, yy, zz, xy, xz, yz].
double
tractionSymTensor(const Dune::FieldVector<double,6> symtensor, Dune::FieldVector<double,3> normal);

// Normal traction of a tensor in Voigt order [xx, yy, zz, yz, xz, xy].
double
tractionVoigt(Dune::FieldVector<double,6> voigt, Dune::FieldVector<double,3> normal);

void assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Opm::FractureMesh& mesh);
void assembleMatrix_fast(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Opm::FractureMesh& mesh);

// Planar lattice the fracture mesh lives on: node (i, j) sits at origin + i*e1 + j*e2.
struct LatticeFrame
{
    Dune::FieldVector<double, 3> origin;
    Dune::FieldVector<double, 3> e1;
    Dune::FieldVector<double, 3> e2;
};

// Kernel values of lattice element pairs, keyed by element shape and centroid offset in lattice units.
// Kept across regrids; cleared when E, nu or the frame change.
struct KernelTable
{
    struct Box
    {
        int d0min = 0, d1min = 0, n0 = 0, n1 = 0;
        std::vector<double> values; // NaN marks an entry not yet computed
    };
    double E = 0.0;
    double nu = 0.0;
    LatticeFrame frame;
    std::map<std::array<int, 4>, int> shapes; // corner offsets from corner 0 -> shape index
    std::vector<Box> boxes; // one per shape
};

// Same matrix as assembleMatrix_fast up to round-off; pairs of off-lattice elements use the direct kernel.
void assembleMatrix_lattice(Dune::DynamicMatrix<double>& matrix, const double E, const double nu,
                            const Opm::FractureMesh& mesh, const LatticeFrame& frame, KernelTable& table);
void applyMatrix_fast(Dune::DynamicVector<double>& result,
                      const Dune::DynamicVector<double>& x,
                      const double E,
                      const double nu,
                      const Opm::FractureMesh& mesh);

Dune::FieldVector<double, 6>
strain(const Dune::FieldVector<double, 3>& obs,
       const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
       const Opm::FractureMesh& mesh,
       const double E,
       const double nu);

Dune::FieldVector<double, 3>
disp(const Dune::FieldVector<double, 3>& obs,
       const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
       const Opm::FractureMesh& mesh,
       const double E,
       const double nu);

#if HAVE_DUNE_FOAMGRID
// FoamGrid forms for tools and tests: build the mesh and forward.
void assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid);
void assembleMatrix_fast(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid);
void assembleMatrix_lattice(Dune::DynamicMatrix<double>& matrix, const double E, const double nu,
                            const Dune::FoamGrid<2, 3>& grid, const LatticeFrame& frame, KernelTable& table);
void applyMatrix_fast(Dune::DynamicVector<double>& result,
                      const Dune::DynamicVector<double>& x,
                      const double E,
                      const double nu,
                      const Dune::FoamGrid<2, 3>& grid);
#endif

} // namespace ddm
