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

#include "config.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <opm/common/TimingMacros.hpp>
#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/FractureMeshFoamGrid.hpp>

#include <cmath>
#include <limits>
#include <utility>
namespace ddm
{
namespace
{
struct GeometryCache
{
    std::vector<Dune::FieldVector<double, 3>> centers;
    std::vector<Dune::FieldVector<double, 3>> normals;
    std::vector<std::array<Real3, 3>> tris;
};

std::array<Real3, 3> toReal3(const Opm::FractureMesh::Triangle& t)
{
    std::array<Real3, 3> tri;
    for (int i = 0; i < 3; ++i) {
        tri[i] = make3(t[i][0], t[i][1], t[i][2]);
    }
    return tri;
}

GeometryCache buildGeometryCache(const Opm::FractureMesh& mesh)
{
    const int nc = mesh.numCells();

    GeometryCache cache;
    cache.centers.resize(nc);
    cache.normals.resize(nc);
    cache.tris.resize(nc);

    OPM_TIMEBLOCK(GeometryOFfoamgrid);
    for (int idx = 0; idx < nc; ++idx) {
        cache.centers[idx] = mesh.cellCentroid(idx);
        cache.normals[idx] = mesh.cellNormal(idx);
        cache.tris[idx] = toReal3(mesh.cellTriangles(idx)[0]);
    }

    return cache;
}

double kernelNormalTraction(const Dune::FieldVector<double, 3>& center,
                            const Dune::FieldVector<double, 3>& normal,
                            const std::array<Real3, 3>& tri,
                            const double slip_value,
                            const double E,
                            const double nu)
{
    Dune::FieldVector<double, 3> slip;
    // The CutDe wrapper uses TDCS ordering with the normal-opening component first.
    slip[0] = slip_value;
    slip[1] = 0.0;
    slip[2] = 0.0;

    const Dune::FieldVector<double, 6> strain = TDStrainFSTri(center, tri, slip, nu);
    const Dune::FieldVector<double, 6> stress = strainToStress(E, nu, strain);
    return tractionSymTensor(stress, normal);
}
} // namespace

double
fractureK1(double dist, double width, double E, double nu)
{
    double K1;
    double mu = E / (2 * (1 + nu)); //??
    K1 = (mu * sqrt(M_PI) / (2 * std::sqrt(dist) * (1.0 - nu))) * width;
    K1 /= 1.834; // calibrated numerically; with the exact near-tip opening the factor is sqrt(2)
    return K1;
}
double
traceSymTensor(const Dune::FieldVector<double, 6> symtensor)
{
    double trace = 0;
    for (int i = 0; i < 3; ++i) {
        trace += symtensor[i];
    }
    return trace;
}

Dune::FieldVector<double, 6>
strainToStress(const double E, const double nu, const Dune::FieldVector<double, 6> strain)
{
    Dune::FieldVector<double, 6> stress;
    double volume_strain = traceSymTensor(strain);
    double mu = E / (2 * (1 + nu)); //??
    double lambda = 2 * mu * nu / (1 - 2 * nu);
    for (int i = 0; i < 3; ++i) {
        stress[i] = 2 * mu * strain[i] + (lambda * volume_strain);
        stress[i + 3] += 2 * mu * strain[i + 3]; // [xy,xz,yz]
    }
    return stress;
}

Dune::FieldMatrix<double, 3, 3>
symTensor2Matrix(const Dune::FieldVector<double, 6> symtensor)
{
    Dune::FieldMatrix<double, 3, 3> mat;
    for (int i = 0; i < 3; ++i) {
        mat[i][i] = symtensor[i];
    }
    mat[2][1] = symtensor[3];
    mat[1][2] = mat[2][1];
    mat[2][0] = symtensor[4];
    mat[0][2] = mat[2][0];
    mat[1][0] = symtensor[5];
    mat[0][1] = mat[1][0];
    return mat;
}


double
tractionSymTensor(const Dune::FieldVector<double, 6> symtensor, Dune::FieldVector<double, 3> normal)
{
    assert(std::abs(normal.two_norm() - 1) < 1e-13);
    double traction = 0.0;
    for (int i = 0; i < 3; ++i) {
        traction += symtensor[i] * normal[i] * normal[i];
    }
    traction += 2 * symtensor[3] * normal[0] * normal[1]; // xy*nx*ny;
    traction += 2 * symtensor[4] * normal[0] * normal[2]; // xz*nx*nz
    traction += 2 * symtensor[5] * normal[1] * normal[2]; // yz*ny*nz
    return traction;
}

double
tractionVoigt(Dune::FieldVector<double, 6> voigt, Dune::FieldVector<double, 3> normal)
{
    std::swap(voigt[3], voigt[5]);
    return tractionSymTensor(voigt, normal);
}

const Dune::FieldVector<double, 6>
TDStrainFSTri(const Dune::FieldVector<double, 3>& obs,
              const std::array<Real3, 3>& tri,
              const Dune::FieldVector<double, 3>& slip,
              double nu)
{
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
Dune::FieldVector<double, 3>
TDDispFSTri(const Dune::FieldVector<double, 3>& obs,
            const std::array<Real3, 3>& tri,
            const Dune::FieldVector<double, 3>& slip,
            double nu)
{
    Real3 obs_tmp = make3(obs[0], obs[1], obs[2]);
    Real3 slip_tmp = make3(slip[0], slip[1], slip[2]);
    Real3 disp_tmp = disp_fs(obs_tmp, tri, slip_tmp, nu);
    Dune::FieldVector<double, 3> disp;
    disp[0] = disp_tmp.x;
    disp[1] = disp_tmp.y;
    disp[2] = disp_tmp.z;
    return disp;
}

void
assembleMatrix(Dune::DynamicMatrix<double>& matrix,
               const double E,
               const double nu,
               const Opm::FractureMesh& mesh)

{
    OPM_TIMEFUNCTION();
    const int nc = mesh.numCells();
    for (int idx1 = 0; idx1 < nc; ++idx1) {
        auto center = mesh.cellCentroid(idx1);
        auto normal = mesh.cellNormal(idx1);
        for (int idx2 = 0; idx2 < nc; ++idx2) {
            // check if this is defined in relative coordinates
            Dune::FieldVector<double, 3> slip; // = make3(1.0,0.0, 0.0);
            slip[0] = 1;
            slip[1] = 0;
            slip[2] = 0;

            // symmetric stress voit notation
            Dune::FieldVector<double, 6> strain
                = TDStrainFSTri(center, toReal3(mesh.cellTriangles(idx2)[0]), slip, nu);
            Dune::FieldVector<double, 6> stress = strainToStress(E, nu, strain);

            double ntraction = tractionSymTensor(stress, normal);
            // matrix relate to pure traction not area weighted
            matrix[idx1][idx2] = ntraction;
        }
    }
}

void
assembleMatrix_fast(Dune::DynamicMatrix<double>& matrix,
                    const double E,
                    const double nu,
                    const Opm::FractureMesh& mesh)

{
    OPM_TIMEFUNCTION();
    const int nc = mesh.numCells();
    const auto cache = buildGeometryCache(mesh);
    {
        OPM_TIMEBLOCK(AssembleDDMMatrix);
#pragma omp parallel for
        for (int idx1 = 0; idx1 < nc; ++idx1) {
            for (int idx2 = 0; idx2 < nc; ++idx2) {
                matrix[idx1][idx2] = kernelNormalTraction(cache.centers[idx1],
                                                          cache.normals[idx1],
                                                          cache.tris[idx2],
                                                          1.0,
                                                          E,
                                                          nu);
            }
        }
    }
}

void
assembleMatrix_lattice(Dune::DynamicMatrix<double>& matrix,
                       const double E,
                       const double nu,
                       const Opm::FractureMesh& mesh,
                       const LatticeFrame& frame,
                       KernelTable& table)
{
    OPM_TIMEFUNCTION();
    using V = Dune::FieldVector<double, 3>;
    const int nc = mesh.numCells();
    const auto cache = buildGeometryCache(mesh);

    const auto sameFrame = [](const LatticeFrame& f, const LatticeFrame& g) {
        return f.origin == g.origin && f.e1 == g.e1 && f.e2 == g.e2;
    };
    if (table.E != E || table.nu != nu || !sameFrame(table.frame, frame)) {
        table = KernelTable {};
        table.E = E;
        table.nu = nu;
        table.frame = frame;
    }

    const V& e1 = frame.e1;
    const V& e2 = frame.e2;
    V normal {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};
    normal /= normal.two_norm();
    const double g11 = e1 * e1, g12 = e1 * e2, g22 = e2 * e2;
    const double det = g11 * g22 - g12 * g12;
    const double h = std::sqrt(std::min(g11, g22));
    constexpr double tol = 1e-6; // in lattice units

    // Shape index (-1: off lattice) and three times the centroid, in lattice coordinates.
    std::vector<int> shape(nc, -1);
    std::vector<std::array<int, 2>> c3(nc);
    for (int i = 0; i < nc; ++i) {
        if (std::abs(cache.normals[i] * normal) < 1.0 - 1e-9) {
            continue;
        }
        std::array<std::array<int, 2>, 3> node;
        bool onLattice = true;
        for (int m = 0; m < 3 && onLattice; ++m) {
            V p {cache.tris[i][m].x, cache.tris[i][m].y, cache.tris[i][m].z};
            p -= frame.origin;
            const double x = p * e1, y = p * e2;
            const double a = (g22 * x - g12 * y) / det;
            const double b = (g11 * y - g12 * x) / det;
            onLattice = std::abs(a - std::round(a)) < tol && std::abs(b - std::round(b)) < tol
                && std::abs(p * normal) < tol * h && std::abs(a) < 1e8 && std::abs(b) < 1e8;
            node[m] = {static_cast<int>(std::lround(a)), static_cast<int>(std::lround(b))};
        }
        if (!onLattice) {
            continue;
        }
        const std::array<int, 4> key {node[1][0] - node[0][0], node[1][1] - node[0][1],
                                      node[2][0] - node[0][0], node[2][1] - node[0][1]};
        const auto [it, inserted] = table.shapes.emplace(key, static_cast<int>(table.boxes.size()));
        if (inserted) {
            table.boxes.emplace_back();
        }
        shape[i] = it->second;
        c3[i] = {node[0][0] + node[1][0] + node[2][0], node[0][1] + node[1][1] + node[2][1]};
    }

    // Grow each shape's box to cover every centroid offset of this grid.
    const int nshape = static_cast<int>(table.boxes.size());
    constexpr int big = std::numeric_limits<int>::max();
    std::array<int, 2> imin {big, big}, imax {-big, -big};
    std::vector<std::array<int, 4>> jrange(nshape, {big, big, -big, -big});
    for (int i = 0; i < nc; ++i) {
        if (shape[i] < 0) {
            continue;
        }
        auto& r = jrange[shape[i]];
        for (int d = 0; d < 2; ++d) {
            imin[d] = std::min(imin[d], c3[i][d]);
            imax[d] = std::max(imax[d], c3[i][d]);
            r[d] = std::min(r[d], c3[i][d]);
            r[2 + d] = std::max(r[2 + d], c3[i][d]);
        }
    }
    long total = 0;
    for (int k = 0; k < nshape; ++k) {
        const auto& r = jrange[k];
        if (r[0] > r[2]) {
            total += static_cast<long>(table.boxes[k].values.size());
            continue;
        }
        auto& box = table.boxes[k];
        const int lo0 = r[0] - imax[0], hi0 = r[2] - imin[0];
        const int lo1 = r[1] - imax[1], hi1 = r[3] - imin[1];
        const bool covered = !box.values.empty() && lo0 >= box.d0min && lo1 >= box.d1min
            && hi0 < box.d0min + box.n0 && hi1 < box.d1min + box.n1;
        if (!covered) {
            // Margin so a growing fracture does not reallocate on every regrid.
            const int m0 = (hi0 - lo0) / 4 + 6, m1 = (hi1 - lo1) / 4 + 6;
            KernelTable::Box grown;
            grown.d0min = (box.values.empty() ? lo0 : std::min(lo0, box.d0min)) - m0;
            grown.d1min = (box.values.empty() ? lo1 : std::min(lo1, box.d1min)) - m1;
            const int top0 = (box.values.empty() ? hi0 : std::max(hi0, box.d0min + box.n0 - 1)) + m0;
            const int top1 = (box.values.empty() ? hi1 : std::max(hi1, box.d1min + box.n1 - 1)) + m1;
            grown.n0 = top0 - grown.d0min + 1;
            grown.n1 = top1 - grown.d1min + 1;
            grown.values.assign(static_cast<std::size_t>(grown.n0) * grown.n1,
                                std::numeric_limits<double>::quiet_NaN());
            for (int d0 = 0; d0 < box.n0; ++d0) {
                for (int d1 = 0; d1 < box.n1; ++d1) {
                    grown.values[static_cast<std::size_t>(d0 + box.d0min - grown.d0min) * grown.n1
                                 + (d1 + box.d1min - grown.d1min)]
                        = box.values[static_cast<std::size_t>(d0) * box.n1 + d1];
                }
            }
            box = std::move(grown);
        }
        total += static_cast<long>(box.values.size());
    }
    if (total > 50'000'000) { // 400 MB: rather compute directly than hold such a table
        table = KernelTable {};
        assembleMatrix_fast(matrix, E, nu, mesh);
        return;
    }

    const auto entry = [&table](const int k, const std::array<int, 2>& ci, const std::array<int, 2>& cj) {
        const auto& box = table.boxes[k];
        return static_cast<std::size_t>(cj[0] - ci[0] - box.d0min) * box.n1 + (cj[1] - ci[1] - box.d1min);
    };

    std::vector<std::vector<unsigned char>> need(nshape);
    for (int k = 0; k < nshape; ++k) {
        need[k].assign(table.boxes[k].values.size(), 0);
    }
#pragma omp parallel for
    for (int i = 0; i < nc; ++i) {
        if (shape[i] < 0) {
            continue;
        }
        for (int j = 0; j < nc; ++j) {
            const int k = shape[j];
            if (k >= 0) {
                const auto e = entry(k, c3[i], c3[j]);
                if (std::isnan(table.boxes[k].values[e])) {
#pragma omp atomic write
                    need[k][e] = 1;
                }
            }
        }
    }
    std::vector<std::pair<int, std::size_t>> missing;
    for (int k = 0; k < nshape; ++k) {
        for (std::size_t e = 0; e < need[k].size(); ++e) {
            if (need[k][e]) {
                missing.emplace_back(k, e);
            }
        }
    }
    std::vector<std::array<int, 4>> offsets(nshape);
    for (const auto& [key, k] : table.shapes) {
        offsets[k] = key;
    }
    {
        OPM_TIMEBLOCK(FillDDMTable);
        const V obs(0.0);
#pragma omp parallel for
        for (std::size_t q = 0; q < missing.size(); ++q) {
            const auto [k, e] = missing[q];
            auto& box = table.boxes[k];
            const int d0 = static_cast<int>(e / box.n1) + box.d0min;
            const int d1 = static_cast<int>(e % box.n1) + box.d1min;
            const auto& o = offsets[k];
            // Corner m of element j relative to the centroid of element i, times three.
            const std::array<std::array<int, 2>, 3> rel {{{d0 - o[0] - o[2], d1 - o[1] - o[3]},
                                                          {d0 + 2 * o[0] - o[2], d1 + 2 * o[1] - o[3]},
                                                          {d0 - o[0] + 2 * o[2], d1 - o[1] + 2 * o[3]}}};
            std::array<Real3, 3> tri;
            for (int m = 0; m < 3; ++m) {
                V p = e1;
                p *= rel[m][0] / 3.0;
                p.axpy(rel[m][1] / 3.0, e2);
                tri[m] = make3(p[0], p[1], p[2]);
            }
            box.values[e] = kernelNormalTraction(obs, normal, tri, 1.0, E, nu);
        }
    }
    {
        OPM_TIMEBLOCK(AssembleDDMMatrix);
#pragma omp parallel for
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nc; ++j) {
                const int k = shape[j];
                matrix[i][j] = (shape[i] >= 0 && k >= 0)
                    ? table.boxes[k].values[entry(k, c3[i], c3[j])]
                    : kernelNormalTraction(cache.centers[i], cache.normals[i], cache.tris[j], 1.0, E, nu);
            }
        }
    }
}

void
applyMatrix_fast(Dune::DynamicVector<double>& result,
                 const Dune::DynamicVector<double>& x,
                 const double E,
                 const double nu,
                 const Opm::FractureMesh& mesh)
{
    OPM_TIMEFUNCTION();
    const int nc = mesh.numCells();
    assert(static_cast<int>(x.size()) == nc);

    result.resize(nc);
    const auto cache = buildGeometryCache(mesh);

    OPM_TIMEBLOCK(ApplyDDMMatrix);
#pragma omp parallel for
    for (int idx1 = 0; idx1 < nc; ++idx1) {
        double value = 0.0;
        for (int idx2 = 0; idx2 < nc; ++idx2) {
            value += kernelNormalTraction(cache.centers[idx1],
                                          cache.normals[idx1],
                                          cache.tris[idx2],
                                          x[idx2],
                                          E,
                                          nu);
        }
        result[idx1] = value;
    }
}




Dune::FieldVector<double, 6>
strain(const Dune::FieldVector<double, 3>& obs,
       const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
       const Opm::FractureMesh& mesh,
       const double /*E*/,
       const double nu)
{
    Dune::FieldVector<double, 6> strain;
    strain = 0.0;
    for (int idx1 = 0; idx1 < mesh.numCells(); ++idx1) {
        strain += TDStrainFSTri(obs, toReal3(mesh.cellTriangles(idx1)[0]), slips[idx1], nu);
    }
    return strain;
}

Dune::FieldVector<double, 3>
disp(const Dune::FieldVector<double, 3>& obs,
     const Dune::BlockVector<Dune::FieldVector<double, 3>>& slips,
     const Opm::FractureMesh& mesh,
     const double /*E*/,
     const double nu)
{
    Dune::FieldVector<double, 3> disp;
    disp = 0.0;
    for (int idx1 = 0; idx1 < mesh.numCells(); ++idx1) {
        disp += TDDispFSTri(obs, toReal3(mesh.cellTriangles(idx1)[0]), slips[idx1], nu);
    }
    return disp;
}

void
assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)
{
    assembleMatrix(matrix, E, nu, Opm::buildFractureMesh(grid));
}

void
assembleMatrix_fast(Dune::DynamicMatrix<double>& matrix,
                    const double E,
                    const double nu,
                    const Dune::FoamGrid<2, 3>& grid)
{
    assembleMatrix_fast(matrix, E, nu, Opm::buildFractureMesh(grid));
}

void
assembleMatrix_lattice(Dune::DynamicMatrix<double>& matrix,
                       const double E,
                       const double nu,
                       const Dune::FoamGrid<2, 3>& grid,
                       const LatticeFrame& frame,
                       KernelTable& table)
{
    assembleMatrix_lattice(matrix, E, nu, Opm::buildFractureMesh(grid), frame, table);
}

void
applyMatrix_fast(Dune::DynamicVector<double>& result,
                 const Dune::DynamicVector<double>& x,
                 const double E,
                 const double nu,
                 const Dune::FoamGrid<2, 3>& grid)
{
    applyMatrix_fast(result, x, E, nu, Opm::buildFractureMesh(grid));
}
} // namespace ddm
