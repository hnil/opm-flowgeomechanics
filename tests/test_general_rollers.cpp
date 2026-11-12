/*
  Copyright 2026 Equinor ASA.

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

#include <config.h>
#include <dune/istl/schwarz.hh>
#include <opm/geomech/VemElasticitySolver.hpp>
#include <opm/geomech/BoundaryUtils.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Vec3 = std::array<double, 3>;
using Opm::Elasticity::NodeConstraint;

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

constexpr int nx = 4, ny = 3, nz = 3;

// Corner-point grid with flat layers; pillar (i, j) at plan position xy(i, j).
Dune::CpGrid makeGrid(const std::function<std::array<double, 2>(int, int)>& xy,
                      std::vector<double>& coord, std::vector<double>& zcorn, std::vector<int>& actnum)
{
    coord.clear();
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const auto p = xy(i, j);
            coord.insert(coord.end(), {p[0], p[1], 0.0, p[0], p[1], 10.0});
        }
    }
    zcorn.assign(8 * nx * ny * nz, 0.0);
    for (int k = 0; k < 2 * nz; ++k) {
        for (int j = 0; j < 2 * ny; ++j) {
            for (int i = 0; i < 2 * nx; ++i) {
                zcorn[k * 4 * nx * ny + j * 2 * nx + i] = 0.5 * ((k + 1) / 2);
            }
        }
    }
    actnum.assign(nx * ny * nz, 1);
    grdecl raw;
    raw.dims[0] = nx;
    raw.dims[1] = ny;
    raw.dims[2] = nz;
    raw.coord = coord.data();
    raw.zcorn = zcorn.data();
    raw.actnum = actnum.data();
    Dune::CpGrid grid;
    grid.processEclipseFormat(raw, false);
    return grid;
}

std::vector<Vec3> vertexPositions(const Dune::CpGrid& grid)
{
    const auto gv = grid.leafGridView();
    std::vector<Vec3> pos(gv.size(3));
    for (const auto& v : vertices(gv)) {
        const auto c = v.geometry().center();
        pos[gv.indexSet().index(v)] = {c[0], c[1], c[2]};
    }
    return pos;
}

// Global node displacements for the cell pressure `load`.
std::vector<Vec3> solve(const Dune::CpGrid& grid, const std::vector<NodeConstraint>& constraints,
                        const std::vector<double>& load)
{
    using namespace Opm::Elasticity;
    const int nc = grid.leafGridView().size(0);
    VemDiscretisation<Dune::CpGrid> disc(grid);
    disc.setMaterial(std::vector<double>(nc, 1e9), std::vector<double>(nc, 0.25));
    disc.setOptions(1, false, false);
    disc.setBodyForce(0.0, std::vector<double>(nc, 0.0));
    disc.fixNodes(constraints);
    disc.assembleMatrix();
    disc.buildOperators();
    Vector p(nc);
    for (int c = 0; c < nc; ++c) {
        p[c] = load[c];
    }
    Vector b;
    disc.rhs(p, b);
    const auto& A = disc.matrix();
    const int n = A.N();
    Dune::DynamicMatrix<double> dense(n, n, 0.0);
    for (auto row = A.begin(); row != A.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            dense[row.index()][col.index()] = (*col)[0][0];
        }
    }
    Dune::DynamicVector<double> x(n), rhs(n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = b[i];
    }
    dense.solve(x, rhs);
    Vector u(n);
    for (int i = 0; i < n; ++i) {
        u[i] = x[i];
    }
    disc.toGlobal(u);
    std::vector<Vec3> result(n / 3);
    for (int k = 0; k < n / 3; ++k) {
        result[k] = {u[3 * k], u[3 * k + 1], u[3 * k + 2]};
    }
    return result;
}

double norm(const Vec3& a)
{
    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// Rollers on all sides of a box rotated about z: the solution is the rotated
// solution of the axis-aligned box.
void rotatedBox()
{
    const double angle = M_PI / 6, c = std::cos(angle), s = std::sin(angle);
    const double dx = 1.0, dy = 1.5;
    std::vector<double> coord0, zcorn0, coord1, zcorn1;
    std::vector<int> act0, act1;
    const auto grid0 = makeGrid([&](int i, int j) { return std::array<double, 2>{dx * i, dy * j}; },
                                coord0, zcorn0, act0);
    const auto grid1 = makeGrid([&](int i, int j) {
        return std::array<double, 2>{c * dx * i - s * dy * j, s * dx * i + c * dy * j};
    }, coord1, zcorn1, act1);

    using Opm::Elasticity::OuterBoundary;
    std::vector<NodeConstraint> con0, con1;
    Opm::Elasticity::nodesAtOuterBoundary(con0, grid0, OuterBoundary::Roller);
    Opm::Elasticity::nodesAtOuterBoundary(con1, grid1, OuterBoundary::Roller);
    for (const auto& con : con0) {
        require(con.axisAligned(), "An axis-aligned box must not rotate any node");
    }
    int rotated = 0;
    for (const auto& con : con1) {
        rotated += !con.axisAligned();
    }
    require(rotated > 0, "A rotated box must rotate its side nodes");

    // pressure linear along the box's own x axis
    const int nc = grid0.leafGridView().size(0);
    std::vector<double> load0(nc), load1(nc);
    for (int e = 0; e < nc; ++e) {
        load0[e] = 1e6 * (1.0 + grid0.cellCentroid(e)[0]);
        const auto x1 = grid1.cellCentroid(e);
        load1[e] = 1e6 * (1.0 + c * x1[0] + s * x1[1]);
    }
    const auto u0 = solve(grid0, con0, load0);
    const auto u1 = solve(grid1, con1, load1);
    const auto pos0 = vertexPositions(grid0);
    const auto pos1 = vertexPositions(grid1);
    double umax = 0.0, err = 0.0;
    for (std::size_t k = 0; k < u0.size(); ++k) {
        umax = std::max(umax, norm(u0[k]));
    }
    for (std::size_t m = 0; m < u1.size(); ++m) {
        const Vec3 back{c * pos1[m][0] + s * pos1[m][1], -s * pos1[m][0] + c * pos1[m][1], pos1[m][2]};
        std::size_t k = 0;
        while (k < pos0.size() && norm({pos0[k][0] - back[0], pos0[k][1] - back[1],
                                        pos0[k][2] - back[2]}) > 1e-9) {
            ++k;
        }
        require(k < pos0.size(), "Rotated vertex has no axis-aligned twin");
        const Vec3 expect{c * u0[k][0] - s * u0[k][1], s * u0[k][0] + c * u0[k][1], u0[k][2]};
        err = std::max(err, norm({u1[m][0] - expect[0], u1[m][1] - expect[1], u1[m][2] - expect[2]}));
    }
    std::cout << "rotated box: max |u| " << umax << ", max error " << err << '\n';
    require(umax > 0.0 && err < 1e-9 * umax, "Rotated box differs from the rotated axis-aligned solution");
}

// Rollers on two vertical planes through the z axis (one at 30 degrees) and on
// z = 0: uniform pressure gives the exact isotropic expansion u = eps x.
double wedgeError(bool snapToAxis)
{
    const double sector = M_PI / 6;
    std::vector<double> coord, zcorn;
    std::vector<int> act;
    const auto grid = makeGrid([&](int i, int j) {
        const double r = 2.0 + i, t = sector * j / ny;
        return std::array<double, 2>{r * std::cos(t), r * std::sin(t)};
    }, coord, zcorn, act);
    const auto gv = grid.leafGridView();
    const Vec3 tilted{-std::sin(sector), std::cos(sector), 0.0};
    std::vector<std::vector<Vec3>> normals(gv.size(3));
    for (const auto& cell : elements(gv)) {
        const int e = gv.indexSet().index(cell);
        for (int f = 0; f < grid.numCellFaces(e); ++f) {
            const int face = grid.cellFace(e, f);
            if (grid.faceCell(face, 0) >= 0 && grid.faceCell(face, 1) >= 0) {
                continue;
            }
            const auto x = grid.faceCentroid(face);
            const bool onRoller = std::abs(x[1]) < 1e-9 || std::abs(x[2]) < 1e-9
                || std::abs(tilted[0] * x[0] + tilted[1] * x[1]) < 1e-9;
            if (!onRoller) {
                continue;
            }
            auto n = grid.faceAreaNormalEcl(face);
            Vec3 normal{n[0], n[1], n[2]};
            if (snapToAxis) {
                int axis = 0;
                for (int d = 1; d < 3; ++d) {
                    if (std::abs(normal[d]) > std::abs(normal[axis])) {
                        axis = d;
                    }
                }
                const double area = norm(normal);
                normal = {0.0, 0.0, 0.0};
                normal[axis] = area;
            }
            for (int v = 0; v < grid.numFaceVertices(face); ++v) {
                normals[grid.faceVertex(face, v)].push_back(normal);
            }
        }
    }
    std::vector<NodeConstraint> constraints;
    for (std::size_t k = 0; k < normals.size(); ++k) {
        if (!normals[k].empty()) {
            constraints.push_back(Opm::Elasticity::rollerConstraint(k, normals[k]));
        }
    }
    const auto u = solve(grid, constraints, std::vector<double>(gv.size(0), 1e6));
    const auto pos = vertexPositions(grid);
    double num = 0.0, den = 0.0;
    for (std::size_t k = 0; k < u.size(); ++k) {
        for (int d = 0; d < 3; ++d) {
            num += u[k][d] * pos[k][d];
            den += pos[k][d] * pos[k][d];
        }
    }
    const double eps = num / den;
    double err = 0.0, umax = 0.0;
    for (std::size_t k = 0; k < u.size(); ++k) {
        umax = std::max(umax, norm(u[k]));
        err = std::max(err, norm({u[k][0] - eps * pos[k][0], u[k][1] - eps * pos[k][1],
                                  u[k][2] - eps * pos[k][2]}));
    }
    std::cout << "wedge" << (snapToAxis ? " (axis snap)" : "") << ": eps " << eps
              << ", max |u| " << umax << ", max error " << err << '\n';
    require(eps != 0.0, "Uniform pressure must expand the wedge");
    return err / umax;
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    try {
        rotatedBox();
        require(wedgeError(false) < 1e-9, "Tilted rollers must reproduce the isotropic expansion");
        require(wedgeError(true) > 1e-3, "Axis-snapped rollers must not reproduce it (test is blind)");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "general rollers: ok\n";
    return 0;
}
