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

// The corner-point FEM_BBAR path must give the same Q1 / B-bar element matrices whether
// they come from Opm::Elasticity (opm-upscaling) or from vem::stiffness_matrix_fem_hex_3D.

#include "config.h"

#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/materials.hh>
#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int nx = 2, ny = 2, nz = 2;

// Pillars inclined by `shear` (top x offset per unit depth); ZCORN perturbed by `warp`
// per corner so that cell faces are not planar.
std::string deckString(double aspect, double shear, double warp)
{
    const double dx = 10.0 * aspect, dy = 10.0, dz = 10.0;
    std::ostringstream os;
    os << "RUNSPEC\nDIMENS\n" << nx << " " << ny << " " << nz << " /\nOIL\nWATER\nMETRIC\n"
       << "START\n1 'JAN' 2000 /\nGRID\nCOORD\n";
    const double ztop = 1000.0, zbot = ztop + nz * dz;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            os << i * dx << " " << j * dy << " " << ztop << " "
               << i * dx + shear * (zbot - ztop) << " " << j * dy << " " << zbot << "\n";
    os << "/\nZCORN\n";
    for (int k = 0; k < nz; ++k)
        for (int s = 0; s < 2; ++s)
            for (int j = 0; j < 2 * ny; ++j)
                for (int i = 0; i < 2 * nx; ++i) {
                    const int layer = k + s;
                    const int pi = (i + 1) / 2, pj = (j + 1) / 2;
                    // interior surfaces move, the outer ones stay flat; same value on both
                    // sides of a shared corner keeps the grid conforming
                    const bool inner = layer > 0 && layer < nz;
                    const double bump = inner ? warp * dz * (((pi + 2 * pj) % 3) - 1.0) : 0.0;
                    os << ztop + layer * dz + bump << "\n";
                }
    os << "/\nPORO\n" << nx * ny * nz << "*0.3 /\nEND\n";
    return os.str();
}

Dune::CpGrid makeGrid(const std::string& deck_str)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_str);
    Opm::EclipseState state(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&state.getInputGrid(), nullptr, false);
    return grid;
}

// The element as assembled by assemble_mech_system_3D_dune before it used vem::.
std::vector<double> legacyElement(const Dune::CpGrid& grid,
                                  const Dune::CpGrid::LeafGridView::Codim<0>::Entity& elem,
                                  double young, double poisson, bool do_bbar)
{
    constexpr int dim = 3, comp = 6, bfunc = 8;
    Opm::Elasticity::Elasticity E(grid);
    Dune::FieldMatrix<double, comp, comp> C;
    Opm::Elasticity::Isotropic material(1, young, poisson, 0.0);
    material.getConstitutiveMatrix(C);
    const auto& rule = Dune::QuadratureRules<double, dim>::rule(elem.type(), 2);
    std::vector<Dune::FieldMatrix<double, comp, dim * bfunc>> Bmats;
    std::vector<double> qweights;
    Dune::FieldVector<double, dim * bfunc> bbar_row(0.0);
    double cell_volume = 0.0;
    for (auto r = rule.begin(); r != rule.end(); ++r) {
        const auto jacInvTra = elem.geometry().jacobianInverseTransposed(r->position());
        const double detJ = elem.geometry().integrationElement(r->position());
        Dune::FieldMatrix<double, comp, dim * bfunc> lB;
        E.getBmatrix(lB, r->position(), jacInvTra);
        const double qw = detJ * r->weight();
        Bmats.push_back(lB);
        qweights.push_back(qw);
        cell_volume += qw;
        if (do_bbar)
            for (int col = 0; col < dim * bfunc; ++col)
                bbar_row[col] += lB[col % dim][col] * qw;
    }
    if (do_bbar)
        bbar_row /= cell_volume;
    Dune::FieldMatrix<double, dim * bfunc, dim * bfunc> K(0.0), Aq;
    for (std::size_t qp = 0; qp < Bmats.size(); ++qp) {
        auto lB = Bmats[qp];
        if (do_bbar)
            for (int col = 0; col < dim * bfunc; ++col) {
                const double corr = (bbar_row[col] - Bmats[qp][col % dim][col]) / dim;
                for (int row = 0; row < dim; ++row)
                    lB[row][col] += corr;
            }
        E.getStiffnessMatrix(Aq, lB, C, qweights[qp]);
        K += Aq;
    }
    std::vector<double> out(24 * 24);
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j)
            out[i * 24 + j] = K[i][j];
    return out;
}

double maxRelDiff(const Dune::CpGrid& grid, bool bbar)
{
    std::vector<double> coords;
    std::vector<int> ncf, nfc, fc;
    vem::getGridVectors(grid, coords, ncf, nfc, fc);
    const auto gv = grid.leafGridView();
    const double young = 1e9, poisson = 0.3;
    double worst = 0.0;
    for (const auto& elem : elements(gv)) {
        std::array<double, 24> hc;
        for (int k = 0; k < 8; ++k) {
            const int node = gv.indexSet().subIndex(elem, k, 3);
            for (int d = 0; d < 3; ++d)
                hc[3 * k + d] = coords[3 * node + d];
        }
        std::vector<double> Kv(24 * 24);
        vem::stiffness_matrix_fem_hex_3D(hc.data(), young, poisson, bbar, Kv.data());
        const auto Kl = legacyElement(grid, elem, young, poisson, bbar);
        double diff = 0.0, scale = 0.0;
        for (int i = 0; i < 24 * 24; ++i) {
            diff = std::max(diff, std::abs(Kl[i] - Kv[i]));
            scale = std::max(scale, std::abs(Kl[i]));
        }
        worst = std::max(worst, diff / scale);
    }
    return worst;
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    struct Case { const char* name; double aspect, shear, warp; };
    const Case cases[] = {{"regular", 1.0, 0.0, 0.0},
                          {"aspect 20", 20.0, 0.0, 0.0},
                          {"sheared 30 deg", 1.0, std::tan(M_PI / 6), 0.0},
                          {"warped", 1.0, 0.0, 0.15},
                          {"sheared + warped, aspect 10", 10.0, 0.3, 0.1}};
    int failures = 0;
    for (const auto& c : cases) {
        const auto grid = makeGrid(deckString(c.aspect, c.shear, c.warp));
        for (const bool bbar : {false, true}) {
            const double d = maxRelDiff(grid, bbar);
            const bool ok = d < 1e-12;
            failures += !ok;
            std::cout << c.name << (bbar ? ", B-bar" : ", Q1") << ": max rel diff " << d
                      << (ok ? "" : "  FAIL") << "\n";
        }
    }
    std::cout << (failures ? "FAILED" : "PASSED") << "\n";
    return failures ? 1 : 0;
}
