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

// Mechanics only, on the simulator's edge-conformal CpGrid of a deck, uniform E and nu.
//   case=bend_x|bend_y : pure bending (sigma_xx resp. sigma_yy = E kappa z), exact displacement on
//                        every boundary node; prints the stress error against the analytic field.
//   case=flex_x|flex_y : plane-strain cantilever (Airy a x z^3 + b x z): moment varies along x (y),
//                        with shear; exact displacement on every boundary node.
//   case=linear        : a linear displacement field (patch test).
//   case=gauss         : isotropic eigenstrain e0*exp(-sum x_i^2/a_i^2) centred at cell src=i,j,k (1-based),
//                        widths ax, ay, az; exact infinite-space displacement on every boundary node.
//                        srcavg=1 loads each cell with its cell average instead of the centroid value.
//   case=load load=<f> : per-cell potential load (lines "global_cell value", Pa), deck boundary
//                        conditions (rollers on the sides, bottom fixed in z, top free).
// out=<f> writes "global_cell cx cy cz s0..s5" with the simulator's delstress (compression > 0).

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/umfpack.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/geomech/vem/vem.hpp>
#include <opm/geomech/vem/vemutils.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
// phi with lap(phi) = exp(-sum x_i^2/a_i^2): gradient and Hessian (xx, yy, zz, xy, yz, xz) by the
// heat-kernel integral phi = -int_0^inf (e^{t lap} g) dt, trapezoid in log t.
void
gauss_potential(const double x[3], const double a[3], double grad[3], double hess[6])
{
    const double amin = std::min({a[0], a[1], a[2]}), amax = std::max({a[0], a[1], a[2]});
    const double s0 = std::log(1e-8 * amin * amin), s1 = std::log(1e9 * amax * amax);
    const int n = 600;
    const double ds = (s1 - s0) / n;
    for (int k = 0; k < 3; ++k)
        grad[k] = 0;
    for (int k = 0; k < 6; ++k)
        hess[k] = 0;
    for (int i = 0; i <= n; ++i) {
        const double t = std::exp(s0 + i * ds), w = (i == 0 || i == n ? 0.5 : 1.0) * ds * t;
        double d[3], ex = 0, pre = 1;
        for (int k = 0; k < 3; ++k) {
            d[k] = a[k] * a[k] + 4 * t;
            pre *= a[k] / std::sqrt(d[k]);
            ex += x[k] * x[k] / d[k];
        }
        if (ex > 700)
            continue;
        const double P = w * pre * std::exp(-ex);
        for (int k = 0; k < 3; ++k) {
            grad[k] += P * 2 * x[k] / d[k];
            hess[k] += P * (2 / d[k] - 4 * x[k] * x[k] / (d[k] * d[k]));
        }
        hess[3] -= P * 4 * x[0] * x[1] / (d[0] * d[1]);
        hess[4] -= P * 4 * x[1] * x[2] / (d[1] * d[2]);
        hess[5] -= P * 4 * x[0] * x[2] / (d[0] * d[2]);
    }
}

// 27-point Gauss rule on the trilinear map of a hex (Dune corner order): points and weights (|J| included)
void
hex_quadrature(const std::array<Dune::FieldVector<double, 3>, 8>& X, std::vector<std::array<double, 3>>& pts,
               std::vector<double>& wts)
{
    const double gp[3] = {0.5 - 0.5 * std::sqrt(0.6), 0.5, 0.5 + 0.5 * std::sqrt(0.6)};
    const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};
    pts.clear();
    wts.clear();
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c) {
                const double r[3] = {gp[a], gp[b], gp[c]};
                std::array<double, 3> p {0, 0, 0};
                double J[3][3] = {{0}};
                for (int k = 0; k < 8; ++k) {
                    double N = 1, dN[3];
                    for (int e = 0; e < 3; ++e)
                        N *= ((k >> e) & 1) ? r[e] : 1 - r[e];
                    for (int e = 0; e < 3; ++e) {
                        dN[e] = ((k >> e) & 1) ? 1.0 : -1.0;
                        for (int f = 0; f < 3; ++f)
                            if (f != e)
                                dN[e] *= ((k >> f) & 1) ? r[f] : 1 - r[f];
                    }
                    for (int i = 0; i < 3; ++i) {
                        p[i] += N * X[k][i];
                        for (int j = 0; j < 3; ++j)
                            J[i][j] += X[k][i] * dN[j];
                    }
                }
                const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                    - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                    + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
                pts.push_back(p);
                wts.push_back(gw[a] * gw[b] * gw[c] * std::abs(det));
            }
}
} // namespace

int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    std::map<std::string, std::string> arg {{"choice", "3"}, {"quad", "triangulated"}, {"case", "bend_x"},
                                            {"E", "14e9"}, {"nu", "0.25"}, {"sigma", "1e7"}};
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        const auto eq = a.find('=');
        if (eq != std::string::npos)
            arg[a.substr(a.rfind('-', 1) == 1 ? 2 : 0, eq - (a.rfind('-', 1) == 1 ? 2 : 0))] = a.substr(eq + 1);
    }
    const double E = std::stod(arg["E"]), nu = std::stod(arg["nu"]);
    const auto sc = vem::StabilityChoice(std::stoi(arg["choice"]));
    vem::set_face_quadrature(arg["quad"] == "kite" ? vem::FaceQuadrature::KITE
                                                   : vem::FaceQuadrature::TRIANGULATED);
    const std::string cs = arg["case"];

    Opm::Parser parser;
    Opm::ParseContext ctx(Opm::InputErrorAction::WARN);
    Opm::ErrorGuard guard;
    const auto deck = parser.parseFile(arg["grid"], ctx, guard);
    Opm::EclipseState state(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&state.getInputGrid(), &state, false, false, false, /*edge_conformal*/ true);
    const auto& gv = grid.leafGridView();

    std::vector<double> coords;
    std::vector<int> ncf, nfc, fc;
    vem::getGridVectors(grid, coords, ncf, nfc, fc);
    const int nn = int(coords.size() / 3), nc = grid.size(0);
    {
        int mism = 0;
        double dmax = 0;
        for (int v = 0; v < nn; ++v) {
            const auto p = grid.vertexPosition(v);
            const double d = std::abs(p[0] - coords[3 * v]) + std::abs(p[1] - coords[3 * v + 1])
                + std::abs(p[2] - coords[3 * v + 2]);
            mism += d > 1e-6;
            dmax = std::max(dmax, d);
        }
        if (mism)
            std::printf("coords vs vertexPosition: %d mismatches (max %.3g)\n", mism, dmax);
    }
    {
        // grid check: every cell must be a closed polyhedron whose divergence volume matches
        int open = 0, nonpos = 0, mism = 0, shortedge = 0;
        std::FILE* bad = arg.count("badcells") ? std::fopen(arg["badcells"].c_str(), "w") : nullptr;
        std::FILE* q = arg.count("cellq") ? std::fopen(arg["cellq"].c_str(), "w") : nullptr;
        for (int c2 = 0, f = 0, k = 0; c2 < nc; ++c2) {
            std::map<std::pair<int, int>, int> dir;
            double vol6 = 0, emin = 1e300;
            const auto cc = grid.cellCentroid(c2);
            for (int ff = 0; ff < ncf[c2]; ++ff, k += nfc[f], ++f) {
                const int m = nfc[f];
                double cen[3] = {0, 0, 0};
                for (int i = 0; i < m; ++i)
                    for (int d = 0; d < 3; ++d)
                        cen[d] += coords[3 * fc[k + i] + d] / m;
                for (int i = 0; i < m; ++i) {
                    const int u = fc[k + i], v = fc[k + (i + 1) % m];
                    dir[{u, v}] += 1;
                    const double* pu = &coords[3 * u];
                    const double* pv = &coords[3 * v];
                    emin = std::min(emin, std::hypot(pu[0] - pv[0], pu[1] - pv[1], pu[2] - pv[2]));
                    double a[3], bq[3], cq[3];
                    for (int d = 0; d < 3; ++d)
                        a[d] = cen[d] - cc[d], bq[d] = pu[d] - cc[d], cq[d] = pv[d] - cc[d];
                    vol6 += a[0] * (bq[1] * cq[2] - bq[2] * cq[1]) - a[1] * (bq[0] * cq[2] - bq[2] * cq[0])
                        + a[2] * (bq[0] * cq[1] - bq[1] * cq[0]);
                }
            }
            int unpaired = 0;
            for (const auto& [e, n] : dir) {
                const auto it = dir.find({e.second, e.first});
                unpaired += it == dir.end() || it->second != n;
            }
            const double vol = grid.cellVolume(c2);
            const bool is_open = unpaired > 0, is_neg = !(vol6 > 0);
            const bool is_mism = std::abs(vol6 / 6 - vol) > 0.05 * std::abs(vol);
            open += is_open, nonpos += is_neg, mism += is_mism, shortedge += emin < 0.01;
            if (bad && (is_open || is_neg || is_mism))
                std::fprintf(bad, "%d\n", grid.globalCell()[c2]);
            if (q) // global cell, shortest edge, volume, faces, nodes
                std::fprintf(q, "%d %.4g %.4g %d %zu\n", grid.globalCell()[c2], emin, vol, ncf[c2], dir.size());
        }
        if (bad)
            std::fclose(bad);
        if (q)
            std::fclose(q);
        std::printf("grid check: %d cells, %d open surfaces, %d non-positive volume, %d volume mismatch > 5%%, "
                    "%d with an edge < 1 cm\n", nc, open, nonpos, mism, shortedge);
        if (arg.count("checkonly"))
            return 0;
    }
    std::vector<double> young(nc, E), poisson(nc, nu), body(3 * nc, 0.0), field(nc, 0.0);
    // gauss: source centre, widths, exact cell-average stress
    const bool gauss = cs == "gauss";
    std::array<double, 3> src {0, 0, 0}, aw {std::stod(arg.count("ax") ? arg["ax"] : "200"),
                                             std::stod(arg.count("ay") ? arg["ay"] : "200"),
                                             std::stod(arg.count("az") ? arg["az"] : "20")};
    const double e0 = std::stod(arg.count("e0") ? arg["e0"] : "1e-3");
    const double mu = E / (2 * (1 + nu)), Kb = E / (3 * (1 - 2 * nu)), Af = (1 + nu) / (1 - nu) * e0;
    std::vector<std::array<double, 6>> exact_stress(nc, {0, 0, 0, 0, 0, 0});
    if (gauss) {
        int si = 0, sj = 0, sk = 0;
        std::sscanf(arg["src"].c_str(), "%d,%d,%d", &si, &sj, &sk);
        const auto dims = grid.logicalCartesianSize();
        const int gsrc = (si - 1) + dims[0] * ((sj - 1) + dims[1] * (sk - 1));
        bool found = false;
        for (int c = 0; c < nc; ++c)
            if (grid.globalCell()[c] == gsrc) {
                const auto cc = grid.cellCentroid(c);
                src = {cc[0], cc[1], cc[2]};
                found = true;
            }
        if (!found)
            throw std::runtime_error("src cell is not active");
        const bool avg = arg.count("srcavg") && arg["srcavg"] == "1";
        std::vector<std::array<double, 3>> qp;
        std::vector<double> qw;
        for (const auto& e : elements(gv)) {
            const int c = gv.indexSet().index(e);
            std::array<Dune::FieldVector<double, 3>, 8> X;
            for (int k = 0; k < 8; ++k)
                X[k] = grid.vertexPosition(gv.indexSet().subIndex(e, k, 3));
            hex_quadrature(X, qp, qw);
            double vol = 0, gavg = 0;
            std::array<double, 6> sav {0, 0, 0, 0, 0, 0};
            for (std::size_t q = 0; q < qp.size(); ++q) {
                double x[3], gr[3], h[6];
                for (int d = 0; d < 3; ++d)
                    x[d] = qp[q][d] - src[d];
                const double g = std::exp(-(x[0] * x[0] / (aw[0] * aw[0]) + x[1] * x[1] / (aw[1] * aw[1])
                                            + x[2] * x[2] / (aw[2] * aw[2])));
                gauss_potential(x, aw.data(), gr, h);
                vol += qw[q];
                gavg += qw[q] * g;
                // total stress 2 mu (grad grad psi - lap psi I), psi = Af phi, lap phi = g
                for (int k = 0; k < 6; ++k)
                    sav[k] += qw[q] * 2 * mu * Af * (h[k] - (k < 3 ? g : 0.0));
            }
            for (int k = 0; k < 6; ++k)
                exact_stress[c][k] = sav[k] / vol;
            const auto cc = grid.cellCentroid(c);
            const double xc[3] = {cc[0] - src[0], cc[1] - src[1], cc[2] - src[2]};
            const double gc = std::exp(-(xc[0] * xc[0] / (aw[0] * aw[0]) + xc[1] * xc[1] / (aw[1] * aw[1])
                                         + xc[2] * xc[2] / (aw[2] * aw[2])));
            field[c] = 3 * Kb * e0 * (avg ? gavg / vol : gc);
        }
    }
    if (cs == "load") {
        std::map<int, double> load;
        std::ifstream in(arg["load"]);
        int g;
        double v;
        while (in >> g >> v)
            load[g] = v;
        for (int c = 0; c < nc; ++c)
            field[c] = load.count(grid.globalCell()[c]) ? load[grid.globalCell()[c]] : 0.0;
    }

    // centre and size of the model, for the analytic fields
    std::array<double, 3> lo {1e300, 1e300, 1e300}, hi {-1e300, -1e300, -1e300};
    for (int n = 0; n < nn; ++n)
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], coords[3 * n + d]);
            hi[d] = std::max(hi[d], coords[3 * n + d]);
        }
    const std::array<double, 3> x0 {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
    const double kappa = std::stod(arg["sigma"]) / (E * 0.5 * (hi[2] - lo[2]));
    // flex: sigma_ll = 6 a L Z, sigma_lz = 3 a (c^2 - Z^2), plane strain across l (l = x or y)
    const double half = 0.5 * (hi[2] - lo[2]);
    const int l = cs == "flex_y" ? 1 : 0;
    const double a = std::stod(arg["sigma"]) / (6 * 0.5 * (hi[l] - lo[l]) * half);
    const double Ep = E / (1 - nu * nu), nup = nu / (1 - nu), G2 = 2 * (1 + nu) / E, bb = -3 * a * half * half;
    const auto exact = [&](const double* p, double u[3]) {
        const double X = p[0] - x0[0], Y = p[1] - x0[1], Z = p[2] - x0[2];
        if (cs == "gauss") {
            const double x[3] = {p[0] - src[0], p[1] - src[1], p[2] - src[2]};
            double gr[3], h[6];
            gauss_potential(x, aw.data(), gr, h);
            for (int d = 0; d < 3; ++d)
                u[d] = Af * gr[d];
        } else if (cs == "flex_x" || cs == "flex_y") {
            const double L = l == 0 ? X : Y;
            const double ul = 3 * a * L * L * Z / Ep + nup * a * Z * Z * Z / Ep - G2 * (a * Z * Z * Z + bb * Z);
            u[l] = ul, u[1 - l] = 0, u[2] = -3 * nup * a * L * Z * Z / Ep - a * L * L * L / Ep;
        } else if (cs == "bend_x") {
            u[0] = kappa * X * Z, u[1] = -nu * kappa * Y * Z, u[2] = -0.5 * kappa * (X * X + nu * (Z * Z - Y * Y));
        } else if (cs == "bend_y") {
            u[1] = kappa * Y * Z, u[0] = -nu * kappa * X * Z, u[2] = -0.5 * kappa * (Y * Y + nu * (Z * Z - X * X));
        } else {
            u[0] = 1e-4 * X + 2e-4 * Y - 1e-4 * Z, u[1] = -3e-4 * Y + 1e-4 * Z, u[2] = 2e-4 * X + 1e-4 * Z;
        }
    };

    // boundary nodes and the outward normal of their boundary faces
    std::vector<std::array<bool, 3>> fix(nn, {false, false, false});
    for (int f = 0; f < grid.numFaces(); ++f) {
        if (grid.faceCell(f, 0) >= 0 && grid.faceCell(f, 1) >= 0)
            continue;
        const int m = grid.numFaceVertices(f);
        if (m < 3)
            continue;
        const int cell = grid.faceCell(f, 0) >= 0 ? grid.faceCell(f, 0) : grid.faceCell(f, 1);
        std::array<double, 3> nrm {0, 0, 0}, fcen {0, 0, 0};
        for (int i = 0; i < m; ++i) {
            const auto a = grid.vertexPosition(grid.faceVertex(f, i));
            const auto b = grid.vertexPosition(grid.faceVertex(f, (i + 1) % m));
            nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
            nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
            nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
            for (int d = 0; d < 3; ++d)
                fcen[d] += a[d] / m;
        }
        const auto cc = grid.cellCentroid(cell);
        double s = 0;
        for (int d = 0; d < 3; ++d)
            s += nrm[d] * (fcen[d] - cc[d]);
        if (s < 0)
            for (auto& x : nrm)
                x = -x;
        const int dmax = std::abs(nrm[0]) >= std::abs(nrm[1]) && std::abs(nrm[0]) >= std::abs(nrm[2]) ? 0
            : std::abs(nrm[1]) >= std::abs(nrm[2])                                               ? 1
                                                                                                 : 2;
        for (int i = 0; i < m; ++i) {
            auto& fx = fix[grid.faceVertex(f, i)];
            if (cs != "load")
                fx = {true, true, true}; // analytic cases: exact displacement everywhere on the boundary
            else if (dmax < 2)
                fx[dmax] = true; // rollers on the sides
            else if (nrm[2] > 0)
                fx[2] = true; // bottom (z grows downwards); top is free
        }
    }
    std::vector<int> fixed;
    std::vector<double> fixval;
    for (int n = 0; n < nn; ++n) {
        double u[3] = {0, 0, 0};
        if (cs != "load")
            exact(&coords[3 * n], u);
        for (int d = 0; d < 3; ++d)
            if (fix[n][d]) {
                fixed.push_back(3 * n + d);
                fixval.push_back(u[d]);
            }
    }

    {
        std::vector<int> used(nn, 0);
        for (int v : fc)
            used[v] = 1;
        int unused = 0;
        for (int n = 0; n < nn; ++n)
            if (!used[n]) {
                ++unused;
                for (int d = 0; d < 3; ++d)
                    fix[n][d] = true;
            }
        if (unused)
            std::printf("%d nodes belong to no cell face: fixed\n", unused);
        fixed.clear();
        fixval.clear();
        for (int n = 0; n < nn; ++n) {
            double u[3] = {0, 0, 0};
            if (cs != "load")
                exact(&coords[3 * n], u);
            for (int d = 0; d < 3; ++d)
                if (fix[n][d]) {
                    fixed.push_back(3 * n + d);
                    fixval.push_back(u[d]);
                }
        }
    }
    std::vector<std::tuple<int, int, double>> A;
    std::vector<double> b;
    const int nfree = vem::assemble_mech_system_3D_dune(grid, &coords[0], nc, &ncf[0], &nfc[0], &fc[0], &young[0],
                                                        &poisson[0], &body[0], int(fixed.size()), &fixed[0],
                                                        &fixval[0], 0, nullptr, nullptr, A, b, sc, true);
    std::vector<int> free_dofs;
    for (int i = 0, k = 0; i < 3 * nn; ++i) {
        if (k < int(fixed.size()) && fixed[k] == i)
            ++k;
        else
            free_dofs.push_back(i);
    }
    if (cs == "load" || gauss) {
        std::vector<double> fgrad;
        std::vector<std::tuple<int, int, double>> div;
        vem::potential_gradient_force_3D(&coords[0], nc, &ncf[0], &nfc[0], &fc[0], &field[0], fgrad, div, false);
        for (int i = 0; i < nfree; ++i)
            b[i] += fgrad[free_dofs[i]];
    }

    {
        int bad = 0;
        std::vector<double> diag(nfree, 0.0);
        for (const auto& [r, c2, v] : A) {
            bad += !std::isfinite(v);
            if (r == c2)
                diag[r] += v;
        }
        int nonpos = 0;
        for (double d : diag)
            nonpos += !(d > 0);
        if (bad || nonpos)
            std::printf("matrix: %d non-finite entries, %d non-positive diagonals\n", bad, nonpos);
    }
    using Mat = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Vec = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    std::sort(A.begin(), A.end());
    std::vector<int> rowlen(nfree, 0);
    for (std::size_t i = 0; i < A.size(); ++i)
        if (i == 0 || std::get<0>(A[i]) != std::get<0>(A[i - 1]) || std::get<1>(A[i]) != std::get<1>(A[i - 1]))
            ++rowlen[std::get<0>(A[i])];
    Mat M(nfree, nfree, Mat::random);
    for (int r = 0; r < nfree; ++r)
        M.setrowsize(r, rowlen[r]);
    M.endrowsizes();
    for (const auto& [r, c, v] : A)
        M.addindex(r, c);
    M.endindices();
    M = 0.0;
    for (const auto& [r, c, v] : A)
        M[r][c] += v;
    if (arg.count("verbose") && cs == "linear") {
        // residual of the exact field at the free dofs: identifies nodes where interior forces do not cancel
        std::vector<double> r(b.begin(), b.begin() + nfree);
        for (const auto& [rr, cc2, v] : A) {
            double u[3];
            exact(&coords[3 * (free_dofs[cc2] / 3)], u);
            r[rr] -= v * u[free_dofs[cc2] % 3];
        }
        std::vector<std::pair<double, int>> big;
        double rmax = 0;
        for (int i = 0; i < nfree; ++i) {
            rmax = std::max(rmax, std::abs(r[i]));
            big.push_back({std::abs(r[i]), free_dofs[i] / 3});
        }
        std::sort(big.rbegin(), big.rend());
        int nbad = 0;
        for (const auto& [v, n] : big)
            nbad += v > 1e-8 * big.front().first + 1e-3;
        std::printf("residual of exact field: max %.3g, %d dofs above 1e-3\n", rmax, nbad);
        for (int w = 0; w < 6 && w < int(big.size()); ++w)
            std::printf("  node %d at %.2f %.2f %.2f residual %.3g\n", big[w].second, coords[3 * big[w].second],
                        coords[3 * big[w].second + 1], coords[3 * big[w].second + 2], big[w].first);
        // the cells around the worst nodes: closed surface? volume? shortest edge?
        std::vector<int> fstart(nc + 1, 0), ffirst(nc + 1, 0);
        for (int c2 = 0, fpos = 0; c2 < nc; ++c2) {
            int cnt = 0;
            for (int f = 0; f < ncf[c2]; ++f)
                cnt += nfc[fpos + f];
            fstart[c2 + 1] = fstart[c2] + cnt;
            ffirst[c2 + 1] = ffirst[c2] + ncf[c2];
            fpos += ncf[c2];
        }
        std::set<int> worstnodes;
        for (int w = 0; w < 20 && w < int(big.size()); ++w)
            worstnodes.insert(big[w].second);
        for (int c2 = 0; c2 < nc; ++c2) {
            bool touch = false;
            for (int k = fstart[c2]; k < fstart[c2 + 1]; ++k)
                touch |= worstnodes.count(fc[k]) > 0;
            if (!touch)
                continue;
            std::map<std::pair<int, int>, int> dir; // directed edge count
            double emin = 1e300, vol6 = 0;
            for (int f = ffirst[c2], k = fstart[c2]; f < ffirst[c2 + 1]; k += nfc[f], ++f) {
                const int m = nfc[f];
                double cen[3] = {0, 0, 0};
                for (int i = 0; i < m; ++i)
                    for (int d = 0; d < 3; ++d)
                        cen[d] += coords[3 * fc[k + i] + d] / m;
                for (int i = 0; i < m; ++i) {
                    const int u = fc[k + i], v = fc[k + (i + 1) % m];
                    dir[{u, v}] += 1;
                    const double* pu = &coords[3 * u];
                    const double* pv = &coords[3 * v];
                    emin = std::min(emin, std::hypot(pu[0] - pv[0], pu[1] - pv[1], pu[2] - pv[2]));
                    // signed volume of the tetra (origin = cell centroid) on the fan triangle
                    const auto cc = grid.cellCentroid(c2);
                    double a[3], bq[3], cq[3];
                    for (int d = 0; d < 3; ++d)
                        a[d] = cen[d] - cc[d], bq[d] = pu[d] - cc[d], cq[d] = pv[d] - cc[d];
                    vol6 += a[0] * (bq[1] * cq[2] - bq[2] * cq[1]) - a[1] * (bq[0] * cq[2] - bq[2] * cq[0])
                        + a[2] * (bq[0] * cq[1] - bq[1] * cq[0]);
                }
            }
            int unpaired = 0;
            for (const auto& [e, n] : dir)
                unpaired += dir.count({e.second, e.first}) == 0 || dir.at({e.second, e.first}) != n;
            const int g = grid.globalCell()[c2];
            std::printf("  cell %d ijk %d %d %d: faces %d, unpaired directed edges %d, volume %.4g (fan %.4g), "
                        "shortest edge %.3g\n",
                        c2, g % 46 + 1, (g / 46) % 112 + 1, g / (46 * 112) + 1, ffirst[c2 + 1] - ffirst[c2], unpaired,
                        grid.cellVolume(c2), vol6 / 6, emin);
        }
    }
    Vec rhs(nfree), x(nfree);
    for (int i = 0; i < nfree; ++i)
        rhs[i] = b[i];
    Dune::UMFPack<Mat> lu(M, arg.count("verbose") ? 2 : 0);
    Dune::InverseOperatorResult res;
    lu.apply(x, rhs, res);

    std::vector<double> disp(3 * nn, 0.0);
    for (std::size_t i = 0; i < fixed.size(); ++i)
        disp[fixed[i]] = fixval[i];
    for (int i = 0; i < nfree; ++i)
        disp[free_dofs[i]] = x[i];
    std::vector<std::array<double, 6>> stress(nc);
    std::vector<std::tuple<int, int, double>> smat;
    vem::compute_stress_3D(&coords[0], nc, &ncf[0], &nfc[0], &fc[0], &young[0], &poisson[0], disp, stress, smat,
                           false, true, sc, false);

    if (gauss) {
        // total stress sigma(u) - field I against the exact cell average; split near / far from the source
        double num = 0, den = 0, numn = 0, denn = 0, emax = 0, smax = 0;
        for (int c = 0; c < nc; ++c)
            for (int k = 0; k < 6; ++k)
                smax = std::max(smax, std::abs(exact_stress[c][k]));
        for (int c = 0; c < nc; ++c) {
            const auto cc = grid.cellCentroid(c);
            double r2 = 0;
            for (int d = 0; d < 3; ++d)
                r2 += (cc[d] - src[d]) * (cc[d] - src[d]) / (aw[d] * aw[d]);
            const double vol = grid.cellVolume(c);
            for (int k = 0; k < 6; ++k) {
                const double d = stress[c][k] - (k < 3 ? field[c] : 0.0) - exact_stress[c][k];
                num += vol * d * d;
                den += vol * exact_stress[c][k] * exact_stress[c][k];
                if (r2 < 4) {
                    numn += vol * d * d;
                    denn += vol * exact_stress[c][k] * exact_stress[c][k];
                }
                emax = std::max(emax, std::abs(d));
            }
        }
        std::printf("case gauss choice %d a %g %g %g srcavg %s: stress rel L2 %.4e (within 2a: %.4e), "
                    "max err / max stress %.4e (max stress %.4g MPa)\n",
                    int(sc), aw[0], aw[1], aw[2], arg.count("srcavg") ? arg["srcavg"].c_str() : "0",
                    std::sqrt(num / den), std::sqrt(numn / denn), emax / smax, smax / 1e6);
    }
    if (cs != "load" && !gauss) {
        // analytic: the bending stress is the only nonzero component, linear in z (exact at the centroid)
        double num = 0, den = 0, emax = 0, umax = 0, uref = 0, smax = std::stod(arg["sigma"]);
        std::vector<std::pair<double, int>> cellerr;
        const int ib = cs == "bend_y" ? 1 : 0; // Voigt order xx, yy, zz, shears
        for (const auto& e : elements(gv)) {
            const int c = gv.indexSet().index(e);
            const auto cen = grid.cellCentroid(c);
            const double vol = grid.cellVolume(c);
            std::array<double, 6> ex {0, 0, 0, 0, 0, 0}; // xx, yy, zz, xy, yz, xz
            const double Z = cen[2] - x0[2];
            if (cs == "bend_x" || cs == "bend_y") {
                ex[ib] = E * kappa * Z;
            } else if (cs == "flex_x" || cs == "flex_y") {
                const double L = cen[l] - x0[l];
                ex[l] = 6 * a * L * Z;
                ex[1 - l] = nu * ex[l];
                ex[l == 0 ? 5 : 4] = 3 * a * (half * half - Z * Z);
            } else { // linear field: constant strain (engineering shears)
                const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
                const double e[6] = {1e-4, -3e-4, 1e-4, 2e-4, 1e-4, 1e-4};
                for (int k = 0; k < 3; ++k)
                    ex[k] = lam * (e[0] + e[1] + e[2]) + 2 * mu * e[k];
                for (int k = 3; k < 6; ++k)
                    ex[k] = mu * e[k];
            }
            double ce = 0;
            for (int k = 0; k < 6; ++k) {
                const double d = stress[c][k] - ex[k];
                num += vol * d * d;
                den += vol * ex[k] * ex[k];
                emax = std::max(emax, std::abs(d));
                ce = std::max(ce, std::abs(d));
            }
            cellerr.push_back({ce / smax, c});
        }
        if (arg.count("verbose")) {
            std::vector<std::pair<double, int>> worst;
            for (int cc2 = 0; cc2 < nc; ++cc2) {
                double e2 = 0;
                for (int k = 0; k < 6; ++k)
                    e2 = std::max(e2, std::abs(stress[cc2][k]));
                worst.push_back({e2, cc2});
            }
            std::sort(worst.rbegin(), worst.rend());
            int fpos = 0;
            std::vector<int> fstart(nc + 1, 0);
            for (int cc2 = 0; cc2 < nc; ++cc2) {
                int cnt = 0;
                for (int f = 0; f < ncf[cc2]; ++f)
                    cnt += nfc[fpos + f];
                fstart[cc2 + 1] = fstart[cc2] + cnt;
                fpos += ncf[cc2];
            }
            for (int w = 0; w < 8; ++w) {
                const int cc2 = worst[w].second;
                std::vector<int> nodes(fc.begin() + fstart[cc2], fc.begin() + fstart[cc2 + 1]);
                std::sort(nodes.begin(), nodes.end());
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
                const int g = grid.globalCell()[cc2];
                std::printf("worst cell %d (ijk %d %d %d): |stress| %.3g faces %d nodes %zu volume %.3g\n", cc2,
                            g % 46 + 1, (g / 46) % 112 + 1, g / (46 * 112) + 1, worst[w].first, ncf[cc2],
                            nodes.size(), grid.cellVolume(cc2));
            }
        }
        for (int n = 0; n < nn; ++n) {
            double u[3];
            exact(&coords[3 * n], u);
            for (int d = 0; d < 3; ++d) {
                umax = std::max(umax, std::abs(disp[3 * n + d] - u[d]));
                uref = std::max(uref, std::abs(u[d]));
            }
        }
        std::sort(cellerr.begin(), cellerr.end());
        const auto pct = [&](double q) { return cellerr[std::min(nc - 1, int(q * nc))].first; };
        std::printf("case %s choice %d cells %d nodes %d: stress rel L2 %.4e, max err %.4e (/sigma %.4e), "
                    "disp max err rel %.4e; cell err/sigma p50 %.3g p90 %.3g p99 %.3g p99.9 %.3g\n",
                    cs.c_str(), int(sc), nc, nn, den > 0 ? std::sqrt(num / den) : std::sqrt(num), emax,
                    emax / smax, umax / std::max(uref, 1e-300), pct(0.5), pct(0.9), pct(0.99), pct(0.999));
        if (arg.count("errout")) {
            std::FILE* fe = std::fopen(arg["errout"].c_str(), "w");
            for (const auto& [e2, c2] : cellerr)
                std::fprintf(fe, "%d %.6g\n", grid.globalCell()[c2], e2);
            std::fclose(fe);
        }
    }
    if (arg.count("out")) {
        std::FILE* fp = std::fopen(arg["out"].c_str(), "w");
        for (int c = 0; c < nc; ++c) {
            const auto cen = grid.cellCentroid(c);
            std::fprintf(fp, "%d %.6f %.6f %.6f", grid.globalCell()[c], cen[0], cen[1], cen[2]);
            for (int k = 0; k < 6; ++k)
                std::fprintf(fp, " %.6e", -stress[c][k] + (k < 3 ? field[c] : 0.0));
            std::fprintf(fp, "\n");
        }
        std::fclose(fp);
    }
    return 0;
}
