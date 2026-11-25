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
#include <opm/geomech/FracturePressureAssemblerAD.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

constexpr int nx = 6, ny = 5, nc = nx * ny;
constexpr double min_width = 1e-3;

// Lattice fracture with widths on both sides of min_width (none at the kink),
// gravity, leak-off and a rate_well row fed by two cells.
FracturePressureInput makeInput(std::mt19937& gen)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    FracturePressureInput in;
    in.num_cells = nc;
    in.min_width = min_width;
    in.max_width = 0.5;
    in.control_type = "rate_well";
    in.num_well_equations = 1;
    in.total_wellindex = 3e-10;
    in.perfinj = {{7, 2e-9}, {8, 1e-9}};
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int c = i + nx * j;
            if (i + 1 < nx) in.htrans.emplace_back(c, c + 1, 0.5 + u(gen), 0.5 + u(gen));
            if (j + 1 < ny) in.htrans.emplace_back(c, c + nx, 0.5 + u(gen), 0.5 + u(gen));
        }
    }
    for (std::size_t k = 0; k < in.htrans.size(); ++k) in.face_gravity.push_back(1e4 * (u(gen) - 0.5));
    for (int c = 0; c < nc; ++c) {
        const double w = (c % 3 == 0) ? 0.3 * min_width * u(gen) : min_width * (1.5 + 4.0 * u(gen));
        in.fracture_width.push_back(w);
        in.fracture_pressure.push_back(2e7 + 1e6 * u(gen));
        in.face_mobility.push_back(1.0 / (1e-3 * (1.0 + u(gen))));
        in.leakof.push_back(1e-11 * u(gen));
        in.density.push_back({1000.0, 0.0});
        in.viscosity.push_back({1.0, 0.0});
    }
    in.fracture_width.push_back(0.0);
    in.fracture_pressure.push_back(2.1e7);
    in.face_mobility.push_back(0.0);
    in.density.push_back({1.0, 0.0});
    in.viscosity.push_back({1.0, 0.0});
    return in;
}

// Max over columns of |FD - assembled| / max|assembled entry in the block|.
double columnError(const BCRSMatrix1x1& J, int col, const BlockVector1& fd, double blockScale)
{
    double err = 0.0;
    for (std::size_t r = 0; r < fd.size(); ++r) {
        const double a = J.exists(r, col) ? J[r][col][0][0] : 0.0;
        err = std::max(err, std::abs(fd[r][0] - a) / blockScale);
    }
    return err;
}

double blockMax(const BCRSMatrix1x1& J)
{
    double m = 0.0;
    for (auto row = J.begin(); row != J.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col) m = std::max(m, std::abs((*col)[0][0]));
    return m;
}

// dR/dp (M) and dR/dw (C) of assemblePressureAD against central differences of its residual.
void flowBlocksMatchFiniteDifferences()
{
    std::mt19937 gen(42);
    const FracturePressureInput in = makeInput(gen);
    const auto base = assemblePressureAD(in);
    const double sM = blockMax(*base.pressure_matrix), sC = blockMax(*base.coupling_matrix);
    double errM = 0.0, errC = 0.0;
    for (int k = 0; k < nc + 1; ++k) {
        FracturePressureInput ip = in, im = in;
        const double h = 10.0; // constant mobility: residual linear in p
        ip.fracture_pressure[k] += h;
        im.fracture_pressure[k] -= h;
        BlockVector1 fd = assemblePressureAD(ip).residual;
        fd -= assemblePressureAD(im).residual;
        fd *= 0.5 / h;
        errM = std::max(errM, columnError(*base.pressure_matrix, k, fd, sM));
    }
    for (int k = 0; k < nc; ++k) {
        FracturePressureInput ip = in, im = in;
        const double h = 1e-4 * std::max(in.fracture_width[k], min_width);
        ip.fracture_width[k] += h;
        im.fracture_width[k] -= h;
        BlockVector1 fd = assemblePressureAD(ip).residual;
        fd -= assemblePressureAD(im).residual;
        fd *= 0.5 / h;
        errC = std::max(errC, columnError(*base.coupling_matrix, k, fd, sC));
    }
    std::cout << "flow rows: max rel error dR/dp " << errM << ", dR/dw " << errC << "\n";
    require(errM < 1e-6, "pressure block matches finite differences");
    require(errC < 1e-6, "coupling block matches finite differences");
    require(sC > 0.0, "coupling block is not empty (some widths above min_width)");
}

// Hand copy of Fracture::addSource's face gravity term T(w) g: rhs - M p with it
// must equal minus the AD residual.
void solverResidualFormMatchesAdResidual()
{
    std::mt19937 gen(7);
    const FracturePressureInput in = makeInput(gen);
    const auto ad = assemblePressureAD(in);
    BlockVector1 rhs(nc + 1);
    rhs = 0;
    for (std::size_t k = 0; k < in.htrans.size(); ++k) {
        const auto [i, j, t1, t2] = in.htrans[k];
        const double h1 = std::clamp(in.fracture_width[i], in.min_width, in.max_width);
        const double h2 = std::clamp(in.fracture_width[j], in.min_width, in.max_width);
        const double T = 0.5 * (in.face_mobility[i] + in.face_mobility[j])
            / (12.0 / (h1 * h1 * h1 * t1) + 12.0 / (h2 * h2 * h2 * t2));
        rhs[i] += T * in.face_gravity[k];
        rhs[j] -= T * in.face_gravity[k];
    }
    BlockVector1 p(nc + 1);
    for (int i = 0; i < nc + 1; ++i) p[i] = in.fracture_pressure[i];
    ad.pressure_matrix->mmv(p, rhs);
    double err = 0.0, scale = 0.0;
    for (int i = 0; i < nc + 1; ++i) {
        err = std::max(err, std::abs(rhs[i][0] + ad.residual[i][0]));
        scale = std::max(scale, std::abs(ad.residual[i][0]));
    }
    std::cout << "solver residual form vs AD residual: rel " << err / scale << "\n";
    require(err <= 1e-9 * scale, "rhs - M p equals minus the AD residual");
}

// Fischer-Burmeister mechanics rows: phi = a + b - sqrt(a^2 + b^2), a = c_i w_i,
// b = t - A w - p. Checks the derivation used in Fracture::assembleMechanicsRows_
// (B = diag(-db) A + diag(da c), I = diag(-db)), re-implemented here.
void fbRowsMatchFiniteDifferences()
{
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int n = 12;
    std::vector<double> A(n * n), t(n), w(n), p(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i * n + j] = (i == j) ? 1e10 * (2.0 + u(gen)) : 1e9 * (u(gen) - 0.5);
    for (int i = 0; i < n; ++i) {
        t[i] = -2.1e7 + 3e6 * (u(gen) - 0.5);
        p[i] = 2e7 + 1e6 * u(gen);
        w[i] = (i % 3 == 0) ? 0.0 : 1e-3 * u(gen);
    }
    const auto phi = [&](const std::vector<double>& ww, const std::vector<double>& pp, int i) {
        double b = t[i] - pp[i];
        for (int j = 0; j < n; ++j) b -= A[i * n + j] * ww[j];
        const double a = A[i * n + i] * ww[i];
        return a + b - std::sqrt(a * a + b * b);
    };
    double errW = 0.0, errP = 0.0, sW = 0.0;
    for (int i = 0; i < n; ++i) {
        double b = t[i] - p[i];
        for (int j = 0; j < n; ++j) b -= A[i * n + j] * w[j];
        const double c = A[i * n + i], a = c * w[i], r = std::sqrt(a * a + b * b);
        const double da = 1.0 - a / r, db = 1.0 - b / r;
        for (int k = 0; k < n; ++k) {
            const double Jik = -db * A[i * n + k] + (i == k ? da * c : 0.0);
            auto wp = w, wm = w;
            const double h = 1e-9;
            wp[k] += h;
            wm[k] -= h;
            errW = std::max(errW, std::abs((phi(wp, p, i) - phi(wm, p, i)) / (2 * h) - Jik));
            sW = std::max(sW, std::abs(Jik));
        }
        auto pp = p, pm = p;
        pp[i] += 1.0;
        pm[i] -= 1.0;
        errP = std::max(errP, std::abs((phi(w, pp, i) - phi(w, pm, i)) / 2.0 + db));
    }
    std::cout << "FB rows: max rel error d/dw " << errW / sW << ", d/dp " << errP << "\n";
    require(errW < 1e-5 * sW, "FB width block matches finite differences");
    require(errP < 1e-5, "FB pressure block matches finite differences");
}
} // namespace

int main()
{
    try {
        flowBlocksMatchFiniteDifferences();
        solverResidualFormMatchesAdResidual();
        fbRowsMatchFiniteDifferences();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
