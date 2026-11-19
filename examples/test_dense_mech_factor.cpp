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

// DenseMechFactor: a Woodbury row update with tol = 0 must reproduce a direct solve.

#include <config.h>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

namespace
{
using Opm::FMatrix;
using Opm::Vector;

FMatrix
block(const FMatrix& A, const std::vector<double>& alpha, const std::vector<double>& beta, double s)
{
    const std::size_t n = A.N();
    FMatrix B(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j)
            B[i][j] = s * alpha[i] * A[i][j];
        B[i][i] += s * beta[i];
    }
    return B;
}

// normwise backward error ||b - Bx|| / (||B|| ||x|| + ||b||)
double
relError(const FMatrix& B, const Opm::DenseMechFactor& f, const Vector& rhs)
{
    Vector x(rhs.size()), r(rhs);
    f.solve(x, rhs);
    B.mmv(x, r);
    return r.infinity_norm() / (B.infinity_norm() * x.infinity_norm() + rhs.infinity_norm());
}

int
check(bool inverse)
{
    const int n = 300;
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    FMatrix A(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = (i == j) ? 1e9 * (5.0 + u(gen)) : 1e8 * u(gen) / (1 + std::abs(i - j));
    std::vector<double> alpha(n, 1.0), beta(n, 0.0);
    for (int i = 0; i < n; i += 7) {
        alpha[i] = 0.0;
        beta[i] = 1.0;
    }
    Vector rhs(n);
    for (int i = 0; i < n; ++i)
        rhs[i] = u(gen);

    Opm::DenseMechFactor::Options opt;
    opt.inverse = inverse;
    opt.tol = 0.0;
    opt.maxRankFraction = 0.2;
    opt.refactorIterations = 0;
    Opm::DenseMechFactor f(opt);
    Opm::MechRows rows {&A, &alpha, &beta, 1};
    int fails = 0;
    const auto expect = [&](bool ok, const char* what) {
        if (!ok) {
            std::cerr << (inverse ? "inverse: " : "lu: ") << what << " failed\n";
            ++fails;
        }
    };

    double s = 1e-9;
    FMatrix B = block(A, alpha, beta, s);
    expect(f.update(B, s, rows), "first update factors");
    expect(relError(B, f, rhs) < 1e-14, "direct solve");

    // toggle contact on some rows, smooth FB-like changes on others, new scale
    for (int i = 3; i < n; i += 29) {
        std::swap(alpha[i], beta[i]);
        beta[i] *= 1e9;
    }
    for (int i = 5; i < n; i += 41) {
        alpha[i] = -0.6;
        beta[i] = 0.3 * std::abs(A[i][i]);
    }
    s = 2.5e-9;
    B = block(A, alpha, beta, s);
    expect(!f.update(B, s, rows), "row change is a rank update");
    expect(f.rank() > 0, "rank > 0");
    expect(relError(B, f, rhs) < 1e-13, "Woodbury solve");

    // more rows on top of the existing update
    for (int i = 11; i < n; i += 53)
        alpha[i] = 0.5;
    B = block(A, alpha, beta, s);
    const int r0 = f.rank();
    expect(!f.update(B, s, rows), "second rank update");
    expect(f.rank() > r0, "rank grows");
    expect(relError(B, f, rhs) < 1e-13, "Woodbury solve after growth");

    // new elasticity matrix: refactor
    rows.generation = 2;
    expect(f.update(B, s, rows), "generation change refactors");
    expect(f.rank() == 0, "rank reset");
    expect(relError(B, f, rhs) < 1e-14, "solve after refactor");

    // too many rows: refactor
    for (int i = 0; i < n; i += 2)
        alpha[i] = 0.9 * alpha[i] + 0.01;
    B = block(A, alpha, beta, s);
    expect(f.update(B, s, rows), "rank cap refactors");
    expect(relError(B, f, rhs) < 1e-14, "solve after cap");
    return fails;
}
} // namespace

int
main()
{
    const int fails = check(true) + check(false);
    std::cout << (fails == 0 ? "test_dense_mech_factor: OK" : "test_dense_mech_factor: FAILED") << "\n";
    return fails == 0 ? 0 : 1;
}
