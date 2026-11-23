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

// NearFieldMechILU: ILU(0) is exact when the pattern holds every nonzero.

#include <config.h>

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include <opm/geomech/FractureMechanicsPreconditioner.hpp>

namespace
{
using Opm::FMatrix;
using Opm::Vector;

// cells on a line with unit spacing; radius in spacings
int
checkNearField()
{
    const int n = 200;
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> xyz(3 * n, 0.0);
    for (int i = 0; i < n; ++i)
        xyz[3 * i] = i;
    std::vector<double> alpha(n, 1.0), beta(n, 0.0);
    Vector rhs(n);
    for (int i = 0; i < n; ++i)
        rhs[i] = u(gen);
    int fails = 0;
    const auto residual = [&](const FMatrix& B, const Opm::NearFieldMechILU& p) {
        Vector x(n), r(rhs);
        p.solve(x, rhs);
        B.mmv(x, r);
        return r.infinity_norm() / (B.infinity_norm() * x.infinity_norm() + rhs.infinity_norm());
    };
    FMatrix A(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = (i == j) ? 5.0 + u(gen) : u(gen) / (1 + std::abs(i - j));
    Opm::MechRows rows {&A, &alpha, &beta, 1, &xyz};

    Opm::NearFieldMechILU full(10.0 * n);
    full.update(A, rows);
    if (residual(A, full) > 1e-13) {
        std::cerr << "nearfield: full pattern is not exact\n";
        ++fails;
    }
    FMatrix T(n, n);
    T = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = std::max(0, i - 1); j <= std::min(n - 1, i + 1); ++j)
            T[i][j] = A[i][j];
    Opm::NearFieldMechILU band(1.5);
    band.update(T, rows);
    if (residual(T, band) > 1e-13 || band.nnzPerRow() > 3.0) {
        std::cerr << "nearfield: tridiagonal pattern is not exact\n";
        ++fails;
    }
    std::vector<double> bad(3 * (n - 1), 0.0);
    Opm::MechRows badRows {&A, &alpha, &beta, 2, &bad};
    try {
        band.update(A, badRows);
        std::cerr << "nearfield: wrong centre count accepted\n";
        ++fails;
    } catch (const std::exception&) {
    }
    return fails;
}
} // namespace

int
main()
{
    const int fails = checkNearField();
    std::cout << (fails == 0 ? "test_nearfield_mech_ilu: OK" : "test_nearfield_mech_ilu: FAILED") << "\n";
    return fails == 0 ? 0 : 1;
}
