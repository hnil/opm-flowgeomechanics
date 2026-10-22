/*
  Regression test for the penny-crack DDM assembly.

  Solves a penny-shaped crack problem on a small 2-layer triangular mesh and
  compares the computed crack-opening displacement (COD) against stored reference
  values.  The reference was generated with the direct solver at layers=2 (10
  cells, E=1e9 Pa, nu=0.25, p=1e6 Pa).

  The test detects silent changes in ddm::assembleMatrix_fast or in the mesh
  generation that would alter the numerical solution.
*/

#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <opm/geomech/DiscreteDisplacement.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

// Stored reference values for layers=2, E=1e9, nu=0.25, p=1e6.
// Generated with the direct (dense LU) solver; update here if the assembly
// intentionally changes.
static const std::vector<double> reference_opening = {
    0.00209105307267588,
    0.00279548562329015,
    0.00209105307267588,
    0.00209105307267588,
    0.00279548562329017,
    0.00307711615538572,
    0.00279548562329018,
    0.00209105307267587,
    0.00209105307267588,
    0.00209105307267586,
};

int main()
{
    const int    layers = 2;
    const double E  = 1.0e9;
    const double nu = 0.25;
    const double p  = 1.0e6;
    const double tol = 1.0e-6; // relative tolerance for each opening value

    Opm::RegularTrimesh trimesh(layers,
                                {0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.5, std::sqrt(3.0) / 2.0, 0.0},
                                {1.0, 1.0});

    auto [grid, cellmap, bcells] = trimesh.createDuneGrid(0, {}, false);
    const int nc = grid->leafGridView().size(0);

    if (nc != static_cast<int>(reference_opening.size())) {
        std::cerr << "FAILED: expected " << reference_opening.size()
                  << " cells, got " << nc << "\n";
        return 1;
    }

    Dune::DynamicVector<double> rhs(nc, -p);
    Dune::DynamicVector<double> opening(nc, 0.0);
    Dune::DynamicMatrix<double> A(nc, nc, 0.0);

    ddm::assembleMatrix_fast(A, E, nu, *grid);
    A.solve(opening, rhs);

    int failures = 0;
    std::cout << std::setprecision(15);
    for (int i = 0; i < nc; ++i) {
        const double ref = reference_opening[static_cast<std::size_t>(i)];
        const double rel_err = std::abs(opening[i] - ref) / std::abs(ref);
        const bool ok = rel_err <= tol;
        if (!ok) {
            std::cerr << "  FAIL  w[" << i << "]  got=" << opening[i]
                      << "  ref=" << ref
                      << "  rel_err=" << rel_err << "\n";
            ++failures;
        }
    }

    if (failures == 0) {
        std::cout << "PASSED: DDM assembly regression (layers=2, " << nc
                  << " cells)\n";
        return 0;
    } else {
        std::cerr << "FAILED: " << failures << " of " << nc
                  << " opening values deviate from reference\n";
        return 1;
    }
}
