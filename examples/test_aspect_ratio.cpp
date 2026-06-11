/*
  Standalone test of the pure-mechanics element discretizations (VEM with the
  various stabilization choices, and Q1 FEM with/without B-bar) with emphasis on
  robustness for high-aspect-ratio cells.

  The test uses only the array-based vem:: API (no Dune grid required) and a small
  dense solver, so it runs in milliseconds and can serve as a regression test.

  Three experiments are run for every StabilityChoice variant:

  1. Patch test (3D and 2D): a linear displacement field imposed on the boundary of
     a small Cartesian grid with stretched cells must be reproduced exactly in the
     interior (first-order consistency).  Hard failure if violated.

  2. Single-element spectrum: a single hexahedron with increasing aspect ratio must
     have exactly 6 (near-)zero eigenvalues (the rigid body modes) and no spurious
     zero-energy modes.  Hard failure for the new variants (ANISO_DIAG,
     ANISO_HARMONIC, FEM, FEM_BBAR); legacy variants only report.

  3. Cantilever bending sweep: a clamped beam with fixed element count and
     increasing element aspect ratio, loaded by an end shear traction.  The tip
     deflection is compared with the analytic Timoshenko solution; the ratio
     (numeric/analytic) is reported per variant and aspect ratio.  This is a
     diagnostic: locking shows up as a ratio dropping far below 1 as the aspect
     ratio grows.  In addition, a near-incompressible (nu=0.499) variant verifies
     that FEM_BBAR suffers less from volumetric locking than plain FEM.
*/

#include <opm/geomech/vem/vem.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace
{

// ----------------------------------------------------------------------------
// dense linear solver: Gaussian elimination with partial pivoting
std::vector<double>
dense_solve(std::vector<double> A, std::vector<double> b)
{
    const int n = int(b.size());
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int row = col + 1; row < n; ++row)
            if (std::fabs(A[row * n + col]) > std::fabs(A[piv * n + col]))
                piv = row;
        if (piv != col) {
            for (int k = 0; k < n; ++k)
                std::swap(A[col * n + k], A[piv * n + k]);
            std::swap(b[col], b[piv]);
        }
        const double d = A[col * n + col];
        for (int row = col + 1; row < n; ++row) {
            const double f = A[row * n + col] / d;
            if (f == 0.0)
                continue;
            for (int k = col; k < n; ++k)
                A[row * n + k] -= f * A[col * n + k];
            b[row] -= f * b[col];
        }
    }
    std::vector<double> x(n);
    for (int row = n - 1; row >= 0; --row) {
        double s = b[row];
        for (int k = row + 1; k < n; ++k)
            s -= A[row * n + k] * x[k];
        x[row] = s / A[row * n + row];
    }
    return x;
}

// ----------------------------------------------------------------------------
// eigenvalues of a symmetric matrix by cyclic Jacobi rotations
std::vector<double>
sym_eigenvalues(std::vector<double> A, const int n)
{
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                off += A[i * n + j] * A[i * n + j];
        double diagnorm = 0.0;
        for (int i = 0; i < n; ++i)
            diagnorm += A[i * n + i] * A[i * n + i];
        if (off < 1e-28 * (diagnorm + 1e-300))
            break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                const double apq = A[p * n + q];
                if (std::fabs(apq) < 1e-300)
                    continue;
                const double theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0)
                    / (std::fabs(theta) + std::sqrt(theta * theta + 1));
                const double c = 1.0 / std::sqrt(t * t + 1);
                const double s = t * c;
                for (int k = 0; k < n; ++k) {
                    const double akp = A[k * n + p];
                    const double akq = A[k * n + q];
                    A[k * n + p] = c * akp - s * akq;
                    A[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const double apk = A[p * n + k];
                    const double aqk = A[q * n + k];
                    A[p * n + k] = c * apk - s * aqk;
                    A[q * n + k] = s * apk + c * aqk;
                }
            }
    }
    std::vector<double> eig(n);
    for (int i = 0; i < n; ++i)
        eig[i] = A[i * n + i];
    std::sort(eig.begin(), eig.end());
    return eig;
}

// ----------------------------------------------------------------------------
// Cartesian hex grid of nx x ny x nz cells with cell sizes dx, dy, dz, described
// in the face-based format expected by vem::assemble_mech_system_3D.  Faces of
// each cell are listed in the order [x-, x+, y-, y+, z-, z+], oriented outwards.
struct HexGrid {
    int nx, ny, nz;
    double dx, dy, dz;
    std::vector<double> points;
    std::vector<int> num_cell_faces, num_face_corners, face_corners;

    HexGrid(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_)
        : nx(nx_), ny(ny_), nz(nz_), dx(dx_), dy(dy_), dz(dz_)
    {
        for (int k = 0; k <= nz; ++k)
            for (int j = 0; j <= ny; ++j)
                for (int i = 0; i <= nx; ++i) {
                    points.push_back(i * dx);
                    points.push_back(j * dy);
                    points.push_back(k * dz);
                }
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int n000 = node(i, j, k), n100 = node(i + 1, j, k);
                    const int n010 = node(i, j + 1, k), n110 = node(i + 1, j + 1, k);
                    const int n001 = node(i, j, k + 1), n101 = node(i + 1, j, k + 1);
                    const int n011 = node(i, j + 1, k + 1), n111 = node(i + 1, j + 1, k + 1);
                    num_cell_faces.push_back(6);
                    const int faces[6][4] = {{n000, n001, n011, n010},  // x-
                                             {n100, n110, n111, n101},  // x+
                                             {n000, n100, n101, n001},  // y-
                                             {n010, n011, n111, n110},  // y+
                                             {n000, n010, n110, n100},  // z-
                                             {n001, n101, n111, n011}}; // z+
                    for (int f = 0; f < 6; ++f) {
                        num_face_corners.push_back(4);
                        for (int c = 0; c < 4; ++c)
                            face_corners.push_back(faces[f][c]);
                    }
                }
    }
    int node(int i, int j, int k) const
    {
        return i + (nx + 1) * (j + (ny + 1) * k);
    }
    int num_nodes() const
    {
        return (nx + 1) * (ny + 1) * (nz + 1);
    }
    int num_cells() const
    {
        return nx * ny * nz;
    }
    bool boundary_node(int i, int j, int k) const
    {
        return i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz;
    }
};

const std::vector<std::pair<vem::StabilityChoice, std::string>> variants
    = {{vem::SIMPLE, "SIMPLE"},
       {vem::HARMONIC, "HARMONIC"},
       {vem::D_RECIPE, "D_RECIPE"},
       {vem::ANISO_DIAG, "ANISO_DIAG"},
       {vem::ANISO_HARMONIC, "ANISO_HARM"},
       {vem::FEM, "FEM"},
       {vem::FEM_BBAR, "FEM_BBAR"}};

// ----------------------------------------------------------------------------
// expand reduced solution to full dof vector
std::vector<double>
expand_solution(const std::vector<double>& reduced,
                const std::vector<int>& fixed_ixs,
                const std::vector<double>& fixed_vals,
                const int num_dofs)
{
    std::vector<double> full(num_dofs);
    std::size_t fi = 0, ri = 0;
    for (int d = 0; d < num_dofs; ++d) {
        if (fi < fixed_ixs.size() && fixed_ixs[fi] == d)
            full[d] = fixed_vals[fi++];
        else
            full[d] = reduced[ri++];
    }
    return full;
}

// ----------------------------------------------------------------------------
// 3D patch test: returns max abs error (relative to field magnitude) at interior nodes
double
patch_test_3D(const vem::StabilityChoice sc, const double aspect, const bool diag_scaling)
{
    const HexGrid g(3, 3, 3, 1.0, 1.0, 1.0 / aspect);
    const double E = 1e9, nu = 0.3;
    std::vector<double> young(g.num_cells(), E), poisson(g.num_cells(), nu);
    std::vector<double> body(3 * g.num_cells(), 0.0);

    // linear displacement field u = A x
    const double A[3][3] = {{1e-2, 5e-3, 3e-3}, {2e-3, -4e-3, 6e-3}, {1e-3, 2e-3, 3e-3}};
    auto exact = [&](const double* x, int d) {
        return A[d][0] * x[0] + A[d][1] * x[1] + A[d][2] * x[2];
    };

    std::vector<int> fixed_ixs;
    std::vector<double> fixed_vals;
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int i = 0; i <= g.nx; ++i) {
                if (!g.boundary_node(i, j, k))
                    continue;
                const int n = g.node(i, j, k);
                for (int d = 0; d < 3; ++d) {
                    fixed_ixs.push_back(3 * n + d);
                    fixed_vals.push_back(exact(&g.points[3 * n], d));
                }
            }

    std::vector<std::tuple<int, int, double>> A_entries;
    std::vector<double> b;
    const int n_unknown = vem::assemble_mech_system_3D(&g.points[0],
                                                       g.num_cells(),
                                                       &g.num_cell_faces[0],
                                                       &g.num_face_corners[0],
                                                       &g.face_corners[0],
                                                       &young[0],
                                                       &poisson[0],
                                                       &body[0],
                                                       int(fixed_ixs.size()),
                                                       &fixed_ixs[0],
                                                       &fixed_vals[0],
                                                       0,
                                                       nullptr,
                                                       nullptr,
                                                       A_entries,
                                                       b,
                                                       sc,
                                                       true);
    std::vector<double> scale;
    if (diag_scaling)
        vem::diagonal_scale_system(A_entries, b, scale);

    auto x = dense_solve(vem::sparse2full(A_entries, n_unknown, n_unknown), b);
    if (diag_scaling)
        for (int i = 0; i < n_unknown; ++i)
            x[i] *= scale[i];
    const auto full = expand_solution(x, fixed_ixs, fixed_vals, 3 * g.num_nodes());

    double maxerr = 0.0, scale_u = 0.0;
    for (int n = 0; n < g.num_nodes(); ++n)
        for (int d = 0; d < 3; ++d) {
            const double ex = exact(&g.points[3 * n], d);
            maxerr = std::max(maxerr, std::fabs(full[3 * n + d] - ex));
            scale_u = std::max(scale_u, std::fabs(ex));
        }
    return maxerr / scale_u;
}

// ----------------------------------------------------------------------------
// 2D patch test on stretched quads
double
patch_test_2D(const vem::StabilityChoice sc, const double aspect)
{
    const int nx = 3, ny = 3;
    const double dx = 1.0, dy = 1.0 / aspect;
    std::vector<double> points;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            points.push_back(i * dx);
            points.push_back(j * dy);
        }
    auto node = [&](int i, int j) { return i + (nx + 1) * j; };
    std::vector<int> num_cell_faces, cell_corners;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            num_cell_faces.push_back(4);
            // counterclockwise
            cell_corners.push_back(node(i, j));
            cell_corners.push_back(node(i + 1, j));
            cell_corners.push_back(node(i + 1, j + 1));
            cell_corners.push_back(node(i, j + 1));
        }
    const int num_cells = nx * ny;
    const int num_nodes = (nx + 1) * (ny + 1);
    std::vector<double> young(num_cells, 1e9), poisson(num_cells, 0.3);
    std::vector<double> body(2 * num_cells, 0.0);

    const double A[2][2] = {{1e-2, 4e-3}, {3e-3, -5e-3}};
    auto exact = [&](const double* x, int d) { return A[d][0] * x[0] + A[d][1] * x[1]; };

    std::vector<int> fixed_ixs;
    std::vector<double> fixed_vals;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            if (i != 0 && i != nx && j != 0 && j != ny)
                continue;
            const int n = node(i, j);
            for (int d = 0; d < 2; ++d) {
                fixed_ixs.push_back(2 * n + d);
                fixed_vals.push_back(exact(&points[2 * n], d));
            }
        }

    std::vector<std::tuple<int, int, double>> A_entries;
    std::vector<double> b;
    const int n_unknown = vem::assemble_mech_system_2D(&points[0],
                                                       num_cells,
                                                       &num_cell_faces[0],
                                                       &cell_corners[0],
                                                       &young[0],
                                                       &poisson[0],
                                                       &body[0],
                                                       int(fixed_ixs.size()),
                                                       &fixed_ixs[0],
                                                       &fixed_vals[0],
                                                       0,
                                                       nullptr,
                                                       nullptr,
                                                       A_entries,
                                                       b,
                                                       sc,
                                                       true);
    const auto x = dense_solve(vem::sparse2full(A_entries, n_unknown, n_unknown), b);
    const auto full = expand_solution(x, fixed_ixs, fixed_vals, 2 * num_nodes);

    double maxerr = 0.0, scale_u = 0.0;
    for (int n = 0; n < num_nodes; ++n)
        for (int d = 0; d < 2; ++d) {
            const double ex = exact(&points[2 * n], d);
            maxerr = std::max(maxerr, std::fabs(full[2 * n + d] - ex));
            scale_u = std::max(scale_u, std::fabs(ex));
        }
    return maxerr / scale_u;
}

// ----------------------------------------------------------------------------
// single-element spectrum: returns number of (near-)zero eigenvalues, and the
// stiffness ratio between smallest nonzero and largest eigenvalue
std::pair<int, double>
single_element_spectrum(const vem::StabilityChoice sc, const double aspect)
{
    const HexGrid g(1, 1, 1, aspect, 1.0, 1.0);
    std::array<double, 3> centroid;
    std::vector<int> indexing;
    std::vector<double> K;
    vem::assemble_stiffness_matrix_3D(&g.points[0],
                                      &g.face_corners[0],
                                      &g.num_face_corners[0],
                                      6,
                                      1e9,
                                      0.3,
                                      sc,
                                      centroid,
                                      indexing,
                                      K);
    const auto eig = sym_eigenvalues(K, 24);
    const double emax = eig.back();
    int nzero = 0;
    double emin_nonzero = emax;
    for (const double e : eig) {
        if (std::fabs(e) < 1e-9 * emax)
            ++nzero;
        else
            emin_nonzero = std::min(emin_nonzero, std::fabs(e));
    }
    return {nzero, emin_nonzero / emax};
}

// ----------------------------------------------------------------------------
// cantilever bending: returns delta_numeric / delta_timoshenko
double
cantilever_ratio(const vem::StabilityChoice sc, const double length, const double nu)
{
    const int nx = 10, ny = 2, nz = 2;
    const double W = 1.0, H = 1.0;
    const HexGrid g(nx, ny, nz, length / nx, W / ny, H / nz);
    const double E = 1e9;
    std::vector<double> young(g.num_cells(), E), poisson(g.num_cells(), nu);
    std::vector<double> body(3 * g.num_cells(), 0.0);

    // clamp x=0 face
    std::vector<int> fixed_ixs;
    std::vector<double> fixed_vals;
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int d = 0; d < 3; ++d) {
                fixed_ixs.push_back(3 * g.node(0, j, k) + d);
                fixed_vals.push_back(0.0);
            }
    std::sort(fixed_ixs.begin(), fixed_ixs.end());

    // end shear traction on the x+ faces of the last column of cells
    const double P = 1e3; // total force
    const double tz = -P / (W * H);
    std::vector<int> neumann_faces;
    std::vector<double> neumann_forces;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
            const int cell = (nx - 1) + nx * (j + ny * k);
            neumann_faces.push_back(6 * cell + 1); // x+ face
            neumann_forces.push_back(0.0);
            neumann_forces.push_back(0.0);
            neumann_forces.push_back(tz);
        }

    std::vector<std::tuple<int, int, double>> A_entries;
    std::vector<double> b;
    const int n_unknown = vem::assemble_mech_system_3D(&g.points[0],
                                                       g.num_cells(),
                                                       &g.num_cell_faces[0],
                                                       &g.num_face_corners[0],
                                                       &g.face_corners[0],
                                                       &young[0],
                                                       &poisson[0],
                                                       &body[0],
                                                       int(fixed_ixs.size()),
                                                       &fixed_ixs[0],
                                                       &fixed_vals[0],
                                                       int(neumann_faces.size()),
                                                       &neumann_faces[0],
                                                       &neumann_forces[0],
                                                       A_entries,
                                                       b,
                                                       sc,
                                                       true);
    const auto x = dense_solve(vem::sparse2full(A_entries, n_unknown, n_unknown), b);
    const auto full = expand_solution(x, fixed_ixs, fixed_vals, 3 * g.num_nodes());

    // average z-deflection of the tip face nodes
    double tip = 0.0;
    int count = 0;
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j) {
            tip += full[3 * g.node(nx, j, k) + 2];
            ++count;
        }
    tip /= count;

    // Timoshenko cantilever with shear correction (rectangular section)
    const double I = W * H * H * H / 12.0;
    const double Ar = W * H;
    const double G = E / (2 * (1 + nu));
    const double kappa = 10.0 * (1 + nu) / (12.0 + 11.0 * nu);
    const double delta = P * std::pow(length, 3) / (3.0 * E * I) + P * length / (kappa * G * Ar);
    return -tip / delta;
}

} // anonymous namespace

// ----------------------------------------------------------------------------
int
main()
{
    int failures = 0;
    const std::vector<double> aspects = {1.0, 10.0, 100.0};

    std::printf("=== 3D patch test (linear field reproduced exactly), rel error ===\n");
    std::printf("%-12s", "variant");
    for (const double a : aspects)
        std::printf("  aspect=%-8g", a);
    std::printf("\n");
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double a : aspects) {
            const double err = patch_test_3D(sc, a, false);
            std::printf("  %-15.3e", err);
            if (err > 1e-6) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }
    // patch test must also survive diagonal scaling of the system
    {
        const double err = patch_test_3D(vem::ANISO_DIAG, 100.0, true);
        std::printf("%-12s  %-15.3e (with diagonal scaling)\n", "ANISO_DIAG", err);
        if (err > 1e-6)
            ++failures;
    }

    std::printf("\n=== 2D patch test, rel error ===\n");
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double a : aspects) {
            const double err = patch_test_2D(sc, a);
            std::printf("  %-15.3e", err);
            if (err > 1e-6) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }

    std::printf("\n=== single element: #zero eigenvalues (must be 6) and min/max stiffness ===\n");
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double a : aspects) {
            const auto [nzero, ratio] = single_element_spectrum(sc, a);
            std::printf("  %d/%-10.2e", nzero, ratio);
            const bool is_new = (sc == vem::ANISO_DIAG || sc == vem::ANISO_HARMONIC
                                 || sc == vem::FEM || sc == vem::FEM_BBAR);
            if (nzero != 6) {
                std::printf(is_new ? "(FAIL)" : "(warn)");
                if (is_new)
                    ++failures;
            }
        }
        std::printf("\n");
    }

    std::printf("\n=== cantilever bending, tip deflection / Timoshenko (nu=0.25) ===\n");
    std::printf("(10x2x2 elements; element aspect ratio = L/5)\n");
    const std::vector<double> lengths = {5.0, 25.0, 125.0, 250.0};
    std::printf("%-12s", "variant");
    for (const double L : lengths)
        std::printf("  aspect=%-6g", L / 5.0);
    std::printf("\n");
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double L : lengths) {
            const double r = cantilever_ratio(sc, L, 0.25);
            std::printf("  %-13.3f", r);
            if (!(r > 0.0 && r < 1.5)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }

    std::printf("\n=== volumetric locking: cantilever at nu=0.4999, ratio to Timoshenko ===\n");
    std::printf("(volumetric locking adds to shear locking as the aspect ratio grows;\n");
    std::printf(" B-bar removes the volumetric part: FEM_BBAR must beat FEM everywhere)\n");
    std::printf("%-12s", "variant");
    for (const double L : lengths)
        std::printf("  aspect=%-6g", L / 5.0);
    std::printf("\n");
    std::map<std::string, std::vector<double>> vol_ratios;
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double L : lengths) {
            const double r = cantilever_ratio(sc, L, 0.4999);
            vol_ratios[name].push_back(r);
            std::printf("  %-13.3f", r);
            if (!(r > 0.0 && r < 1.5)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }
    // hard checks: B-bar must reduce the volumetric locking of plain Q1 at every
    // aspect ratio, and must essentially cure it where shear locking is absent
    for (std::size_t i = 0; i < lengths.size(); ++i) {
        if (!(vol_ratios["FEM_BBAR"][i] > vol_ratios["FEM"][i])) {
            std::printf("FAIL: FEM_BBAR (%.3f) not softer than FEM (%.3f) at aspect %g\n",
                        vol_ratios["FEM_BBAR"][i], vol_ratios["FEM"][i], lengths[i] / 5.0);
            ++failures;
        }
    }
    if (!(vol_ratios["FEM_BBAR"][0] > 0.7)) {
        std::printf("FAIL: FEM_BBAR should cure volumetric locking at aspect 1 (got %.3f)\n",
                    vol_ratios["FEM_BBAR"][0]);
        ++failures;
    }

    std::printf("\n%s (%d failures)\n", failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED", failures);
    return failures == 0 ? 0 : 1;
}
