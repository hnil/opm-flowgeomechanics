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
     SMOOTH_HOURGLASS must stay near 1 at every aspect, also on sheared cells.
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
       {vem::FEM_BBAR, "FEM_BBAR"},
       {vem::FACE_D_RECIPE, "FACE_D_REC"},
       {vem::FACE_ANISO_HARMONIC, "FACE_ANISO"},
       {vem::FACE_SIMPLE, "FACE_SIMPLE"},
       {vem::TRACE, "TRACE"},
       {vem::SMOOTH_HOURGLASS, "SMOOTH_HG"}};

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
// 3D patch test with every node displaced pseudo-randomly by up to `warp` cell sizes,
// so that almost no face is planar
double
patch_test_3D_warped(const vem::StabilityChoice sc, const double warp, const vem::FaceQuadrature fq)
{
    HexGrid g(3, 3, 3, 1.0, 1.0, 1.0);
    unsigned seed = 12345u;
    for (auto& x : g.points) {
        seed = 1664525u * seed + 1013904223u;
        x += warp * (double(seed >> 8) / double(1u << 24) - 0.5);
    }
    const double E = 1e9, nu = 0.3;
    std::vector<double> young(g.num_cells(), E), poisson(g.num_cells(), nu);
    std::vector<double> body(3 * g.num_cells(), 0.0);
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
    const auto saved = vem::face_quadrature();
    vem::set_face_quadrature(fq);
    std::vector<std::tuple<int, int, double>> A_entries;
    std::vector<double> b;
    const int n_unknown = vem::assemble_mech_system_3D(&g.points[0], g.num_cells(), &g.num_cell_faces[0],
                                                       &g.num_face_corners[0], &g.face_corners[0], &young[0],
                                                       &poisson[0], &body[0], int(fixed_ixs.size()),
                                                       &fixed_ixs[0], &fixed_vals[0], 0, nullptr, nullptr,
                                                       A_entries, b, sc, true);
    vem::set_face_quadrature(saved);
    const auto x = dense_solve(vem::sparse2full(A_entries, n_unknown, n_unknown), b);
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
// Unit cube whose x+ face is split at height s (a fault throw: neighbour corners on the two
// x+ pillar edges).  Returns the energy of a unit x-displacement of the fault node alone.
double
fault_node_energy(const vem::StabilityChoice sc, const double s)
{
    std::vector<double> pts;
    for (int k = 0; k < 8; ++k) {
        pts.push_back(k & 1);
        pts.push_back((k >> 1) & 1);
        pts.push_back((k >> 2) & 1);
    }
    pts.insert(pts.end(), {1.0, 0.0, s, 1.0, 1.0, s}); // nodes 8, 9
    const std::vector<std::vector<int>> faces = {{0, 4, 6, 2},     {0, 1, 8, 5, 4}, {2, 6, 7, 9, 3},
                                                 {0, 2, 3, 1},     {4, 5, 7, 6},    {1, 3, 9, 8},
                                                 {8, 9, 7, 5}};
    std::vector<int> fc, nfc;
    for (const auto& f : faces) {
        nfc.push_back(int(f.size()));
        fc.insert(fc.end(), f.begin(), f.end());
    }
    std::array<double, 3> centroid;
    std::vector<int> indexing;
    std::vector<double> K;
    vem::assemble_stiffness_matrix_3D(&pts[0], &fc[0], &nfc[0], int(faces.size()), 1e9, 0.3, sc, centroid,
                                      indexing, K);
    const int n = 3 * int(indexing.size());
    const int loc = int(std::find(indexing.begin(), indexing.end(), 8) - indexing.begin());
    return K[(3 * loc) * n + 3 * loc];
}

// ----------------------------------------------------------------------------
// Flat fault cell (a x a x 1, fault nodes at s = 1/2 on the two x+ pillars), top face twisted by
// twist_deg about the vertical axis. Returns the stiffness and the bending mode
// u_x = (2xi-1)(2zeta-1) (extra nodes interpolated); hex_corners = Dune-ordered corners or none.
std::vector<double>
fault_cell_K(const vem::StabilityChoice sc, const double a, const double twist_deg, const bool corners,
             std::vector<double>& ubend)
{
    std::vector<double> pts;
    const double t = twist_deg * M_PI / 180.0;
    for (int k = 0; k < 8; ++k) {
        double x = a * ((k & 1) - 0.5), y = a * (((k >> 1) & 1) - 0.5);
        const double z = (k >> 2) & 1;
        if (z > 0) {
            const double xr = std::cos(t) * x - std::sin(t) * y, yr = std::sin(t) * x + std::cos(t) * y;
            x = xr;
            y = yr;
        }
        pts.insert(pts.end(), {x, y, z});
    }
    for (const int c : {1, 3}) // nodes 8, 9: mid-points of the x+ pillars
        for (int d = 0; d < 3; ++d)
            pts.push_back(0.5 * (pts[3 * c + d] + pts[3 * (c + 4) + d]));
    const std::vector<std::vector<int>> faces = {{0, 4, 6, 2},     {0, 1, 8, 5, 4}, {2, 6, 7, 9, 3},
                                                 {0, 2, 3, 1},     {4, 5, 7, 6},    {1, 3, 9, 8},
                                                 {8, 9, 7, 5}};
    std::vector<int> fc, nfc;
    for (const auto& f : faces) {
        nfc.push_back(int(f.size()));
        fc.insert(fc.end(), f.begin(), f.end());
    }
    const int hc[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    std::array<double, 3> centroid;
    std::vector<int> indexing;
    std::vector<double> K;
    vem::assemble_stiffness_matrix_3D(&pts[0], &fc[0], &nfc[0], int(faces.size()), 1e9, 0.3, sc, centroid,
                                      indexing, K, corners ? hc : nullptr);
    ubend.assign(3 * indexing.size(), 0.0);
    for (std::size_t i = 0; i < indexing.size(); ++i) {
        const int v = indexing[i];
        const double xi = v < 8 ? (v & 1) : 1.0, zeta = v < 8 ? ((v >> 2) & 1) : 0.5;
        ubend[3 * i] = (2 * xi - 1) * (2 * zeta - 1);
    }
    return K;
}

// ----------------------------------------------------------------------------
// Cantilever (length L, 1 x 1 section, one cell through the thickness) made of nx columns, each
// either one hex or two wedges given as collapsed hexes (corner ids repeat). Tip load as nodal
// forces; returns tip deflection / Timoshenko.
double
wedge_cantilever_ratio(const vem::StabilityChoice sc, const double L, const int nx, const bool wedges,
                       const bool corners, const int only_column = -1)
{
    const auto node = [](int i, int j, int k) { return 4 * i + 2 * k + j; };
    std::vector<double> pts;
    for (int i = 0; i <= nx; ++i)
        for (int k = 0; k <= 1; ++k)
            for (int j = 0; j <= 1; ++j)
                pts.insert(pts.end(), {L * i / nx, double(j), double(k)});
    std::vector<std::array<int, 8>> cells; // Dune-ordered corners, repeated where collapsed
    for (int i = 0; i < nx; ++i) {
        std::array<int, 8> h;
        for (int c = 0; c < 8; ++c)
            h[c] = node(i + (c & 1), (c >> 1) & 1, (c >> 2) & 1);
        if (!wedges || (only_column >= 0 && i != only_column)) {
            cells.push_back(h);
            continue;
        }
        auto a = h, b = h;
        a[5] = h[1], a[7] = h[3]; // x+ top edge down: triangle (0,0) (a,0) (0,1)
        b[0] = h[4], b[2] = h[6]; // x- bottom edge up: triangle (a,0) (a,1) (0,1)
        cells.push_back(a);
        cells.push_back(b);
    }
    const int N = 3 * int(pts.size() / 3);
    std::vector<double> K(std::size_t(N) * N, 0.0), f(N, 0.0);
    for (const auto& h : cells) {
        // hex faces from the corners, dropping collapsed edges and degenerate faces; outward winding
        const int quad[6][4] = {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
        std::vector<int> fc, nfc;
        double cc[3] = {0, 0, 0};
        for (int c = 0; c < 8; ++c)
            for (int d = 0; d < 3; ++d)
                cc[d] += pts[3 * h[c] + d] / 8;
        for (const auto& q : quad) {
            std::vector<int> f2;
            for (int v = 0; v < 4; ++v)
                if (f2.empty() || (h[q[v]] != f2.back() && !(v == 3 && h[q[v]] == f2.front())))
                    f2.push_back(h[q[v]]);
            if (f2.size() < 3)
                continue;
            double nrm[3] = {0, 0, 0}, fcen[3] = {0, 0, 0};
            for (std::size_t v = 0; v < f2.size(); ++v) {
                const double* p = &pts[3 * f2[v]];
                const double* r = &pts[3 * f2[(v + 1) % f2.size()]];
                nrm[0] += (p[1] - r[1]) * (p[2] + r[2]);
                nrm[1] += (p[2] - r[2]) * (p[0] + r[0]);
                nrm[2] += (p[0] - r[0]) * (p[1] + r[1]);
                for (int d = 0; d < 3; ++d)
                    fcen[d] += p[d] / f2.size();
            }
            if (nrm[0] * (fcen[0] - cc[0]) + nrm[1] * (fcen[1] - cc[1]) + nrm[2] * (fcen[2] - cc[2]) < 0)
                std::reverse(f2.begin(), f2.end());
            nfc.push_back(int(f2.size()));
            fc.insert(fc.end(), f2.begin(), f2.end());
        }
        std::array<double, 3> centroid;
        std::vector<int> idx;
        std::vector<double> Ke;
        vem::assemble_stiffness_matrix_3D(&pts[0], &fc[0], &nfc[0], int(nfc.size()), 1e9, 0.25, sc, centroid,
                                          idx, Ke, corners ? &h[0] : nullptr);
        const int n = 3 * int(idx.size());
        for (int a2 = 0; a2 < n; ++a2)
            for (int b2 = 0; b2 < n; ++b2)
                K[std::size_t(3 * idx[a2 / 3] + a2 % 3) * N + 3 * idx[b2 / 3] + b2 % 3] += Ke[a2 * n + b2];
    }
    const double P = 1e3;
    for (int j = 0; j <= 1; ++j)
        for (int k = 0; k <= 1; ++k)
            f[3 * node(nx, j, k) + 2] = -P / 4;
    std::vector<int> fixed;
    for (int j = 0; j <= 1; ++j)
        for (int k = 0; k <= 1; ++k)
            for (int d = 0; d < 3; ++d)
                fixed.push_back(3 * node(0, j, k) + d);
    for (const int I : fixed) {
        for (int J = 0; J < N; ++J)
            K[std::size_t(I) * N + J] = K[std::size_t(J) * N + I] = 0.0;
        K[std::size_t(I) * N + I] = 1.0;
        f[I] = 0.0;
    }
    const auto u = dense_solve(K, f);
    double tip = 0;
    for (int j = 0; j <= 1; ++j)
        for (int k = 0; k <= 1; ++k)
            tip += u[3 * node(nx, j, k) + 2] / 4;
    const double E = 1e9, nu = 0.25, I = 1.0 / 12, G = E / (2 * (1 + nu)), kap = 10 * (1 + nu) / (12 + 11 * nu);
    return -tip / (P * L * L * L / (3 * E * I) + P * L / (kap * G));
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
// energies of the bending u_x = (2xi-1)(2zeta-1) and twist u_x = (2eta-1)(2zeta-1) hourglass
// modes on one hex of size aspect x 1 x 1 (same nodal magnitudes, true energies differ ~aspect^2)
std::pair<double, double>
hourglass_energies(const vem::StabilityChoice sc, const double aspect)
{
    const HexGrid g(1, 1, 1, aspect, 1.0, 1.0);
    std::array<double, 3> centroid;
    std::vector<int> indexing;
    std::vector<double> K;
    vem::assemble_stiffness_matrix_3D(&g.points[0], &g.face_corners[0], &g.num_face_corners[0], 6, 1e9, 0.3,
                                      sc, centroid, indexing, K);
    const int n = 3 * int(indexing.size());
    std::vector<double> ub(n, 0.0), ut(n, 0.0);
    for (int a = 0; a < int(indexing.size()); ++a) {
        const double* x = &g.points[3 * indexing[a]];
        const double xi = x[0] / aspect, eta = x[1], zeta = x[2];
        ub[3 * a] = (2 * xi - 1) * (2 * zeta - 1);
        ut[3 * a] = (2 * eta - 1) * (2 * zeta - 1);
    }
    double eb = 0, et = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            eb += ub[i] * K[i * n + j] * ub[j];
            et += ut[i] * K[i * n + j] * ut[j];
        }
    return {eb, et};
}

// ----------------------------------------------------------------------------
// cantilever bending: returns delta_numeric / delta_timoshenko.  If a rotation
// matrix is given (row-major 3x3), the whole problem (grid, traction, deflection
// direction) is rigidly rotated; the result should be unchanged for any frame-
// invariant discretization.  tan_pillar inclines the cell pillars in the x-z plane.
double
cantilever_ratio(const vem::StabilityChoice sc,
                 const double length,
                 const double nu,
                 const std::array<double, 9>* rot = nullptr,
                 const double tan_pillar = 0.0)
{
    const int nx = 10, ny = 2, nz = 2;
    const double W = 1.0, H = 1.0;
    HexGrid g(nx, ny, nz, length / nx, W / ny, H / nz);
    for (std::size_t p = 0; p < g.points.size(); p += 3)
        g.points[p] += tan_pillar * (g.points[p + 2] - 0.5 * H);
    if (rot) {
        const auto& R = *rot;
        for (std::size_t p = 0; p < g.points.size(); p += 3) {
            const double x = g.points[p], y = g.points[p + 1], z = g.points[p + 2];
            for (int d = 0; d < 3; ++d)
                g.points[p + d] = R[3 * d] * x + R[3 * d + 1] * y + R[3 * d + 2] * z;
        }
    }
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
    const double tz = -P / (W * H) / std::sqrt(1.0 + tan_pillar * tan_pillar); // end face area grows with the incline
    double tvec[3] = {0.0, 0.0, tz};
    double ez[3] = {0.0, 0.0, 1.0};
    if (rot) {
        const auto& R = *rot;
        for (int d = 0; d < 3; ++d) {
            tvec[d] = R[3 * d + 2] * tz;
            ez[d] = R[3 * d + 2];
        }
    }
    std::vector<int> neumann_faces;
    std::vector<double> neumann_forces;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
            const int cell = (nx - 1) + nx * (j + ny * k);
            neumann_faces.push_back(6 * cell + 1); // x+ face
            neumann_forces.push_back(tvec[0]);
            neumann_forces.push_back(tvec[1]);
            neumann_forces.push_back(tvec[2]);
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

    // average deflection (along the rotated z direction) of the tip face nodes
    double tip = 0.0;
    int count = 0;
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j) {
            for (int d = 0; d < 3; ++d)
                tip += ez[d] * full[3 * g.node(nx, j, k) + d];
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

    // warped (non-planar) faces: only the triangulated face quadrature is exact
    std::printf("\n=== 3D patch test, warped faces: kite / triangulated face quadrature ===\n");
    for (const double w : {0.1, 0.3}) {
        for (const auto& [sc, name] : variants) {
            if (vem::is_fem_choice(sc))
                continue;
            const double ek = patch_test_3D_warped(sc, w, vem::FaceQuadrature::KITE);
            const double et = patch_test_3D_warped(sc, w, vem::FaceQuadrature::TRIANGULATED);
            std::printf("warp=%-4g %-12s  %-12.3e %-12.3e", w, name.c_str(), ek, et);
            if (et > 1e-6) {
                ++failures;
                std::printf("(FAIL)");
            }
            std::printf("\n");
        }
    }

    // moving a node on a short edge must cost more as the edge shrinks (gradient ~ 1/s); a
    // stabilization scaled by the cell instead of the faces leaves a near-free mode there
    std::printf("\n=== fault cell split at s: energy of moving the fault node, s=1/32 relative to s=1/2 ===\n");
    for (const auto sc : {vem::SIMPLE, vem::D_RECIPE, vem::ANISO_HARMONIC, vem::TRACE, vem::SMOOTH_HOURGLASS}) {
        const double r = fault_node_energy(sc, 1.0 / 32) / fault_node_energy(sc, 0.5);
        std::printf("%-12s  %.2f", sc == vem::TRACE ? "TRACE" : sc == vem::SIMPLE ? "SIMPLE"
                                   : sc == vem::D_RECIPE ? "D_RECIPE" : sc == vem::ANISO_HARMONIC ? "ANISO_HARM"
                                   : "SMOOTH_HG", r);
        if ((sc == vem::TRACE || sc == vem::SMOOTH_HOURGLASS) && r < 4.0) {
            ++failures;
            std::printf(" (FAIL)");
        }
        std::printf("\n");
    }

    // the grid's corner-point corners replace the face-merge guess; it must agree where the guess
    // works and keep the smooth-field term where it fails (strongly bent faces)
    std::printf("\n=== SMOOTH_HG with given cell corners: flat fault cell 9x9x1 ===\n");
    {
        std::vector<double> u;
        const auto energy = [&u](const std::vector<double>& K) {
            const std::size_t n = u.size();
            double e = 0;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    e += u[i] * K[i * n + j] * u[j];
            return e;
        };
        const auto Kg = fault_cell_K(vem::SMOOTH_HOURGLASS, 9.0, 0.0, true, u);
        const auto Km = fault_cell_K(vem::SMOOTH_HOURGLASS, 9.0, 0.0, false, u);
        double dmax = 0, kmax = 0;
        for (std::size_t i = 0; i < Kg.size(); ++i) {
            dmax = std::max(dmax, std::fabs(Kg[i] - Km[i]));
            kmax = std::max(kmax, std::fabs(Kg[i]));
        }
        std::printf("untwisted: given vs merged corners, max |dK|/max|K| = %.2e", dmax / kmax);
        if (!(dmax < 1e-10 * kmax)) {
            ++failures;
            std::printf(" (FAIL)");
        }
        const auto Kt = fault_cell_K(vem::TRACE, 9.0, 30.0, false, u);
        const auto Kn = fault_cell_K(vem::SMOOTH_HOURGLASS, 9.0, 30.0, false, u);
        const auto Kc = fault_cell_K(vem::SMOOTH_HOURGLASS, 9.0, 30.0, true, u);
        const int nzero = [&] {
            const auto eig = sym_eigenvalues(Kc, int(u.size()));
            int z = 0;
            for (const double e : eig)
                z += std::fabs(e) < 1e-9 * eig.back();
            return z;
        }();
        std::printf("\ntwisted 30 deg (faces bent past the merge angle): bending energy / TRACE: merged %.3f, given corners %.4f; zero modes %d",
                    energy(Kn) / energy(Kt), energy(Kc) / energy(Kt), nzero);
        if (!(std::fabs(energy(Kn) / energy(Kt) - 1) < 1e-12 && energy(Kc) < 0.75 * energy(Kt) && nzero == 6)) {
            ++failures;
            std::printf(" (FAIL)");
        }
        std::printf("\n");
    }

    // collapsed hexes (wedges) lock in their constant-strain consistency term whatever S is; the
    // smooth-field term must keep the hexes around one wedge column accurate
    std::printf("\n=== cantilever of wedges (collapsed hexes), aspect 1 / 5 / 25: all wedges; one wedge column ===\n");
    for (const auto sc : {vem::TRACE, vem::SMOOTH_HOURGLASS}) {
        std::printf("%-12s", sc == vem::TRACE ? "TRACE" : "SMOOTH_HG");
        for (const double L : {10.0, 50.0, 250.0})
            std::printf("  %-7.3f", wedge_cantilever_ratio(sc, L, 10, true, true));
        std::printf(" |");
        for (const double L : {10.0, 50.0, 250.0}) {
            const double r = wedge_cantilever_ratio(sc, L, 10, true, true, 5);
            std::printf("  %-7.3f", r);
            if (sc == vem::SMOOTH_HOURGLASS && !(r > 0.85 && r < 1.1)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
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
                                 || sc == vem::FEM || sc == vem::FEM_BBAR || sc == vem::TRACE
                                 || sc == vem::SMOOTH_HOURGLASS);
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
    std::map<std::string, std::vector<double>> bend_ratios;
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double L : lengths) {
            const double r = cantilever_ratio(sc, L, 0.25);
            bend_ratios[name].push_back(r);
            std::printf("  %-13.3f", r);
            if (!(r > 0.0 && r < 1.5)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }

    // the smooth hourglass term must not lock in bending: beat ANISO_HARMONIC at aspect 25 and 50
    for (std::size_t i = 2; i < lengths.size(); ++i)
        if (!(bend_ratios["SMOOTH_HG"][i] > bend_ratios["ANISO_HARM"][i] && bend_ratios["SMOOTH_HG"][i] > 0.8)) {
            std::printf("FAIL: SMOOTH_HG (%.3f) locks at aspect %g (ANISO_HARM %.3f)\n",
                        bend_ratios["SMOOTH_HG"][i], lengths[i] / 5.0, bend_ratios["ANISO_HARM"][i]);
            ++failures;
        }

    // twist must keep its stiffness on thin cells (ANISO_HARM/SIMPLE: ~1/aspect^2, spurious
    // hourglassing); bending must not stiffen like Q1 (~aspect^2). Energies relative to Q1.
    std::printf("\n=== hourglass modes on one hex, energy / Q1 energy: bending, twist ===\n");
    std::printf("%-12s", "variant");
    for (const double a : aspects)
        std::printf("  aspect=%-14g", a);
    std::printf("\n");
    for (const auto& [sc, name] : variants) {
        if (vem::is_fem_choice(sc))
            continue;
        std::printf("%-12s", name.c_str());
        for (const double a : aspects) {
            const auto [eb, et] = hourglass_energies(sc, a);
            const auto [qb, qt] = hourglass_energies(vem::FEM, a);
            std::printf("  %-9.3g %-9.3g", eb / qb, et / qt);
            if (sc == vem::SMOOTH_HOURGLASS && (et / qt < 0.9 || et / qt > 1.1 || eb / qb > 1.1)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }

    // inclined pillars (sheared cells): Q1-like terms lock further, the smooth field must not
    std::printf("\n=== cantilever on sheared cells, aspect 25: pillar angle 10 / 20 / 30 / 45 deg ===\n");
    for (const auto& [sc, name] : variants) {
        std::printf("%-12s", name.c_str());
        for (const double deg : {10.0, 20.0, 30.0, 45.0}) {
            const double r = cantilever_ratio(sc, 125.0, 0.25, nullptr, std::tan(deg * M_PI / 180.0));
            std::printf("  %-8.3f", r);
            if (sc == vem::SMOOTH_HOURGLASS && !(r > 0.9 && r < 1.1)) {
                ++failures;
                std::printf("(FAIL)");
            }
        }
        std::printf("\n");
    }

    std::printf("\n=== frame invariance: rotated cantilever (aspect 25, nu=0.25) ===\n");
    std::printf("(grid + load rigidly rotated by Ry(30)*Rz(20); difference in deflection ratio)\n");
    {
        // R = Ry(30 deg) * Rz(20 deg), row-major
        const double a = 30.0 * M_PI / 180.0, b = 20.0 * M_PI / 180.0;
        const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b);
        const std::array<double, 9> R = {ca * cb, -ca * sb, sa,
                                         sb,      cb,       0.0,
                                         -sa * cb, sa * sb, ca};
        const double L = 125.0;
        for (const auto& [sc, name] : variants) {
            const double r0 = cantilever_ratio(sc, L, 0.25);
            const double r1 = cantilever_ratio(sc, L, 0.25, &R);
            const double diff = std::fabs(r1 - r0);
            std::printf("%-12s  %.6f -> %.6f   |diff| = %.2e", name.c_str(), r0, r1, diff);
            // hard requirement for the frame-invariant-by-construction variants;
            // the tolerance is limited by the conditioning of the solved system
            // (~1e-7 at this aspect ratio), not by the element formulation
            const bool must_be_invariant
                = (sc == vem::ANISO_HARMONIC || sc == vem::FEM || sc == vem::FEM_BBAR || sc == vem::TRACE
                   || sc == vem::SMOOTH_HOURGLASS);
            if (must_be_invariant && diff > 1e-6) {
                ++failures;
                std::printf("  (FAIL)");
            }
            std::printf("\n");
        }
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
