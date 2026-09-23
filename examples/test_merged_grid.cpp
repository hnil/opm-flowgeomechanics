// Check the merged grid: volumes, and whether each cell's oriented faces close.
#include "config.h"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

Opm::Coarsening::Grdecl tensorGrid(int nx, int ny, int nz, double dx, double dz)
{
    Opm::Coarsening::Grdecl g;
    g.dims = {nx, ny, nz};
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = i*dx, y = j*dx;
            for (double z : {0.0, nz*dz}) {
                g.coord.insert(g.coord.end(), {x, y, z});
            }
        }
    }
    g.zcorn.assign(8ull*nx*ny*nz, 0.0);
    for (int k = 0; k < nz; ++k) {
        for (int dk = 0; dk < 2; ++dk) {
            for (int j = 0; j < 2*ny; ++j) {
                for (int i = 0; i < 2*nx; ++i) {
                    g.zcorn[i + 2ull*nx*j + 4ull*nx*ny*(2*k + dk)] = (k + dk)*dz;
                }
            }
        }
    }
    g.actnum.assign(1ull*nx*ny*nz, 1);
    return g;
}

/// Minimal grdecl reader: SPECGRID/COORD/ZCORN/ACTNUM.
Opm::Coarsening::Grdecl readGrdecl(const std::string& path)
{
    std::ifstream is(path);
    Opm::Coarsening::Grdecl g;
    std::string token;
    const auto read = [&is](auto& out) {
        std::string t;
        while (is >> t) {
            if (t == "/") { return; }
            if (t.rfind("--", 0) == 0) { std::getline(is, t); continue; }
            if (t == "F" || t == "T") { continue; }
            const auto star = t.find('*');
            if (star != std::string::npos) {
                const int n = std::stoi(t.substr(0, star));
                out.insert(out.end(), n,
                           static_cast<typename std::decay_t<decltype(out)>::value_type>(
                               std::stod(t.substr(star + 1))));
            } else {
                out.push_back(static_cast<typename std::decay_t<decltype(out)>::value_type>(
                                  std::stod(t)));
            }
        }
    };
    while (is >> token) {
        if (token.rfind("--", 0) == 0) { std::getline(is, token); }
        else if (token == "SPECGRID" || token == "DIMENS") {
            std::vector<int> v; read(v); g.dims = {v.at(0), v.at(1), v.at(2)};
        }
        else if (token == "COORD")  { read(g.coord); }
        else if (token == "ZCORN")  { read(g.zcorn); }
        else if (token == "ACTNUM") { read(g.actnum); }
    }
    return g;
}

std::vector<Opm::Coarsening::CoarsenRequest> readRecords(const std::string& path)
{
    std::ifstream is(path);
    std::vector<Opm::Coarsening::CoarsenRequest> out;
    std::string line;
    while (std::getline(is, line)) {
        const auto comment = line.find("--");
        if (comment != std::string::npos) { line.erase(comment); }
        std::istringstream ls(line);
        std::vector<int> v;
        for (int x; ls >> x;) { v.push_back(x); }
        if (v.size() != 9) { continue; }
        Opm::Coarsening::CoarsenRequest r;
        for (int d = 0; d < 3; ++d) {
            r.startIJK[d] = v[2*d] - 1;
            r.endIJK[d] = v[2*d + 1];
            r.cellsPerDim[d] = v[6 + d];
        }
        out.push_back(r);
    }
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    // With two arguments: a grdecl include file and a file of COARSEN records.
    // Without: a small built-in case, the top two layers coarsened 2x2, which
    // cannot be written as a grdecl.
    Opm::Coarsening::Grdecl fine;
    std::vector<Opm::Coarsening::CoarsenRequest> requests;
    if (argc > 2) {
        fine = readGrdecl(argv[1]);
        requests = readRecords(argv[2]);
    } else {
        fine = tensorGrid(4, 4, 4, 100.0, 10.0);
        Opm::Coarsening::CoarsenRequest r;
        r.startIJK = {0, 0, 0};
        r.endIJK = {4, 4, 2};
        r.cellsPerDim = {2, 2, 1};
        requests = {r};
    }
    const int nx = fine.dims[0], ny = fine.dims[1], nz = fine.dims[2];
    if (fine.actnum.empty()) { fine.actnum.assign(1ull*nx*ny*nz, 1); }
    const auto layout = Opm::Coarsening::blockLayout(fine.dims, requests);
    std::cout << "blocks: " << layout.boxes.size() << '\n';

    grdecl input{};
    input.dims[0] = nx; input.dims[1] = ny; input.dims[2] = nz;
    input.coord = fine.coord.data();
    input.zcorn = fine.zcorn.data();
    input.actnum = fine.actnum.data();

    Dune::CpGrid grid;
    grid.processEclipseFormatCoarsened(input, layout.blockOfCartesian, layout.boxes, true);

    const auto& gv = grid.leafGridView();
    std::cout << "cells: " << gv.size(0) << "  nodes: " << gv.size(3) << '\n';

    double total = 0.0;
    int openCells = 0, worstFaces = 0;
    for (const auto& cell : elements(gv)) {
        total += cell.geometry().volume();
        std::array<double,3> sum{0.0, 0.0, 0.0};
        int faces = 0;
        for (const auto& is : intersections(gv, cell)) {
            const auto n = is.centerUnitOuterNormal();
            const double a = is.geometry().volume();
            for (int d = 0; d < 3; ++d) {
                sum[d] += n[d]*a;
            }
            ++faces;
        }
        worstFaces = std::max(worstFaces, faces);
        const double scale = std::pow(cell.geometry().volume(), 2.0/3.0);
        const double close = std::sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2])/scale;
        if (close > 1e-8) {
            if (openCells < 5) {
                std::cout << "cell " << gv.indexSet().index(cell) << " does not close: "
                          << close << " with " << faces << " faces, volume "
                          << cell.geometry().volume() << '\n';
            }
            ++openCells;
        }
    }
    // Is the average of a cell's nodes inside it? VEM starts its star-point
    // search there.
    int outside = 0;
    for (const auto& cell : elements(gv)) {
        std::array<double,3> avg{0.0, 0.0, 0.0};
        int n = 0;
        for (int c = 0; c < cell.geometry().corners(); ++c) {
            const auto p = cell.geometry().corner(c);
            for (int d = 0; d < 3; ++d) { avg[d] += p[d]; }
            ++n;
        }
        for (int d = 0; d < 3; ++d) { avg[d] /= n; }
        for (const auto& is : intersections(gv, cell)) {
            const auto nrm = is.centerUnitOuterNormal();
            const auto ctr = is.geometry().center();
            double dot = 0.0;
            for (int d = 0; d < 3; ++d) { dot += (avg[d] - ctr[d])*nrm[d]; }
            if (dot > 1e-9*std::pow(cell.geometry().volume(), 1.0/3.0)) {
                if (outside < 5) {
                    std::cout << "cell " << gv.indexSet().index(cell)
                              << ": corner average is outside, past a face by " << dot << '\n';
                }
                ++outside;
                break;
            }
        }
    }
    std::cout << "cells whose corner average is outside: " << outside << '\n';

    // Replicate the node ordering the VEM assembly builds (vemutils), and
    // check that each face's polygon normal points away from the cell.
    int wrongWay = 0;
    for (const auto& cell : elements(gv)) {
        const int c = gv.indexSet().index(cell);
        const auto centre = cell.geometry().center();
        const int nf = grid.numCellFaces(c);
        for (int f = 0; f < nf; ++f) {
            const int face = grid.cellFace(c, f);
            const int nv = grid.numFaceVertices(face);
            if (nv < 3) { continue; }
            std::vector<std::array<double,3>> poly;
            const bool forward = grid.faceCell(face, 1) != c;
            for (int v = 0; v < nv; ++v) {
                const int idx = forward ? v : nv - 1 - v;
                const auto p = grid.vertexPosition(grid.faceVertex(face, idx));
                poly.push_back({p[0], p[1], p[2]});
            }
            std::array<double,3> mid{0,0,0};
            for (const auto& p : poly) {
                for (int d = 0; d < 3; ++d) { mid[d] += p[d]/poly.size(); }
            }
            std::array<double,3> nrm{0,0,0};
            for (std::size_t v = 0; v < poly.size(); ++v) {
                const auto& a = poly[v];
                const auto& b = poly[(v + 1) % poly.size()];
                const double u[3] = {a[0]-mid[0], a[1]-mid[1], a[2]-mid[2]};
                const double w[3] = {b[0]-mid[0], b[1]-mid[1], b[2]-mid[2]};
                nrm[0] += u[1]*w[2] - u[2]*w[1];
                nrm[1] += u[2]*w[0] - u[0]*w[2];
                nrm[2] += u[0]*w[1] - u[1]*w[0];
            }
            double dot = 0.0;
            for (int d = 0; d < 3; ++d) { dot += nrm[d]*(mid[d] - centre[d]); }
            if (dot < 0.0) {
                if (wrongWay < 5) {
                    std::cout << "cell " << c << " face " << face
                              << ": nodes ordered the wrong way round\n";
                }
                ++wrongWay;
            }
        }
    }
    std::cout << "faces whose node order points inwards: " << wrongWay << '\n';

    // Is cell_to_face consistent with face_to_cell?
    int notMine = 0, bothSides = 0;
    for (const auto& cell : elements(gv)) {
        const int c = gv.indexSet().index(cell);
        for (int f = 0; f < grid.numCellFaces(c); ++f) {
            const int face = grid.cellFace(c, f);
            const int from = grid.faceCell(face, 0), to = grid.faceCell(face, 1);
            if (from != c && to != c) { ++notMine; }
            if (from == c && to == c) { ++bothSides; }
        }
    }
    std::cout << "cell-face entries whose face does not know the cell: " << notMine
              << ", faces with the cell on both sides: " << bothSides << '\n';
    std::cout << "total volume: " << total << '\n';
    std::cout << "cells whose faces do not close: " << openCells
              << ", most faces on a cell: " << worstFaces << '\n';
    return openCells == 0 ? 0 : 1;
}
