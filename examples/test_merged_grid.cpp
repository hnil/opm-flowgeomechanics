// Check the merged grid: volumes, and whether each cell's oriented faces close.
#include "config.h"

#include <opm/grid/CpGrid.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
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

/// Write each face as a polygon, so a coarse cell's subdivided faces and the
/// nodes on them are visible. Cell data: how many nodes the face has, and the
/// two cells it separates.
void writeFaceVtu(const Dune::CpGrid& grid, const std::string& path)
{
    const auto& gv = grid.leafGridView();

    // Faces reachable from the cells, in order.
    std::vector<int> faces;
    std::vector<char> seen;
    for (const auto& cell : elements(gv)) {
        const int c = gv.indexSet().index(cell);
        for (int f = 0; f < grid.numCellFaces(c); ++f) {
            const int face = grid.cellFace(c, f);
            if (face >= static_cast<int>(seen.size())) {
                seen.resize(face + 1, 0);
            }
            if (!seen[face]) {
                seen[face] = 1;
                faces.push_back(face);
            }
        }
    }

    std::ofstream os(path);
    os << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
       << "  <UnstructuredGrid>\n";

    std::size_t nodes = 0;
    for (const int face : faces) {
        nodes += grid.numFaceVertices(face);
    }
    os << "    <Piece NumberOfPoints=\"" << nodes << "\" NumberOfCells=\"" << faces.size()
       << "\">\n";

    // One point per face corner: the polygons then draw their own nodes, and a
    // node shared by a coarse and a fine face is not silently welded.
    os << "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
          "format=\"ascii\">\n";
    for (const int face : faces) {
        for (int v = 0; v < grid.numFaceVertices(face); ++v) {
            const auto& p = grid.vertexPosition(grid.faceVertex(face, v));
            os << "          " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
    }
    os << "        </DataArray>\n      </Points>\n";

    os << "      <Cells>\n        <DataArray type=\"Int64\" Name=\"connectivity\" "
          "format=\"ascii\">\n";
    std::size_t next = 0;
    std::vector<std::size_t> offsets;
    for (const int face : faces) {
        os << "         ";
        for (int v = 0; v < grid.numFaceVertices(face); ++v) {
            os << ' ' << next++;
        }
        os << '\n';
        offsets.push_back(next);
    }
    os << "        </DataArray>\n";
    os << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n         ";
    for (const auto o : offsets) {
        os << ' ' << o;
    }
    os << "\n        </DataArray>\n";
    os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         ";
    for (std::size_t f = 0; f < faces.size(); ++f) {
        os << " 7";                       // VTK_POLYGON
    }
    os << "\n        </DataArray>\n      </Cells>\n";

    os << "      <CellData Scalars=\"nodes\">\n";
    const auto array = [&os, &faces](const char* name, auto value) {
        os << "        <DataArray type=\"Int32\" Name=\"" << name << "\" format=\"ascii\">\n"
           << "         ";
        for (const int face : faces) {
            os << ' ' << value(face);
        }
        os << "\n        </DataArray>\n";
    };
    array("nodes", [&grid](int f) { return grid.numFaceVertices(f); });
    array("cell0", [&grid](int f) { return grid.faceCell(f, 0); });
    array("cell1", [&grid](int f) { return grid.faceCell(f, 1); });
    // How many faces the cell on each side has: a merged cell that meets
    // finer neighbours has many more than the six of a hexahedron.
    array("facesOnCell0", [&grid](int f) {
        const int c = grid.faceCell(f, 0);
        return (c < 0) ? 0 : grid.numCellFaces(c);
    });
    array("facesOnCell1", [&grid](int f) {
        const int c = grid.faceCell(f, 1);
        return (c < 0) ? 0 : grid.numCellFaces(c);
    });
    array("boundary", [&grid](int f) {
        return (grid.faceCell(f, 0) < 0 || grid.faceCell(f, 1) < 0) ? 1 : 0;
    });
    os << "      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";

    std::cout << "wrote " << path << " (" << faces.size() << " faces)\n";
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    // --collapse: one face between two coarse cells, as an LGR would have it.
    std::vector<char*> args(argv, argv + argc);
    const auto flag = std::find(args.begin(), args.end(), std::string("--collapse"));
    const bool collapse = flag != args.end();
    if (collapse) {
        args.erase(flag);
    }
    argc = static_cast<int>(args.size());
    argv = args.data();

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
    const std::string vtkPrefix = (argc > 3) ? argv[3] : std::string{};
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
    grid.processEclipseFormatCoarsened(input, layout.blockOfCartesian, layout.boxes, true, collapse);

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

    // Edge conformality: no node may sit inside another face's edge. A node on
    // an edge has to be a corner of every face that edge belongs to, or a cell
    // meeting it does not carry that degree of freedom.
    {
        std::vector<std::array<double,3>> pts;
        pts.reserve(gv.size(3));
        for (int v = 0; v < gv.size(3); ++v) {
            const auto& p = grid.vertexPosition(v);
            pts.push_back({p[0], p[1], p[2]});
        }
        // bucket the nodes so each edge only looks at what is near it
        double lo[3] = {1e30, 1e30, 1e30}, hi[3] = {-1e30, -1e30, -1e30};
        for (const auto& p : pts) {
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], p[d]);
                hi[d] = std::max(hi[d], p[d]);
            }
        }
        const int N = 64;
        const auto bucket = [&](const std::array<double,3>& p) {
            int b = 0;
            for (int d = 0; d < 3; ++d) {
                const double t = (p[d] - lo[d])/std::max(1e-30, hi[d] - lo[d]);
                const int i = std::min(N - 1, std::max(0, int(t*N)));
                b = b*N + i;
            }
            return b;
        };
        std::map<int, std::vector<int>> cells;
        for (std::size_t v = 0; v < pts.size(); ++v) {
            cells[bucket(pts[v])].push_back(int(v));
        }

        std::size_t onEdge = 0;
        for (int c = 0; c < gv.size(0); ++c) {
            for (int f = 0; f < grid.numCellFaces(c); ++f) {
                const int face = grid.cellFace(c, f);
                const int nv = grid.numFaceVertices(face);
                for (int e = 0; e < nv; ++e) {
                    const int a = grid.faceVertex(face, e);
                    const int b = grid.faceVertex(face, (e + 1) % nv);
                    const auto& pa = pts[a];
                    const auto& pb = pts[b];
                    double len2 = 0.0;
                    for (int d = 0; d < 3; ++d) {
                        len2 += (pb[d]-pa[d])*(pb[d]-pa[d]);
                    }
                    if (len2 < 1e-18) { continue; }
                    // candidates: nodes bucketed near the midpoint
                    std::array<double,3> mid{0.5*(pa[0]+pb[0]), 0.5*(pa[1]+pb[1]),
                                             0.5*(pa[2]+pb[2])};
                    const int mb = bucket(mid);
                    for (const int v : cells[mb]) {
                        if (v == a || v == b) { continue; }
                        const auto& p = pts[v];
                        double t = 0.0;
                        for (int d = 0; d < 3; ++d) {
                            t += (p[d]-pa[d])*(pb[d]-pa[d]);
                        }
                        t /= len2;
                        if (t <= 1e-9 || t >= 1 - 1e-9) { continue; }
                        double off = 0.0;
                        for (int d = 0; d < 3; ++d) {
                            const double proj = pa[d] + t*(pb[d]-pa[d]);
                            off += (p[d]-proj)*(p[d]-proj);
                        }
                        if (off < 1e-12*len2) {
                            if (onEdge < 5) {
                                std::cout << "node " << v << " sits inside an edge of face "
                                          << face << " (cell " << c << ")\n";
                            }
                            ++onEdge;
                        }
                    }
                }
            }
        }
        std::cout << "nodes sitting inside another face's edge: " << onEdge << '\n';
    }
    std::cout << "total volume: " << total << '\n';
    std::cout << "cells whose faces do not close: " << openCells
              << ", most faces on a cell: " << worstFaces << '\n';

    if (!vtkPrefix.empty()) {
        // The cells are written as hexahedra through their eight corners, so
        // the picture shows each coarse cell's extent, not the faces it is
        // subdivided into.
        const auto write = [&vtkPrefix](const Dune::CpGrid& g, const std::string& name) {
            Dune::VTKWriter<Dune::CpGrid::LeafGridView> writer(g.leafGridView());
            writer.write(vtkPrefix + name);
            std::cout << "wrote " << vtkPrefix + name << ".vtu\n";
        };
        write(grid, "_merged");
        writeFaceVtu(grid, vtkPrefix + "_merged_faces.vtu");

        Dune::CpGrid fineGrid;
        grdecl asIs{};
        asIs.dims[0] = nx; asIs.dims[1] = ny; asIs.dims[2] = nz;
        asIs.coord = fine.coord.data();
        asIs.zcorn = fine.zcorn.data();
        asIs.actnum = fine.actnum.data();
        fineGrid.processEclipseFormat(asIs, false, false, true);
        write(fineGrid, "_fine");
        writeFaceVtu(fineGrid, vtkPrefix + "_fine_faces.vtu");

        // The same coarsening as a corner-point description, where possible.
        try {
            Opm::Coarsening::Options options;
            options.activity = Opm::Coarsening::Activity::FillHoles;
            const auto coarse = Opm::Coarsening::coarsenCornerPoint(fine, requests, options);
            grdecl cg{};
            cg.dims[0] = coarse.grid.dims[0];
            cg.dims[1] = coarse.grid.dims[1];
            cg.dims[2] = coarse.grid.dims[2];
            cg.coord = coarse.grid.coord.data();
            cg.zcorn = coarse.grid.zcorn.data();
            cg.actnum = coarse.grid.actnum.data();
            Dune::CpGrid grdeclGrid;
            grdeclGrid.processEclipseFormat(cg, false, false, true);
            write(grdeclGrid, "_grdecl");
            writeFaceVtu(grdeclGrid, vtkPrefix + "_grdecl_faces.vtu");
        } catch (const std::invalid_argument& e) {
            std::cout << "no corner-point version of this coarsening: " << e.what() << '\n';
        }
    }
    return openCells == 0 ? 0 : 1;
}
