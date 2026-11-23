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

#include <opm/geomech/FractureMesh.hpp>
#include <opm/geomech/FractureMeshFoamGrid.hpp>
#include <opm/geomech/FractureVtkWriter.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <zlib.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <span>
#include <string>
#include <vector>

using namespace Opm;

namespace
{
int failures = 0;

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
bool same(const FractureMesh::Vec3& a, const FractureMesh::Vec3& b)
{
    return same(a[0], b[0]) && same(a[1], b[1]) && same(a[2], b[2]);
}
bool same(int a, int b) { return a == b; }
bool same(const FractureMesh::Triangle& a, const FractureMesh::Triangle& b)
{
    return same(a[0], b[0]) && same(a[1], b[1]) && same(a[2], b[2]);
}
template <class T>
bool same(std::span<const T> a, std::span<const T> b)
{
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same(a[i], b[i]))
            return false;
    return true;
}

void check(bool ok, const std::string& what, int i)
{
    if (!ok && failures++ < 20)
        std::cerr << "differs: " << what << " at " << i << '\n';
}

// Both builders must give the same numbering and the same bits.
void compare(const Triangulation& tri, const std::string& name)
{
    const auto grid = RegularTrimesh::createDuneGrid(tri);
    const auto g = buildFractureMesh(*grid);
    const auto t = buildFractureMesh(tri.nodes, tri.triangles);
    const int before = failures;
    check(g.numCells() == t.numCells() && g.numFaces() == t.numFaces() && g.numNodes() == t.numNodes(),
          name + " sizes", 0);
    if (failures > before)
        return;
    for (int c = 0; c < g.numCells(); ++c) {
        check(same(g.cellArea(c), t.cellArea(c)), name + " cellArea", c);
        check(same(g.cellCentroid(c), t.cellCentroid(c)), name + " cellCentroid", c);
        check(same(g.cellNormal(c), t.cellNormal(c)), name + " cellNormal", c);
        check(g.cellLevel(c) == t.cellLevel(c), name + " cellLevel", c);
        check(same(g.cellTriangles(c), t.cellTriangles(c)), name + " cellTriangles", c);
        check(same(g.cellNodes(c), t.cellNodes(c)), name + " cellNodes", c);
        check(same(g.cellCorners(c), t.cellCorners(c)), name + " cellCorners", c);
        check(same(g.cellFaces(c), t.cellFaces(c)), name + " cellFaces", c);
        check(same(g.cellNeighbours(c), t.cellNeighbours(c)), name + " cellNeighbours", c);
        for (int j = 0; j < static_cast<int>(g.cellFaces(c).size()); ++j) {
            check(same(g.cellFaceCentroid(c, j), t.cellFaceCentroid(c, j)), name + " cellFaceCentroid", c);
            check(same(g.cellFaceLength(c, j), t.cellFaceLength(c, j)), name + " cellFaceLength", c);
            check(same(g.cellFaceNormal(c, j), t.cellFaceNormal(c, j)), name + " cellFaceNormal", c);
            check(same(g.cellFaceCorners(c, j)[0], t.cellFaceCorners(c, j)[0])
                      && same(g.cellFaceCorners(c, j)[1], t.cellFaceCorners(c, j)[1]),
                  name + " cellFaceCorners", c);
        }
    }
    for (int f = 0; f < g.numFaces(); ++f) {
        check(same(g.faceCells(f), t.faceCells(f)), name + " faceCells", f);
        check(g.faceNodes(f) == t.faceNodes(f), name + " faceNodes", f);
        check(same(g.faceCentroid(f), t.faceCentroid(f)), name + " faceCentroid", f);
        check(same(g.faceLength(f), t.faceLength(f)), name + " faceLength", f);
    }
    for (int n = 0; n < g.numNodes(); ++n)
        check(same(g.nodeCoord(n), t.nodeCoord(n)), name + " nodeCoord", n);
    std::cout << name << ": " << g.numCells() << " cells, " << g.numFaces() << " faces, "
              << (failures == before ? "identical" : "DIFFER") << '\n';
}

std::string slurp(const std::string& file)
{
    std::ifstream in(file, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string unbase64(const std::string& t)
{
    static const std::string abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    unsigned buf = 0;
    int bits = 0;
    for (char ch : t) {
        const auto k = abc.find(ch);
        if (k == std::string::npos)
            continue;
        buf = (buf << 6) | static_cast<unsigned>(k);
        if ((bits += 6) >= 8) {
            bits -= 8;
            out.push_back(static_cast<char>((buf >> bits) & 0xff));
        }
    }
    return out;
}

// Payloads of all inline DataArrays, in file order.
std::vector<std::string> payloads(const std::string& text)
{
    std::vector<std::string> p;
    for (std::size_t at = 0; (at = text.find("<DataArray", at)) != std::string::npos;) {
        const auto b = text.find('>', at) + 1;
        const auto e = text.find("</DataArray>", b);
        p.push_back(text.substr(b, e - b));
        at = e;
    }
    return p;
}

std::string noSpace(std::string t)
{
    t.erase(std::remove_if(t.begin(), t.end(), [](char ch) { return std::isspace(ch); }), t.end());
    return t;
}

// The 4-byte size header is base64-encoded (and padded) on its own: 8 characters.
std::string decodeBinary(const std::string& payload)
{
    return unbase64(noSpace(payload).substr(8));
}

std::string decodeZlib(const std::string& text)
{
    const std::string payload = noSpace(text);
    std::uint32_t h[3];
    std::memcpy(h, unbase64(payload.substr(0, 16)).data(), sizeof(h));
    const std::size_t hbytes = 4 * (3 + h[0]);
    const std::size_t hchars = 4 * ((hbytes + 2) / 3);
    const std::string header = unbase64(payload.substr(0, hchars));
    const std::string data = unbase64(payload.substr(hchars));
    std::string out;
    std::size_t at = 0;
    for (std::uint32_t b = 0; b < h[0]; ++b) {
        std::uint32_t clen;
        std::memcpy(&clen, header.data() + 4 * (3 + b), 4);
        uLongf len = (b + 1 == h[0] && h[2] != 0) ? h[2] : h[1];
        std::string block(len, '\0');
        if (uncompress(reinterpret_cast<Bytef*>(block.data()), &len,
                       reinterpret_cast<const Bytef*>(data.data() + at), clen) != Z_OK)
            return "inflate failed";
        out += block.substr(0, len);
        at += clen;
    }
    return out;
}

// The mesh writer must give Dune::VTKWriter's bytes on the FoamGrid of the same triangulation.
void compareVtk(const Triangulation& tri, const std::string& name)
{
    using Format = FractureVtkWriter::Format;
    const auto grid = RegularTrimesh::createDuneGrid(tri);
    const auto mesh = buildFractureMesh(tri.nodes, tri.triangles);
    std::vector<double> a(mesh.numCells()), b(mesh.numCells());
    for (int c = 0; c < mesh.numCells(); ++c) {
        a[c] = mesh.cellCentroid(c)[2] + 1e-7 * c;
        b[c] = mesh.cellArea(c);
    }
    const std::string dir = std::filesystem::temp_directory_path().string();
    const int before = failures;
    for (const bool conforming : {true, false}) {
        const auto prec = conforming ? Dune::VTK::Precision::float64 : Dune::VTK::Precision::float32;
        const auto tag = name + (conforming ? " conforming" : " nonconforming");
        FractureVtkWriter mine;
        mine.addCellData(a, "A");
        mine.addCellData(b, "B");
        Dune::VTKWriter dune(grid->leafGridView(), conforming ? Dune::VTK::conforming : Dune::VTK::nonconforming,
                             prec);
        dune.addCellData(a, "A");
        dune.addCellData(b, "B");
        const std::pair<Format, Dune::VTK::OutputType> formats[]
            = {{Format::ascii, Dune::VTK::ascii}, {Format::binary, Dune::VTK::base64},
               {Format::appended, Dune::VTK::appendedraw}};
        for (const auto& [mf, df] : formats) {
            const auto m = slurp(mine.write(dir + "/fvw_mine", mesh, conforming, prec, mf));
            const auto d = slurp(dune.write(dir + "/fvw_dune", df));
            check(m == d, tag + " format " + std::to_string(static_cast<int>(mf)), 0);
        }
        const auto bin = payloads(slurp(mine.write(dir + "/fvw_mine", mesh, conforming, prec, Format::binary)));
        const auto zlib = slurp(mine.write(dir + "/fvw_mine", mesh, conforming, prec, Format::zlib));
        check(zlib.find("compressor=\"vtkZLibDataCompressor\"") != std::string::npos, tag + " zlib attribute", 0);
        const auto zp = payloads(zlib);
        check(zp.size() == bin.size(), tag + " zlib array count", 0);
        for (std::size_t i = 0; i < std::min(zp.size(), bin.size()); ++i)
            check(decodeZlib(zp[i]) == decodeBinary(bin[i]), tag + " zlib array", static_cast<int>(i));
    }
    std::cout << name << " vtk: " << (failures == before ? "identical" : "DIFFER") << '\n';
}
} // namespace

int main()
{
    // large offsets and an oblique plane, so round-off is not trivially exact
    const std::array<double, 3> origin {456789.123, 7123456.789, 2534.5};
    const double s = 1.0 / std::sqrt(2.0);
    const std::array<double, 3> ax1 {s, s, 0.0};
    const std::array<double, 3> ax2 {0.5 * s, 0.5 * s, std::sqrt(3.0) / 2};
    RegularTrimesh mesh(60.0, origin, ax1, ax2, {7.3, 7.3});
    mesh.removeSawtooths();
    mesh.expandGrid(CellRef {9, 0, 0});
    mesh.expandGrid(CellRef {-3, 8, 1});
    mesh.removeSawtooths();
    const auto ring = RegularTrimesh::inner_ring_cells();

    compare(mesh.createTriangulation(1, ring), "seed level 1");
    compare(mesh.createTriangulation(0, {}, true), "fine, smoothing");
    compare(mesh.createTriangulation(3, ring, false, 10, 1), "coarsen 3");
    compare(mesh.createTriangulation(3, ring, true, 4, 2), "coarsen 3, smoothing, band 2");
    compareVtk(mesh.createTriangulation(1, ring), "seed level 1");
    compareVtk(mesh.createTriangulation(3, ring, true, 4, 2), "coarsen 3, smoothing, band 2");

    // large enough for several zlib blocks per array
    Triangulation big;
    const unsigned n = 60;
    for (unsigned j = 0; j <= n; ++j)
        for (unsigned i = 0; i <= n; ++i)
            big.nodes.push_back({origin[0] + 3.1 * i, origin[1] + 0.7 * j, origin[2] + 2.9 * j});
    for (unsigned j = 0; j < n; ++j)
        for (unsigned i = 0; i < n; ++i) {
            const unsigned v = j * (n + 1) + i;
            big.triangles.push_back({v, v + 1, v + n + 2});
            big.triangles.push_back({v, v + n + 2, v + n + 1});
        }
    compareVtk(big, "60x60 squares");

    std::cout << (failures == 0 ? "PASSED" : "FAILED") << '\n';
    return failures == 0 ? 0 : 1;
}
