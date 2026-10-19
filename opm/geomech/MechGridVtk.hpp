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
#ifndef OPM_MECH_GRID_VTK_HPP
#define OPM_MECH_GRID_VTK_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/grid/CpGrid.hpp>

#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

/// VTK output of a mechanics grid of its own, laid out as the flow grid's geomechanical VTK
/// (ASCII, hexahedra, "disp" first among the point data, tensors as rows "stress[0..2]"),
/// so ResInsight opens it as a geomechanical case.
namespace Opm::MechGridVtk
{

struct Field
{
    std::string name;
    int components = 1;
    std::vector<double> values;    // size = entities * components
};

inline void writeArray(std::ofstream& os, const Field& f)
{
    os << "<DataArray type=\"Float64\" Name=\"" << f.name << "\" NumberOfComponents=\""
       << f.components << "\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < f.values.size(); ++i) {
        os << f.values[i] << (((i + 1) % (3 * f.components) == 0) ? "\n" : " ");
    }
    os << "\n</DataArray>\n";
}

inline void writeVtu(const std::string& file, const Dune::CpGrid& grid,
                     const std::vector<Field>& cellData, const std::vector<Field>& pointData)
{
    std::ofstream os(file);
    if (!os) {
        OPM_THROW(std::runtime_error, "Cannot write " + file);
    }
    os.precision(10);
    const auto& gv = grid.leafGridView();
    const auto& is = gv.indexSet();
    const std::size_t nc = gv.size(0), nv = gv.size(3);

    os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n<Piece NumberOfCells=\""
       << nc << "\" NumberOfPoints=\"" << nv << "\">\n";

    os << "<PointData>\n";
    for (const auto& f : pointData) {
        writeArray(os, f);
    }
    os << "</PointData>\n<CellData>\n";
    for (const auto& f : cellData) {
        writeArray(os, f);
    }
    os << "</CellData>\n";

    std::vector<double> xyz(3 * nv);
    for (const auto& v : vertices(gv)) {
        const auto x = v.geometry().center();
        for (int d = 0; d < 3; ++d) {
            xyz[3 * is.index(v) + d] = x[d];
        }
    }
    os << "<Points>\n";
    writeArray(os, Field{"Coordinates", 3, xyz});
    os << "</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    // Dune numbers hexahedron corners lexicographically, VTK walks each face round.
    static constexpr int vtkOrder[8] = {0, 1, 3, 2, 4, 5, 7, 6};
    for (const auto& e : elements(gv)) {
        for (const int c : vtkOrder) {
            os << is.index(e.template subEntity<3>(c)) << " ";
        }
        os << "\n";
    }
    os << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t c = 1; c <= nc; ++c) {
        os << 8 * c << ((c % 12 == 0) ? "\n" : " ");
    }
    os << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 1; c <= nc; ++c) {
        os << 12 << ((c % 24 == 0) ? "\n" : " ");
    }
    os << "\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

inline void writePvd(const std::string& file, const std::vector<std::pair<double, std::string>>& steps)
{
    std::ofstream os(file);
    os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for (const auto& [time, vtu] : steps) {
        os << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << vtu << "\"/>\n";
    }
    os << "</Collection>\n</VTKFile>\n";
}

} // namespace Opm::MechGridVtk

#endif // OPM_MECH_GRID_VTK_HPP
