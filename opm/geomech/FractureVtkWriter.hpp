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

#pragma once

#include <dune/grid/io/file/vtk/common.hh>

#include <array>
#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

class FractureMesh;

// Cells of one VTK type given by node lists, e.g. a well trajectory as line segments.
struct VtkCellList {
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::vector<unsigned>> cells;
    Dune::VTK::GeometryType type = Dune::VTK::line;
};

// VTU output of a FractureMesh with cell data, laid out as Dune::VTKWriter lays out the same triangles.
class FractureVtkWriter
{
public:
    // binary = inline base64, appended = raw appended, zlib = inline base64 of vtkZLibDataCompressor blocks
    enum class Format { ascii, binary, appended, zlib };
    static Format formatFromString(const std::string& name);

    template <class Vector>
    void addCellData(const Vector& v, const std::string& name)
    {
        std::vector<double> values(v.size());
        for (std::size_t i = 0; i < v.size(); ++i) {
            values[i] = static_cast<double>(v[i]);
        }
        cell_data_.emplace_back(name, std::move(values));
    }
    void clear() { cell_data_.clear(); }

    // Writes <name>.vtu (.vtp for line cells) and returns that file name.
    std::string write(const std::string& name,
                      const FractureMesh& mesh,
                      bool conforming,
                      Dune::VTK::Precision coord_prec,
                      Format format) const;
    std::string write(const std::string& name,
                      const VtkCellList& cells,
                      bool conforming,
                      Dune::VTK::Precision coord_prec,
                      Format format) const;

private:
    std::vector<std::pair<std::string, std::vector<double>>> cell_data_;
};

// <dir>/<sim>-NNNNN.vtu per call plus <dir>/<sim>.pvd, named and formatted as Opm::VtkMultiWriter does.
class FractureVtkSeriesWriter
{
public:
    FractureVtkSeriesWriter(std::string output_dir, std::string sim_name, FractureVtkWriter::Format format);

    FractureVtkWriter& beginWrite(double time);
    void endWrite(const FractureMesh& mesh);
    void endWrite(const VtkCellList& cells);

private:
    template <class Cells>
    void endWrite_(const Cells& cells);
    void finishCollection_();

    std::string output_dir_;
    std::string sim_name_;
    FractureVtkWriter::Format format_;
    std::ofstream pvd_;
    FractureVtkWriter writer_;
    double time_ = 0.0;
    int num_written_ = 0;
};

} // namespace Opm
