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

#include <opm/geomech/FractureVtkWriter.hpp>

#include <opm/geomech/FractureMesh.hpp>

#include <dune/common/indent.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
#include <dune/grid/io/file/vtk/streams.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace Opm
{

namespace
{

class ZlibDataArrayWriter final : public Dune::VTK::DataArrayWriter
{
public:
    ZlibDataArrayWriter(std::ostream& s, const std::string& name, unsigned ncomps, Dune::VTK::Precision prec)
        : DataArrayWriter(prec), s_(s), name_(name), ncomps_(ncomps), prec_(prec)
    {
    }

    ~ZlibDataArrayWriter() override
    {
        // Blocks and header as vtkZLibDataCompressor with UInt32 headers (the VTK file format 0.1 default).
        constexpr std::size_t block = 1 << 15;
        const std::size_t nblocks = (raw_.size() + block - 1) / block;
        std::vector<std::uint32_t> header {static_cast<std::uint32_t>(nblocks),
                                           static_cast<std::uint32_t>(block),
                                           static_cast<std::uint32_t>(raw_.size() % block)};
        std::vector<Bytef> packed;
        for (std::size_t b = 0; b < nblocks; ++b) {
            const std::size_t len = std::min(block, raw_.size() - b * block);
            uLongf clen = compressBound(len);
            const std::size_t at = packed.size();
            packed.resize(at + clen);
            if (compress2(packed.data() + at, &clen, raw_.data() + b * block, len, Z_DEFAULT_COMPRESSION) != Z_OK) {
                std::terminate(); // only Z_MEM_ERROR is possible with a compressBound-sized buffer
            }
            packed.resize(at + clen);
            header.push_back(static_cast<std::uint32_t>(clen));
        }

        const Dune::Indent indent("  ", 4);
        s_ << indent << "<DataArray type=\"" << Dune::VTK::toString(prec_) << "\" Name=\"" << name_
           << "\" NumberOfComponents=\"" << ncomps_ << "\" format=\"binary\">\n";
        s_ << indent + 1;
        Dune::Base64Stream h(s_);
        for (auto& x : header) {
            h.write(x);
        }
        h.flush();
        Dune::Base64Stream d(s_);
        for (auto& x : packed) {
            d.write(x);
        }
        d.flush();
        s_ << "\n" << indent << "</DataArray>\n";
    }

private:
    template <class T>
    void put_(T v)
    {
        const auto* p = reinterpret_cast<const Bytef*>(&v);
        raw_.insert(raw_.end(), p, p + sizeof(T));
    }
    void writeFloat64(double v) override { put_(v); }
    void writeFloat32(float v) override { put_(v); }
    void writeInt32(std::int32_t v) override { put_(v); }
    void writeUInt32(std::uint32_t v) override { put_(v); }
    void writeUInt8(std::uint8_t v) override { put_(v); }

    std::ostream& s_;
    std::string name_;
    unsigned ncomps_;
    Dune::VTK::Precision prec_;
    std::vector<Bytef> raw_;
};

// Dune's VTUWriter with zlib arrays; the caller adds the compressor attribute to the VTKFile line.
class ZlibVTUWriter : public Dune::VTK::VTUWriter
{
public:
    explicit ZlibVTUWriter(std::ostream& s)
        : VTUWriter(s, Dune::VTK::ascii, Dune::VTK::unstructuredGrid)
    {
    }
    Dune::VTK::DataArrayWriter*
    makeArrayWriter(const std::string& name, unsigned ncomps, unsigned, Dune::VTK::Precision prec)
    {
        return new ZlibDataArrayWriter(stream, name, ncomps, prec);
    }
};

struct Piece {
    std::vector<FractureMesh::Vec3> points;
    std::vector<int> connectivity;
};

// Point order and numbering of Dune::VTKWriter: first visit over cells and their corners.
Piece
makePiece(const FractureMesh& mesh, bool conforming)
{
    Piece p;
    std::vector<int> number(conforming ? mesh.numNodes() : 0, -1);
    for (int c = 0; c < mesh.numCells(); ++c) {
        const auto nodes = mesh.cellNodes(c);
        const auto corners = mesh.cellCorners(c);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (!conforming) {
                p.connectivity.push_back(static_cast<int>(p.points.size()));
                p.points.push_back(corners[i]);
                continue;
            }
            if (number[nodes[i]] < 0) {
                number[nodes[i]] = static_cast<int>(p.points.size());
                p.points.push_back(corners[i]);
            }
            p.connectivity.push_back(number[nodes[i]]);
        }
    }
    return p;
}

template <class Writer>
void
writeArrays(Writer& w,
            const FractureMesh& mesh,
            const Piece& p,
            const std::vector<std::pair<std::string, std::vector<double>>>& cell_data,
            Dune::VTK::Precision coord_prec)
{
    using Dune::VTK::Precision;
    const unsigned ncells = mesh.numCells();
    auto array = [&w](const std::string& name, unsigned ncomps, unsigned n, Precision prec) {
        return std::unique_ptr<Dune::VTK::DataArrayWriter>(w.makeArrayWriter(name, ncomps, n, prec));
    };

    if (!cell_data.empty()) {
        w.beginCellData(cell_data.front().first, "");
        for (const auto& [name, values] : cell_data) {
            auto a = array(name, 1, ncells, Precision::float32);
            if (!a->writeIsNoop()) {
                for (unsigned c = 0; c < ncells; ++c) {
                    a->write(values[c]);
                }
            }
        }
        w.endCellData();
    }

    w.beginPoints();
    {
        auto a = array("Coordinates", 3, p.points.size(), coord_prec);
        if (!a->writeIsNoop()) {
            for (const auto& x : p.points) {
                for (int d = 0; d < 3; ++d) {
                    a->write(x[d]);
                }
            }
        }
    }
    w.endPoints();

    w.beginCells();
    {
        auto a = array("connectivity", 1, p.connectivity.size(), Precision::int32);
        if (!a->writeIsNoop()) {
            for (int v : p.connectivity) {
                a->write(v);
            }
        }
    }
    {
        auto a = array("offsets", 1, ncells, Precision::int32);
        if (!a->writeIsNoop()) {
            int offset = 0;
            for (unsigned c = 0; c < ncells; ++c) {
                offset += static_cast<int>(mesh.cellNodes(c).size());
                a->write(offset);
            }
        }
    }
    {
        auto a = array("types", 1, ncells, Precision::uint8);
        if (!a->writeIsNoop()) {
            for (unsigned c = 0; c < ncells; ++c) {
                a->write(static_cast<int>(Dune::VTK::triangle));
            }
        }
    }
    w.endCells();
}

template <class Writer>
void
writePiece(Writer& w,
           const FractureMesh& mesh,
           const std::vector<std::pair<std::string, std::vector<double>>>& cell_data,
           Dune::VTK::Precision coord_prec,
           bool conforming)
{
    const Piece p = makePiece(mesh, conforming);
    w.beginMain(mesh.numCells(), p.points.size());
    writeArrays(w, mesh, p, cell_data, coord_prec);
    w.endMain();
    if (w.beginAppended()) {
        writeArrays(w, mesh, p, cell_data, coord_prec);
    }
    w.endAppended();
}

} // namespace

FractureVtkWriter::Format
FractureVtkWriter::formatFromString(const std::string& name)
{
    if (name == "ascii") {
        return Format::ascii;
    }
    if (name == "binary") {
        return Format::binary;
    }
    if (name == "appended") {
        return Format::appended;
    }
    if (name == "zlib") {
        return Format::zlib;
    }
    throw std::invalid_argument("Unknown vtk_format '" + name + "' (ascii, binary, appended or zlib)");
}

std::string
FractureVtkWriter::write(const std::string& name,
                         const FractureMesh& mesh,
                         bool conforming,
                         Dune::VTK::Precision coord_prec,
                         Format format) const
{
    const std::string file_name = name + ".vtu";
    std::ofstream file;
    file.exceptions(std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);
    file.open(file_name, std::ios::binary);

    if (format == Format::zlib) {
        std::ostringstream os;
        {
            ZlibVTUWriter w(os);
            writePiece(w, mesh, cell_data_, coord_prec, conforming);
        }
        std::string text = os.str();
        const std::string key = "byte_order=\"" + Dune::VTK::getEndiannessString() + "\"";
        text.insert(text.find(key) + key.size(), " compressor=\"vtkZLibDataCompressor\"");
        file << text;
    } else {
        const auto type = format == Format::binary ? Dune::VTK::base64
            : format == Format::appended           ? Dune::VTK::appendedraw
                                                   : Dune::VTK::ascii;
        Dune::VTK::VTUWriter w(file, type, Dune::VTK::unstructuredGrid);
        writePiece(w, mesh, cell_data_, coord_prec, conforming);
    }
    file.close();
    return file_name;
}

FractureVtkSeriesWriter::FractureVtkSeriesWriter(std::string output_dir,
                                                 std::string sim_name,
                                                 FractureVtkWriter::Format format)
    : output_dir_(output_dir.empty() ? "." : std::move(output_dir))
    , sim_name_(sim_name.empty() ? "sim" : std::move(sim_name))
    , format_(format)
{
}

FractureVtkWriter&
FractureVtkSeriesWriter::beginWrite(double time)
{
    if (!pvd_.is_open()) {
        pvd_.open(output_dir_ + "/" + sim_name_ + ".pvd");
        pvd_ << "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"Collection\"\n"
                "         version=\"0.1\"\n"
                "         byte_order=\"LittleEndian\"\n"
                "         compressor=\"vtkZLibDataCompressor\">\n"
                " <Collection>\n";
    }
    time_ = time;
    writer_.clear();
    return writer_;
}

void
FractureVtkSeriesWriter::endWrite(const FractureMesh& mesh)
{
    std::ostringstream base;
    base << sim_name_ << "-" << std::setw(5) << std::setfill('0') << num_written_++;
    const std::string file = writer_.write(output_dir_ + "/" + base.str(), mesh, true,
                                           Dune::VTK::Precision::float64, format_);
    writer_.clear();
    pvd_.precision(16);
    pvd_ << "   <DataSet timestep=\"" << time_ << "\" file=\""
         << std::filesystem::path(file).filename().string() << "\"/>\n";
    finishCollection_();
}

// Keeps the .pvd valid after every step: closing tags are written, then overwritten by the next entry.
void
FractureVtkSeriesWriter::finishCollection_()
{
    const auto pos = pvd_.tellp();
    pvd_ << " </Collection>\n"
            "</VTKFile>\n";
    pvd_.seekp(pos);
    pvd_.flush();
}

} // namespace Opm
