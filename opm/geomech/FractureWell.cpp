/*
  Copyright 2025 Equinor ASA.

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

#include "config.h"

#include <opm/geomech/FractureWell.hpp>

#include <dune/common/fvector.hh>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <stddef.h>

namespace Opm
{
FractureWell::FractureWell(const std::string& outputprefix,
                           const std::string& casename,
                           const std::string& name,
                           const std::vector<Connection>& conns,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments)
    : outputprefix_ {outputprefix}
    , casename_ {casename}
    , name_ {name}
    , conns_ {conns}
    , points_ {points}
    , segments_ {segments}
{
    this->init(points, segments);

    reservoir_stress_.resize(this->conns_.size());
    std::fill(reservoir_stress_.begin(), reservoir_stress_.end(),
              100e5); // random initialization
}

void
FractureWell::init(const std::vector<Point3D>& points, const std::vector<Segment>& segments)
{
    points_ = points;
    segments_ = segments;
    this->resetWriters();
}

VtkCellList
FractureWell::vtkCells_() const
{
    VtkCellList cells;
    for (const auto& p : points_) {
        cells.nodes.push_back({p[0], p[1], p[2]});
    }
    for (const auto& seg : segments_) {
        cells.cells.push_back({seg[0], seg[1]});
    }
    return cells;
}

void
FractureWell::write() const
{
    FractureVtkWriter().write(outputprefix_ + "/" + this->name(), vtkCells_(), /*conforming*/ false,
                              Dune::VTK::Precision::float32, FractureVtkWriter::Format::ascii);
}

void
FractureWell::resetWriters()
{
    vtkmultiwriter_ = std::make_unique<FractureVtkSeriesWriter>(
        outputprefix_, casename_ + "_" + this->name(), FractureVtkWriter::Format::ascii);
}

void
FractureWell::writemulti(double time) const
{
    auto& vtkwriter = vtkmultiwriter_->beginWrite(time);
    if (perf_pressure_.size() > 0) {
        vtkwriter.addCellData(perf_pressure_, "PerfPressure");
    }
    if (reservoir_pressure_.size() > 0) {
        vtkwriter.addCellData(reservoir_pressure_, "ReservoirPressure");
    }
    if (reservoir_temperature_.size() > 0) {
        vtkwriter.addCellData(reservoir_temperature_, "ReservoirTemperature");
    }
    vtkmultiwriter_->endWrite(vtkCells_());
}
} // namespace Opm
