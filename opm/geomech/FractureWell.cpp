#include "config.h"

#include <opm/geomech/FractureWell.hpp>

#include <dune/geometry/referenceelements.hh>

namespace Opm {
    void FractureWell::init(const std::vector<Point3D>& points,
                            const std::vector<Segment>& segments)
    {
        Dune::GridFactory<Grid> factory{};

        for (const auto& vertex : points) {
            factory.insertVertex(vertex);
        }

        for (const auto& segment : segments) {
            std::vector<unsigned int> seg(2,0);
            seg[0] = segment[0];
            seg[1] = segment[1];

            factory.insertElement(Dune::GeometryTypes::line, seg); // std::move(mapping));
        }

        grid_ = factory.createGrid();

        this->resetWriters();
    }

    void FractureWell::write() const
    {
        const std::string filename = outputprefix_ + "/" + this->name();
        vtkwriter_->write(filename);
    }

    void FractureWell::resetWriters()
    {
        // nead to be reseat if grid is changed ??
        vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>
            (grid_->leafGridView(), Dune::VTK::nonconforming);

        const std::string outputdir = outputprefix_;
        const std::string simName = casename_ + "_" + this->name();
        const std::string multiFileName =  "";

        vtkmultiwriter_ = std::make_unique<VtkMultiWriter<Grid::LeafGridView, VTKFormat>>
            (/*async*/ false,
             grid_->leafGridView(),
             outputdir,
             simName,
             multiFileName);
    }

    void FractureWell::writemulti(const double time) const
    {
        std::vector<double> reservoir_pressure = reservoir_pressure_;
        std::vector<double> reservoir_temperature = reservoir_temperature_;
        //std::vector<double> reservoir_cells = reservoir_cells_;
        std::vector<double> perf_pressure = perf_pressure_;
        vtkmultiwriter_->beginWrite(time);
        if (perf_pressure.size() > 0) {
            vtkmultiwriter_->attachScalarElementData(perf_pressure, "PerfPressure");
        }
        if (reservoir_pressure.size() > 0) {
            vtkmultiwriter_->attachScalarElementData(reservoir_pressure, "ReservoirPressure");
        }
        if (reservoir_temperature.size() > 0) {
            vtkmultiwriter_->attachScalarElementData(reservoir_temperature, "ReservoirTemperature");
        }
        vtkmultiwriter_->endWrite();
    }

    FractureWell::FractureWell(const std::string& outputprefix,
                               const std::string& casename,
                               const std::string& name,
                               const std::vector<Point3D>& points,
                               const std::vector<Segment>& segments,
                               const std::vector<int>& res_cells)
    {
        outputprefix_ = outputprefix;
        casename_ = casename;
        name_ = name;

        this->init(points, segments);

        vtkwriter_ =
            std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>
            (grid_->leafGridView(), Dune::VTK::nonconforming);

        reservoir_cells_ = res_cells;
        reservoir_stress_.resize(res_cells.size());
        std::fill(reservoir_stress_.begin(), reservoir_stress_.end(), 100.0e5); // random initialization
    }
}
