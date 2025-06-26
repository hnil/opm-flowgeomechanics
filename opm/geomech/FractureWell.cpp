#include "config.h"

#include <opm/geomech/FractureWell.hpp>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <opm/models/io/vtkmultiwriter.hh>

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
        : outputprefix_ { outputprefix }
        , casename_     { casename }
        , name_         { name }
        , conns_        { conns }
    {
        this->init(points, segments);

        this->vtkwriter_ =
            std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>
            (grid_->leafGridView(), Dune::VTK::nonconforming);

        reservoir_stress_.resize(this->conns_.size());
        std::fill(reservoir_stress_.begin(),
                  reservoir_stress_.end(),
                  100e5);// random initialization
    }

    void FractureWell::init(const std::vector<Point3D>& points,
                            const std::vector<Segment>& segments)
    {
        Dune::GridFactory<Grid> factory;
        for(size_t i=0; i < points.size(); ++i){
            factory.insertVertex(points[i]);
        }
        for (size_t i = 0; i < segments.size(); ++i) {
            std::vector<unsigned int> seg(2,0);
            seg[0] = segments[i][0];seg[1] = segments[i][1];
            factory.insertElement(Dune::GeometryTypes::line, seg); // std::move(mapping));
        }
        grid_ = factory.createGrid();
        this->resetWriters();
    }

    void FractureWell::write() const{
        std::string filename = outputprefix_ + "/" + this->name();
        vtkwriter_->write(filename.c_str());
    }

    void FractureWell::resetWriters(){
        // nead to be reseat if grid is changed ??
        vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(), Dune::VTK::nonconforming);
        std::string outputdir = outputprefix_;
        std::string simName = casename_ + "_" + this->name();
        std::string multiFileName =  "";
        vtkmultiwriter_ = std::make_unique< Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat > >(/*async*/ false,
                                                                                              grid_->leafGridView(),
                                                                                              outputdir,
                                                                                              simName,
                                                                                              multiFileName
        );
    }

    void FractureWell::writemulti(double time) const{
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
        vtkmultiwriter_->endWrite(false);
    }
}
