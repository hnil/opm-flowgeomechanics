#include "config.h"
#include <opm/geomech/FractureWell.hpp>
#include <dune/geometry/referenceelements.hh>
namespace Opm{
    void FractureWell::init(const std::vector<Point3D>& points,
                   const std::vector<Segment>& segments){
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

    FractureWell::FractureWell(std::string outputprefix, 
                           std::string name,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments,
                           const std::vector<int>& res_cells){
                            outputprefix_ = outputprefix;
        name_ = name;
        this->init(points, segments);
        vtkwriter_ =
            std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>
            (grid_->leafGridView(), Dune::VTK::nonconforming);
        reservoir_cells_ = res_cells;
        reservoir_stress_.resize(res_cells.size());
        std::fill(reservoir_stress_.begin(),reservoir_stress_.end(),100e5);// random initialization    
    }
}
