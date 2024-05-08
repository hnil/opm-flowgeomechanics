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
    }

    FractureWell::FractureWell(std::string name,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments){
        name_ = name;
        this->init(points, segments);
        vtkwriter_ =
            std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>
            (grid_->leafGridView(), Dune::VTK::nonconforming);
    }
}
