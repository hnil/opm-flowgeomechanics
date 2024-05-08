#pragma once
#include <opm/geomech/Fracture.hpp>
namespace Opm{
class FractureWell
{
    using Grid = Dune::FoamGrid<1,3>;
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned,2>;
private:
    std::string name_;
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
public:
    FractureWell(std::string name,
                 const std::vector<Point3D>& points,
                 const std::vector<Segment>& segments);
    void init(const std::vector<Point3D>& points,
              const std::vector<Segment>& segments);
    template<class Grid,class Well>
    FractureWell(std::string name,Grid grid,const Well& well){
        name_ = name;
    }
    const Grid& grid() const{return *grid_;}
    std::string name() const{return name_;};
    void write() const{
        vtkwriter_->write(this->name().c_str());
    }
};
}
