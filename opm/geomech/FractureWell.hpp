#pragma once
#include <opm/geomech/Fracture.hpp>
namespace Opm{
class FractureWell
{
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned,2>;
public:
    using Grid = Dune::FoamGrid<1,3>;
    FractureWell(std::string outputprefix, 
                 std::string name,
                 const std::vector<Point3D>& points,
                 const std::vector<Segment>& segments,
                 const std::vector<int>& res_cells);
    void init(const std::vector<Point3D>& points,
              const std::vector<Segment>& segments);

    template<class Problem>
    void updateReservoirStress(const Problem& problem){
        assert(grid_->leafGridView().size(0) == reservoir_cells_.size());
        reservoir_stress_.resize(reservoir_cells_.size());
        for(size_t i=0;i < reservoir_cells_.size(); ++i){
            reservoir_stress_[i] = problem.stress(reservoir_cells_[i]);
        }
    }

    const Grid& grid() const{return *grid_;}
    std::string name() const{return name_;};
    void write() const{
        std::string filename = outputprefix_ + "/" + this->name();
        vtkwriter_->write(filename.c_str());
    }
    int reservoirCell(int wellcell) const {return reservoir_cells_[wellcell];};
    Dune::FieldVector<double,6> reservoirStress(int wellcell) const{return reservoir_stress_[wellcell];};
private: 
    std::string outputprefix_;
    std::string name_;
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    // should probably be separated for easier testing
    std::vector<int> reservoir_cells_;  
    std::vector<Dune::FieldVector<double,6>> reservoir_stress_;    
};
}
