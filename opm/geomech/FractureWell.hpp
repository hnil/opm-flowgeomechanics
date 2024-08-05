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
                 std::string casename,
                 std::string name,
                 const std::vector<Point3D>& points,
                 const std::vector<Segment>& segments,
                 const std::vector<int>& res_cells);

    void init(const std::vector<Point3D>& points,
              const std::vector<Segment>& segments);

    template<class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator)
    {
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        assert(grid_->leafGridView().size(0) == reservoir_cells_.size());
        const auto& problem = simulator.problem();
        perf_pressure_.resize(reservoir_cells_.size());
        reservoir_stress_.resize(reservoir_cells_.size());
        reservoir_pressure_.resize(reservoir_cells_.size());
        reservoir_temperature_.resize(reservoir_cells_.size());
        int water_index = 0;
        for(size_t i=0;i < reservoir_cells_.size(); ++i){
            reservoir_stress_[i] = problem.stress(reservoir_cells_[i]);
            {
                const auto& intQuants =
                    simulator.model().intensiveQuantities(reservoir_cells_[i], /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                auto val = fs.pressure(FluidSystem::waterPhaseIdx);
                reservoir_pressure_[i] = Opm::getValue(val);
                unsigned dummy = 0;
                reservoir_temperature_[i] = Opm::getValue(fs.temperature(dummy));
            }

        }
    }
    void setPerfPressure(int perf_index, double pressure){perf_pressure_[perf_index] = pressure;}
    const Grid& grid() const{return *grid_;}
    std::string name() const{return name_;};
    void write() const;
    void writemulti(double time) const;
    void resetWriters();
    int reservoirCell(int wellcell) const {return reservoir_cells_[wellcell];};
    Dune::FieldVector<double,6> reservoirStress(int wellcell) const{return reservoir_stress_[wellcell];};
private:
    std::string outputprefix_;
    std::string casename_;
    std::string name_;
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    // should probably be separated for easier testing
    std::vector<int> reservoir_cells_;
    // reservoir data
    std::vector<Dune::FieldVector<double,6>> reservoir_stress_;
    std::vector<double> reservoir_pressure_;
    std::vector<double> reservoir_temperature_;
    // well data
    std::vector<double> perf_pressure_;

    static constexpr int VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_;
};
}
