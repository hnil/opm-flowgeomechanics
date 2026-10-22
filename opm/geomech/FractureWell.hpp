#pragma once

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>

#include <array>
#include <cassert>
#include <optional>
#include <utility>
#include <memory>
#include <vector>

namespace Opm::Properties {
    template <class TypeTag, class MyType>
    struct FluidSystem;
} // namespace Opm::Properties

namespace Opm {

class FractureWell
{
    using Point3D = Dune::FieldVector<double,3>;
    using Segment = std::array<unsigned,2>;

public:
    using Grid = Dune::FoamGrid<1,3>;

    struct Connection
    {
        int cell{};
        int segment{};
        std::optional<std::pair<double, double>> perf_range{};
    };

    FractureWell(const std::string& outputprefix,
                 const std::string& casename,
                 const std::string& name,
                 const std::vector<Connection>& conns,
                 const std::vector<Point3D>& points,
                 const std::vector<Segment>& segments);

    void init(const std::vector<Point3D>& points,
              const std::vector<Segment>& segments);

    template<class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator)
    {
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        assert(grid_->leafGridView().size(0) == this->conns_.size());

        this->perf_pressure_.resize(this->conns_.size());
        this->perf_active_.resize(this->conns_.size(), false);
        this->reservoir_stress_.resize(this->conns_.size());
        this->reservoir_pressure_.resize(this->conns_.size());
        this->reservoir_temperature_.resize(this->conns_.size());

        for (auto i = 0*this->conns_.size(); i < this->conns_.size(); ++i) {
            reservoir_stress_[i] = simulator.problem().stress(this->reservoirCell(i));

            const auto& intQuants = simulator.model()
                .intensiveQuantities(this->reservoirCell(i), /*timeIdx*/ 0);

            const auto& fs = intQuants.fluidState();

            reservoir_pressure_[i] = fs.pressure(FluidSystem::waterPhaseIdx).value();

            const unsigned dummy = 0u;
            reservoir_temperature_[i] = fs.temperature(dummy).value();
        }
    }
    void setPerfPressure(int perf_index, double pressure){perf_pressure_[perf_index] = pressure;}
    const Grid& grid() const{return *grid_;}
    const std::string& name() const{return name_;}
    void write() const;
    void writemulti(double time) const;
    void resetWriters();
    int reservoirCell(int wellcell) const { return this->conns_[wellcell].cell; }
    int segment(const int wellcell) const { return this->conns_[wellcell].segment; }
    decltype(auto) perfRange(const int wellcell) const { return this->conns_[wellcell].perf_range; }
    const Dune::FieldVector<double,6>& reservoirStress(int wellcell) const { return this->reservoir_stress_[wellcell]; }
    bool isActive() const {return active_;};
    void setActive(bool active) {active_ = active;};
    void setPerfActive(int perf_index, bool active)
    {
        assert(perf_index>= 0);
        assert(perf_index < static_cast<int>(perf_active_.size()));
        perf_active_[perf_index] = active;
    }

private:
    std::string outputprefix_;
    std::string casename_;
    std::string name_;

    // should probably be separated for easier testing
    std::vector<Connection> conns_;

    std::unique_ptr<Grid> grid_{};
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_{};

    // reservoir data
    std::vector<Dune::FieldVector<double,6>> reservoir_stress_{};
    std::vector<double> reservoir_pressure_{};
    std::vector<double> reservoir_temperature_{};

    // well data
    std::vector<double> perf_pressure_{};
    std::vector<double> perf_active_{};
    bool active_ = false; // if the well is active
  //
  // welldata : state, WI, ?
    static constexpr auto VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_{};
};

} // namespace Opm
