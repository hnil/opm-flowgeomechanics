/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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
#ifndef OPM_FRACTURE_MODEL_HH
#define OPM_FRACTURE_MODEL_HH
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <opm/grid/CpGrid.hpp>
// from pdelab
#include "GeometryHelpers.hpp"
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>

// for linear solve
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
// #include <opm/models/utils/parametersystem.hh>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <dune/istl/matrixmarket.hh>

#include <opm/geomech/GridStretcher.hpp>

namespace Opm
{
struct WellInfo {
    std::string name;
    int perf;
    int well_cell;
};

/// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
class Fracture
{
public:
    using Grid = Dune::FoamGrid<2, 3>;
    using Point3D = Dune::FieldVector<double, 3>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using DynamicMatrix = Dune::DynamicMatrix<double>;
  
    void init(std::string well, int perf, int well_cell, Point3D origo,
              Point3D normal, Opm::PropertyTree prm);
    void grow(int layers, int method);
    std::string name() const;
    void write(int reportStep = -1) const;
    void writemulti(double time) const;
    template <class Grid3D>
    void updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                              const Grid3D& grid3D);

    // solver related
    void updateReservoirProperties();
    void removeCells();

    template <class TypeTag, class Simulator>
    void updateReservoirProperties(const Simulator& simulator, bool init_constant_vals=false)
    {
        // if `init_contant_vals` is true, the fields that should normally not
        // change, i.e.  E_, nu_ and reservoir_perm_ will be updated. This should
        // normally only be needed the first time.
        if (init_constant_vals) {
            initReservoirProperties<TypeTag, Simulator>(simulator);
        }
      
        std::cout << "updateReservoirProperties (simulator)" << std::endl;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto& problem = simulator.problem();
        // NB burde truleg interpolere
        // NB reservoir dist not calculated
        size_t ncf = reservoir_cells_.size();
        assert(ncf > 0); 
        reservoir_pressure_.resize(ncf);
        reservoir_stress_.resize(ncf);
        reservoir_mobility_.resize(ncf);
        reservoir_perm_.resize(ncf);
        
        // should be calculated
        double dist = prm_.get<double>("reservoir.dist");
        reservoir_dist_.resize(ncf, dist);
        
        double numax = 0, Emax = 0;
        
        for (size_t i = 0; i < ncf; ++i) {
            size_t cell = reservoir_cells_[i];
            if (!(cell < 0)) {
                //auto normal = this->cell_normals_[i];
                const auto& intQuants = simulator.model().intensiveQuantities(cell, /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                {
                  auto val = fs.pressure(FluidSystem::waterPhaseIdx);
                  reservoir_pressure_[i] = val.value();
                }
                enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                  if (FluidSystem::phaseIsActive(phaseIdx)) {
                    // assume sum should only be water;
                    auto val = intQuants.mobility(phaseIdx);
                    reservoir_mobility_[i] += val.value();
                  }
                }
                for (int dim = 0; dim < 3; ++dim) {
                  reservoir_stress_[i] = problem.stress(cell);
                }
            } else {
                // probably outside reservoir set all to zero
                double stressval = 0;
                reservoir_stress_[i][0] = stressval;
                reservoir_stress_[i][1] = stressval;
                reservoir_stress_[i][2] = stressval; //???
                reservoir_pressure_[i] = 0.0;
                reservoir_mobility_[i] = 0.0;
                reservoir_perm_[i] = 0.0;
            }
            // assume reservoir distance is calculated
        }
    }
  
    void initFracturePressureFromReservoir();
    void initFractureStates();
    void initFractureWidth();
    void solveFractureWidth();
    void solvePressure();
    void solve();

    void printPressureMatrix() const; // debug purposes
    void printMechMatrix() const; // debug purposes
    void writeFractureSystem()  const;
    void writePressureSystem()  const;
    void setFractureGrid(std::unique_ptr<Fracture::Grid> gptr = nullptr); // a hack to allow use of another grid
    std::vector<std::tuple<int, double, double>> wellIndices() const;
    WellInfo& wellInfo(){return wellinfo_;}
    std::vector<double> leakOfRate() const;
    double injectionPressure() const;
    void setPerfPressure(double perfpressure){perf_pressure_ = perfpressure;}
private:

    template <class TypeTag, class Simulator>
    void initReservoirProperties(const Simulator& simulator)
    {
      std::cout << "initReservoirProperties (simulator)" << std::endl;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto& problem = simulator.problem();
        // NB burde truleg interpolere
        // NB reservoir dist not calculated
        size_t ncf = reservoir_cells_.size();
        reservoir_perm_.resize(ncf);
        // should be calcualted
        double dist = prm_.get<double>("reservoir.dist");
        reservoir_dist_.resize(ncf, dist);
        double numax = -1e99;
        double Emax = -1e99;
        assert(ncf>0);
        for (size_t i = 0; i < ncf; ++i) {
            size_t cell = reservoir_cells_[i];
            auto normal = this->cell_normals_[i];
            {
                auto permmat = problem.intrinsicPermeability(cell);
                auto np = normal;
                permmat.mv(normal, np);
                double value = np.dot(normal);
                reservoir_perm_[i] = value;
            }
            Emax = std::max(Emax,problem.yModule(cell));
            numax = std::max(numax,problem.pRatio(cell));
        }
        E_ = Emax;
        nu_ = numax;
        assert(E_ > 0);
        assert(nu_>0 && nu_ < 1);
    }

  void resetWriters();
    // helpers for growing grid
    void insertLinear(const std::vector<unsigned int>& inner_indices);
    void insertExp(const std::vector<unsigned int>& inner_indices);
    void initFracture(); // create a new fracture grid from scratch
    Point3D surfaceMap(double x, double y);

    // should be called every time a grid is updated/changed
    //void updateGridDiscretizations() {assembleFractureMatrix(); initPressureMatrix();}

    std::unique_ptr<GridStretcher> grid_stretcher_;
  
    std::unique_ptr<Grid> grid_;
    Point3D origo_;
    std::array<Point3D, 3> axis_;
    WellInfo wellinfo_;
    std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
    static constexpr int VTKFormat = Dune::VTK::ascii;
    std::unique_ptr<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>> vtkmultiwriter_;
    std::vector<unsigned int> out_indices_;

    // should probably not be needed
    int layers_;
    int nlinear_;

    // help function for solving
    void assemblePressure();
    void setSource();
    void initPressureMatrix();
    void setupPressureSolver();
    void updateFractureRHS();
    void updateLeakoff();
    void updateCellNormals();
    void normalFractureTraction(Dune::BlockVector<Dune::FieldVector<double, 1>>& traction) const;
    double normalFractureTraction(size_t ix) const;

    // one nonlinear iteration of fully coupled system.  Returns 'true' if converged
    bool fullSystemIteration(const double tol);

    void assembleFractureMatrix() const;
    std::vector<double> stressIntensityK1() const;


    //double well_pressure_;// for now using prm object for definition
    std::vector<int> well_source_;
    // for reservoir
    std::vector<int> reservoir_cells_;
    // std::vector< Dune::FieldMatrix<double, 3, 3> > reservoir_perm_;
    std::vector<double> reservoir_perm_;
    std::vector<double> reservoir_mobility_;
    std::vector<double> reservoir_dist_;
    std::vector<double> reservoir_pressure_;
    std::vector<Dune::FieldVector<double, 6>> reservoir_stress_;

    // only for radom access need to be updater after trid change
    Dune::BlockVector<Dune::FieldVector<double, 3>> cell_normals_;

    // solution variables (only to avoid memory allocation, do not trust their state)
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_width_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> fracture_pressure_;
    mutable Dune::BlockVector<Dune::FieldVector<double, 1>> rhs_pressure_;

    // transmissibilities
    using Htrans = std::tuple<size_t, size_t, double, double>;
    std::vector<Htrans> htrans_;
    std::vector<std::tuple<int,double>> perfinj_;
    double perf_pressure_;
    std::vector<double> leakof_;
    //
    Opm::PropertyTree prmpressure_;
    using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
    mutable std::unique_ptr<Matrix> pressure_matrix_;
    mutable std::unique_ptr<PressureOperatorType> pressure_operator_;
    mutable std::unique_ptr<FlexibleSolverType> pressure_solver_;
    mutable std::unique_ptr<Matrix> coupling_matrix_; // will be updated by `fullSystemIteration`
  
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    // using DynamicMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
    mutable std::unique_ptr<DynamicMatrix> fracture_matrix_; 

    // function ensuring that the fracture matrix exists, and returning a reference to it
    DynamicMatrix& fractureMatrix() const {
        if (fracture_matrix_ == nullptr)
            assembleFractureMatrix();
        return *fracture_matrix_;
    }
  
    double E_;
    double nu_;
    Opm::PropertyTree prm_;
};
} // namespace Opm
#endif
