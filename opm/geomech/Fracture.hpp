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
#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <functional>
#include <string>
#include <iostream>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <opm/grid/CpGrid.hpp>
#include <dune/foamgrid/foamgrid.hh>
// from pdelab
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>
#include "GeometryHelpers.hpp"

// for linear solve
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/blackoil/blackoilmodel.hh>
//#include <opm/models/utils/parametersystem.hh>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
namespace Opm {
    struct WellInfo{
        std::string name;
        int perf;
        int well_cell;

    };


    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    class Fracture
    {
    public:
        using Grid =  Dune::FoamGrid<2, 3>;
        using Point3D =  Dune::FieldVector<double, 3>;
        void init(std::string well,
                  int perf,
                  int well_cell,
                  Point3D origo,
                  Point3D normal,
                  Opm::PropertyTree prm
            );
        void grow(int layers,int method);
        std::string name() const;
        void write(int reportStep = -1) const;
        void writemulti(double time) const;
        template<class Grid3D>
        void updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                                  const Grid3D& grid3D);

        // solver related
        void updateReservoirProperties();
        template<class TypeTag, class Simulator>
        void updateReservoirProperties(const Simulator& simulator){
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            const auto& problem = simulator.problem();
            //NB burde truleg interpolere
            size_t ncf = reservoir_cells_.size();
            reservoir_perm_.resize(ncf);
            reservoir_pressure_.resize(ncf);
            reservoir_stress_.resize(ncf);
            for(size_t i=0; i < ncf; ++i){
                int cell = reservoir_cells_[i];
                auto normal = this->cell_normals_[i];
                {

                    auto permmat = problem.intrinsicPermeability(cell);
                    auto np = permmat.mv(normal);
                    double value = np.dot(normal);
                    reservoir_perm_[i]  = value;
                } 
                const auto &intQuants = simulator.model().intensiveQuantities(cell, /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();
                reservoir_pressure_[i]  = fs.pressure(FluidSystem::waterPhaseIdx);
                enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)){
                        // assume sum should only be water;
                        reservoir_mobility_[i]  += fs.mobility(phaseIdx);   
                    }
                }
                for(int dim=0; dim < 3; ++ dim){
                    reservoir_stress_[i][dim] = problem.model().stress(cell, dim);
                }
                // assume reservoir distance is calculated
            }
        }
        void initFractureWidth();
        void solveFractureWidth();
        void solvePressure();
        void solve();

      void printPressureMatrix() const; // debug purposes
      void printMechMatrix() const; // debug purposes
      void setFractureGrid(std::unique_ptr<Fracture::Grid>& gptr); // a hack to allow use of another grid
    private:
        // helpers for growing grid
        void insertLinear(const std::vector<unsigned int>& inner_indices);
        void insertExp(const std::vector<unsigned int>& inner_indices);
        void initFracture();
        Point3D surfaceMap(double x, double y);

        std::unique_ptr<Grid> grid_;
        Point3D origo_;
        std::array<Point3D,3> axis_;
        WellInfo wellinfo_;
        std::unique_ptr<Dune::VTKWriter<Grid::LeafGridView>> vtkwriter_;
        static constexpr int VTKFormat = Dune::VTK::ascii;
        std::unique_ptr<Opm::VtkMultiWriter<Grid::LeafGridView,VTKFormat>> vtkmultiwriter_;
        std::vector<unsigned int> out_indices_;

        // should probably not be needed
        int layers_;
        int nlinear_;



        // help function for solving
        void assemblePressure();
        void setSource();
        void initPressureMatrix();
        void setupPressureSolver();


        void assembleFracture();
        std::vector<double> stressIntensityK1() const;


        std::vector<int> well_source_;
        // for reservoir
        std::vector<int> reservoir_cells_;
        //std::vector< Dune::FieldMatrix<double, 3, 3> > reservoir_perm_;
        std::vector< double > reservoir_perm_;
        std::vector<double> reservoir_mobility_;
        std::vector<double> reservoir_dist_;
        std::vector<double> reservoir_pressure_;
        std::vector<Dune::FieldVector<double,6> > reservoir_stress_;

        // only for radom access need to be updater after trid change
        Dune::BlockVector<Dune::FieldVector<double,3>> cell_normals_;

        // solution variables
        Dune::BlockVector<Dune::FieldVector<double,1>> fracture_width_;
        Dune::BlockVector<Dune::FieldVector<double,1>> rhs_width_;
        Dune::BlockVector<Dune::FieldVector<double,1>> fracture_pressure_;
        Dune::BlockVector<Dune::FieldVector<double,1>> rhs_pressure_;

        //transmissibilities
        using Htrans = std::tuple<size_t,size_t, double, double>;
        std::vector<Htrans> htrans_ ;
        std::vector<double> leakof_ ;

        //
        Opm::PropertyTree prmpressure_;
        using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
        using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
        using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
        std::unique_ptr<Matrix> pressure_matrix_;
        std::unique_ptr<PressureOperatorType> pressure_operator_;
        std::unique_ptr<FlexibleSolverType> pressure_solver_;
        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
        //
        //using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
        //using DenseMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
        //using DynamicMatrix = Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>;
         using DynamicMatrix = Dune::DynamicMatrix<double>;
        std::unique_ptr<DynamicMatrix> fracture_matrix_;

        double E_;
        double nu_;
        Opm::PropertyTree prm_;
    };
}
#endif
