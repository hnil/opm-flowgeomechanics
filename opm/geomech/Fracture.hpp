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
                  Point3D normal);
        void grow(int layers,int method);
        std::string name() const;
        void write() const;
        template<class Grid3D>
        void updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                                  Grid3D& grid3D);

        // solver related
        void updateReservoirProperties();

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


        std::vector<int> well_source_;
        // for reservoir
        std::vector<int> reservoir_cells_;
        std::vector<double> reservoir_perm_;
        std::vector<double> reservoir_dist_;
        std::vector<double> reservoir_pressure_;

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
    };
}
#endif
