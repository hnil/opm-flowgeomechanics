// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::VtkEclTracerModule
 */
#ifndef OPM_VTK_GEOMECH_MODULE_HPP
#define OPM_VTK_GEOMECH_MODULE_HPP


#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/baseoutputmodule.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>


#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <cstdio>
#include <type_traits>
#include <vector>

namespace Opm::Parameters {

// set default values for what quantities to output

struct VtkWriteGeoMech {
    static constexpr bool value = true;
};

struct VtkWriteGeoMechNodal {
    static constexpr bool value = false;
};

} // namespace Opm::Properties

namespace Opm {
    /*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the tracer model's parameters.
 */
    template <class TypeTag>
    class VtkGeoMechModule : public BaseOutputModule<TypeTag>
    {
        using ParentType = BaseOutputModule<TypeTag>;

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
        using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;


        using ScalarBuffer = typename ParentType::ScalarBuffer;
        using VectorBuffer = typename ParentType::VectorBuffer;
        using TensorBuffer = typename ParentType::TensorBuffer;
        using Tensor = Dune::DynamicMatrix<double>;
        using SymTensor = Dune::FieldVector<double,6>;
    public:
        VtkGeoMechModule(const Simulator& simulator)
            : ParentType(simulator)
        { }

        /*!
     * \brief Register all run-time parameters for the tracer VTK output
     * module.
     */
        static void registerParameters()
        {
            Parameters::Register<Parameters::VtkWriteGeoMech>(
                                 "Include geomech quentities "
                                 "in the VTK output files");
            Parameters::Register<Parameters::VtkWriteGeoMechNodal>(
                                 "Also write stress, strain, pressure and temperature "
                                 "as point data, recovered from the cell values by a linear least-squares "
                                 "fit over the cells touching each node, and the cells' logical (0-based) ijk");
        }

        /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
        void allocBuffers()
        {
            if (geoMechOutput_()){
                this->resizeScalarBuffer_(pressDiff_, ParentType::BufferType::Element);
                this->resizeVectorBuffer_(disp_, ParentType::BufferType::Vertex);
                this->resizeTensorBuffer_(stress_, ParentType::BufferType::Element);
                this->resizeTensorBuffer_(delstress_, ParentType::BufferType::Element);
                this->resizeTensorBuffer_(strain_, ParentType::BufferType::Element);


                this->resizeScalarBuffer_(pressure_, ParentType::BufferType::Element);
                this->resizeScalarBuffer_(biot_, ParentType::BufferType::Element);
                this->resizeScalarBuffer_(ymodule_, ParentType::BufferType::Element);
                this->resizeScalarBuffer_(pratio_, ParentType::BufferType::Element);
                this->resizeScalarBuffer_(poro_, ParentType::BufferType::Element);

                if (nodalOutput_()) {
                    this->resizeTensorBuffer_(nodalStress_, ParentType::BufferType::Vertex);
                    this->resizeTensorBuffer_(nodalStrain_, ParentType::BufferType::Vertex);
                    this->resizeScalarBuffer_(nodalPressure_, ParentType::BufferType::Vertex);
                    this->resizeScalarBuffer_(nodalTemperature_, ParentType::BufferType::Vertex);
                    this->resizeScalarBuffer_(cellTemperature_, ParentType::BufferType::Element);
                    this->resizeVectorBuffer_(cellIjk_, ParentType::BufferType::Element);
                }


            }

        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void setTensor(Tensor& tensor,const SymTensor& symtensor) const{
                for(int i=0; i< 3; ++i){
                    tensor[i][i] = symtensor[i];
                }
                // voit notation converion
                tensor[0][1] = symtensor[5];//xy
                tensor[0][2] = symtensor[4];//xz
                tensor[1][2] = symtensor[3];//yz
                // fix symmetry of tensor
                for(int i=0; i< 3; ++i){
                    for(int j=0; j < 3; ++j){
                        if(i > j){
                            tensor[i][j] = tensor[j][i];
                        }
                    }
                }
        }

        void processElement(const ElementContext& elemCtx)
        {
            if (!Parameters::Get<Parameters::EnableVtkOutput>())
                return;

            const auto& problem = elemCtx.problem();
            const auto& geoMechModel = problem.geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                // get the global index of the dof
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                pressDiff_[globalDofIdx] = geoMechModel.pressureDiff(globalDofIdx);
                poro_[globalDofIdx] = problem.porosity(globalDofIdx, /*timeIdx=*/0);
                biot_[globalDofIdx] = problem.biotCoef(globalDofIdx);
                pratio_[globalDofIdx] = problem.pRatio(globalDofIdx);
                ymodule_[globalDofIdx] = problem.yModule(globalDofIdx);
                // give the pressure for the geomech model
                pressure_[globalDofIdx] = geoMechModel.pressure(globalDofIdx);
                {
                    const SymTensor& symtensor = geoMechModel.stress(globalDofIdx);
                    Tensor& stress = stress_[globalDofIdx];
                    this->setTensor(stress, symtensor);
                }
                {
                    const SymTensor& symtensor = geoMechModel.delstress(globalDofIdx);
                    Tensor& delstress = delstress_[globalDofIdx];
                    this->setTensor(delstress, symtensor);
                }
                {
                    const SymTensor& symtensor = geoMechModel.strain(globalDofIdx);
                    Tensor& strain = strain_[globalDofIdx];
                    this->setTensor(strain, symtensor);
                }
                if (nodalOutput_()) {
                    std::array<int, 3> ijk;
                    elemCtx.simulator().vanguard().cartesianIndexMapper().cartesianCoordinate(globalDofIdx, ijk);
                    for (int d = 0; d < 3; ++d) {
                        cellIjk_[globalDofIdx][d] = ijk[d];
                    }
                    cellTemperature_[globalDofIdx] =
                        getValue(elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().temperature(/*phaseIdx=*/0));
                }
            }
            // all vertices proably do it to many times for now
            auto gv  = elemCtx.gridView();
            auto elem = elemCtx.element();
            static constexpr int Dim = 3;
            for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<Dim>{})){
                auto index = gv.indexSet().index(vertex);
                disp_[index] = geoMechModel.displacement(index);
            }


        }

        /*!
     * \brief Add all buffers to the VTK output writer.
     */
        void commitBuffers(BaseOutputWriter& baseWriter)
        {
            VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
            if (!vtkWriter)
                return;

            if (geoMechOutput_()){
                {
                    const std::string tmp = "pressureDiff";
                    this->commitScalarBuffer_(baseWriter,tmp.c_str(),
                                              pressDiff_,
                                              ParentType::BufferType::Element);
                }
                {
                    const std::string tmp = "disp";
                    this->commitVectorBuffer_(baseWriter,tmp.c_str(),
                                              disp_, ParentType::BufferType::Vertex);
                }
                {
                    const std::string tmp = "stress";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              stress_,
                                              ParentType::BufferType::Element);
                }
                {
                    const std::string tmp = "delstress";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              delstress_,
                                              ParentType::BufferType::Element);
                }
                {
                    const std::string tmp = "strain";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              strain_,
                                              ParentType::BufferType::Element);
                }
                if (nodalOutput_()) {
                    recoverToNodes_();
                    this->commitTensorBuffer_(baseWriter, "stress", nodalStress_, ParentType::BufferType::Vertex);
                    this->commitTensorBuffer_(baseWriter, "strain", nodalStrain_, ParentType::BufferType::Vertex);
                    this->commitScalarBuffer_(baseWriter, "pressure", nodalPressure_, ParentType::BufferType::Vertex);
                    this->commitScalarBuffer_(baseWriter, "temperature", nodalTemperature_, ParentType::BufferType::Vertex);
                    this->commitVectorBuffer_(baseWriter, "ijk", cellIjk_, ParentType::BufferType::Element);
                }

            }
        }





    private:
        static bool geoMechOutput_(){
            static bool val = Parameters::Get<Parameters::VtkWriteGeoMech>();
            return val;
        }
        static bool nodalOutput_(){
            static bool val = geoMechOutput_() && Parameters::Get<Parameters::VtkWriteGeoMechNodal>();
            return val;
        }

        // Faces of the cell, so hanging nodes of edge-conformal corner-point
        // grids also get the cells they lie on, not only those they are a corner of.
        std::vector<int> cellNodes_(const Grid& grid, const typename GridView::template Codim<0>::Entity& elem,
                                    int cellIdx) const
        {
            std::vector<int> nodes;
            if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {
                for (int f = 0; f < grid.numCellFaces(cellIdx); ++f) {
                    const int face = grid.cellFace(cellIdx, f);
                    if (grid.numFaceVertices(face) < 3) {
                        continue; // NNC pseudo-face
                    }
                    for (int v = 0; v < grid.numFaceVertices(face); ++v) {
                        nodes.push_back(grid.faceVertex(face, v));
                    }
                }
                std::ranges::sort(nodes);
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            } else {
                const auto& gv = this->simulator_.gridView();
                for (const auto& vertex : Dune::subEntities(elem, Dune::Codim<Grid::dimension>{})) {
                    nodes.push_back(gv.indexSet().index(vertex));
                }
            }
            return nodes;
        }

        // First step of patch recovery: a linear least-squares fit through the
        // centroids of the cells touching the node, evaluated at the node.
        void recoverToNodes_()
        {
            const auto& grid = this->simulator_.vanguard().grid();
            const auto& gv = this->simulator_.gridView();
            std::vector<Dune::FieldVector<double, 3>> nodePosition(gv.size(Grid::dimension));
            for (const auto& vertex : Dune::vertices(gv)) {
                nodePosition[gv.indexSet().index(vertex)] = vertex.geometry().center();
            }
            std::vector<std::vector<std::pair<int, Dune::FieldVector<double, 3>>>> patch(nodePosition.size());
            for (const auto& elem : Dune::elements(gv)) {
                const int cell = gv.indexSet().index(elem);
                for (const int node : cellNodes_(grid, elem, cell)) {
                    patch[node].emplace_back(cell, elem.geometry().center());
                }
            }
            for (std::size_t node = 0; node < patch.size(); ++node) {
                if (patch[node].empty()) {
                    continue;
                }
                const auto weight = recoveryWeights_(nodePosition[node], patch[node]);
                for (std::size_t k = 0; k < weight.size(); ++k) {
                    const int cell = patch[node][k].first;
                    Tensor stress = stress_[cell];
                    Tensor strain = strain_[cell];
                    stress *= weight[k];
                    strain *= weight[k];
                    nodalStress_[node] += stress;
                    nodalStrain_[node] += strain;
                    nodalPressure_[node] += weight[k] * pressure_[cell];
                    nodalTemperature_[node] += weight[k] * cellTemperature_[cell];
                }
            }
        }

        // Weights w_c with v(x) = sum w_c v_c for the fit v = a + b.(c - x). The small
        // penalty on b keeps directions the centroids do not span (boundaries, one or
        // two cells) flat instead of extrapolating.
        static std::vector<double>
        recoveryWeights_(const Dune::FieldVector<double, 3>& x,
                         const std::vector<std::pair<int, Dune::FieldVector<double, 3>>>& patch)
        {
            Dune::FieldMatrix<double, 4, 4> normal(0.0);
            std::vector<Dune::FieldVector<double, 4>> rows;
            double spread = 0.0;
            for (const auto& [cell, centroid] : patch) {
                Dune::FieldVector<double, 4> q{1.0, centroid[0] - x[0], centroid[1] - x[1], centroid[2] - x[2]};
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        normal[i][j] += q[i] * q[j];
                    }
                }
                spread += q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
                rows.push_back(q);
            }
            const double penalty = 1e-3 * spread / patch.size();
            for (int i = 1; i < 4; ++i) {
                normal[i][i] += penalty;
            }
            Dune::FieldVector<double, 4> r;
            normal.solve(r, Dune::FieldVector<double, 4>{1.0, 0.0, 0.0, 0.0});
            std::vector<double> weight;
            for (const auto& q : rows) {
                weight.push_back(r * q);
            }
            return weight;
        }

        ScalarBuffer pressDiff_;
        VectorBuffer disp_;
        VectorBuffer symstress_;
        TensorBuffer stress_;
        TensorBuffer delstress_;
        TensorBuffer strain_;

        ScalarBuffer pressure_;
        ScalarBuffer biot_;
        ScalarBuffer ymodule_;
        ScalarBuffer pratio_;
        ScalarBuffer poro_;

        TensorBuffer nodalStress_;
        TensorBuffer nodalStrain_;
        ScalarBuffer nodalPressure_;
        ScalarBuffer nodalTemperature_;
        ScalarBuffer cellTemperature_;
        VectorBuffer cellIjk_;
    };
} // namespace Opm

#endif
