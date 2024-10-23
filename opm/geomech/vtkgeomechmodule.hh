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
#ifndef VTK_GEOMECH_MODULE_HH
#define VTK_GEOMECH_MODULE_HH


#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/baseoutputmodule.hh>

#include <opm/models/utils/propertysystem.hh>
//#include <opm/models/utils/parametersystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>


#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm::Parameters {

// create new type tag for the VTK tracer output
// create the property tags needed for the tracer model
// set default values for what quantities to output

struct VtkWriteGeoMech {
    static constexpr bool value = true;
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

        static constexpr int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
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
        }

        /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
        void allocBuffers()
        {
            if (geoMechOutput_()){
                //pressDiff_.resize(geoMechModel.numCells);
                this->resizeScalarBuffer_(pressDiff_);
                this->resizeVectorBuffer_(disp_,ParentType::BufferType::VertexBuffer);
                this->resizeTensorBuffer_(stress_);
                this->resizeTensorBuffer_(delstress_);
                this->resizeTensorBuffer_(strain_);
                //this->resizeVectorBuffer_(symstress_,ParentType::BufferType::ElementBuffer, 6);
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

            const auto& geoMechModel = elemCtx.problem().geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                pressDiff_[globalDofIdx] = geoMechModel.pressureDiff(globalDofIdx);
                //
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
                                              pressDiff_);
                }
                {
                    const std::string tmp = "disp";
                    this->commitVectorBuffer_(baseWriter,tmp.c_str(),
                                              disp_, ParentType::BufferType::VertexBuffer);
                }
                {
                    const std::string tmp = "stress";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              stress_,
                                              ParentType::BufferType::ElementBuffer);
                }
                {
                    const std::string tmp = "delstress";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              delstress_,
                                              ParentType::BufferType::ElementBuffer);
                }
                {
                    const std::string tmp = "strain";
                    this->commitTensorBuffer_(baseWriter,tmp.c_str(),
                                              strain_,
                                              ParentType::BufferType::ElementBuffer);
                }

            }
        }





    private:
        static bool geoMechOutput_(){
            static bool val = Parameters::Get<Parameters::VtkWriteGeoMech>();
            return val;
        }
        ScalarBuffer pressDiff_;
        VectorBuffer disp_;
        VectorBuffer symstress_;
        TensorBuffer stress_;
        TensorBuffer delstress_;
        TensorBuffer strain_;
    };
} // namespace Opm

#endif
