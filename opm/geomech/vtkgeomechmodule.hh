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
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>


#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm::Properties {

// create new type tag for the VTK tracer output
namespace TTag {
struct VtkGeoMech {};
}

// create the property tags needed for the tracer model
template<class TypeTag, class MyTypeTag>
struct VtkWriteGeoMech {
    using type = UndefinedProperty;
};

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteGeoMech<TypeTag, TTag::VtkGeoMech> {
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
            Parameters::registerParam<TypeTag, Properties::VtkWriteGeoMech>(
                                 "Include the tracer concentration "
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
                //this->resizeVectorBuffer_(symstress_,ParentType::BufferType::ElementBuffer, 6);
            }

        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void processElement(const ElementContext& elemCtx)
        {
            if (!Parameters::get<TypeTag, Properties::EnableVtkOutput>())
                return;

            const auto& geoMechModel = elemCtx.problem().geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                pressDiff_[globalDofIdx] = geoMechModel.pressureDiff(globalDofIdx);
                //
                auto& stress = stress_[globalDofIdx];
                for(int i=0; i< 3; ++i){
                    stress[i][i] =geoMechModel.stress(globalDofIdx,i);
                }
                // voit notation converion
                stress[0][1] = geoMechModel.stress(globalDofIdx,5); //xy
                stress[0][2] = geoMechModel.stress(globalDofIdx,4); //xz
                stress[1][2] = geoMechModel.stress(globalDofIdx,3); //yz
                // fix symmetry of tensor
                for(int i=0; i< 3; ++i){
                    for(int j=0; j < 3; ++j){
                        if(i > j){
                            stress[i][j] = stress[j][i];
                        }
                    }
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
            }
        }



    

    private:
        static bool geoMechOutput_(){
            static bool val = Parameters::get<TypeTag, Properties::VtkWriteGeoMech>();
            return val;
        }
        ScalarBuffer pressDiff_;
        VectorBuffer disp_;
        VectorBuffer symstress_;
        TensorBuffer stress_;
        
    };
} // namespace Opm

#endif
