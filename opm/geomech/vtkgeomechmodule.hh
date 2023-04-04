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
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        static constexpr int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
        using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;


        using ScalarBuffer = typename ParentType::ScalarBuffer;

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
            EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteGeoMech,
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
                const auto& geoMechModel = this->simulator_.problem().geoMechModel();
                //pressDiff_.resize(geoMechModel.numCells);                
                this->resizeScalarBuffer_(pressDiff_);
            }

        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void processElement(const ElementContext& elemCtx)
        {
            if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
                return;

            const auto& geoMechModel = elemCtx.problem().geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                pressDiff_[globalDofIdx] = geoMechModel.pressureDiff(globalDofIdx); 
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
                const auto& geoMechModel = this->simulator_.problem().geoMechModel();
                const std::string tmp = "pressureDiff"; 
                this->commitScalarBuffer_(baseWriter,tmp.c_str(),
                                          pressDiff_);
            }
        }



    

    private:
        static bool geoMechOutput_(){
            static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteGeoMech);
            return val;
        }
        ScalarBuffer pressDiff_;
    };
} // namespace Opm

#endif
