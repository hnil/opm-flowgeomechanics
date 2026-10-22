/*
  Copyright 2025 Equinor ASA.

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

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>
#include <opm/geomech/NonlinearSystemBlackOilReservoirGeoMech.hpp>
namespace Opm
{
    namespace Properties
    {
        namespace TTag
        {
            struct FlowProblemMech {
	      using InheritsFrom = std::tuple<FlowProblem> ;//?? should it be blackoil
            };
        }

        // Set the problem class
        template <class TypeTag>
        struct Problem<TypeTag, TTag::FlowProblemMech> {
            using type = FlowProblemGeoMech<TypeTag>;
        };

        template <class TypeTag>
        struct NonlinearSystem<TypeTag, TTag::FlowProblemMech> {
              using type = NonlinearSystemBlackOilReservoirGeoMech<TypeTag>;
        };

        template <class TypeTag>
        struct EnableMech<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = true;
        };

        template <class TypeTag>
        struct EnableEnergy<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = true;
        };

        template<class TypeTag>
        struct EnergyModuleType<TypeTag, TTag::FlowProblemMech>
        { static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

        // Per-DOF assembly: the element-driven linearizer skips DOFs without a grid element
        // (numerical aquifers, auxiliary fracture cells) and leaves their rows empty.
        template<class TypeTag>
        struct Linearizer<TypeTag, TTag::FlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

        template<class TypeTag>
        struct LocalResidual<TypeTag, TTag::FlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

        template <class TypeTag>
        struct EnableDiffusion<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = false;
        };

        template <class TypeTag>
        struct EnableDisgasInWater<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = false;
        };


        template <class TypeTag>
        struct WellModel<TypeTag, TTag::FlowProblemMech> {
            using type = BlackoilGeoMechWellModel<TypeTag>;
        };

        template <class TypeTag>
        struct Simulator<TypeTag, TTag::FlowProblemMech> {
            using type = Opm::Simulator<TypeTag>;
        };

    }
}
namespace Opm {
namespace Parameters {
}
}
