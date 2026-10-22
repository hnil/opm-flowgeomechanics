/*
  Copyright 2026 Equinor ASA.

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

// Type tag for the fracture-enabled TPSA geomechanics simulator: the
// upstream TPSA flow/mechanics coupling (mirrors flow/flow_blackoil_tpsa.cpp)
// plus the fracture problem, nonlinear system and well model from
// opm-flowgeomechanics.
#pragma once

#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPSA.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/geomech/FlowProblemGeoMechTPSA.hpp>
#include <opm/geomech/NonlinearSystemBlackOilReservoirGeoMechTPSA.hpp>
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>

namespace Opm
{
    namespace Properties
    {
        namespace TTag
        {
            struct FlowProblemMechTpsa {
                using InheritsFrom = std::tuple<FlowProblem, FlowProblemTpsa>;
            };
        }

        // Flow side: TPFA assembly, as in the upstream TPSA type tag
        template<class TypeTag>
        struct Linearizer<TypeTag, TTag::FlowProblemMechTpsa>
        { using type = TpfaLinearizer<TypeTag>; };

        template<class TypeTag>
        struct LocalResidual<TypeTag, TTag::FlowProblemMechTpsa>
        { using type = BlackOilLocalResidualTPFA<TypeTag>; };

        template<class TypeTag>
        struct EnableDiffusion<TypeTag, TTag::FlowProblemMechTpsa>
        { static constexpr bool value = false; };

        template<class TypeTag>
        struct EnableDisgasInWater<TypeTag, TTag::FlowProblemMechTpsa>
        { static constexpr bool value = false; };

        // NOTE: upstream flow_blackoil_tpsa sets AvoidElementContext=true,
        // but the fully implicit energy module is not compatible with it;
        // keep the element-context path like the upstream energy targets.

        // Thermal (fully implicit energy)
        template<class TypeTag>
        struct EnableEnergy<TypeTag, TTag::FlowProblemMechTpsa>
        { static constexpr bool value = true; };

        template<class TypeTag>
        struct EnergyModuleType<TypeTag, TTag::FlowProblemMechTpsa>
        { static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

        // Mechanics + fracture
        template<class TypeTag>
        struct EnableMech<TypeTag, TTag::FlowProblemMechTpsa>
        { static constexpr bool value = true; };

        template<class TypeTag>
        struct Problem<TypeTag, TTag::FlowProblemMechTpsa>
        { using type = FlowProblemGeoMechTPSA<TypeTag>; };

        template<class TypeTag>
        struct NonlinearSystem<TypeTag, TTag::FlowProblemMechTpsa>
        { using type = NonlinearSystemBlackOilReservoirGeoMechTPSA<TypeTag>; };

        template<class TypeTag>
        struct WellModel<TypeTag, TTag::FlowProblemMechTpsa>
        { using type = BlackoilGeoMechWellModel<TypeTag>; };
    }
}
