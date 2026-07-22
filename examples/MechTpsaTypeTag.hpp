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
