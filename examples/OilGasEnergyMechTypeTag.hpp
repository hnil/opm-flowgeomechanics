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

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>
#include <opm/geomech/NonlinearSystemBlackOilReservoirGeoMech.hpp>

#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {
namespace Properties {
namespace TTag {

struct FlowProblemOilGasEnergyMech {
    using InheritsFrom = std::tuple<FlowProblem>;
};

} // namespace TTag

template <class TypeTag>
struct Problem<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = FlowProblemGeoMech<TypeTag>;
};

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = NonlinearSystemBlackOilReservoirGeoMech<TypeTag>;
};

template <class TypeTag>
struct EnableMech<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowProblemOilGasEnergyMech>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

template<class TypeTag>
struct Indices<TypeTag, TTag::FlowProblemOilGasEnergyMech>
{
private:
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr int numEnergyVars = energyModuleType == EnergyModules::FullyImplicitThermal;

public:
    using type = BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         numEnergyVars,
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::waterCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = TpfaLinearizer<TypeTag>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = BlackOilLocalResidualTPFA<TypeTag>;
};

template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct WellModel<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = BlackoilGeoMechWellModel<TypeTag>;
};

template <class TypeTag>
struct Simulator<TypeTag, TTag::FlowProblemOilGasEnergyMech> {
    using type = Opm::Simulator<TypeTag>;
};

} // namespace Properties
} // namespace Opm