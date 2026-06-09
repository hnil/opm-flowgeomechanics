#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>
#include <opm/geomech/BlackoilGeomechWellModel.hpp>
#include <opm/geomech/BlackoilModelGeomech.hpp>

#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {
namespace Properties {
namespace TTag {

struct EclFlowProblemOilGasEnergyMech {
    using InheritsFrom = std::tuple<FlowProblem>;
};

} // namespace TTag

template <class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = EclProblemGeoMech<TypeTag>;
};

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = BlackoilModelGeomech<TypeTag>;
};

template <class TypeTag>
struct EnableMech<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::EclFlowProblemOilGasEnergyMech>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

template<class TypeTag>
struct Indices<TypeTag, TTag::EclFlowProblemOilGasEnergyMech>
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
struct Linearizer<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = TpfaLinearizer<TypeTag>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = BlackOilLocalResidualTPFA<TypeTag>;
};

template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct WellModel<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = BlackoilGeomechWellModel<TypeTag>;
};

template <class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemOilGasEnergyMech> {
    using type = Opm::Simulator<TypeTag>;
};

} // namespace Properties
} // namespace Opm