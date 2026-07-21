#pragma once
#include <opm/simulators/wells/BlackoilWellModel.hpp>
namespace Opm
{
template <typename TypeTag>
class BlackoilGeoMechWellModel : public BlackoilWellModel<TypeTag>
{
    using Parent = BlackoilWellModel<TypeTag>;
    using Simulator = typename Parent::Simulator;
public:
    BlackoilGeoMechWellModel(Simulator& simulator):
    Parent(simulator)
    {

    };
    //using BlackoilWellModel::BlackoilWellModel;

    using NeighborSet = typename Parent::NeighborSet;
    void addNeighbors(std::vector<NeighborSet>& /*neighbors*/) const
    {
        if (!this->param_.matrix_add_well_contributions_) {
            return;
        }
        OPM_THROW(std::runtime_error, "Not implemented");
    };
    void createWellContainer(const int reportStepIdx)
    {
        Parent::createWellContainer(reportStepIdx);
        // only add effect of fracture after one report step
        // NB everything is not explicit and ministeps are not considered
        if (reportStepIdx > 0) {
            const auto& problem = this->simulator_.problem();
            const auto& geoMechModel = problem.geoMechModel();
            if (problem.hasFractures() && geoMechModel.fractureModelActive()) {
                for (auto& wellPtr : this->well_container_) {
                    auto wellName = wellPtr->name();
                    const auto& fracturemodel = geoMechModel.fractureModel();
                    auto wellIndices = fracturemodel.getExtraWellIndices(wellName);
                    wellPtr->addFracturePerforations(wellIndices);
                }
            }
        }
    };
};
} // namespace Opm
