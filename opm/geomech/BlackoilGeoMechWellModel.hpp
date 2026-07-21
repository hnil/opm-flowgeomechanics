#pragma once
#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
namespace Opm
{
template <typename TypeTag>
class BlackoilGeoMechWellModel : public BlackoilWellModel<TypeTag>
{
    using Parent = BlackoilWellModel<TypeTag>;
    using Simulator = typename Parent::Simulator;
public:
    BlackoilGeoMechWellModel(Simulator& simulator,
                             const NewtonIterationContext& iter_ctx):
    Parent(simulator, iter_ctx)
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
    // Fracture-created connections are added to the schedule at runtime and
    // may not (yet) exist on every rank that hosts part of the well; log and
    // continue instead of aborting the run.
    bool continueOnMissingWellConnections() const override
    {
        return true;
    }

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
