#ifndef OPM_FLOW_PROBLEM_GEOMECH_HPP
#define OPM_FLOW_PROBLEM_GEOMECH_HPP

#include <algorithm>
#include <opm/common/ErrorMacros.hpp>
#include <fmt/format.h>

#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/geomech/FlowGeoMechLinearSolverParameters.hpp>
#include <opm/geomech/FlowProblemMech.hpp>
#include <opm/geomech/FractureAuxCells.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/geomech/BoundaryUtils.hpp>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/VtkGeoMechModule.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>


#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/elasticity/material.hh>
#include <opm/elasticity/materials.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm{

    template<typename TypeTag>
    class FlowProblemGeoMech: public FlowProblemMech<TypeTag, FlowProblemBlackoil<TypeTag>, FlowProblemGeoMech<TypeTag>>{
    public:
        using MechParent = FlowProblemMech<TypeTag, FlowProblemBlackoil<TypeTag>, FlowProblemGeoMech<TypeTag>>;
        using Parent = FlowProblemBlackoil<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper = AdaptiveTimeStepping<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
        enum { dim = GridView::dimension };
        enum { dimWorld = GridView::dimensionworld };
        using Toolbox = MathToolbox<Evaluation>;
        using SymTensor = Dune::FieldVector<double,6>;

      //using CellSeedType = typename GridView::template Codim<0>::EntitySeed;
        FlowProblemGeoMech(Simulator& simulator):
            MechParent(simulator),
            geoMechModel_(simulator)
        {
            if(this->simulator().vanguard().eclState().runspec().mech()){
              this->model().addOutputModule(std::make_unique<VtkGeoMechModule<TypeTag>>(simulator));
            }
        }

    private:
        //! Owned by the base problem; this is the typed handle onto it.
        FractureAuxCells<TypeTag>* fractureAuxCells_ = nullptr;
        double embeddedCouplingChange_ = 0.0;
        bool embeddedStatic_ = false;
        // Opt-in (solver.embedded_satnum, 1-based): saturation-function region for
        // the fracture cells instead of the partner's rock table (a fracture has
        // straight-line kr and no capillary pressure). Resolved lazily to the first
        // grid cell of that region, which the material-law lookup is redirected to.
        int embeddedSatnum_ = 0;
        int embeddedCellDump_ = 0; // solver.embedded_cell_dump: worst-N aux cells per report
        mutable int embeddedSatProxyCell_ = -1;
        std::vector<std::size_t> lastSeenFractureLayout_ {}; // see bindFractureAuxCells
        bool embeddedLeakoffReport_ = false;

    public:

        /*!
         * \brief Claim degrees of freedom for the fracture cells, if they are to have any.
         *
         * The model calls this before it sizes anything, which is the only moment at
         * which degrees of freedom can still be added -- and long before any fracture
         * exists, since the fracture model is built at the first report step that seeds
         * one.  So a fixed number is reserved here and handed out as cells appear.
         */
        using MaterialLawParams = typename Parent::MaterialLawParams;

        template <class Context>
        const MaterialLawParams& materialLawParams(const Context& context,
                                                   unsigned spaceIdx, unsigned timeIdx) const
        { return this->materialLawParams(context.globalSpaceIndex(spaceIdx, timeIdx)); }

        const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
        {
            const int proxy = embeddedSatProxy_(globalDofIdx);
            return (proxy >= 0) ? this->materialLawManager()->materialLawParams(proxy)
                                : Parent::materialLawParams(globalDofIdx);
        }

        const MaterialLawParams& materialLawParams(unsigned globalDofIdx, FaceDir::DirEnum facedir) const
        {
            const int proxy = embeddedSatProxy_(globalDofIdx);
            return (proxy >= 0) ? this->materialLawManager()->materialLawParams(proxy, facedir)
                                : Parent::materialLawParams(globalDofIdx, facedir);
        }

        // relperms go through this rather than materialLawParams(), so the
        // fracture-cell redirect has to be applied here as well
        template <class FluidState, class... Args>
        void updateRelperms(std::array<GetPropType<TypeTag, Properties::Evaluation>, GetPropType<TypeTag, Properties::FluidSystem>::numPhases>& mobility,
                            typename Parent::DirectionalMobilityPtr& dirMob,
                            FluidState& fluidState,
                            unsigned globalSpaceIdx) const
        {
            const int proxy = embeddedSatProxy_(globalSpaceIdx);
            if (proxy < 0) {
                Parent::template updateRelperms<FluidState, Args...>(mobility, dirMob, fluidState, globalSpaceIdx);
                return;
            }
            using ContainerT = std::array<GetPropType<TypeTag, Properties::Evaluation>, GetPropType<TypeTag, Properties::FluidSystem>::numPhases>;
            GetPropType<TypeTag, Properties::MaterialLaw>::template relativePermeabilities<ContainerT, FluidState, Args...>
                (mobility, this->materialLawManager()->materialLawParams(proxy), fluidState);
        }

        void registerAuxiliaryCellModules()
        {
            MechParent::registerAuxiliaryCellModules();

            if (!this->hasFractures()) {
                return;
            }

            const Opm::PropertyTree prm = this->getFractureParam();
            if (prm.get<std::string>("solver.fracture_flow_mode", std::string {"wi_upscaling"})
                != "embedded")
            {
                return;
            }

            const auto capacity = prm.get<int>("solver.embedded_capacity", 5000);

            // How the well's perforations of the fracture cells get their well index.
            // "fracture" carries over the factor the fracture's own pressure solve uses
            // -- the current modelling, a constant factor on the cells around the
            // wellbore -- and stays the default so the old behaviour is preserved.
            // "estimate" forms a radial-flow index through a fixed prescribed aperture;
            // making that width follow the solved aperture is the later, dynamic step.
            const auto perfWiModeName =
                prm.get<std::string>("solver.embedded_perf_wi_mode", std::string {"fracture"});
            const auto perfWiMode = (perfWiModeName == "estimate")
                ? FractureAuxCells<TypeTag>::PerfWiMode::Estimate
                : FractureAuxCells<TypeTag>::PerfWiMode::Fracture;
            const auto perfWidth = prm.get<double>("solver.embedded_perf_width", 5e-3);
            const auto perfRw = prm.get<double>("solver.embedded_perf_rw", 0.1);

            // The static gate: bind the fracture into the flow problem once, when it
            // first appears, and hold that description for the rest of the run.  The
            // fracture model itself keeps solving -- its state is simply no longer
            // re-read -- so the flow sees a fixed, fully formed fracture.  This is the
            // validation configuration: growth feedback cannot confound a comparison of
            // the conductances themselves.
            embeddedStatic_ = prm.get<bool>("solver.embedded_static", false);
            embeddedSatnum_ = prm.get<int>("solver.embedded_satnum", 0);
            embeddedCellDump_ = prm.get<int>("solver.embedded_cell_dump", 0);
            embeddedLeakoffReport_ = prm.get<bool>("solver.embedded_leakoff_report", false);

            // The floor under the aperture used for the cells' volume and cubic-law
            // transmissibility.  Deliberately NOT config.min_width: the fracture's own
            // solver may run with that at zero, handling closure through the contact
            // machinery instead -- but a flow cell with zero volume has no storage and,
            // when it holds no oil and sees no pressure difference, an identically zero
            // oil row.  The floor keeps every open cell a well-posed, if small, volume.
            //
            // The default matches the fracture solver's own historical min_width.  Going
            // much below it runs into the absolute singularity threshold of the block
            // inverter (|det| < 1e-40, matrixblock.hh): a fracture row's entries scale
            // with aperture^3 through the cubic law, and at 1e-4 the determinant of a
            // perfectly well-conditioned block falls under the cutoff.
            const auto minWidth = prm.get<double>("solver.embedded_min_aperture", 1e-3);

            fractureAuxCells_ = &this->registerAuxCellModule_
                (std::make_unique<FractureAuxCells<TypeTag>>(this->simulator(),
                                                             static_cast<unsigned>(capacity),
                                                             static_cast<Scalar>(minWidth),
                                                             perfWiMode,
                                                             static_cast<Scalar>(perfWidth),
                                                             static_cast<Scalar>(perfRw)));

            if (this->simulator().gridView().comm().rank() == 0) {
                OpmLog::info(fmt::format("Embedded fracture flow: {} degrees of freedom "
                                         "reserved for fracture cells", capacity));
            }
        }

        /*!
         * \brief Hand the fracture's cells their degrees of freedom.
         *
         * Called once the fracture model has been built or has moved, so that what the
         * reservoir sees matches what the fracture is.
         */
        // grid cell whose saturation functions an auxiliary DOF borrows, -1 = partner
        int embeddedSatProxy_(unsigned globalDofIdx) const
        {
            if ((embeddedSatnum_ <= 0) || (globalDofIdx < this->model().numGridDof())) {
                return -1;
            }
            if (embeddedSatProxyCell_ < 0) {
                const unsigned want = static_cast<unsigned>(embeddedSatnum_ - 1);
                const unsigned n = this->model().numGridDof();
                for (unsigned c = 0; c < n; ++c) {
                    if (this->satnumRegionIndex(c) == want) { embeddedSatProxyCell_ = static_cast<int>(c); break; }
                }
                if (embeddedSatProxyCell_ < 0) {
                    OPM_THROW(std::runtime_error, "solver.embedded_satnum=" + std::to_string(embeddedSatnum_)
                              + ": no grid cell carries that SATNUM region");
                }
            }
            return embeddedSatProxyCell_;
        }

        /*!
         * \brief Hand the fracture's cells their degrees of freedom.
         *
         * \param allowTopologyChange whether the set of cells may change here.
         * \param requireStableLayout only restructure if the fracture asked for
         *        the same shape as at the previous call.  A fracture solve
         *        re-grids while it searches for its propagation front, so a bind
         *        inside the step that followed every one of those would chase a
         *        moving target and the coupling residual would never settle.
         */
        void bindFractureAuxCells(const bool allowTopologyChange = true,
                                  const bool requireStableLayout = false)
        {
            if ((fractureAuxCells_ == nullptr) || !this->geoMechModel().fractureModelActive()) {
                return;
            }

            // In the static configuration the first successful bind is also the last:
            // the flow keeps the fracture exactly as it first appeared.
            if (embeddedStatic_ && (fractureAuxCells_->numActive() > 0)) {
                embeddedCouplingChange_ = 0.0;
                return;
            }

            // A propagation attempt the fracture rolled back leaves its leak-off sized
            // for the grid that was tried; recompute it so real growth is not mistaken
            // for an inconsistent state and silently refused, step boundary after step
            // boundary.
            if (allowTopologyChange) {
                this->geoMechModel().fractureModel().ensureFlowDescriptionCurrent();
            }

            // Inside a time step the topology stays what it was: a fracture solve that
            // grew or reset its grid changes the shape of the flow system, and a Newton
            // iteration already under way cannot converge on a moving target.  Value
            // changes -- apertures, transmissibilities -- pass through; shape changes
            // wait for the step boundary, which is the sequentially implicit contract.
            // Three distinct things can be asked of the binding, at very
            // different cost:
            //
            //  - the shape is unchanged and only the apertures have moved, which
            //    is every iteration of the coupled mechanics-pressure solve:
            //    refresh the pore volumes, the cubic-law transmissibilities and
            //    the connection factors over the binding that exists.  No
            //    sparsity change and no matrix rebuild;
            //  - the shape may change and restructuring is allowed -- a step
            //    boundary, or right after a fracture solve that grew: rebind;
            //  - the shape changed and restructuring is not allowed: nothing
            //    per-cell is well defined, because a regrid renumbers the
            //    trimesh and a cell index stops meaning the same cell, so the
            //    old binding stands until someone may rebind.
            auto& fractureModel = this->geoMechModel().fractureModel();

            // What shape is the fracture asking for now, and is it the same one
            // it asked for last time?
            std::vector<std::size_t> layoutNow;
            for (const auto& wellFractures : fractureModel.wellFractures()) {
                for (const auto& fracture : wellFractures) {
                    layoutNow.push_back(fracture.numCells());
                }
            }
            const bool layoutStable = (layoutNow == lastSeenFractureLayout_);
            lastSeenFractureLayout_ = layoutNow;

            const bool mayRestructure
                = allowTopologyChange && (layoutStable || !requireStableLayout);

            if (!mayRestructure) {
                if (fractureAuxCells_->updateValues(fractureModel)) {
                    embeddedCouplingChange_ = fractureAuxCells_->lastBindChange();
                    this->refreshAuxCellModules_(/*topologyChanged=*/false);
                    // both time levels: the pore volume moved with the aperture,
                    // and the start-of-step state is stored as a volume too
                    this->model().updateAuxiliaryIntQuants(/*timeIdx=*/0);
                    this->model().updateAuxiliaryIntQuants(/*timeIdx=*/1);
                    this->checkFractureCouplingIfRequested_();
                    return;
                }
                if (embeddedCellDump_ > 0) {
                    OpmLog::info("Embedded fracture flow: the fracture changed shape "
                                 "mid-step; keeping the previous binding until it may "
                                 "be rebuilt");
                }
                embeddedCouplingChange_ = 0.0;
                return;
            }

            // A topology change makes the flow problem copy the current state of
            // every auxiliary degree of freedom into the previous-time state.
            // That is right for a cell that has just appeared -- it has no
            // history, so it has moved no mass by coming into existence -- and
            // harmless at a step boundary, where the two are equal anyway.  In
            // the middle of a step it would also erase the start-of-step state
            // of every cell that was already there, which is the reference its
            // accumulation term is measured against.  Keep theirs.
            auto& previous = this->model().solution(/*timeIdx=*/1);
            const auto firstAux = this->model().numGridDof();
            const auto numTotalDof = this->model().numTotalDof();
            std::vector<typename std::decay_t<decltype(previous)>::block_type> previousAux;
            previousAux.reserve(numTotalDof - firstAux);
            for (unsigned dof = firstAux; dof < numTotalDof; ++dof) {
                previousAux.push_back(previous[dof]);
            }

            const bool topologyChanged = fractureAuxCells_->bind(fractureModel);

            this->refreshAuxCellModules_(topologyChanged);

            if (topologyChanged) {
                const auto& newborn = fractureAuxCells_->newbornDofs();
                for (unsigned dof = firstAux; dof < numTotalDof; ++dof) {
                    const bool isNewborn
                        = std::find(newborn.begin(), newborn.end(), dof) != newborn.end();
                    if (!isNewborn) {
                        previous[dof] = previousAux[dof - firstAux];
                    }
                }
            }
            // newborn cells were assigned at both time levels; refresh their cached
            // intensive quantities so the first linearization sees that state
            this->model().updateAuxiliaryIntQuants(/*timeIdx=*/0);
            this->model().updateAuxiliaryIntQuants(/*timeIdx=*/1);
            embeddedCouplingChange_ = fractureAuxCells_->lastBindChange();
            if (topologyChanged && requireStableLayout) {
                // A restructure inside the step handed the fracture different
                // degrees of freedom; the well's perforations of them are stale
                // and must be re-registered, or the well silently loses its
                // fracture (zero connection factor, no rate through it).
                this->addFracturePerforationsToWells();
            }
            if (mayRestructure) {
                fractureAuxCells_->cellDump(this->geoMechModel().fractureModel(), "after-bind", embeddedCellDump_);
            }
            this->checkFractureCouplingIfRequested_();
        }

        /*!
         * \brief Verify the fracture <-> mechanics coupling blocks against finite
         *        differences (opt-in, fractureparam.solver.check_coupling_fd).
         *
         * The blocks themselves are built with AD; this differentiates the same
         * residual kernels numerically and compares. check_coupling_fd_columns
         * keeps the cost bounded by checking only that many columns per
         * fracture, spread over the matrix; -1 checks every one.
         */
        void checkFractureCouplingIfRequested_()
        {
            if ((fractureAuxCells_ == nullptr) || !this->geoMechModel().fractureModelActive()) {
                return;
            }
            const PropertyTree prm = this->getFractureParam();
            CouplingCheckOptions opt;
            opt.enabled = prm.get<bool>("solver.check_coupling_fd", false);
            if (!opt.enabled) {
                return;
            }
            opt.max_columns = prm.get<int>("solver.check_coupling_fd_columns", 4);
            opt.tolerance = prm.get<double>("solver.check_coupling_fd_tolerance", 1e-5);
            opt.perturbation = prm.get<double>("solver.check_coupling_fd_perturbation", 1e-6);
            opt.verbosity = prm.get<int>("solver.check_coupling_fd_verbosity", 0);
            const bool ok = fractureAuxCells_->checkCoupling(
                this->geoMechModel().fractureModel(), opt, this->simulator().timeStepSize());
            if (!ok) {
                OpmLog::warning("Fracture coupling matrices disagree with finite differences; "
                                "see the per-fracture reports above");
            }
        }

        /*!
         * \brief The coupling residual of the embedded representation.
         *
         * The relative change, over the last rebind, of what the flow is actually fed:
         * total fracture-to-reservoir conductance and total fracture pore volume.  The
         * outer loop watches this in embedded mode instead of the well-index change
         * list, which is computed but never applied there.
         */
        double embeddedCouplingChange() const
        { return embeddedCouplingChange_; }

        const FractureAuxCells<TypeTag>* fractureAuxCells() const
        { return fractureAuxCells_; }

        /*!
         * \brief Hand the flow's fracture pressures to the fractures (opt-in
         *        fractureparam.solver.pressure_from_flow).
         *
         * Each bound fracture is switched to external-pressure mode and given the
         * water pressure of its aux cells plus the well's BHP for its well DOF, so
         * its next solve does mechanics, contact and propagation at the pressure
         * the flow (and the well) actually hold. Unbound fractures (seed phase, or
         * a grid the binding has not caught up with) keep solving their own pressure.
         */
        void pushAuxPressuresToFractures()
        {
            if ((fractureAuxCells_ == nullptr) || !this->geoMechModel().fractureModelActive()) {
                return;
            }
            const PropertyTree prm = this->getFractureParam();
            if (!prm.get<bool>("solver.pressure_from_flow", false)) {
                return;
            }
            const double minWidthFactor =
                prm.get<double>("solver.pressure_from_flow_min_width_factor", 2.0);
            auto& fractures = this->geoMechModel().fractureModel();
            const auto& wellState = this->wellModel().wellState();
            std::size_t fidx = 0;
            for (auto& wellFractures : fractures.wellFractures()) {
                for (auto& fracture : wellFractures) {
                    const auto p = fractureAuxCells_->cellPressures(fidx);
                    ++fidx;
                    if (p.size() != fracture.numCells()) {
                        continue;
                    }
                    double bhp = -1.0;
                    if (const auto wi = wellState.index(fracture.wellInfo().name); wi.has_value()) {
                        bhp = wellState.well(*wi).bhp;
                    }
                    if (bhp <= 0.0 && !p.empty()) {
                        bhp = p.front();
                    }
                    // Only hand over the pressure of a fracture that is already
                    // open: while it is establishing, its aperture is at the
                    // cubic-law floor and the flow's pressure is the pressure of a
                    // closed fracture, so pinning it there removes the
                    // width-pressure feedback that opens and propagates it.
                    const auto& w = fracture.fractureWidth();
                    double wmax = 0.0;
                    for (std::size_t i = 0; i < w.size(); ++i) {
                        wmax = std::max(wmax, w[i][0]);
                    }
                    const bool established = wmax > minWidthFactor * fracture.cubicLawMinWidth();
                    if (established && fracture.setExternalPressure(p, bhp)) {
                        fracture.setExternalPressureMode(true);
                    } else {
                        fracture.setExternalPressureMode(false);
                    }
                }
            }
        }

        //! Whether the fracture flows through degrees of freedom of its own.
        bool fractureFlowIsEmbedded() const
        { return fractureAuxCells_ != nullptr; }

        /*!
         * \brief Perforate the fracture's own cells from the well.
         *
         * The upscaled representation gives the well an extra index on each reservoir
         * cell the fracture reaches, which is how the fluid got from the well into the
         * formation without the fracture being part of the flow problem.  Here it is, so
         * the well connects to the fracture and the fracture connects to the formation --
         * each conductance appearing once, and none of them an upscaled q/dp.
         */
        void addFracturePerforationsToWells()
        {
            if (fractureAuxCells_ == nullptr) {
                return;
            }

            // Registered as real perforations, sized into the wells' state and
            // equations by the well model's own dynamic-structure rebuild -- the same
            // path an ACTIONX-driven COMPDAT takes.  Registering identical lists is a
            // no-op, so calling this every step costs nothing when nothing moved.
            using PerfData = PerforationData<Scalar>;

            for (const auto& wname : this->wellModel().schedule().wellNames(this->episodeIndex())) {
                const auto runtime = fractureAuxCells_->wellPerforations(wname);

                std::vector<PerfData> perfs;
                perfs.reserve(runtime.size());
                for (const auto& rp : runtime) {
                    auto& pd = perfs.emplace_back();
                    pd.cell_index = rp.cell;
                    pd.connection_transmissibility_factor = static_cast<Scalar>(rp.ctf);
                }

                this->wellModel().setAuxiliaryPerforations(wname, std::move(perfs));
            }
        }

        static void registerParameters(){
            MechParent::registerParameters();
            VtkGeoMechModule<TypeTag>::registerParameters();
            FlowLinearSolverParametersGeoMech::registerParameters<TypeTag>();
	    Parameters::Register<Parameters::MechPorosityCoupling>
	        ("Feed the geomechanical pore-volume change back into the flow "
	         "equations");
	    Opm::Parameters::SetDefault<Opm::Parameters::MechPorosityCoupling>(false);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableOpmRstFile>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableVtkOutput>(true);
	    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncVtkOutput>(false);
	    Opm::Parameters::SetDefault<Opm::Parameters::EnableAsyncEclOutput>(false);
        }

        void finishInit(){
            OPM_TIMEBLOCK(finishInit);
            MechParent::finishInit();
            const auto& simulator = this->simulator();
            const auto& eclState = simulator.vanguard().eclState();
            if(eclState.runspec().mech()){
                const auto& initconfig = eclState.getInitConfig();
                geoMechModel_.init(initconfig.restartRequested());
                for(size_t i=0; i < this->ymodule_.size(); ++i){
                    using IsoMat = Opm::Elasticity::Isotropic;
                    elasticparams_.push_back(std::make_shared<IsoMat>(i,this->ymodule_[i],this->pratio_[i]));
                }
                // read mechanical boundary conditions
                const auto& vanguard = simulator.vanguard();
                const auto& bcconfigs = vanguard.eclState().getSimulationConfig().bcconfig();
                const auto& bcprops = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
                const auto& gv = this->gridView();
                const auto& cartesianIndexMapper = vanguard.cartesianIndexMapper();
                Opm::Elasticity::nodesAtBoundary(bc_nodes_,
                                                 bcconfigs,
                                                 bcprops,
                                                 gv,
                                                 cartesianIndexMapper);

                bool is_ok = checkBcConfig(bc_nodes_);
                if(!is_ok){
                  // this need to be fixed for parallel runs
                  std::cout << "Error in boundary condition specification not proper for mechanical problem" << std::endl;
                }
            }
        }

        // ///
        // Backend hooks for the shared initial-stress handling
        // ///
        void prepareMechForInit() override{
            this->geoMechModel_.setMaterial(this->ymodule_, this->pratio_);
            this->geoMechModel_.updatePotentialForces();//Neede only of output
        }

        void initializeStressFromMechSolve() override{
            geoMechModel_.solveGeomechanics(/*use_body_force*/ true, /*relative_solve*/ false);
            for (size_t i = 0; i < this->initstress_.size(); ++i) {
                this->initstress_[i] = geoMechModel_.stress(i);
            }
            // stress in output on first step maybe wrong i.e. 2*stress;
            this->geoMechModel_.setFirstSolveTrue();// to do full rebuild next time step
        }

        void applyInitialOutputStress() override{
            this->geoMechModel_.setOutputPutStress(this->initstress_);
        }

        void timeIntegration()
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start TimeIntegration-------------------\n"
                << std::flush;
            }
            Parent::timeIntegration();
        }
        // The sequential mech/fracture coupling re-solves flow inside a step with
        // the Newton iteration context reset (SetupIterationContextGuard), so the
        // parent's first-iteration storage recycling would rebase the start-of-step
        // storage onto a mid-step iterate and silently destroy mass.  Same reason
        // the parent already disables it for TPSA.
        bool recycleFirstIterationStorage() const
        {
            if (this->simulator().vanguard().eclState().runspec().mech()) {
                return false;
            }
            return Parent::recycleFirstIterationStorage();
        }

        void beginTimeStep() override{
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start beginTimeStep-------------------\n"
                << std::flush;
            }
            Parent::beginTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures()){
                  if(!(int(this->cstress_.size()) == int(this->gridView().size(0)))){
                        OPM_THROW(std::runtime_error,"CSTRESS not set but fractures exists");
                    }
                }
                geoMechModel_.beginTimeStep();
                if(this->hasFractures()){
                    if (this->fractureFlowIsEmbedded()) {
                        // The fracture is part of the flow problem in its own right, so
                        // the reservoir has to be told what it now looks like.
                        this->bindFractureAuxCells();

                        // The well perforates the fracture, not the reservoir cells the
                        // fracture leaks into.  Registered before the well model starts
                        // the step, so its dynamic-structure rebuild sizes the wells,
                        // their state and their equations around the new perforations.
                        this->addFracturePerforationsToWells();
                        this->wellModel().beginTimeStep();

                        // After the well model has rebuilt around them: what the wells
                        // now carry at the fracture's degrees of freedom.  Reported here
                        // as well as at the end of the step so that a step which never
                        // converges still says whether the fracture was connected.
                        if (embeddedLeakoffReport_) {
                            fractureAuxCells_->perforationReport();
                        }
                    }
                    else {
                        this->wellModel().beginTimeStep();// just to be sure well conteiner is reinitialized
                        this->addConnectionsToWell(); // modify wells WI wiht fracture well 
                    }
                }
                
            }
            OPM_END_PARALLEL_TRY_CATCH("Begin time step geomech failed:",this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
        }
        void endTimeStep() override{
            // The state the step just converged on is still bound; compare the two
            // representations' view of the leak-off before anything moves.
            if (embeddedLeakoffReport_ && (fractureAuxCells_ != nullptr)
                && this->geoMechModel().fractureModelActive())
            {
                fractureAuxCells_->leakoffReport(this->geoMechModel().fractureModel());
            }
            if ((fractureAuxCells_ != nullptr) && this->geoMechModel().fractureModelActive()) {
                fractureAuxCells_->cellDump(this->geoMechModel().fractureModel(), "step-end", embeddedCellDump_);
            }

            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;
            }
            //Parent::FlowProblemType::endTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            this->wellModel().traceWellSources();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                geoMechModel_.endTimeStep();
                if(this->hasFractures() && this->geoMechModel().fractureModelActive()){
                    // method for handling extra connections from fractures
                    // it is options for not including them in fractures i.e. addconnections
                    //if(addPerfsToSchedule_){
                        if (!this->fractureFlowIsEmbedded()) {
                            // In the embedded representation the schedule owns no fracture
                            // completions: the well perforates the fracture's degrees of
                            // freedom directly, refreshed each time the fracture is bound.
                            this->addConnectionsToSchedual();
                        }
                        this->gridView().comm().barrier();
                    //}else{
                    // not not working ... more work...
                    // will only work if structure is ok
                    //    assert(false);
                    //    this->addConnectionsToWell();
                    //}
                
                    this->geoMechModel_.fractureModel().
                        assignGeoMechWellState(this->wellModel_.wellState());
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("End time step geomech failed: ", this->simulator().vanguard().grid().comm());
            this->emptyFractureLogger();
            Parent::endTimeStep();
            if(this->simulator().vanguard().eclState().runspec().mech()){
                if(this->hasFractures() ){
                    // need to be here ?? to have updated values
                    this->geoMechModel_.updateFilterCakePropertiesOnFractures();
                }
            }
            //
            // if(first_fracture_solve_){
            //     first_fracture_solve_ = false;
            //     geoMechModel_.writeFractureSolutionFirst();
            // }
            if(Parameters::Get<Parameters::EnableWriteAllSolutions>()){
              //OPM_BEGIN_PARALLEL_TRY_CATCH();
                geoMechModel_.writeFractureSolution();
                //  OPM_END_PARALLEL_TRY_CATCH("Writing fracture solution failed: ", this->simulator().vanguard().grid().comm());
            }
            //
            //Parent::FlowProblemType::endTimeStep();
            //Parent::endStepApplyAction();
        }

        const GeoMechModel<TypeTag>& geoMechModel() const
        { return geoMechModel_; }

        //! The fracture-driving host used by the shared problem layer and
        //! the coupling loop (for the VEM backend this is the GeoMechModel).
        GeoMechModel<TypeTag>& fractureHost()
        { return geoMechModel_; }

        const GeoMechModel<TypeTag>& fractureHost() const
        { return geoMechModel_; }

        GeoMechModel<TypeTag>& geoMechModel()
        { return geoMechModel_; }

        //! Mechanical pore-volume change fed back into the flow equations.
        //! Off by default for the VEM backend (historic behaviour); when the
        //! common MechPorosityCoupling parameter is enabled the feedback is
        //! biot * tr(eps) from the last mechanics solve, mirroring the TPSA
        //! backend's native coupling.
        Scalar rockMechPoroChange(unsigned elementIdx, unsigned /*timeIdx*/) const
        {
            if (!Parameters::Get<Parameters::MechPorosityCoupling>()) {
                return 0.0;
            }
            const auto& eps = geoMechModel_.strain(elementIdx);
            return this->biotCoef(elementIdx) * (eps[0] + eps[1] + eps[2]);
        }

        const std::vector<std::tuple<size_t,MechBCValue>>& bcNodes() const{
            return bc_nodes_;
        }

        Dune::FieldVector<double,6> stress(size_t globalIdx) const{
            return geoMechModel_.stress(globalIdx);
        }
        // double getFieldProps(const std::string& field, unsigned globalIdx) const{
        //     const auto& eclState = this->simulator().vanguard().eclState();
        //     const auto& fp = eclState.fieldProps();
        //     const auto& myvec = fp.get_double(field);
        //     return myvec[globalIdx];
        // }
      //const std::vector< GridView::Codim<0>::EntitySeed >& elementEntitySeed(){return entitity_seed_;}
      //const std::vector< CellSeedType >& elementEntitySeed(){return entity_seed_;}
        void advanceTimeLevel(){
            Parent::advanceTimeLevel();
            //this->simulator_.problem().geoMechModel().advanceTimeLevel();// this is done in begin timestep
            //wellModel_.serialize(res);
            //aquiferModel_.serialize(res);
        }
 
    private:
        // ----------------------------------------------------------------------------
        // Heuristic rotational constraint check based only on fixed-direction
        // masks in bc_nodes.
        //
        // Since coordinates are not available here, this is a structural check:
        // to constrain rotation around an axis, we require fixed DOFs in the two
        // transverse directions on at least two distinct nodes.
        // ----------------------------------------------------------------------------
        static bool checkRotationsConstrained(
            const std::vector<std::tuple<size_t, MechBCValue>>& bc_nodes)
        {
            auto constrainedAroundAxis = [&](int d1,
                                             int d2,
                                             const std::string& axisName) -> bool
            {
                std::set<size_t> nodes_d1;
                std::set<size_t> nodes_d2;

                for (const auto& [node_idx, bc] : bc_nodes) {
                    if (bc.fixeddir[d1]) {
                        nodes_d1.insert(node_idx);
                    }
                    if (bc.fixeddir[d2]) {
                        nodes_d2.insert(node_idx);
                    }
                }

                if (nodes_d1.empty() || nodes_d2.empty()) {
                    OpmLog::warning(
                        "Mechanical BC: rotation around " + axisName
                        + " may be unconstrained (missing fixed DOFs in one or "
                          "both transverse directions).");
                    return false;
                }

                std::set<size_t> union_nodes = nodes_d1;
                union_nodes.insert(nodes_d2.begin(), nodes_d2.end());

                if (union_nodes.size() < 2) {
                    OpmLog::warning(
                        "Mechanical BC: rotation around " + axisName
                        + " may be unconstrained (transverse constraints are "
                          "applied at a single node only).");
                    return false;
                }

                return true;
            };

            const bool rx = constrainedAroundAxis(1, 2, "X");
            const bool ry = constrainedAroundAxis(0, 2, "Y");
            const bool rz = constrainedAroundAxis(0, 1, "Z");

            return rx && ry && rz;
        }

        // ----------------------------------------------------------------------------
        // Check that the mechanical boundary conditions are self-consistent and
        // cover enough DOFs to prevent rigid-body translations.
        //
        // Returns true if the configuration is valid, false otherwise.
        // An empty bc_nodes list is accepted (e.g. stress-only BCs on all faces).
        // A warning is emitted for any spatial direction that has no fixed node,
        // since the resulting system may be singular.
        // ----------------------------------------------------------------------------
        static bool checkBcConfig(
            const std::vector<std::tuple<size_t, MechBCValue>>& bc_nodes)
        {
            // Count how many nodes fix each direction, and check that every
            // listed node actually constrains something.
            std::array<int, 3> fixed_count = {0, 0, 0};
            const std::array<std::string, 3> dir_name = {"X", "Y", "Z"};

            for (const auto& [node_idx, bc] : bc_nodes) {
                bool node_fixes_anything = false;
                for (int d = 0; d < 3; ++d) {
                    if (bc.fixeddir[d]) {
                        ++fixed_count[d];
                        node_fixes_anything = true;
                    }
                }
                if (!node_fixes_anything) {
                    OpmLog::warning("BC node " + std::to_string(node_idx)
                        + " is listed in bc_nodes but does not fix any"
                          " displacement direction – check BCMECH input.");
                    //return false;
                }
            }

            // Warn (but do not fail) when a direction has no fixed nodes.
            // The system may still be well-posed via stress BCs or symmetry.
            if(bc_nodes.empty()){
                OpmLog::warning(
                    "Mechanical BC configuration: no nodes are fixed.  The system may be singular");
                return false;
            }
            if (!bc_nodes.empty()) {
                for (int d = 0; d < 3; ++d) {
                    if (fixed_count[d] == 0) {
                        OpmLog::warning(
                            "Mechanical BC configuration: no nodes are fixed in "
                            + dir_name[d]
                            + "-direction.  The system may is singular");
                        return false;      
                    };
                    //  else {
                    //     OpmLog::info(
                    //         "Mechanical BC: " + std::to_string(fixed_count[d])
                    //         + " node(s) fixed in " + dir_name[d] + "-direction.");
                    // }
                }
            }

            if (!checkRotationsConstrained(bc_nodes)) {
                OpmLog::warning(
                    "Mechanical BC configuration: rotational modes are not "
                    "sufficiently constrained based on bc_nodes.");
                return false;
            }
            return true;
        }

        GeoMechModel<TypeTag> geoMechModel_;

        std::vector<std::tuple<size_t,MechBCValue>> bc_nodes_;
        //std::vector<Opm::Elasticity::Material> elasticparams_;
        std::vector<std::shared_ptr<Opm::Elasticity::Material>> elasticparams_;
      //std::vector< CellSeedType > entity_seed_;
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };
}
#endif
