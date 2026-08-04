#ifndef OPM_FLOW_PROBLEM_GEOMECH_HPP
#define OPM_FLOW_PROBLEM_GEOMECH_HPP

#include <opm/common/ErrorMacros.hpp>

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
        void bindFractureAuxCells(const bool allowTopologyChange = true)
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
            if (!allowTopologyChange
                && !fractureAuxCells_->layoutMatches(this->geoMechModel().fractureModel()))
            {
                // The fracture wants a different shape; that waits for the step
                // boundary, and what the flow sees is unchanged.
                embeddedCouplingChange_ = 0.0;
                return;
            }

            const bool topologyChanged =
                fractureAuxCells_->bind(this->geoMechModel().fractureModel());

            this->refreshAuxCellModules_(topologyChanged);
            embeddedCouplingChange_ = fractureAuxCells_->lastBindChange();
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

            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start endTimeStep-------------------\n"
                << std::flush;
            }
            //Parent::FlowProblemType::endTimeStep();
            OPM_BEGIN_PARALLEL_TRY_CATCH();
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
