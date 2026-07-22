#ifndef OPM_FRACTURE_MECH_HOST_HPP
#define OPM_FRACTURE_MECH_HOST_HPP

// Fracture-driving host, extracted verbatim from the VEM GeoMechModel so
// that any mechanics backend (VEM, TPSA, ...) can own and drive the same
// fracture machinery.  Everything the fracture solve needs it reaches via
// simulator.problem(); the host itself only owns the FractureModel.

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/geomech/FractureModel.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm{
    template<typename TypeTag>
    class FractureMechHost
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    public:
        explicit FractureMechHost(Simulator& simulator)
            : simulator_(simulator)
        {}

        //! Called at the start of each time step.
        void beginTimeStep(){
            if(fracturemodel_){
                fracturemodel_->moveForwardInTime(simulator_.timeStepSize());
            }
        }

        //! Whether the fracture model adds contributions to disp/strain/
        //! stress output (include_fracture_contributions parameter; set when
        //! the fracture model is created).
        bool includeFractureContributions() const{
            return include_fracture_contributions_;
        }

        //! Raw access for output-side null checks; may be null before the
        //! first fracture solve.
        const FractureModel* fractureModelPtr() const{
            return fracturemodel_.get();
        }

       void solveFractures(){
            OPM_TIMEBLOCK(solveFractures);
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            int reportStepIdx = simulator_.episodeIndex();
            const auto& schedule =  this->simulator_.vanguard().schedule();
            int end_step = schedule.size() - 1;
            bool no_seeds = schedule[end_step].wseed().empty();
            const auto& prm = this->simulator_.problem().getFractureParam();
            int verbosity = prm.get("fractureparam.verbosity",1);
            if (simulator_.gridView().comm().rank() == 0) {
                std::stringstream os;
                if (!no_seeds) {
                    // std::cout << "No fracture seeds found, on this step " << reportStepIdx <<
                    // std::endl;
                    if(verbosity > 1){
                        os << "Fracture seeds found, on this step " << std::endl;
                    }
                } else {
                    if(verbosity > 1){
                        os << "No fracture seeds found, on this step " << std::endl;
                    }
                }
                if (fracturemodel_ && (verbosity > 1)) {
                    os << "Fracture model already initialized, solving fractures using previous "
                          "fractures"
                       << std::endl;
                }
                OpmLog::info(os.str());
            }
            if(!no_seeds && !fracturemodel_){
                    if(simulator_.gridView().comm().rank() == 0 && verbosity > 1){
                        std::stringstream os;    
                        os << "Fracture model not initialized, initializing now. report step" <<  reportStepIdx << std::endl;
                        OpmLog::info(os.str());
                    }
                    const auto& problem = simulator_.problem();
                    //NB could probably be moved to some initialization
                    // let fracture contain all wells
                    Opm::PropertyTree param = problem.getFractureParam();
                    include_fracture_contributions_ = param.get<bool>("include_fracture_contributions");
                    //param.read("fractureparam.json");
                    //const auto& schedule =  this->simulator_.vanguard().schedule();
                    //int reportStepIdx = simulator_.episodeIndex();
                    // take all wells and perforations
                    //int end_step = schedule.size() - 1;
                    const std::vector<Opm::Well>& wells = problem.wellModel().getLocalWells(end_step);
                    //const std::vector<Opm::Well>& wells = schedule.getWells(reportStepIdx);
                    //const Opm::EclipseGrid& eclgrid = simulator_.vanguard().eclState().getInputGrid();
                    const auto& grid = simulator_.vanguard().grid();
                    std::string outputDir = Parameters::Get<Parameters::OutputDir>();
                    std::string caseName  = simulator_.vanguard().caseName();
                    param.put("outputdir", outputDir);
                    param.put("casename", caseName);
                    //
                    

                    // LGR-aware connection-to-cell resolution: COMPDATL
                    // connections carry LGR-local global indices which must
                    // be mapped through the vanguard's LGR machinery, not the
                    // level-zero Cartesian mapper.
                    const auto& eclGrid = simulator_.vanguard().eclState().getInputGrid();
                    FractureModel::ConnCellResolver connCellResolver =
                        [this, &eclGrid](const Opm::Well&, const Opm::Connection& conn) -> int {
                            if (conn.get_lgr_level() > 0) {
                                const auto& tag = eclGrid.get_lgr_labels_by_number(conn.get_lgr_level());
                                return simulator_.vanguard().compressedIndexForInteriorLGR(tag, conn);
                            }
                            return simulator_.vanguard().compressedIndexForInterior(conn.global_index());
                        };
                    try{
                        fracturemodel_ = std::make_unique<FractureModel>(grid,
                                                                     wells,
                                                                     param,
                                                                     connCellResolver
                        );
                    }catch (...){
                       fracturemodel_ = nullptr;
                       OPM_THROW(std::runtime_error,"Error in initialising fracture model");
                    }
                    // not to get the reservoir properties along the well before initialising the well
                    // most important stress
                    fracturemodel_->updateReservoirWellProperties<TypeTag,Simulator>(simulator_);
                    // add fractures along the wells
                    //fracturemodel_->addFractures(schedule[reportStepIdx]);
                    fracturemodel_->addFractures(
                        schedule[end_step],
                        &simulator_.vanguard().eclState().getInputGrid());

                    fracturemodel_->updateFractureReservoirCells(grid);
                    fracturemodel_->initReservoirProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->updateReservoirAndWellProperties<TypeTag,Simulator>(simulator_);
                    fracturemodel_->updateActive(schedule[reportStepIdx].wseed);
                    fracturemodel_->initFractureStates();
                    this->writeFractureSolutionFirst();
            }
            // get reservoir properties on fractures
            // simulator need
            
            if(fracturemodel_ && !schedule[reportStepIdx].wseed().empty()){
                const auto& current_wseed = schedule[reportStepIdx].wseed;    
                if(simulator_.gridView().comm().rank() == 0){
                    if(verbosity > 1){
                        std::ostringstream os;
                        os << "Frac modelfound, updating reservoir properties and solving fractures";// << std::endl;
                        OpmLog::info(os.str());
                    }
                }
                fracturemodel_->updateReservoirAndWellProperties<TypeTag,Simulator>(simulator_);// set all fractures active if well is active
                fracturemodel_->updateActive(current_wseed);//only set fracture active if seed is active
                fracturemodel_->solve<TypeTag, Simulator>(simulator_);
                if(simulator_.gridView().comm().rank() == 0){
                    const auto& last_stats = fracturemodel_->lastSolveStats();
                    const auto& total_stats = fracturemodel_->totalSolveStats();
                    std::ostringstream os;
                    os << "Fracture solve stats: fractures_solved=" << last_stats.fractures_solved
                       << " (total " << total_stats.fractures_solved << ")"
                       << ", nonlinear_iterations=" << last_stats.nonlinear_iterations
                       << " (total " << total_stats.nonlinear_iterations << ")"
                       << ", linear_solves=" << last_stats.linear_solves
                       << " (total " << total_stats.linear_solves << ")"
                       << ", linear_iterations=" << last_stats.linear_iterations
                       << " (total " << total_stats.linear_iterations << ")"
                       << ", closed_cell_toggles=" << last_stats.closed_cell_toggles
                       << " (total " << total_stats.closed_cell_toggles << ")"
                       << ", linear_solve_failures=" << last_stats.linear_solve_failures
                       << " (total " << total_stats.linear_solve_failures << ")"
                       << ", ladder_rescues=" << last_stats.ladder_rescues
                       << " (total " << total_stats.ladder_rescues << ")"
                       << ", solve_time_s=" << last_stats.solve_time_seconds
                       << " (total " << total_stats.solve_time_seconds << ")"
                       << ", converged=" << (last_stats.converged ? "true" : "false");
                    OpmLog::info(os.str());
                }
            }else{
                if(simulator_.gridView().comm().rank() == 0){
                    std::ostringstream os;
                    os << "Fracture model not initialized, not solving fractures";// << std::endl;
                    OpmLog::info(os.str());
                }
            }
                // copy from apply action
            OPM_END_PARALLEL_TRY_CATCH("Solving fracture failed: ", simulator_.vanguard().grid().comm());  
       }
       

        void writeFractureSolutionFirst(){
            const auto& problem = simulator_.problem();
            if(problem.hasFractures() && fracturemodel_){
                // write first solution in standard format
                // this may ad some extra output of static variables
                //int reportStepIdx = simulator_.episodeIndex();
                    //fracturemodel_->write(reportStepIdx);
                    // hack to get correct number of fracture output
                    fracturemodel_->writemulti(0.0);
            }
        }
        
        void writeFractureSolution(){
            const auto& problem = simulator_.problem();
            if(problem.hasFractures() && fracturemodel_){
                // write first solution in standard format
                // this may ad some extra output of static variables
                //int reportStepIdx = simulator_.episodeIndex();
                double time = simulator_.time();
                fracturemodel_->writemulti(time);
            }

        }


        std::vector<RuntimePerforation> getExtraWellIndices(const std::string& wellname){
            if(fracturemodel_){
                return fracturemodel_->getExtraWellIndices(wellname);
            }else{
                return std::vector<RuntimePerforation>();
            }
        }
      
        bool fractureModelActive() const{
            if(!fracturemodel_){                
                return false;
            }else{
                return true;
            }           
        }

      void updateFilterCakePropertiesOnFractures(){
        if(this->fractureModelActive()){
          fracturemodel_->updateFilterCakeProperties<TypeTag,Simulator>(simulator_);
        }
      }

      FractureModel& fractureModel() {
            assert(fracturemodel_);
            return *fracturemodel_;
      }
      const FractureModel& fractureModel() const{
            if(!fracturemodel_){
                std::cout << "Fracture model not initialized, returning nullptr" << std::endl;
                throw std::runtime_error("Fracture model not initialized");
            }
            return *fracturemodel_;
        }
        void resetFractureModel(){
            if(fracturemodel_){
                fracturemodel_->resetFractures<TypeTag, Simulator>(simulator_);
            }
        }

    private:
        Simulator& simulator_;
        bool include_fracture_contributions_{false};
        std::unique_ptr<FractureModel> fracturemodel_;
    };
}
#endif
