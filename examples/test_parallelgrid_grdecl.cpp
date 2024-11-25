//==============================================================================
//!
//! \file upscale_elasticity.cpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscaling on cornerpoint grids
//!
//==============================================================================
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
# define DISABLE_ALUGRID_SFC_ORDERING 1
#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
//#include <dune/grid/io/file/vtk/subsampleingvtkwriter.hh>

#include <dune/istl/matrixmarket.hh>
#include <dune/common/fmatrix.hh>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#if HAVE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/common/fromtogridfactory.hh>
#endif
#include <opm/grid/utility/StopWatch.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>
#include <opm/elasticity/matrixops.hpp>
#include <dune/grid/utility/globalindexset.hh>

#include <cstring>
#include <iostream>
#include "vectorfunctions.hh"
#include <unistd.h>
#include <opm/geomech/boundaryutils.hh>
#include <unistd.h>

#include <opm/geomech/dune_utilities.hpp>
#include <opm/geomech/DuneCommunicationHelpers.hpp>

using namespace Opm::Elasticity;


//! \brief Display the available command line parameters
void syntax(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " gridfilename=filename.grdecl [method=]" << std::endl
            << "\t[xmax=] [ymax=] [zmax=] [xmin=] [ymin=] [zmin=] [linsolver_type=]" << std::endl
            <<" \t[Emin=] [ctol=] [ltol=] [rock_list=] [vtufilename=] [output=] [verbose=]" << std::endl << std::endl
            << "\t gridfilename             - the grid file. can be 'uniform'" << std::endl
            << "\t vtufilename              - save results to vtu file" << std::endl;
}



//! \brief Structure holding parameters configurable from command line
struct Params {
  std::string file;
  //! \brief Rocklist overriding material params in .grdecl
  std::string vtufile;
  //! \brief Text output file
  std::string output;
};

//! \brief Parse the command line arguments
void parseCommandLine(int argc, char** argv, Params& p)
{
  Opm::ParameterGroup param(argc, argv);
  p.file     = param.get<std::string>("gridfilename");
  p.vtufile  = param.getDefault<std::string>("vtufilename","test_grid");
}

//! \brief Write a log of the simulation to a text file
void writeOutput(const Params& p, Opm::time::StopWatch& watch, int cells)
{
  // get current time
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  // get hostname
  char hostname[1024];
  gethostname(hostname,1024);

  // write log
  std::ofstream f;
  f.open(p.output.c_str());
  f << "######################################################################" << std::endl
    << "# Results from upscaling elastic moduli." << std::endl
    << "#" << std::endl
    << "# Finished: " << asctime(timeinfo)
    << "# Hostname: " << hostname << std::endl
    << "#" << std::endl
    << "# Upscaling time: " << watch.secsSinceStart() << " secs" << std::endl
    << "#" << std::endl;

    f  << "# Eclipse file: " << p.file << std::endl
       << "#\t cells: " << cells << std::endl;

  f << "#" << std::endl;
}
using PolyGrid = Dune::PolyhedralGrid<3, 3>;
#if HAVE_ALUGRID
using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming >;
void createGrids(std::unique_ptr<AluGrid3D>& grid ,Opm::EclipseState& eclState,std::vector<unsigned int>& ordering){
    Dune::CpGrid cpgrid;
    const auto& input_grid = eclState.getInputGrid();
    const auto& global_porv = eclState.fieldProps().porv(true);
    cpgrid.processEclipseFormat(&input_grid,
                                &eclState,
                                /*isPeriodic=*/false,
                                /*flipNormals=*/false,
                                /*clipZ=*/false);

    //auto  cartesianCellId = cpgrid.globalCell();
    //std::vector<std::size_t>  nums = cpgrid.processEclipseFormat(eclState.getInputGrid(), nullptr, false);
    Dune::FromToGridFactory<AluGrid3D> factory;
    //std::vector<unsigned int> ordering;
    auto cartesianCellIndx = cpgrid.globalCell();
    grid = factory.convert(cpgrid, cartesianCellIndx, ordering);
}
#endif
void createGrids(std::unique_ptr<PolyGrid>& grid ,const Opm::EclipseState& eclState,std::vector<unsigned int>& /*ordering*/)
{
    grid  = std::make_unique<PolyGrid>(eclState.getInputGrid(), eclState.fieldProps().porv(true));
}

void createGrids(std::unique_ptr<Dune::CpGrid>& grid ,const Opm::EclipseState& eclState,std::vector<unsigned int>& /*ordering*/){
    grid = std::make_unique<Dune::CpGrid>();
    std::vector<std::size_t>  nums = grid->processEclipseFormat(&eclState.getInputGrid(), nullptr, false);
}

//! \brief Main solution loop. Allows templating over the AMG type
template<class GridType>
int run(Params& p, const std::string& name, Dune::MPIHelper& mpihelper)
{    
    Dune::Communication<MPI_Comm> world_comm = Dune::MPIHelper::getCommunication();
    constexpr int dim = 3;
    if (0 == mpihelper.rank()){
        std::cout << "Using " << mpihelper.size() << " Processes." << std::endl;
    }
  try {
    static constexpr int dim = GridType::dimension;
    //static constexpr int dimensionworld = GridType::dimensionworld;
    //static const int dim = 3;

    Opm::time::StopWatch watch;
    watch.start();


    std::unique_ptr<GridType> grid_ptr;
    Opm::Parser parser;
    // process grid
    auto deck = parser.parseFile(p.file);
    Opm::EclipseState eclState(deck);
    Opm::EclipseGrid inputGrid(deck);
    // create grids depeing on grid type
    std::vector<unsigned int> ordering;
    createGrids(grid_ptr, eclState, ordering);
    GridType& grid = *grid_ptr;
    auto python = std::make_shared<Opm::Python>();
    //const Opm::Schedule schedule(deck, eclState);
    //
    //Dune::CpGrid cpGrid;
    //std::vector<std::size_t>  nums = cpGrid.processEclipseFormat(&eclState.getInputGrid(), nullptr, false);
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;//maybe wrong
    if(mpihelper.rank() == 0){
        std::cout << "Before Loadbalance" << std::endl;
    }
    world_comm.barrier();

    for(int rank=0; rank<mpihelper.size(); ++rank){
        if(rank == mpihelper.rank()){
            std::cout << "Grid size: " << grid.size(0) << " on rank " << mpihelper.rank() <<std::endl;
        }
        usleep(100);
        world_comm.barrier();
    }
    
    int overlapLayers=1;
    int partitionMethod = Dune::PartitionMethod::zoltan;
    //int partitionMethod = Dune::PartitionMethod::simple;
    double imbalanceTol = 1.1;
    bool addCornerCells = false;
    bool output_vertex = false;
    usleep(100);
    if(mpihelper.rank() == 0){
        std::cout << "Start Loadbalance" << " cells." << std::endl;
    }
    world_comm.barrier();
    grid.loadBalance(overlapLayers,partitionMethod,imbalanceTol, addCornerCells);
    grid.switchToDistributedView();
    for(int rank=0; rank<mpihelper.size(); ++rank){
        if(rank == mpihelper.rank()){
            std::cout << "Grid size: " << grid.size(0) << " on rank " << mpihelper.rank() <<std::endl;
        }
        usleep(100);
        world_comm.barrier();
    }
   
    auto gv = grid.leafGridView();
    auto indexset = gv.indexSet();
    auto gidSet = gv.grid().globalIdSet();
    int codim = 3;
    Dune::Codim<3> mycodim;
    using Vector = std::vector<int>;
    Vector myrank(gv.size(codim), world_comm.rank());
    using GridView =  typename GridType::LeafGridView;
    Dune::MaxEntityVectorVectorDataHandle<GridView,Vector> datahandle(myrank, gv, codim);
    //gv.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    gv.communicate(datahandle, Dune::All_All_Interface, Dune::ForwardCommunication);
    auto allranks = datahandle.all_data();
    std::vector<int> maxrank(gv.size(codim), world_comm.rank());
    for(size_t j=0; j<allranks.size(); ++j){
         const auto& all = allranks[j];  
        for(size_t i=0; i<all.size(); ++i){
            maxrank[i] = std::max(maxrank[i], all[i]);
        }
    }


    for (int rank = 0; rank < mpihelper.size(); ++rank) {
          if (rank == mpihelper.rank()) {
              std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;
              for (auto& elem : Dune::entities(gv, mycodim)) {
                  auto lindex = elem.index();
                  // auto gindex = gindexset.index(elem);
                  auto gid = gidSet.id(elem);
                  std::cout << "Entity<" << codim << "> global: id " << gid << " local " << lindex << " type ";
                  std::cout << elem.partitionType() ; //<< " owner rank " << myrank[lindex];
                  //std::cout << " num pros " << numpros[lindex];
                  std::cout << " all ranks ";
                  for (auto& r : allranks[lindex]) {
                       // if(r!=1000){
                       std::cout << r << " ";
                       //}
                  }
                  std::cout << " max rank " << maxrank[lindex];
                  std::cout << std::endl;
              }
          }
          // for bliding communicator see ParallelIstlInformation. (or Dune::OwnerOverlapCopyCommunication)
          usleep(100);
          world_comm.barrier();
    }
    
    bool do_all = false;
    if(do_all){
    std::cout<<"cell-cell communication"<<std::endl;
    Opm::entityEntityCommunication<0>(grid, mpihelper);
    std::cout<<"cell-cell communication"<<std::endl;
    Opm::entityEntityCommunication<3>(grid, mpihelper);
    }

    auto cell_parallindex = Opm::makeEntityEntityCommunication<0>(grid,false);
    auto vertex_parallindex = Opm::makeEntityEntityCommunication<3>(grid,false);
    auto dof_parallel_index = Opm::entityToDofIndexSet(vertex_parallindex, 3);
    //Opm::makeEntityEntityCommunication<3>(grid, mpihelper);
    if (!p.output.empty()){
        writeOutput(p, watch, grid.size(0));
    }
    return 0;
  }
  catch (Dune::Exception &e) {
      throw e;
  }
  catch (...) {
      throw;
  }
  return 1;
}


//! \brief Main driver
int
main(int argc, char** argv)
{
    try {
        if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-?") == 0) {
            syntax(argv);
            exit(1);
        }

        Dune::MPIHelper& mpi = Dune::MPIHelper::instance(argc, argv);
        const int size = mpi.size();
        //if (size != 1) {
         //   std::cerr << "This program does not support MPI parallelization" << std::endl;
         //   return 2;
        //}

        Params p;
        parseCommandLine(argc, argv, p);
        // using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming >;
        using GridType = Dune::CpGrid;
        // using GridType = AluGrid3D;
        // using GridType = PolyGrid;
        int ok;
        ok = run<GridType>(p, "cpgrid", mpi);
        // ok = run<AluGrid3D>(p);
        // ok = run<PolyGrid>(p);
        return ok;
    } catch (Dune::Exception& e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (const std::exception& e) {
        throw e;
    } catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
    return 1;
}
