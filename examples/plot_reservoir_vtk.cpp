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
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/readDeck.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
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

#include <cstring>
#include <iostream>
#include "vectorfunctions.hh"
#include <unistd.h>
#include <opm/geomech/boundaryutils.hh>
using namespace Opm::Elasticity;

#include <opm/geomech/FractureModel.hpp>

//! \brief Display the available command line parameters
void syntax(char** argv)
{
    std::cerr << "Usage: " << argv[0] << " gridfilename=filename.grdecl [option=]" << std::endl
            << "options ctol Emin vtufilename resultfilename output verbose inpect with_gravity with_pressure" << std::endl
            << "\t gridfilename             - the grid file. can be 'uniform'" << std::endl
            << "\t vtufilename              - save results to vtu file" << std::endl;
}



//! \brief Structure holding parameters configurable from command line
struct Params {
    double ctol;
    double Emin;
  //! \brief The eclipse grid file
  std::string file;
  //! \brief Rocklist overriding material params in .grdecl
  std::string vtufile;
  //! \brief Text output file
  std::string output;
  //! \brief verbose output
  bool verbose;
  //! \brief Run a inspection only, currently 'mesh, results, load'
  std::string inspect;
  //! \brief Result template filename (input/output)
  std::string resultfilename;
  bool with_gravity;
  bool with_pressure;
};

//! \brief Parse the command line arguments
void parseCommandLine(int argc, char** argv, Params& p)
{
  Opm::ParameterGroup param(argc, argv);
  p.ctol     = param.getDefault<double>("ctol",1.e-6);
  p.Emin     = param.getDefault<double>("Emin",1);
  p.file     = param.get<std::string>("gridfilename");
  p.vtufile  = param.getDefault<std::string>("vtufilename","test_elast");
  p.resultfilename  = param.getDefault<std::string>("resultfilename","");
  p.output   = param.getDefault<std::string>("output","");
  p.verbose  = param.getDefault<bool>("verbose",false);
  p.inspect  = param.getDefault<std::string>("inspect","");
  p.with_gravity  = param.getDefault<bool>("with_gravity",false);
  p.with_pressure  = param.getDefault<bool>("with_pressure",false);
  size_t i;
  if ((i=p.vtufile.find(".vtu")) != std::string::npos){
    p.vtufile = p.vtufile.substr(0,i);
  }
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
int run(Params& p, const std::string& name)
{
    static constexpr int dim = GridType::dimension;
    //static constexpr int dimensionworld = GridType::dimensionworld;
    //static const int dim = 3;

    Opm::time::StopWatch watch;
    watch.start();
    std::unique_ptr<GridType> grid_ptr;
    Opm::Parser parser;
    // process grid
    std::unique_ptr<Opm::ParseContext> parsercontext = Opm::setupParseContext(false);
    auto deck = parser.parseFile(p.file, *parsercontext.get());
    Opm::EclipseState eclState(deck);
    Opm::EclipseGrid inputGrid(deck);
    // create grids depeing on grid type
    std::vector<unsigned int> ordering;
    createGrids(grid_ptr, eclState, ordering);
    const GridType& grid = *grid_ptr;
    auto python = std::make_shared<Opm::Python>();
    const Opm::Schedule schedule(deck, eclState);
    std::vector<Opm::Well> wells = schedule.getWells(0);

    //const auto& wells = schedule.getWells(/*reportstep*/0);
    //
    const Opm::EclipseGrid& eclgrid = eclState.getInputGrid();
    Dune::CpGrid cpGrid;
    std::vector<std::size_t>  nums = cpGrid.processEclipseFormat(&eclgrid, nullptr, false);
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;

    Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView(), Dune::VTK::nonconforming);
    std::string outputfile = name + "_" + p.vtufile;
    std::cout << "Writing output to " << outputfile << std::endl;
    vtkwriter.write(outputfile);

    // look at fracture model
    Opm::PropertyTree prm;
    Opm::PropertyTree prmtmp;
    prm.put("hasfracture",true);
    prm.put("fractureparam.radius",10);
    Opm::FractureModel fracturemodel(grid, wells, eclgrid, prm,/*default_fractures*/ true);
    fracturemodel.updateReservoirProperties();
    fracturemodel.solve();
    fracturemodel.write();

    if (!p.output.empty()){
        writeOutput(p, watch, grid.size(0));
    }
    return 1;
}


//! \brief Main driver
int main(int argc, char** argv)
try
{
  try {
    if (argc < 2 || strcmp(argv[1],"-h") == 0
                 || strcmp(argv[1],"--help") == 0
                 || strcmp(argv[1],"-?") == 0) {
      syntax(argv);
      exit(1);
    }

    Dune::MPIHelper& mpi=Dune::MPIHelper::instance(argc, argv);
    const int size = mpi.size();
    if (size != 1) {
        std::cerr << "This program does not support MPI parallelization" << std::endl;
        return 2;
    }

    Params p;
    parseCommandLine(argc,argv,p);
    //using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming >;
    //using GridType = Dune::CpGrid;
    //using GridType = AluGrid3D;
    //using GridType = PolyGrid;
    int ok;
    //ok = run<AluGrid3D>(p);
    //ok = run<PolyGrid>(p);
    using GridType = Dune::CpGrid;
    ok = run<Dune::CpGrid>(p, "cpgrid");
    return ok;
  } catch (Dune::Exception &e) {
      std::cerr << "Dune reported error: " << e << std::endl;
  } catch (const std::exception &e) {
      throw e;
  }
  catch (...) {
      std::cerr << "Unknown exception thrown!" << std::endl;
  }
  return 1;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
