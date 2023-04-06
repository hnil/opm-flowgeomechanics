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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/common/fmatrix.hh>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/utility/StopWatch.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/elasticity/matrixops.hpp>

#include <cstring>
#include <iostream>

#include <unistd.h>

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
};

//! \brief Parse the command line arguments
void parseCommandLine(int argc, char** argv, Params& p)
{
  Opm::ParameterGroup param(argc, argv);
  p.ctol     = param.getDefault<double>("ctol",1.e-6);
  p.Emin     = param.getDefault<double>("Emin",1);
  p.file     = param.get<std::string>("gridfilename");
  p.vtufile  = param.getDefault<std::string>("vtufilename","");
  p.resultfilename  = param.getDefault<std::string>("resultfilename","");
  p.output   = param.getDefault<std::string>("output","");
  p.verbose  = param.getDefault<bool>("verbose",false);
  p.inspect  = param.getDefault<std::string>("inspect","");
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

//! \brief Main solution loop. Allows templating over the AMG type
template<class GridType>
int run(Params& p)
{
  try {
    static const int dim = 3;

    Opm::time::StopWatch watch;
    watch.start();

    GridType grid;
    Opm::Parser parser;
    auto deck = parser.parseFile(p.file);
    Opm::EclipseGrid inputGrid(deck);
    grid.processEclipseFormat(&inputGrid, nullptr, false);

    ElasticitySolver<GridType> esolver(grid, p.ctol, p.Emin, p.verbose);
    std::vector<std::shared_ptr<Opm::Elasticity::Material>> materials;
    Opm::EclipseState eclState(deck);
    //const auto& initconfig = eclState.getInitConfig();
    const auto& fp = eclState.fieldProps();            
    std::vector<double> ymodule = fp.get_double("YMODULE");
    std::vector<double> pratio = fp.get_double("PRATIO");
    std::vector<double> biotcoef = fp.get_double("BIOTCOEF");
    for(size_t i=0; i < ymodule.size(); ++i){
        using IsoMat = Opm::Elasticity::Isotropic;
        if(pratio[i]>0.5 || pratio[i] < 0){
            OPM_THROW(std::runtime_error,"Pratio not valid");
        }
        materials.push_back(std::make_shared<IsoMat>(i,ymodule[i],pratio[i]));
    }    
    esolver.setMaterial(materials);
    
    std::cout << "logical dimension: " << grid.logicalCartesianSize()[0]
              << "x"                   << grid.logicalCartesianSize()[1]
              << "x"                   << grid.logicalCartesianSize()[2]
              << std::endl;

    if (p.inspect == "mesh")
      return 0;


    Opm::Elasticity::Vector pressforce;
    pressforce.resize(grid.size(0));
    pressforce = 1.0;
    Dune::loadMatrixMarket(pressforce,"pressforce.mtx");
    //   upscale.fixCorners(p.min, p.max);
    bool do_matrix = true;//assemble matrix
    bool do_vector = true;//assemble matrix
    esolver.A.initForAssembly();
    esolver.assemble(pressforce, do_matrix, do_vector);
    Opm::PropertyTree prm("mechsolver.json");
    esolver.setupSolver(prm);


     esolver.A.printOperator();
     esolver.A.printLoadVector();
     esolver.solve();
     Opm::Elasticity::Vector field;
     esolver.A.expandSolution(field,esolver.u);
     Dune::storeMatrixMarket(esolver.A.getOperator(), "A.mtx");
     Dune::storeMatrixMarket(esolver.A.getLoadVector(), "b.mtx");
     Dune::storeMatrixMarket(esolver.u, "u.mtx");
     Dune::storeMatrixMarket(field, "field.mtx");
     Dune::FieldMatrix<double,6,6> C;
    

    if (!p.vtufile.empty()) {
      Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView());

      for (int i=0;i<6;++i) {
        std::stringstream str;
        str << "sol " << i+1;
        vtkwriter.addVertexData(field[i], str.str().c_str(), dim);
      }
      vtkwriter.write(p.vtufile);
    }

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
    return run<Dune::CpGrid>(p);
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
