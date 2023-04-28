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
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCMECHConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <dune/alugrid/grid.hh>
#include <opm/grid/utility/StopWatch.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/geomech/vem_elasticity_solver.hpp>
#include <opm/elasticity/matrixops.hpp>
#include <dune/alugrid/common/fromtogridfactory.hh>

#include <cstring>
#include <iostream>
#include "vectorfunctions.hh"
#include <unistd.h>
#include <opm/geomech/boundaryutils.hh>
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
  p.vtufile  = param.getDefault<std::string>("vtufilename","test_elast");
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
using PolyGrid = Dune::PolyhedralGrid<3, 3>;
using AluGrid3D = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming >;

void createGrids(std::unique_ptr<PolyGrid>& grid ,const Opm::EclipseState& eclState,std::vector<unsigned int>& /*ordering*/)
{
    grid  = std::make_unique<PolyGrid>(eclState.getInputGrid(), eclState.fieldProps().porv(true));
}

void createGrids(std::unique_ptr<Dune::CpGrid>& grid ,const Opm::EclipseState& eclState,std::vector<unsigned int>& /*ordering*/){
    grid = std::make_unique<Dune::CpGrid>();
    std::vector<std::size_t>  nums = grid->processEclipseFormat(&eclState.getInputGrid(), nullptr, false);
}

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


//! \brief Main solution loop. Allows templating over the AMG type
template<class GridType,class ElasticitySolverType>
int run(Params& p, bool with_pressure, bool with_gravity, std::string name)
{    
  try {
    static constexpr int dim = GridType::dimension;
    //static constexpr int dimensionworld = GridType::dimensionworld;  
    //static const int dim = 3;

    Opm::time::StopWatch watch;
    watch.start();
    
    
    std::unique_ptr<GridType> grid_ptr;
    using GridView = typename GridType::LeafGridView;//Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
    Opm::Parser parser;
    // process grid
    auto deck = parser.parseFile(p.file);
    Opm::EclipseState eclState(deck);    
    Opm::EclipseGrid inputGrid(deck);
    // create grids depeing on grid type
    std::vector<unsigned int> ordering;
    createGrids(grid_ptr, eclState, ordering);
    const GridType& grid = *grid_ptr;
    Dune::CpGrid cpGrid;
    std::vector<std::size_t>  nums = cpGrid.processEclipseFormat(&eclState.getInputGrid(), nullptr, false);
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;//maybe wrong
    CartesianIndexMapper cartesianIndexMapper(cpGrid);
    ElasticitySolverType esolver(grid);// p.ctol, p.Emin, p.verbose);
    std::vector<std::shared_ptr<Opm::Elasticity::Material>> materials;
    
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
    //esolver.setMaterial(materials);
    esolver.setMaterial(ymodule,pratio);
    std::vector<size_t> fixed_nodes;
    const auto& bcconfig = eclState.getSimulationConfig().bcconfig();
    const auto& gv = grid.leafGridView();
    Opm::Elasticity::fixNodesAtBoundary(fixed_nodes,
                                        bcconfig,
                                        gv,
                                        cartesianIndexMapper
        );
    
    std::cout << "Effected nodes" << std::endl;
    for(int i:fixed_nodes){
        std::cout << i << std::endl;
    }
    
    // std::cout << "logical dimension: " << grid.logicalCartesianSize()[0]
    //           << "x"                   << grid.logicalCartesianSize()[1]
    //           << "x"                   << grid.logicalCartesianSize()[2]
    //           << std::endl;

    if (p.inspect == "mesh")
      return 0;


    Opm::Elasticity::Vector pressforce;
    pressforce.resize(grid.size(0));
    pressforce = 0.0;
    if(with_pressure){
        Dune::loadMatrixMarket(pressforce,"pressforce.mtx");
    }
    if(with_gravity){
        esolver.setBodyForce(9.8);
    }else{
        esolver.setBodyForce(0.0);
    }
    bool do_matrix = true;//assemble matrix
    bool do_vector = true;//assemble matrix
    esolver.fixNodes(fixed_nodes);
    esolver.initForAssembly();
    esolver.assemble(pressforce, do_matrix, do_vector);
    Opm::PropertyTree prm("mechsolver.json");
    esolver.setupSolver(prm);


     esolver.solve();
     std::cout << "\tsolution norm: " << esolver.u.two_norm() << std::endl;
     Opm::Elasticity::Vector field;
     field.resize(grid.size(dim)*dim);
     esolver.expandSolution(field,esolver.u);
     Dune::storeMatrixMarket(esolver.A.getOperator(), "A.mtx");
     Dune::storeMatrixMarket(esolver.A.getLoadVector(), "b.mtx");
     Dune::storeMatrixMarket(esolver.u, "u.mtx");
     Dune::storeMatrixMarket(field, "field.mtx");
     Dune::FieldMatrix<double,6,6> C;
     Dune::BlockVector<Dune::FieldVector<double,3>> disp;
     disp.resize(pressforce.size());
     for(size_t i = 0; i < pressforce.size(); ++i){
         for(size_t k = 0; k < dim; ++k){
             disp[i][k] = field[i*dim+k];
         }
     }

    if (!p.vtufile.empty()) {
        Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView(), Dune::VTK::nonconforming);
        
        //for (int i=0;i<6;++i) {
        std::stringstream str;
        
        str << "sol ";
        vtkwriter.addVertexData(field, str.str().c_str(), dim);
        using Container = Dune::BlockVector<Dune::FieldVector<double,1>>;
        using Function =  Dune::P1VTKFunctionVector<typename GridType::LeafGridView, Container>;
        //Dune::VTK::FieldInfo dispinfo("disp",Dune::VTK::FieldInfo::Type::vector,3);
        Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
        std::string vtkname("disp");
        Dune::VTKFunction<typename GridType::LeafGridView>* fp = new Function(grid.leafGridView(),
                                                       field,
                                                       vtkname,
                                                       3, 0, prec);
        vtkwriter.addVertexData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(fp));
        //vtkwriter.addVertexData(disp, dispinfo);
        vtkwriter.addCellData(pressforce, "pressforce");
        //vtkwriter.addVertexData(disp, stress);
        //}
        std::string outputfile = name + "_" + p.vtufile;
        vtkwriter.write(outputfile);
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
    static const int dim = 3;
    //using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming >;
    //using GridType = Dune::CpGrid;
    //using GridType = AluGrid3D;
    //using GridType = PolyGrid;
    int ok;
    //ok = run<AluGrid3D>(p);
    //ok = run<PolyGrid>(p);
    using GridType = Dune::CpGrid;
    using ElasticitySolverTypeVem = Opm::Elasticity::VemElasticitySolver<GridType>;
    using ElasticitySolverTypeFem = Opm::Elasticity::ElasticitySolver<GridType>;
    bool with_pressure =true;
    bool with_gravity = false;
    std::string name = "vem";
    ok = run<Dune::CpGrid, ElasticitySolverTypeVem>(p, with_pressure, with_gravity, std::string("vem"));
    ok = run<Dune::CpGrid, ElasticitySolverTypeFem>(p, with_pressure, with_gravity, std::string("fem"));
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
