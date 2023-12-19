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
  p.with_pressure  = param.getDefault<bool>("with_pressure",true);
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
    auto python = std::make_shared<Opm::Python>();
    const Opm::Schedule schedule(deck, eclState);
    //
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
        if(biotcoef[i]>1.0 || biotcoef[i] < 0.0){
            OPM_THROW(std::runtime_error,"BIOTCOEF not valid");
        }
        materials.push_back(std::make_shared<IsoMat>(i,ymodule[i],pratio[i]));
    }
    //esolver.setMaterial(materials);
    esolver.setMaterial(ymodule,pratio);
    std::vector< std::tuple<size_t, Opm::MechBCValue> > bc_nodes;
    const auto& bcconfigs = eclState.getSimulationConfig().bcconfig();
    const auto& bcprops  = schedule[0].bcprop;
    const auto& gv = grid.leafGridView();
    Opm::Elasticity::nodesAtBoundary(bc_nodes,
                                     bcconfigs,
                                     bcprops,
                                     gv,
                                     cartesianIndexMapper
        );

    // std::cout << "Effected nodes" << std::endl;
    // for(int i:fixed_nodes){
    //     std::cout << i << std::endl;
    // }

    // std::cout << "logical dimension: " << grid.logicalCartesianSize()[0]
    //           << "x"                   << grid.logicalCartesianSize()[1]
    //           << "x"                   << grid.logicalCartesianSize()[2]
    //           << std::endl;

    if (p.inspect == "mesh")
      return 0;


    Opm::Elasticity::Vector pressforce;
    size_t num_cells = grid.size(0);
    pressforce.resize(num_cells);
    pressforce = 1.0;
    if(p.with_pressure){
        Dune::loadMatrixMarket(pressforce,"pressforce.mtx");
    }
    if(p.with_gravity){
        esolver.setBodyForce(9.8);
    }else{
        esolver.setBodyForce(0.0);
    }
    bool do_matrix = true;//assemble matrix
    bool do_vector = true;//assemble matrix
    esolver.fixNodes(bc_nodes);
    esolver.initForAssembly();
    esolver.assemble(pressforce, do_matrix, do_vector);
    esolver.updateRhsWithGrad(pressforce);
    Opm::PropertyTree prm("mechsolver.json");
    esolver.setupSolver(prm);


     esolver.solve();
     std::cout << "\tsolution norm: " << esolver.u.two_norm() << std::endl;
     Opm::Elasticity::Vector field;
     field.resize(grid.size(dim)*dim);
     esolver.expandSolution(field,esolver.u);
     Dune::storeMatrixMarket(esolver.A.getOperator(), name + "_" + std::string("A.mtx"));
     Dune::storeMatrixMarket(esolver.A.getLoadVector(), name + "_" + std::string("b.mtx"));
     Dune::storeMatrixMarket(esolver.u, name + "_" + std::string("u.mtx"));
     Dune::storeMatrixMarket(field, name + "_" + std::string("field.mtx"));

     Dune::BlockVector<Dune::FieldVector<double,3>> disp;
     Opm::Elasticity::Vector lindisp;
     disp.resize(pressforce.size());
     lindisp.resize(pressforce.size()*dim);
     disp = 0.0;
     for (const auto& cell: elements(gv)){
         auto cellindex = gv.indexSet().index(cell);
         const auto& vertices = Dune::subEntities(cell, Dune::Codim<GridType::dimension>{});
         const auto numvert = vertices.size();
         for (const auto& vertex : vertices){
             auto nodeidex = gv.indexSet().index(vertex);
             for(int k=0; k < dim; ++k){
                 disp[cellindex][k] += field[nodeidex*dim+k] ;
                 lindisp[cellindex*dim+k] += field[nodeidex*dim+k]/numvert;
             }
         }
         disp[cellindex] /= numvert;
     }
     std::vector<double> stresslin(num_cells*6,0.0);
     std::vector<double> strainlin(num_cells*6,0.0);

     esolver.calculateStress(true);
     esolver.calculateStrain(true);

     const Dune::BlockVector<Dune::FieldVector<double,6>>& stress = esolver.stress();
     const Dune::BlockVector<Dune::FieldVector<double,6>>& strain = esolver.strain();

     for(int i=0; i < grid.size(0); ++i){
         for(int k=0; k < 6; ++k){
             stresslin[i*6+k] = stress[i][k];
             strainlin[i*6+k] = strain[i][k];
         }
     }

    if (!p.vtufile.empty()) {
        Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView(), Dune::VTK::nonconforming);

        //for (int i=0;i<6;++i) {
        std::stringstream str;

        str << "sol ";
        vtkwriter.addVertexData(field, str.str().c_str(), dim);
        vtkwriter.addCellData(lindisp, "celldisp_", dim);
        {
            using Container = Dune::BlockVector<Dune::FieldVector<double,1>>;
            using Function =  Dune::P1VTKFunctionVectorLin<typename GridType::LeafGridView, Container>;
            //Dune::VTK::FieldInfo dispinfo("disp",Dune::VTK::FieldInfo::Type::vector,3);
            Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
            std::string vtkname("disp");
            Dune::VTKFunction<typename GridType::LeafGridView>* Fp = new Function(grid.leafGridView(),
                                                                                  field,
                                                                                  vtkname,
                                                                                  3, prec);
            vtkwriter.addVertexData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(Fp));
        }
        vtkwriter.addCellData(pressforce, "pressforce");
        //vtkwriter.addVertexData(disp, dispinfo);

        {
            using Container = Dune::BlockVector<Dune::FieldVector<double,3>>;
            using Function =  Dune::P0VTKFunctionVector<typename GridType::LeafGridView, Container>;
            Dune::VTK::FieldInfo dispinfo("veccelldisp",Dune::VTK::FieldInfo::Type::vector,3);
            Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
            std::string vtkname("veccelldisp");
            Dune::VTKFunction<typename GridType::LeafGridView>* Fp = new Function(grid.leafGridView(),
                                                                                  disp,
                                                                                  vtkname,
                                                                                  3, prec);
            auto ptr = std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(Fp);
            vtkwriter.addCellData(ptr);
            //vtkwriter.addCellData(std::make_unique<Dune::VTKFunctionWrapper>(ptr), dispinfo);
            //vtkwriter.addCellData(ptr, dispinfo, 1);

        }

        vtkwriter.addCellData(stresslin, "stresscomp", 6);
        {
            using Container = std::vector<double>;//Dune::BlockVector<Dune::FieldVector<double,1>>;
            using Function =  Dune::P0VTKFunctionVectorLin<typename GridType::LeafGridView, Container>;
            Dune::VTK::FieldInfo dispinfo("disp",Dune::VTK::FieldInfo::Type::vector,3);
            Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
            std::string vtkname("stresstensor");
            Dune::VTKFunction<typename GridType::LeafGridView>* Fp = new Function(grid.leafGridView(),
                                                                                  stresslin,
                                                                                  vtkname,
                                                                                  6, prec);
            vtkwriter.addCellData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(Fp));
            //vtkwriter.addCellData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(fp));
        }
        {
            using Container = std::vector<double>;//Dune::BlockVector<Dune::FieldVector<double,1>>;
            using Function =  Dune::P0VTKFunctionVectorLin<typename GridType::LeafGridView, Container>;
            Dune::VTK::FieldInfo dispinfo("disp",Dune::VTK::FieldInfo::Type::vector,3);
            Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
            std::string vtkname("straintensor");
            Dune::VTKFunction<typename GridType::LeafGridView>* Fp = new Function(grid.leafGridView(),
                                                                                  strainlin,
                                                                                  vtkname,
                                                                                  6, prec);
            vtkwriter.addCellData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(Fp));
            //vtkwriter.addCellData(std::shared_ptr< const Dune::VTKFunction<typename GridType::LeafGridView> >(fp));
        }
        //vtkwriter.addTensorData(stress, "stress");
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
    //ok = run<Dune::CpGrid, ElasticitySolverTypeFem>(p, with_pressure, with_gravity, std::string("fem"));
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
