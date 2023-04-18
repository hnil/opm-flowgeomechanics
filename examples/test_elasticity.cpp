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
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCMECHConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/utility/StopWatch.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/geomech/elasticity_solver.hpp>
#include <opm/elasticity/matrixops.hpp>

#include <cstring>
#include <iostream>
#include "vectorfunctions.hh"
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

template<int dimension>
unsigned cartesianIndex(const std::array<int,dimension>& coords,
                        const std::array<int, dimension>& cartesianDimensions
    )
{
    unsigned cartIndex = coords[0];
    int factor = cartesianDimensions[0];
    for (unsigned i = 1; i < dimension; ++i) {
        cartIndex += coords[i]*factor;
        factor *= cartesianDimensions[i];
    }

    return cartIndex;
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
    using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
    Opm::Parser parser;
    auto deck = parser.parseFile(p.file);
    Opm::EclipseGrid inputGrid(deck);
    grid.processEclipseFormat(&inputGrid, nullptr, false);
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
    CartesianIndexMapper cartesianIndexMapper(grid);
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
    std::vector<size_t> fixed_nodes;
    const auto& bcconfig = eclState.getSimulationConfig().bcconfig();
    if (bcconfig.size() > 0) {
        //nonTrivialBoundaryConditions_ = true;
        auto gv = grid.leafGridView();
        size_t numCartDof = cartesianIndexMapper.cartesianSize();
        unsigned numElems = gv.size(/*codim=*/0);
        std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx){
            cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
        }
        std::array<int, 3> cartdim = cartesianIndexMapper.cartesianDimensions();    
        for (const auto& bcface : bcconfig) {
            const auto& type = bcface.bcmechtype;
            if((bcface.i1 < 0) || (bcface.j1<0) || (bcface.k1<0)){
                throw std::logic_error("Lower range of BC wrong");
            }
            if( (bcface.i2 > cartdim[0]) || (bcface.j2> cartdim[1]) || (bcface.k2 > cartdim[2])){
                throw std::logic_error("Upper range of BC wrong");
            }
            if (type == Opm::BCMECHType::FREE) {
                // do nothing
            }else if (type == Opm::BCMECHType::FIXED) {
                std::set<size_t> effected_cells;
                for (int i = bcface.i1; i <= bcface.i2; ++i) {
                    for (int j = bcface.j1; j <= bcface.j2; ++j) {
                        for (int k = bcface.k1; k <= bcface.k2; ++k) {
                            
                            
                            std::array<int, 3> tmp = {i,j,k};
                            int cartindex =
                                cartesianIndex<3>(tmp,cartdim);
                            auto elemIdx = cartesianToCompressedElemIdx[cartindex];
                            if (elemIdx>-1){
                                effected_cells.insert(elemIdx);
                            }
                        }
                    }
                }
                std::cout << "Effected cells" << std::endl;
                for(int i:effected_cells){
                    std::cout << i << std::endl;
                }
                const auto& gv = grid.leafGridView();
                for(const auto& cell:elements(gv)){
                    auto index = gv.indexSet().index(cell);
                    auto it = effected_cells.find(index);
                    if(!(it == effected_cells.end())){
                        // fix all noted for now
                        for (const auto& vertex : Dune::subEntities(cell, Dune::Codim<dim>{})){
                            fixed_nodes.push_back(gv.indexSet().index(vertex));
                        }
                    }
                }                    
            } else {    
                throw std::logic_error("invalid type for BC. Use FREE or RATE");
            }
        }
    }
    std::sort(fixed_nodes.begin(), fixed_nodes.end()); // {1 1 2 3 4 4 5}
    auto last = std::unique(fixed_nodes.begin(), fixed_nodes.end());
    // v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
    fixed_nodes.erase(last, fixed_nodes.end());
    std::cout << "Effected nodes" << std::endl;
    for(int i:fixed_nodes){
        std::cout << i << std::endl;
    }
    
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

    {
        //auto gv = grid.leafGridView();
        // using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
        // Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(), Dune::mcmgVertexLayout());        
        // for(const auto& vert : Dune::vertices(gv,Dune::Partitions::border)){
        //     int indexi = 0 ;//mapper.index(vert);
        //     Dune::FieldVector<double,3> value;
        //     value = 0;
        //     esolver.A.updateFixedNode(indexi,std::make_pair(Opm::Elasticity::XYZ,value));
        // }
        // for (const auto& cell: Dune::elements(gv)){
        //     for (const auto& is: Dune::intersections(gv,cell)){
        //         if(is.boundary()){
        //             auto normal = is.centerUnitOuterNormal();
        //             //for (const auto& vert: Dune::edges(gv,is)){
        //             //}
        //         }
        //     }
        // }
    }

//   upscale.fixCorners(p.min, p.max);
    bool do_matrix = true;//assemble matrix
    bool do_vector = true;//assemble matrix
    esolver.fixNodes(fixed_nodes); 
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
     Dune::BlockVector<Dune::FieldVector<double,3>> disp;
     disp.resize(pressforce.size());
     for(size_t i = 0; i < pressforce.size(); ++i){
         for(size_t k = 0; k < dim; ++k){
             disp[i][k] = field[i*dim+k];
         }
     }

    if (!p.vtufile.empty()) {
        Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView());
        
        //for (int i=0;i<6;++i) {
        std::stringstream str;
        
        str << "sol ";
        vtkwriter.addVertexData(field, str.str().c_str(), dim);
        using Container = Dune::BlockVector<Dune::FieldVector<double,1>>;
        using Function =  Dune::P1VTKFunctionVector<GridView, Container>;
        //Dune::VTK::FieldInfo dispinfo("disp",Dune::VTK::FieldInfo::Type::vector,3);
        Dune::VTK::Precision prec = Dune::VTK::Precision::float32;
        std::string vtkname("disp");
        Dune::VTKFunction<GridView>* fp = new Function(grid.leafGridView(),
                                                       field,
                                                       vtkname,
                                                       3, 0, prec);
        vtkwriter.addVertexData(std::shared_ptr< const Dune::VTKFunction<GridView> >(fp));
        //vtkwriter.addVertexData(disp, dispinfo);
        vtkwriter.addCellData(pressforce, "pressforce");
        //vtkwriter.addVertexData(disp, stress);
        //}
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
