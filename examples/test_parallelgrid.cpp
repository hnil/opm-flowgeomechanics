#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <set>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/owneroverlapcopy.hh>
//#include <dune/grid/utility/parmetisgridpartitioner.hh>

template<class vec>
bool my_equal(const vec& s1, const vec& s2)
{
    auto tmp = s1-s2;
    return tmp.two_norm()< 1e-10;
}

int main(int argc, char** argv){
    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc,argv);
    Dune::Communication<MPI_Comm> world_comm = Dune::MPIHelper::getCommunication();
    constexpr int dim = 3;
    if (0 == mpihelper.rank()){
        std::cout << "Using " << mpihelper.size() << " Processes." << std::endl;
    }
    typedef std::size_t GlobalId;
	typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
	Communication comm(world_comm);
    using Grid = Dune::YaspGrid<dim>;
    using GridView = Grid::LeafGridView;
    using Geometry = Grid::template Codim<0>::Geometry;
    using CoordType =  Geometry::GlobalCoordinate;
    std::bitset<dim> periodic("000");
    Dune::FieldVector<double,dim> L = {4.0, 4.0, 4.0};
    std::array<int,dim> s = {4,4,4};
    Grid grid( L , s,periodic,1,MPI_COMM_WORLD);
    
    {
    auto gv = grid.leafGridView();
    for(int rank=0; rank<mpihelper.size(); ++rank){
        if(rank == mpihelper.rank()){
            std::cout << "Rank " << rank << " has " << gv.size(0) << " cells." << std::endl;
        }
        world_comm.barrier();
        //Dune::Communication<MPI_Comm>::barrier();
        //mpihelper.barrier();
    }
    }

    //         partition = Dune::graphPartition(MatrixGraph(A),
	// 				 comm,
	// 				 static_cast<int>(world_comm.size()),
	// 				 false);
	// nparts=static_cast<int>(world_comm.size());

    world_comm.barrier();
    {
    auto gv = grid.leafGridView();
    //std::vector<unsigned int> part = Dune::ParMetisGridPartitioner<GridView>::partition(gv, mpiHelper);
    std::vector<unsigned int> part(gv.size(0),0);
    grid.loadBalance(part);
    }
    world_comm.barrier();
    {
    auto gv = grid.leafGridView();
    for(int rank=0; rank<mpihelper.size(); ++rank){
        if(rank == mpihelper.rank()){
            std::cout << "Rank " << rank << " has " << gv.size(0) << " cells." << std::endl;
        }
        world_comm.barrier();
        //mpihelper.barrier();
    }
    }
    //auto set = gv.indexSet();
    //asmhandler_impl.hpp
    //elasticity_upscale_impl.hpp
    // if(false){
    // for(const auto& cell:elements(gv)){
    //     //auto gm = cell->geometry();
    //     //int index2 = set.subIndex(*cell,j,dim);
    //     // coordinates of subindices
    //     std::vector<CoordType> bcoords;
    //     for (const auto& itsec: Dune::intersections(gv,cell)){            
    //         if(itsec.boundary()){
    //             auto boundary = itsec.geometry();
    //             int numc = boundary.corners();
    //             for(int i = 0; i<numc; ++i){
    //                 bcoords.push_back(boundary.corner(i));
    //             }
    //             //std::cout << boundary << std::endl;
    //             //bcoords.push_back();                 
    //         }
    //     }

    //     std::vector<CoordType> coords;
    //     std::vector<size_t> gnodes;
    //     for (const auto& vertex : Dune::subEntities(cell, Dune::Codim<Grid::dimension>{})){
    //         auto point = vertex.geometry().center();
    //         auto is_equal = [point](CoordType v1){return my_equal(v1,point);};
    //         auto it = std::find_if(std::begin(bcoords),std::end(bcoords), is_equal);
    //         bool is_boundary = !(it==bcoords.end());
    //         if(is_boundary){
    //             gnodes.push_back(gv.indexSet().index(vertex));
    //             coords.push_back(point);
    //             std::cout << gv.indexSet().index(vertex) << " " << vertex.geometry().center() << std::endl;
    //         }
    //     }
        
    //     //get coordinates of boundary
    //     //for(auto it = itsec.begin(); it!=itsec.end();++it){
    //     //it->index();
    //     //}
        
    //     //for(size_t i = 0; 
    // }
        
    // }
    return 0;
}
