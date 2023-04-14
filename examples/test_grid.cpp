#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <set>
#include <dune/grid/yaspgrid.hh>

template<class vec>
bool my_equal(const vec& s1, const vec& s2)
{
    auto tmp = s1-s2;
    return tmp.two_norm()< 1e-10;
}

int main(int argc, char** argv){
    Dune::MPIHelper::instance(argc,argv); 
    constexpr int dim = 3;
    using Grid = Dune::YaspGrid<dim>;
    using Geometry = Grid::template Codim<0>::Geometry;
    using CoordType =  Geometry::GlobalCoordinate;
    Dune::FieldVector<double,dim> L = {4.0, 4.0, 4.0};
    std::array<int,dim> s = {2,2,2};
    Grid grid( L , s);
    auto gv = grid.leafGridView();
    //auto set = gv.indexSet();
    //asmhandler_impl.hpp
    //elasticity_upscale_impl.hpp
    for(const auto& cell:elements(gv)){
        //auto gm = cell->geometry();
        //int index2 = set.subIndex(*cell,j,dim);
        // coordinates of subindices
        std::vector<CoordType> bcoords;
        for (const auto& itsec: Dune::intersections(gv,cell)){            
            if(itsec.boundary()){
                auto boundary = itsec.geometry();
                int numc = boundary.corners();
                for(int i = 0; i<numc; ++i){
                    bcoords.push_back(boundary.corner(i));
                }
                //std::cout << boundary << std::endl;
                //bcoords.push_back();                 
            }
        }

        std::vector<CoordType> coords;
        std::vector<size_t> gnodes;
        for (const auto& vertex : Dune::subEntities(cell, Dune::Codim<Grid::dimension>{})){
            auto point = vertex.geometry().center();
            auto is_equal = [point](CoordType v1){return my_equal(v1,point);};
            auto it = std::find_if(std::begin(bcoords),std::end(bcoords), is_equal);
            bool is_boundary = !(it==bcoords.end());
            if(is_boundary){
                gnodes.push_back(gv.indexSet().index(vertex));
                coords.push_back(point);
                std::cout << gv.indexSet().index(vertex) << " " << vertex.geometry().center() << std::endl;
            }
        }
        
        //get coordinates of boundary
        //for(auto it = itsec.begin(); it!=itsec.end();++it){
        //it->index();
        //}
        
        //for(size_t i = 0; 
        
    }
}
