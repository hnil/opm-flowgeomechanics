// SPDX-License-Identifier: GPL-3.0-or-later
#include "config.h"
#include <opm/geomech/GeometryHelpers.hpp>
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv)
{
    Dune::MPIHelper::instance(argc,argv);
    Dune::CpGrid grid;
    const int variant=argc>1?std::stoi(argv[1]):0;
    const bool refined=argc>2&&std::string(argv[2])=="lgr";
    std::string error;
    try {
        if(refined) {
            grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
            Opm::EclipseGrid input(4,4,4);
            grid.processEclipseFormat(&input,nullptr,false,false,false,true,0.);
            // This synthetic input has no EclipseState/CARFIN. Supply the
            // refinement builder explicitly rather than relying on retention
            // enabled only for deck-declared refinements.
            grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
                std::array<int,3>{4,4,4}, input.getCOORD(), input.getZCORN(), input.getACTNUM()));
            grid.addLgrsUpdateLeafView({{2,2,2}},{{0,0,0}},{{1,1,1}},{"FINE"});
        }else grid.createCartesian({4,4,4},{1.,1.,1.});
        if(grid.comm().size()>1) {
            const Dune::CartesianIndexMapper<Dune::CpGrid> rootMapper(grid);
            std::vector<int> parts(grid.size(0));
            for(int c=0;c<grid.size(0);++c) {
                const int root=rootMapper.cartesianIndex(c);
                parts[c]=(variant==0?root%4:(root/4)%4)%grid.comm().size();
            }
            grid.loadBalance(parts,false,true,1);
        }
        const auto gv=grid.leafGridView();
        const Dune::CartesianIndexMapper<Dune::CpGrid> mapper(grid);
        using Seed=Dune::CpGrid::Codim<0>::Entity::EntitySeed;
        std::vector<Seed> seeds(gv.size(0));
        int lower=-1,upper=-1;
        for(const auto& cell:elements(gv)) {
            const int c=gv.indexSet().index(cell);seeds[c]=cell.seed();
            if(cell.partitionType()!=Dune::InteriorEntity)continue;
            if(refined) {
                const auto q=cell.geometry().center();
                if((q-Dune::FieldVector<double,3>{.25,.25,.25}).two_norm()<1e-12)lower=c;
                if((q-Dune::FieldVector<double,3>{.25,.25,.75}).two_norm()<1e-12)upper=c;
            }else {
                if(mapper.cartesianIndex(c)==0)lower=c;
                if(mapper.cartesianIndex(c)==16)upper=c;
            }
        }
        if(grid.comm().rank()==0) {
            if(lower<0||upper<0)throw std::runtime_error("Test partition lost its two owned face neighbours");
            // The upper cell is the known winner on every tested partition:
            // Cartesian global IDs 16 > 0; conforming LGR IDs 193 > 189.
            const double scale=refined?.5:1.;
            for(bool reverse:{false,true}) {
                std::vector<std::size_t> ids{std::size_t(lower),std::size_t(upper)};
                if(reverse)std::reverse(ids.begin(),ids.end());
                // Identical conservative boxes deliberately change broad-phase
                // order while preserving both legitimate nearest candidates.
                external::cvf::BoundingBox box;
                box.add(external::cvf::Vec3d(0,0,0));box.add(external::cvf::Vec3d(scale,scale,2*scale));
                std::vector<external::cvf::BoundingBox> boxes(2,box);
                external::cvf::ref<external::cvf::BoundingBoxTree> tree(new external::cvf::BoundingBoxTree);
                tree->buildTreeFromBoundingBoxes(boxes,&ids);
                for(double z:{.25,1.,1.75}) {
                    const int selected=external::cellOfPoint(tree,grid,seeds,external::cvf::Vec3d(.5*scale,.5*scale,z*scale));
                    const int expected=z<1?lower:upper;
                    if(selected!=expected) {
                        std::cerr << "z=" << z << " selected=" << selected << " expected=" << expected << " lowerId=" << grid.globalIdSet().id(grid.entity(seeds[lower])) << " upperId=" << grid.globalIdSet().id(grid.entity(seeds[upper])) << "\n";
                        throw std::runtime_error("Reservoir-face selection depends on BVH order or local cell numbering");
                    }
                }
                if(external::cellOfPoint(tree,grid,seeds,external::cvf::Vec3d(4,4,4))!=-1)
                    throw std::runtime_error("Outside point acquired a reservoir cell");
            }
        }
    }catch(const std::exception& e){error=e.what();}
    catch(...){error="Unexpected exception during cell-search regression";}
    const int failed=grid.comm().max(int(!error.empty()));
    if(!error.empty())std::cerr<<error<<'\n';
    if(!failed&&grid.comm().rank()==0)std::cout<<"Canonical reservoir-face selection passed\n";
    return failed;
}
