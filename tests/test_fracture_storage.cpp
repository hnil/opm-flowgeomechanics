#include <config.h>
#include <opm/geomech/FractureStorageCheckpoint.hpp>
#include <opm/geomech/RegularTrimesh.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
using Storage = Opm::FractureStorageCheckpoint;
using Triangle = Storage::Triangle;
using Point = Storage::Point;
namespace {
void close(double a, double b, const char* message) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b) > 1e-11*std::max({1.0,std::abs(a),std::abs(b)}))
        throw std::runtime_error(std::string(message)+": "+std::to_string(a)+" != "+std::to_string(b));
}
double total(const std::vector<double>& values) {return std::accumulate(values.begin(),values.end(),0.0);}
struct Mesh {std::vector<Triangle> triangles; std::vector<double> widths, volumes;};
Mesh geometry(const Opm::RegularTrimesh& mesh, int coarsen=0, bool smoothing=false) {
    auto [grid,map,boundary]=mesh.createDuneGrid(coarsen,{},smoothing,10);
    Mesh result;result.triangles.resize(map.size());result.widths.resize(map.size());result.volumes.resize(map.size());
    for (const auto& e:Dune::elements(grid->leafGridView())) {
        const auto i=grid->leafGridView().indexSet().index(e); const auto g=e.geometry();
        for (int c=0;c<3;++c)for(int d=0;d<3;++d)result.triangles[i][c][d]=g.corner(c)[d];
        result.widths[i]=0.001+0.00001*(g.center()[0]+2*g.center()[1]+30);
        result.volumes[i]=g.volume()*result.widths[i];
    }
    return result;
}
void analytic() {
    const Triangle t={Point{0,0,0},Point{2,0,0},Point{0,2,0}};
    const Triangle other={Point{2,2,0},Point{2,0,0},Point{0,2,0}};
    const Triangle half={Point{0,0,0},Point{1,0,0},Point{0,2,0}};
    Storage s; close(total(s.volumes({t})),0,"new fracture has no old fluid");
    s.capture({t},{0.003});close(total(s.volumes({t})),0.006,"identity");
    close(total(s.volumes({half})),0.003,"partial domain overlap");
    close(total(s.volumes({t,other})),0.006,"new area has zero old volume");
    const Triangle c1={Point{0,0,0},Point{1,0,0},Point{0,1,0}};
    const Triangle c2={Point{1,0,0},Point{2,0,0},Point{1,1,0}};
    const Triangle c3={Point{0,1,0},Point{1,1,0},Point{0,2,0}};
    const Triangle c4={Point{1,0,0},Point{1,1,0},Point{0,1,0}};
    for (double v:s.volumes({c1,c2,c3,c4}))close(v,.0015,"four-way refinement");
    s.capture({c1,c2,c3,c4},{.001,.002,.003,.004});
    close(total(s.volumes({t})),.005,"nonuniform fine-to-coarse volume");
    close(total(s.volumes({c4,c3,c2,c1})),.005,"repeated remap starts at immutable snapshot");
    s.capture({t,other},{.003,.007}); const auto a=s.volumes({other,t});
    close(a[0],.014,"permutation first cell");close(a[1],.006,"permutation second cell");
    s.volumes({half}); // Trial mesh; rollback must recover the unchanged checkpoint.
    const auto restored=s.volumes({t,other});close(restored[0],.006,"rollback first cell");close(restored[1],.014,"rollback second cell");
    auto offplane=t;for(auto& p:offplane)p[2]=1;
    bool rejected=false;try{s.volumes({offplane,other});}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("nonplanar transfer accepted");
    close(total(s.volumes({t,other})),.02,"failed transfer leaves cache usable");
    // Rotation, large reservoir-coordinate translation and reversed orientation.
    auto transform=[](Triangle q){for(auto& p:q){auto x=p[0],y=p[1];p={1e6+x,2e6+.6*y,3e6+.8*y};}std::swap(q[1],q[2]);return q;};
    s.capture({transform(t)},{.003});close(total(s.volumes({transform(half)})),.003,"rotated translated fracture");
    s.capture({t},{0});close(total(s.volumes({t})),0,"no artificial aperture floor");
}
void production() {
    Opm::RegularTrimesh mesh(12);auto old=geometry(mesh);
    Storage s;s.capture(old.triangles,old.widths);
    const double original=total(old.volumes);
    auto check=[&](const char* name,const std::vector<Triangle>& triangles){
        const auto start=std::chrono::steady_clock::now();auto values=s.volumes(triangles);
        close(total(values),original,name);
        std::cout<<name<<" cells="<<triangles.size()<<" old_volume="<<total(values)<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    };
    auto reversed=old.triangles;std::reverse(reversed.begin(),reversed.end());check("permutation",reversed);
    const auto perm=s.volumes(reversed);for(size_t i=0;i<perm.size();++i)close(perm[i],old.volumes[perm.size()-1-i],"permuted cell volume");
    auto grown=mesh;grown.expandGrid();check("growth",geometry(grown).triangles);
    check("refinement",geometry(mesh.refine()).triangles);
    auto coarse=mesh.coarsen(false);coarse.expandGrid();coarse.expandGrid();check("coarsening",geometry(coarse).triangles);
    check("multiresolution",geometry(mesh,2,false).triangles);
    // Smoothing may change the boundary domain; using an expanded target isolates
    // the remapping from that separate domain change.
    check("smoothed_multiresolution_growth",geometry(grown,2,true).triangles);
    check("retry_original",old.triangles);
    const auto restored=s.volumes(old.triangles);for(size_t i=0;i<restored.size();++i)close(restored[i],old.volumes[i],"rollback cell volume");
    double indexed=0;const auto target=geometry(coarse);
    for(size_t i=0;i<target.volumes.size()&&i<old.widths.size();++i)indexed+=target.volumes[i]/target.widths[i]*old.widths[i];
    if(std::abs(indexed-original)<1e-5)throw std::runtime_error("fixture does not expose old indexed-storage error");
    std::cout<<"regression old_volume="<<original<<" direct_index_volume="<<indexed<<" transferred_volume="<<total(s.volumes(target.triangles))<<'\n';
}
}
int main(int argc,char** argv) {
    Dune::MPIHelper::instance(argc,argv);
    std::cout.precision(17);
    try {analytic();production();std::cout<<"PASS fracture storage transfer\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
