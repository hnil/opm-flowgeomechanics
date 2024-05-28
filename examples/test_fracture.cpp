/** construct a grid with vertices on the unit circle and element parametrization */

#include <config.h>
#include <iostream>
#include <cmath>
#include <memory>
#include <functional>
#include <string>
#include <iostream>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/foamgrid/foamgrid.hh>

// from pdelab
#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#define ALBERTA_DIM 2
#include <dune/pdelab.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/grid/yaspgrid/partitioning.hh>
/**
 * \brief Mapping class mapping from a secant on the unit circle onto the circle
 */
template<typename ctype, int dimgrid, int dimworld>
class UnitCircleMapping
{
  double fromAngle_;
  double toAngle_;

public:
  UnitCircleMapping(double fromAngle, double toAngle) : fromAngle_(fromAngle), toAngle_(toAngle) {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  Dune::FieldVector<ctype,dimworld> operator() (const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    const double angle = fromAngle_ + x[0]*(toAngle_ - fromAngle_);
    return {std::cos(angle), std::sin(angle)};
  }
};




template<typename ctype, int dimgrid, int dimworld>
class IdentityMapping
{
  const std::array<Dune::FieldVector<double, 3>, 3> vertices_;

public:
  IdentityMapping(const std::array<Dune::FieldVector<double, 3>, 3>& vertices) : vertices_(vertices) {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  Dune::FieldVector<ctype,dimworld> operator() (const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    // calculate global coordinate
    Dune::FieldVector<double, 3> shapeFunctions = evaluateShapeFunctions_(x);
    Dune::FieldVector<ctype,dimworld> y{0,0,0};
    for(size_t i = 0; i < y.size(); i++)
      for(size_t j = 0; j < 3; j++)
        y[j] += vertices_[i][j]*shapeFunctions[i];
    // project it on the unit sphere
    //y /= y.two_norm();
    return y;
  }
private:
  Dune::FieldVector<double, 3> evaluateShapeFunctions_(const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    Dune::FieldVector<double, 3> out;
    out[0] = 1.0;
    for (size_t i=0; i<2; i++)
    {
      out[0]  -= x[i];
      out[i+1] = x[i];
    }
    return out;
  }

};
/**
 * \brief Mapping class mapping from a triangle with points on the unit sphere onto the sphere
  *       with theta in [0, pi] and phi in [0, 2*pi)
 */
template<typename ctype, int dimgrid, int dimworld>
class UnitSphereMapping
{
  const std::array<Dune::FieldVector<double, 3>, 3> vertices_;

public:
  UnitSphereMapping(const std::array<Dune::FieldVector<double, 3>, 3>& vertices) : vertices_(vertices) {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  Dune::FieldVector<ctype,dimworld> operator() (const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    // calculate global coordinate
    Dune::FieldVector<double, 3> shapeFunctions = evaluateShapeFunctions_(x);
    Dune::FieldVector<ctype,dimworld> y{0,0,0};
    for(size_t i = 0; i < y.size(); i++)
      for(size_t j = 0; j < 3; j++)
        y[j] += vertices_[i][j]*shapeFunctions[i];
    // project it on the unit sphere
    y /= y.two_norm();
    return y;
  }
private:
  Dune::FieldVector<double, 3> evaluateShapeFunctions_(const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    Dune::FieldVector<double, 3> out;
    out[0] = 1.0;
    for (size_t i=0; i<2; i++)
    {
      out[0]  -= x[i];
      out[i+1] = x[i];
    }
    return out;
  }

};

/**
 * \brief Method to calculate the vector update for a single time step advance
 */
template<class GridView, class Mapper>
void evolve (const GridView& gridView,
           const Mapper& mapper,
           std::vector<double>& temperature,
           const double lambda,
           double& dt)
{
  // allocate a temporary vector for the update
  std::vector<double> update(temperature.size());
  std::fill(update.begin(), update.end(), 0.0);

  // initialize dt very large
  dt = std::numeric_limits<double>::max();
  double h = std::numeric_limits<double>::max();

  for (auto&& element : elements(gridView))
  {
    int eIdx = mapper.index(element);
    auto eCenter = element.geometry().center();
    // iterator over all intersections
    for (auto&& is : intersections(gridView, element))
    {
      // index of the neighbour
        if(!is.boundary()){
            int nIdx = mapper.index(is.outside());


            // calculate distance between the midpoints

            auto nCenter = is.outside().geometry().center();
            auto isCenter = is.geometry().center();
            eCenter -= isCenter; nCenter -= isCenter;
            double dist = eCenter.two_norm() + nCenter.two_norm();

            //approximate h as the distance to the neihgbour center
            h = std::min(h, dist);

            // approximate gradient
            double gradTn = (temperature[nIdx] - temperature[eIdx])/dist;

            // add to update
            update[eIdx] += lambda*gradTn;

        }

    }
    if(eCenter.two_norm() < 3){
        update[eIdx] += (10 - temperature[eIdx]);
    }
  }

  // CFL criterion
  dt = std::min(dt, h*h/2.0/lambda);

  // scale dt with safety factor
  dt *= 0.99;

  // update the concentration vector
  for (unsigned int i=0; i<temperature.size(); ++i)
    temperature[i] += dt*update[i];
}

struct RestrictedValue
{
  double value;
  int count;
  RestrictedValue ()
  {
    value = 0;
    count = 0;
  }
};

template<class Grid, class Mapper>
bool finitevolumeadapt (Grid& grid, Mapper& mapper, std::vector<double>& temperature, int lmin, int lmax)
{
  // tol value for refinement strategy
  const double refinetol  = 0.05;
  const double coarsentol = 0.001;

  // grid view types
  typedef typename Grid::LeafGridView LeafGridView;

  // get grid view on leaf grid
  LeafGridView leafGridView = grid.leafGridView();

  // compute cell indicators
  std::vector<double> indicator(temperature.size(), std::numeric_limits<double>::lowest());
  double globalmax = std::numeric_limits<double>::lowest();
  double globalmin = std::numeric_limits<double>::max();
  for (auto&& element : elements(leafGridView))
  {
    // element index
    int eIdx = mapper.index(element);

    // global min/max
    globalmax = std::max(globalmax,temperature[eIdx]);
    globalmin = std::min(globalmin,temperature[eIdx]);

    for (auto&& intersection : intersections(leafGridView, element))
    {
      if( !intersection.neighbor() )
        continue;

      // get the neighbor
      const auto& outside = intersection.outside();
      int nIdx = mapper.index(outside);

      // handle face from one side only
      if ( element.level() > outside.level() ||
           (element.level() == outside.level() && eIdx < nIdx) )
      {
        double localdelta = std::abs(temperature[nIdx]-temperature[eIdx]);
        indicator[eIdx] = std::max(indicator[eIdx],localdelta);
        indicator[nIdx] = std::max(indicator[nIdx],localdelta);
      }
    }
  }

  // mark cells for refinement/coarsening
  double globaldelta = globalmax-globalmin;
  int marked=0;
  for (auto&& element : elements(leafGridView))
  {
    if (indicator[mapper.index(element)]>refinetol*globaldelta
        && (element.level() < lmax || !element.isRegular()))
    {
      grid.mark( 1, element );
      ++marked;
tracessy      for(auto&& intersection : intersections(leafGridView, element))
      {
        if( !intersection.neighbor() )
          continue;

        const auto& outside = intersection.outside();
        if( (outside.level() < lmax) || !outside.isRegular() )
          grid.mark( 1, outside );
      }
    }
    if (indicator[mapper.index(element)] < coarsentol*globaldelta && element.level() > lmin)
    {
      grid.mark( -1, element);
      ++marked;
    }
  }
  if( marked==0 )
    return false;

  grid.preAdapt();

  typedef Dune::PersistentContainer<Grid,RestrictedValue> RestrictionMap;
  RestrictionMap restrictionmap(grid,0); // restricted temperature
  for (int level=grid.maxLevel(); level>=0; level--)
  {
    // get grid view on level grid
    for (auto&& element : elements(grid.levelGridView(level)))
    {
      // get your map entry
      RestrictedValue& rv = restrictionmap[element];
      // put your value in the map
      if (element.isLeaf())
      {
        int eIdx = mapper.index(element);
        rv.value = temperature[eIdx];
        rv.count = 1;
      }

      // average in father
      if (element.level() > 0)
      {
        RestrictedValue& rvf = restrictionmap[element.father()];
        rvf.value += rv.value/rv.count;
        rvf.count += 1;
      }
    }
  }

  // adapt mesh and mapper
  bool refined=grid.adapt();
  mapper.update(leafGridView);
  restrictionmap.resize();
  temperature.resize(mapper.size());

  // interpolate new cells, restrict coarsened cells
  for (int level=0; level<=grid.maxLevel(); level++)
  {
    for (auto&& element : elements(grid.levelGridView(level)))
    {
      // check map entry
      if (!element.isNew() )
      {
        // entry is in map, write in leaf
        if (element.isLeaf())
        {
          RestrictedValue& rv = restrictionmap[element];
          int eIdx = mapper.index(element);
          temperature[eIdx] = rv.value/rv.count;
        }
      }
      else
      {
        // value is not in map, interpolate from father element
        assert (element.level() > 0);
        RestrictedValue& rvf = restrictionmap[element.father()];
        if (element.isLeaf())
        {
          int eIdx = mapper.index(element);
          temperature[eIdx] = rvf.value/rvf.count;
        }
        else
        {
          // create new entry
          RestrictedValue& rv = restrictionmap[element];
          rv.value = rvf.value/rvf.count;
          rv.count = 1;
        }
      }
    }
  }
  grid.postAdapt();

  return refined;
}


void twoDimensionalTest ()
{
  typedef Dune::FoamGrid<2, 3> Grid;
  typedef Grid::ctype ctype;
  const int dimgrid = Grid::dimension;
  const int dimworld = Grid::dimensionworld;

  // Start grid creation
  Dune::GridFactory<Grid> factory;

  // Create an icosahedron
  // The list of grid vertex positions
  //double tao = (1.0 + std::sqrt(5))/2.0; // golden ratio
  unsigned N = 6;
  double radius = 1;
  std::vector<Dune::FieldVector<double, dimworld> > vertices;
  {
      Dune::FieldVector<double, dimworld> vertex({0,0, 0});
      vertices.push_back(Dune::FieldVector<double, dimworld>({0, 0, 0}));
  }
  std::vector<unsigned> inner_indices;
  for(unsigned i=0;i < N; ++i){
      inner_indices.push_back(i+1);
      double theta = (i*2*M_PI)/N;
      double x = radius*cos(theta);
      double y = radius*sin(theta);
      double z = x*x+y*y;
      Dune::FieldVector<double, dimworld> vertex({x,y, z});
      vertices.push_back(vertex);
  }


  //std::vector<Dune::FieldVector<double, dimworld> > vertices({{1,tao,0},{-1,tao,0},{1,-tao,0},{-1,-tao,0},
  //                                                            {0,1,tao},{0,-1,tao},{0,1,-tao},{0,-1,-tao},
  //                                                            {tao,0,1},{-tao,0,1},{tao,0,-1},{-tao,0,-1}});

  // Create the grid vertices

  for (size_t i=0; i<vertices.size(); i++)
  {
    // make vertices live on unit sphere
    // insert vertices
    factory.insertVertex(vertices[i]);
  }
  std::vector<std::vector<unsigned int> > cornerIDs;
  for(unsigned i=0;i < N; ++i){
      unsigned next = inner_indices[(i+1)%N];
      std::vector<unsigned int> cornerID({0,inner_indices[i],next});
      cornerIDs.push_back(cornerID);
  }
  for(unsigned i=0;i < unsigned(N); ++i){
      //std::vector<unsigned int> cornerID({0,i+1,i+2});
      // const auto mapping = IdentityMapping<ctype, dimgrid, dimworld>(
      //    {vertices[cornerIDs[i][0]],vertices[cornerIDs[i][1]],vertices[cornerIDs[i][2]]});
      // const auto mapping = IdentityMapping<ctype, dimgrid, dimworld>(
      //     {vertices[cornerID[0]],vertices[cornerID[1]],vertices[cornerID[2]]});
      factory.insertElement(Dune::GeometryTypes::simplex(2), cornerIDs[i]);// std::move(mapping));
  }

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView(), Dune::VTK::nonconforming);
  writer.write("disk_0");

  // refine the grid
  //grid->globalRefine(3);
  int rfac=2;
  int nlinear = 3;
  int layers = 0;
  std::vector<unsigned> out_indices;
  while(layers < 6){
      out_indices.resize(0);
      N= inner_indices.size();
      if(layers < nlinear){
          for(unsigned i = 0; i < N*rfac; ++i){
              double theta = (i*2*M_PI)/(N*rfac);
              double out_radius = radius + layers+1;
              double x = out_radius*cos(theta);
              double y = out_radius*sin(theta);
              double z = (x*x+y*y)*0.1;
              Dune::FieldVector<double, dimworld> vertex({x, y, z});
              int new_ind = grid->insertVertex(vertex);
              out_indices.push_back(new_ind);
          }


          for(unsigned i = 0; i < N; ++i){
              {
                  std::vector<unsigned int> cornerID({inner_indices[i],out_indices[(2*i)%(N*2)],out_indices[(2*i+1)%(N*2)]});
                  grid->insertElement(Dune::GeometryTypes::simplex(2),cornerID);
              }
              {
                  std::vector<unsigned int> cornerID({inner_indices[i],out_indices[(2*i+1)%(N*2)],inner_indices[(i+1)%N]});
                  grid->insertElement(Dune::GeometryTypes::simplex(2),cornerID);
              }
              {
                  std::vector<unsigned int> cornerID({inner_indices[(i+1)%N],out_indices[(2*i+1)%(N*2)],out_indices[(2*i+2)%(N*2)]});
                  grid->insertElement(Dune::GeometryTypes::simplex(2),cornerID);
              }


          }
      }else{
          for(unsigned i = 0; i < N; ++i){

              double theta = (i*2*M_PI)/(N);
              theta+=(layers-nlinear)*0.5*(2*M_PI/N);
              // if((layers%2)==1){
              //     theta+=0.5*(2*M_PI/N);
              // }else{
              //     theta-=0.5*(2*M_PI/N);
              // }
              double out_radius = radius + layers+1;
              double x = out_radius*cos(theta);
              double y = out_radius*sin(theta);
              double z = std::atan(((x*x+y*y)*0.1)/10);
              Dune::FieldVector<double, dimworld> vertex({x, y, z});
              int new_ind = grid->insertVertex(vertex);
              out_indices.push_back(new_ind);
          }
          for(unsigned i = 0; i < N; ++i){
              {
                  std::vector<unsigned int> cornerID({inner_indices[i%N],out_indices[(i)%(N)],inner_indices[(i+1)%(N)]});
                  grid->insertElement(Dune::GeometryTypes::simplex(2),cornerID);
              }

              {
                  std::vector<unsigned int> cornerID({inner_indices[(i+1)%N],out_indices[(i)%(N)],out_indices[(i+1)%(N)]});
                   grid->insertElement(Dune::GeometryTypes::simplex(2),cornerID);
              }

          }
      }
      ++layers;
      inner_indices = out_indices;
      grid->grow();
      grid->postGrow();
      std::string filename = "disk_" + std::to_string(layers);
      //std::stringstream ss("disk_");
      //ss << layers;
      std::cout << "Write file " << filename << std::endl;
      writer.write(filename.c_str());
  }

  // global refine
  grid->globalRefine(1);
  {
      std::string filename = "disk_" + std::to_string(7);
      writer.write(filename.c_str());
  }
  //move grid
  for(const auto& vertex: Dune::vertices(grid->leafGridView())){
          auto vgeo = vertex.geometry();
          auto center = vgeo.center();

          double x = center[0];
          double y = center[1];
          double z = center[2];
          x=x*2;
          y=y*1.5;
          Dune::FieldVector<double, dimworld> vertexPos({x, y, z});
          grid->setPosition(vertex,vertexPos);
  }
  {
      std::string filename = "disk_" + std::to_string(8);
      writer.write(filename.c_str());
  }


  // start try/catch block to get error messages from dune
  try {
      // make a mapper for codim 0 entities in the leaf grid
      Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>
          mapper(grid->leafGridView(), Dune::mcmgElementLayout());

      // the primary variable vector
      std::vector<double> temperature(mapper.size());

      // initial conditions
      temperature[0] = 1.0;

      // the time
      double t = 0.0;
      double dt;
      double tEnd = 10.0;
      int timestep = 0;

      // write output only every nth timestep
      int episode = 10;

      // the heat conductivity
      double lambda = 1.0;

      // Write pvd header
      Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView(), "temperature_2d", ".", "");
      vtkWriter.addCellData(temperature, "celldata");
      vtkWriter.write(t);

      // do the time integration
      while(t <= tEnd)
      {
          // do adaptation
          // int lmax = 4;
          // for(int i = 0; i < lmax; i++)
          //     finitevolumeadapt(*grid, mapper, temperature, 1, lmax);

          // apply finite volume scheme
          evolve(grid->leafGridView(), mapper, temperature, lambda, dt);

          //one time step forward
          t += dt;

          //write a vtk
          if(!(timestep%episode))
              vtkWriter.write(t);

          // Output some infos
          std::cout << "Time step " << timestep << " done, t = " << t << ", dt = " << dt << std::endl;
          // Increment time step counter
          ++timestep;
      }

      {
          Dune::FieldVector<double,3> L(10);
           Dune::FieldVector<double,3> L2(-10);
          std::array<int,3> dims(Dune::filledArray<3,int>(2));
          using GRID3D = Dune::YaspGrid<3>;
          //Dune::YLoadBalance< 3 >* part;
          //Dune::Communication<MPI_Comm> comm;
          //Dune::CollectiveCommunication comm()
          //Dune::YLoadBalance<dim>* lb = Dune::defaultLoadbalancer()
          int overlap =1;
          std::bitset<std::size_t{3}> periodic(0);
          //GRID3D grid3D(L2,L,dims, periodic, overlap, comm, lb);// comm, part );
          GRID3D grid3D(L,dims);//, periodic, overlap, comm, lb);// comm, part );
          Dune::VTKSequenceWriter<GRID3D::LeafGridView> vtkWriter3D(grid3D.leafGridView(),"grid3D",".","");
          vtkWriter3D.write(0);
          grid3D.globalRefine(2);
          vtkWriter3D.write(1);
      }

  }
  catch (std::exception & e) {
     std::cout << "Exception: " << e.what() << std::endl;
     throw;
  }
  // catch (...) {
  //   std::cout << "Unknown ERROR" << std::endl;
  //   return 1;
  // }



}


int main (int argc, char *argv[]) try
{
    Dune::MPIHelper::instance(argc, argv);
    std::cout << std::endl << "Running example for FoamGrid<2, 3>" << std::endl;
    twoDimensionalTest();

}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
  std::cout << e << std::endl;
  return 1;
}
