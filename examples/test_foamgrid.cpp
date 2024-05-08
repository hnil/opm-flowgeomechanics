/** construct a grid with vertices on the unit circle and element parametrization */

#include <config.h>
#include <iostream>
#include <cmath>
#include <memory>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/foamgrid/foamgrid.hh>

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

/**
 * \brief Mapping class mapping from a triangle with points on the unit sphere onto the sphere
  *       with theta in [0, pi] and phi in [0, 2*pi)
 */
template<typename ctype, int dimgrid, int dimworld>
class IdenityMapping
{
  const std::array<Dune::FieldVector<double, 3>, 3> vertices_;

public:
  IdenityMapping(const std::array<Dune::FieldVector<double, 3>, 3>& vertices) : vertices_(vertices) {}

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
    Dune::FieldVector<ctype,dimworld> y{0,0,-1};
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
    // iterator over all intersections
    for (auto&& is : intersections(gridView, element))
    {
      // index of the neighbour
      int nIdx = mapper.index(is.outside());

      // calculate distance between the midpoints
      auto eCenter = element.geometry().center();
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
      for(auto&& intersection : intersections(leafGridView, element))
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

void oneDimensionalTest ()
{
  typedef Dune::FoamGrid<1, 2> Grid;
  typedef Grid::ctype ctype;
  const int dimgrid = Grid::dimension;
  const int dimworld = Grid::dimensionworld;

  // Start grid creation
  Dune::GridFactory<Grid> factory;

  // The list of grid vertex positions
  int numVertices = 3;
  std::vector<Dune::FieldVector<double, dimworld> > vertices({{0.0, 1.0},
                                                              {-0.5*std::sqrt(3), -0.5},
                                                              {0.5*std::sqrt(3), -0.5}});

  // Create the grid vertices
  for (int i=0; i<numVertices; i++)
    factory.insertVertex(vertices[i]);

  // Create the element geometries
  const int numElements = 3;
  std::vector<std::vector<unsigned int> > cornerIDs = {{0,1},{1,2},{2,0}};

  double angle = M_PI/2;
  for(int i = 0; i < numElements; i++)
  {
    const auto mapping = UnitCircleMapping<ctype, dimgrid, dimworld>{angle, angle+2.0/3.0*M_PI};
    factory.insertElement(Dune::GeometryTypes::line, cornerIDs[i], std::move(mapping));
    angle += 2.0/3.0*M_PI;
  }

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView(), Dune::VTK::nonconforming);
  writer.write("initial");

  // refine the grid
  grid->globalRefine(1);

  // output VTK
  writer.write("refine-1");

  // coarden and then refine four times
  grid->globalRefine(-1);
  grid->globalRefine(4);

  // output VT
  writer.write("refine-4");

  // Solve the heat equation on the refined grid
  // dT/dt + div(lambda*grad(T)) = 0
  // using a finite volume method with an explicit Euler time discretization

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
  double tEnd = 1.0;
  int timestep = 0;

  // write output only every nth timestep
  int episode = 10;

  // the heat conductivity
  double lambda = 1.0;

  // Write pvd header
  Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView(), "temperature_1d", ".", "");
  vtkWriter.addCellData(temperature, "celldata");
  vtkWriter.write(t);

  // do the time integration
  while(t <= tEnd)
  {
    // do adaptation
    int lmax = 6;
    for(int i = 0; i < lmax; i++)
      finitevolumeadapt(*grid, mapper, temperature, 1, lmax);


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
  double tao = (1.0 + std::sqrt(5))/2.0; // golden ratio
  std::vector<Dune::FieldVector<double, dimworld> > vertices({{1,tao,0},{-1,tao,0},{1,-tao,0},{-1,-tao,0},
                                                              {0,1,tao},{0,-1,tao},{0,1,-tao},{0,-1,-tao},
                                                              {tao,0,1},{-tao,0,1},{tao,0,-1},{-tao,0,-1}});

  // Create the grid vertices
  for (size_t i=0; i<vertices.size(); i++)
  {
    // make vertices live on unit sphere
    vertices[i] /= vertices[i].two_norm();
    // insert vertices
    factory.insertVertex(vertices[i]);
  }

  // Create the element geometries
  std::vector<std::vector<unsigned int> > cornerIDs({{1,9,11},{3,9,11},{2,3,5},{2,3,7},
                                                     {0,1,4},{0,1,6},{4,5,8},{4,5,9},
                                                     {0,8,10},{2,8,10},{6,7,10},{6,7,11},
                                                     {4,0,8},{4,1,9},{6,10,0},{6,11,1},
                                                     {5,2,8},{5,3,9},{7,2,10},{7,3,11}});
  for(size_t i = 0; i < cornerIDs.size(); i++)
  {
    const auto mapping = UnitSphereMapping<ctype, dimgrid, dimworld>(
        {vertices[cornerIDs[i][0]],vertices[cornerIDs[i][1]],vertices[cornerIDs[i][2]]});

    factory.insertElement(Dune::GeometryTypes::simplex(2), cornerIDs[i], std::move(mapping));
  }

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView(), Dune::VTK::nonconforming);
  writer.write("sphere_0");

  // refine the grid
  //grid->globalRefine(3);
  writer.write("sphere_1");

  // Solve the heat equation on the refined grid
  // dT/dt + div(lambda*grad(T)) = 0
  // using an h-adaptive finite volume method with an explicit Euler time discretization

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
  double tEnd = 1.0;
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
    int lmax = 4;
    for(int i = 0; i < lmax; i++)
        finitevolumeadapt(*grid, mapper, temperature, 1, lmax);

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
}


int main (int argc, char *argv[]) try
{
  std::cout << std::endl << "Running example for FoamGrid<1, 2>" << std::endl;
  oneDimensionalTest();

  // TODO: This currently fails (but we still compile it): fix adaptivity in 2d (reported in issue #6)
  if (true)
  {
    std::cout << std::endl << "Running example for FoamGrid<2, 3>" << std::endl;
    twoDimensionalTest();
  }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
  std::cout << e << std::endl;
  return 1;
}
