#include <dune/foamgrid/foamgrid.hh>
#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/gmsh4/gridcreators/lagrangegridcreator.hh>

int main()
{
  using namespace Dune;
  using Grid = FoamGrid<2,3>;
  GridFactory<Grid> factory;

  // The creator is responsible for filling a GridFactory from the data
  // found in the file. Additionally, it might provide a grid-function representing
  // a geometry parametrization
  Gmsh4::LagrangeGridCreator creator{factory};
  Gmsh4Reader reader{creator};

  reader.read("filename.msh");
  auto grid = factory.createGrid();
}
