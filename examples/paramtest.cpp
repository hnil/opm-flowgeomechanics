#include <config.h>

#include "opm/geomech/param_interior.hpp"
#include "opm/geomech/GridStretcher.hpp"

#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace Opm;

using Grid = Dune::FoamGrid<2, 3>;

// ----------------------------------------------------------------------------
int simpletest(const string& fname)
// ----------------------------------------------------------------------------
{
  ifstream is(fname);
  if (!is.is_open()) {
    cout << "Failed to open file" << endl;
    return 1;
  }
  
  int num_p;
  is >> num_p;

  vector<double> points(num_p * 2);
  for (int i = 0; i != num_p * 2; ++i)
    is >> points[i];
  
  int num_q;
  is >> num_q;

  vector<double> q(num_q * 2);
  for (int i = 0; i != num_q * 2; ++i)
    is >> q[i];
  
  // -------------------------- compute parametrization --------------------------
  vector<double> result;

  parametrize_interior_2D(points, q, result);

  
  return 0;
};

double dist2D(const Dune::FieldVector<double, 3>& p1,
              const Dune::FieldVector<double, 3>& p2)
{
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  return sqrt(dx*dx + dy*dy);
}

// ----------------------------------------------------------------------------
vector<int> boundary_node_indices(const Grid& g)
// ----------------------------------------------------------------------------
{
  //set<int> bnodes;
  auto view = g.leafGridView();

  //using FVec = Dune::FieldVector<double, 3>;
  using CoordType = Dune::FieldVector<double, 3>;

  // register all node coordinates
  vector<CoordType> vcoords;
  for (auto v : vertices(view)) 
    vcoords.push_back(v.geometry().corner(0));

  // determine boundary edges
  vector<int> count(view.size(1), 0); // number of edges (codim 1)

  for (const auto& elem : elements(view))
    for (int i = 0; i != elem.subEntities(1); ++i)
      count[view.indexSet().index(elem.subEntity<1>(i))] += 1;

  // @@ This is silly - there must be a better way to identify boundary nodes
  // than by geometric comparison.
  int ix = 0;
  set<int> bix;
  const double TOL = 1e-3;
  for (auto ep = view.begin<1>(); ep != view.end<1>(); ++ep, ++ix)
    if (count[ix] < 2)
      for (int i = 0; i != 2; ++i) {
        auto pt = ep->geometry().corner(i);
        bool found = false;
        for (int j = 0; j != vcoords.size(); ++j && !found)
          if (dist2D(pt, vcoords[j]) < TOL) {
            found = true;
            bix.insert(j);
          }
      }

  return vector(bix.begin(), bix.end());
  
  // for (auto it = bix.begin(); it != bix.end(); ++it)
  //   cout << *it << endl;

  // for (int i = 0; i != boundary_ind.size(); ++i)
  //   cout << boundary_ind[i] << endl;

  // for (int i = 0; i != count.size(); ++i)
  //   cout << count[i] << endl;
      
  
  // vector<int> count(view.size(2), 0); // number of vertices (codim 2)
  
  // for (const auto& elem : elements(view))
  //   if (elem.hasBoundaryInter
    
  //   for (int i = 0; i != elem.subEntities(2); ++i)
  //     count[view.indexSet().index(elem.subEntity<2>(i))] += 1;

  
    // cout << elem.subEntities(2) << endl;
    // cout << view.indexSet().index(elem.subEntity<2>(0)) << endl;
    // cout << view.indexSet().index(elem.subEntity<2>(0)) << endl;
    // cout << view.indexSet().index(elem.subEntity<2>(0)) << endl;
    //cout << elem.subEntity<2>(1);
    //cout << elem.subEntity<2>(2);
  //}
  
  // for (const auto& elem : elements(view)) {
  //   for (const auto& is : Dune::intersections(view, elem)) { 
  //     if (is.boundary()) {
  //       const auto c1 = is.geometry().corner(0);
  //       const auto c2 = is.geometry().corner(1);

  //       //view.indexSet().index(c1);
  //       //bnodes.insert(view.indexSet().index(c1));
  //       //bnodes.insert(view.indexSet().index(c2));
  //     }
  //   }
  // }
  // return vector<int>(bnodes.begin(), bnodes.end());
  return vector<int>();
}
  
// ----------------------------------------------------------------------------
int meshtest(const string& fname)
// ----------------------------------------------------------------------------
{
  std::unique_ptr grid = Dune::GmshReader<Grid>::read(fname);

  
  // extract coordinates

  auto view = grid->leafGridView();
  
  for (int i = 0; i != 3; ++i)
    cout << view.size(i) << endl;

  //view.begin<0>();
  //for (auto it = view.begin<0>(); it != view.end<0>(); ++it);

  // extract all node coordinates (2D only)
  vector<double> coords;
  for (const auto& vertex : vertices(view)) {
    coords.push_back(vertex.geometry().corner(0)[0]);
    coords.push_back(vertex.geometry().corner(0)[1]);
  }
  
  // identify boundary nodes
  const auto bindices = boundary_node_indices(*grid);

  // identifying interior nodes
  vector<int> tmp(coords.size()/2, 0);
  for (auto b : bindices)
    tmp[b] = 1;
  vector<int> iindices;
  for (int ix = 0; ix != tmp.size(); ++ix)
    if (tmp[ix] == 0)
      iindices.push_back(ix);

  // sort node coordinates into boundary and interior
  vector<double> bcoords, icoords;
  for (int i = 0; i != tmp.size(); ++i) {
    if (tmp[i] == 0) {
      // interior
      icoords.push_back(coords[2*i]);
      icoords.push_back(coords[2*i+1]);
    } else {
      // boundary
      bcoords.push_back(coords[2*i]);
      bcoords.push_back(coords[2*i+1]);
    }
  }
  
  // call parametrization function
  vector<double> result;
  parametrize_interior_2D(bcoords, icoords, result);

  cout << "Num boundary vertices: " << bcoords.size() / 2 << endl;
  cout << "Num interior vertices: " << icoords.size() / 2 << endl;
  cout << "Num params: " << result.size() << endl;

  
  for (int i = 0; i != result.size(); ++i)
    cout << result[i] << endl;


  
  // write result to disk
  auto vtkwriter =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid->leafGridView(),
                                                          Dune::VTK::nonconforming);
  vtkwriter->write("undeformed");


  // modify grid coordinates
  const int Nb = bcoords.size()/2; // num boundary points
  for (int i = 0; i != Nb; ++i) {
    bcoords[2*i  ] += 3; // translate in x direction
    if (i < Nb/3)
      bcoords[2*i+1] *= 2; // stretch in y direction
  }

  const int Ni = icoords.size()/2; // num interior points
  fill(icoords.begin(), icoords.end(), 0); // reset icoords
  for (int i = 0; i != Ni; ++i)
    for (int b = 0; b != Nb; ++b)
      for (int c = 0; c != 2; ++c)
        icoords[2*i + c] += bcoords[2*b + c] * result[Nb * i + b];

  const int N = coords.size() / 2;
  int bc = 0;
  int ic = 0;
  for (int i = 0; i != N; ++i) {
    for (int c = 0; c != 2; ++c) {
      coords[2*i + c] = (tmp[i] == 0) ? icoords[2*ic + c] : bcoords[2*bc+c];
    }
    if (tmp[i] == 0)
      ++ic;
    else
      ++bc;
  }
  
  int vcount = 0;
  for (auto& vertex : vertices(view)) 
    grid->setPosition(vertex,
                      Dune::FieldVector<double, 3> {coords[2*vcount],
                                                    coords[2*vcount++ + 1],
                                                    0});
    //    const VertexHandle& vh = vertex.impl().hostEntity();
  //    vertex.geometry().corner(0)[0] = 0;//coords[2*vcount];
  //    vertex.geometry().corner(0)[1] = 0;//coords[2*vcount + 1];
  //   vcount++;
  

  for (auto& vertex : vertices(view)) {
    cout << vertex.geometry().corner(0)[0] << ", ";
    cout << vertex.geometry().corner(0)[1] << endl;
  }

  
  vtkwriter->write("deformed");  

  // repeat deformation, this time using a grid stretcher
  std::unique_ptr grid2 = Dune::GmshReader<Grid>::read(fname);

  GridStretcher gs(*grid2);

  vector<double> coords2;
  for (const auto& vertex : vertices(grid2->leafGridView())) {
    coords2.push_back(vertex.geometry().corner(0)[0]);
    coords2.push_back(vertex.geometry().corner(0)[1]);
  }

  
  std::vector<GridStretcher::CoordType> disp(Nb, {0, 0, 0});
  for (int i = 0; i != Nb; ++i) {
    disp[i][0] += 6; // translate in x direction
    if (i < Nb/3)
      disp[i][1] = coords2[2*i+1];//bcoords[2*i+1]/2; // stretch in y direction
  }
  
  gs.applyBoundaryNodeDisplacements(disp);
  
  auto vtkwriter2 =
    std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid2->leafGridView(),
                                                          Dune::VTK::nonconforming);
  vtkwriter2->write("deformed2");

  return 0;
}

// ----------------------------------------------------------------------------
int projectiontest(const string& fname, const string& ftarget)
  // ----------------------------------------------------------------------------
{
  // open file and read 3D points to a std::vector<double>
  ifstream is(fname);
  if (!is.is_open()) {
    cout << "Failed to open file" << endl;
    return 1;
  }
  vector<double> points;
  while (!is.eof()) {
    double x, y, z;
    is >> x >> y >> z;
    if (is.eof())
      break;
    points.push_back(x);
    points.push_back(y);
    points.push_back(z);
  }

  // call projection function
  const int N = points.size() / 3;
  vector<double> result(N * 2);
  const Axis3D axis = project_to_2D(points, result);

  // write axis, then result to file
  ofstream os(ftarget);
  os.precision(16);
  for (int i = 0; i != 3; ++i) {
    for (int dim = 0; dim != 3; ++dim)  {
      os << axis[i][dim] << " ";
    }
    os << endl;
  }
  for (int i = 0; i != N; ++i) {
    os << result[2*i] << " " << result[2*i+1] << endl;
  }
  
  return 0;
}

//----------------------------------------------------------------------------
int lifttest(const string& fname, const string& target)
//----------------------------------------------------------------------------
{
  // reading axis and 2D points from file
  ifstream is(fname);
  if (!is.is_open()) {
    cout << "Failed to open file" << endl;
    return 1;
  }
  Axis3D axis;
  for (int i = 0; i != 3; ++i)
    for (int dim = 0; dim != 3; ++dim)
      is >> axis[i][dim];

  vector<double> points;
  while (!is.eof()) {
    double x, y;
    is >> x >> y;
    if (is.eof())
      break;
    points.push_back(x);
    points.push_back(y);
  }

  // lift the points to 3D, using the axis
  std::vector<double> result;
  lift_to_3D(points, axis, result);

  // write result to file, using high precision
  ofstream os(target);
  os.precision(16);
  for (int i = 0; i != result.size() / 3; ++i) {
    os << result[3*i] << " " << result[3*i+1] << " " << result[3*i+2] << endl;
  }
        
  return 0;
}

//----------------------------------------------------------------------------
int loop_repar_test(const string& fname, const string& target)
// ----------------------------------------------------------------------------
{
  // reading 2D points from dfile
  ifstream is(fname);
  if (!is.is_open()) {
    cout << "Failed to open file" << endl;
    return 1;
  }
  vector<double> points;
  while (!is.eof()) {
    double x, y;
    is >> x >> y;
    if (is.eof())
      break;
    points.push_back(x);
    points.push_back(y);
  }
  
  // redistribute points
  vector<double> result;
  redistribute_2D(points, result);

  // write result to file
  ofstream os(target);
  for (int i = 0; i != result.size() / 2; ++i) {
    os << result[2*i] << " " << result[2*i+1] << endl;
  }
          
  return 0;
}

// ============================================================================
int main(int varnum, char** vararg)
// ============================================================================
{
  string fname(vararg[1]);
  string ftype(fname.substr(fname.size()-3));

  if (ftype == "txt")
    return simpletest(fname);
  else if (ftype == "msh")
    return meshtest(fname);
  else if (ftype == "p3d")
    return projectiontest(fname, string(vararg[2]));
  else if (ftype == "p2d")
    return lifttest(fname, string(vararg[2]));
  else if (ftype == "lop")
    return loop_repar_test(fname, string(vararg[2]));

  return 0;
}


