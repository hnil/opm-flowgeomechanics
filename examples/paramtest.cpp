#include "opm/geomech/param_interior.hpp"

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Opm;

int main(int varnum, char** vararg) {

  // --------------------------- read points from file ---------------------------
  
  ifstream is(vararg[1]);
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
  
  // const vector<double> points = {
  //   0, 0,
  //   1, 0,
  //   1, 1,
  //   0, 1
  // };

  // const vector q = {0.25, 0.5};
  
  // -------------------------- compute parametrization --------------------------
  vector<double> result;

  parametrize_interior_2D(points, q, result);

  for (auto it = result.begin(); it != result.end(); ++it)
    cout << *it << "\n";
  
  return 0;
};
  
