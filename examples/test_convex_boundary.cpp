#include "convex_boundary.hpp"
#include <iostream>
#include <iterator>
#include <random>

using namespace Opm;
using namespace std;

int main(int argnum, char** argstr)
  
{
  const int choice = atoi(argstr[1]);
  vector<double> pts;
  if (choice == -1) {
  pts = vector<double> { 2, 1,
                         4, 1,
                         5, 3,
                         3, 2,
                         1, 3};
  } else if (choice == 0) {
    pts = vector<double> { 2, 1,
                           4, 2,
                           3, 5,
                           6, 5,
                           4, 7,
                           2, 8,
                           7, 8,
                           9, 5,
                           8, 2,
                           6, 1};
  } else {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 5.0);
    const int N = choice;
    for (int i = 0; i < N*2; ++i) 
      pts.push_back(dist(gen));
  }
  copy (pts.begin(), pts.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;
  const vector<size_t> hull_ixs = convex_hull(pts);
  
  copy(hull_ixs.begin(), hull_ixs.end(), ostream_iterator<size_t>(cout, " "));

  cout << endl;


  cout << endl;
};
