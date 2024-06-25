#include "coupledsolver.hpp"

#include <opm/common/ErrorMacros.hpp>

namespace Opm
{

void solve_fully_coupled(ResVector& pvec, // output: fracture pressure
                         ResVector& hvec, // output: aperture
                         const SMatrix& pmat, // pressure matrix (sparse)
                         const FMatrix& amat, // aperture matrix (full)
                         const std::vector<Htrans>& htrans,
                         const double rate,
                         const std::vector<size_t> ratecells,
                         const double bhp,
                         const std::vector<size_t> bhpcells)
{
  // We consider the nonlinear system:
  // | A          -I |  | h |   | 0 |   (mechanics equation)
  // |               |  |   | = |   |
  // | C(h, p)  M(h) |  | p |   | q |   ( flow equation)


  
  // @@ Implement me
  OPM_THROW(std::runtime_error,"if: fully implicit not implemented");
}
  
}; // end namespace Opm
