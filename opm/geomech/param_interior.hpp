#pragma once

#include <vector>

namespace Opm{

  void parametrize_interior_2D(const std::vector<double>& bpoints,
                               const std::vector<double>& ipoints,
                               std::vector<double>& result);
  
}; // end namespace Opm
