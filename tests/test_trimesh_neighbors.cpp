/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/geomech/RegularTrimesh.hpp>

#include <algorithm>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace Opm;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

int sharedNodes(const CellRef& a, const CellRef& b)
{
    const auto na = RegularTrimesh::cellNodes(a);
    const auto nb = RegularTrimesh::cellNodes(b);
    int n = 0;
    for (const auto& x : na) {
        n += static_cast<int>(std::count(nb.begin(), nb.end(), x));
    }
    return n;
}
} // namespace

// removeCells and the front veto rely on cellNeighbors giving the three edge neighbours.
int main()
{
    try {
        int checked = 0;
        for (int i = -4; i < 4; ++i) {
            for (int j = -4; j < 4; ++j) {
                for (int k = 0; k < 2; ++k) {
                    const CellRef c {i, j, k};
                    const auto nb = RegularTrimesh::cellNeighbors(c);
                    require(std::set<CellRef>(nb.begin(), nb.end()).size() == 3,
                            "three distinct neighbours");
                    for (const auto& n : nb) {
                        require(n[2] != c[2], "a neighbour has the other orientation");
                        require(sharedNodes(c, n) == 2, "a neighbour shares one edge");
                        const auto back = RegularTrimesh::cellNeighbors(n);
                        require(std::find(back.begin(), back.end(), c) != back.end(),
                                "the neighbour relation is symmetric");
                    }
                    ++checked;
                }
            }
        }
        std::cout << "neighbours: " << checked << " cells checked\n";
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
