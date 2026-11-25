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
#include <opm/geomech/FracturePieces.hpp>
#include <opm/geomech/RegularTrimesh.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm;
using Neighbours = std::vector<std::vector<int>>;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::vector<int> sorted(std::vector<int> v)
{
    std::sort(v.begin(), v.end());
    return v;
}

// Two chains 0-1-2 and 3-4; only the sink decides which piece is returned.
void sinkDecides()
{
    const Neighbours nb {{1}, {0, 2}, {1}, {4}, {3}};
    const std::vector<double> none(5, 0.0);
    require(sorted(cellsInPiecesWithoutSink(nb, {true, false, false, false, false}, none))
                == std::vector<int> {3, 4},
            "a fed cell keeps its piece out");
    std::vector<double> leak(5, 0.0);
    leak[4] = 1e-20;
    require(sorted(cellsInPiecesWithoutSink(nb, std::vector<bool>(5, false), leak))
                == std::vector<int> {0, 1, 2},
            "any positive leak-off is a sink");
    leak[4] = 0.0;
    leak[2] = -1.0;
    require(cellsInPiecesWithoutSink(nb, std::vector<bool>(5, false), leak).size() == 5,
            "zero or negative leak-off is no sink");
    require(cellsInPiecesWithoutSink(nb, {false, true, false, false, true}, none).empty(),
            "every piece fed: nothing returned");
    require(cellsInPiecesWithoutSink({}, {}, {}).empty(), "empty mesh");
    std::cout << "sink decides: ok\n";
}

// On a trimesh, cells joined to a fed block only at a vertex form a piece without sink.
void vertexJoinedTrimeshPiece()
{
    std::vector<CellRef> cells;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cells.push_back({i, j, 0});
            cells.push_back({i, j, 1});
        }
    }
    cells.push_back({4, 4, 0});
    cells.push_back({4, 4, 1});
    std::map<CellRef, int> index;
    for (size_t c = 0; c < cells.size(); ++c) {
        index[cells[c]] = static_cast<int>(c);
    }
    Neighbours nb(cells.size());
    for (size_t c = 0; c < cells.size(); ++c) {
        for (const auto& n : RegularTrimesh::cellNeighbors(cells[c])) {
            const auto it = index.find(n);
            if (it != index.end()) {
                nb[c].push_back(it->second);
            }
        }
    }
    std::vector<bool> fed(cells.size(), false);
    fed[index.at({1, 1, 0})] = true;
    const auto lone = sorted(cellsInPiecesWithoutSink(nb, fed, std::vector<double>(cells.size(), 0.0)));
    require(lone == sorted({index.at({4, 4, 0}), index.at({4, 4, 1})}),
            "the vertex-joined pair is returned, the fed block is not");
    std::cout << "vertex-joined pair: " << lone.size() << " cells without sink\n";
}
} // namespace

int main()
{
    try {
        sinkDecides();
        vertexJoinedTrimeshPiece();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
