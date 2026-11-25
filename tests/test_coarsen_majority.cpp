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

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm;
using Map = std::vector<std::vector<CellRef>>;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// Coarse cells (level 1) with 1, 2, 3 or 4 of their children present at level 0.
std::map<CellRef, int> childCounts()
{
    std::map<CellRef, int> n;
    int k = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int t = 0; t < 2; ++t) {
                n[{3 * i, 3 * j, t}] = 1 + (k++ % 4);
            }
        }
    }
    return n;
}

RegularTrimesh partialMesh(const std::map<CellRef, int>& counts)
{
    std::vector<CellRef> cells;
    for (const auto& [parent, k] : counts) {
        const auto children = RegularTrimesh::coarse_to_fine(parent, 1);
        cells.insert(cells.end(), children.begin(), children.begin() + k);
    }
    return RegularTrimesh(cells.begin(), cells.end());
}

// A parent exists iff at least two of its four children do.
void parentNeedsTwoChildren()
{
    const auto counts = childCounts();
    const auto coarse = partialMesh(counts).coarsen(false, true).cellIndices();
    const std::set<CellRef> present(coarse.begin(), coarse.end());
    for (const auto& [parent, k] : counts) {
        require(present.count(parent) == (k >= 2 ? 1u : 0u), "majority keeps a parent iff >= 2 children");
    }
    require(present.size() == 12, "12 of 16 parents have two or more children");
    // Negative control: the default rule keeps a parent only through its anchor child.
    const auto plain = partialMesh(counts).coarsen(false, false).cellIndices();
    require(plain != coarse, "negative control: anchor rule gives another coarse mesh");
    std::cout << "majority: " << present.size() << " of " << counts.size() << " parents kept, default rule "
              << plain.size() << "\n";
}

// Volume moved onto the majority mesh sums to the fine total, dropped parents included.
void remapConservesVolume()
{
    const auto counts = childCounts();
    const RegularTrimesh fine = partialMesh(counts);
    const RegularTrimesh coarse = fine.coarsen(false, true);
    Map fmap, cmap;
    RegularTrimesh::CellVolumes src;
    double total = 0.0;
    for (const auto& c : fine.cellIndices()) {
        fmap.push_back({c});
        const double v = 1.0 + 0.1 * static_cast<double>(fmap.size());
        src.volume.push_back(v);
        src.group.push_back(0);
        src.centre.push_back({0.0, 0.0, 0.0});
        total += v;
    }
    std::vector<double> area;
    std::vector<int> group;
    std::vector<Coord3D> centre;
    for (const auto& c : coarse.cellIndices()) {
        cmap.push_back({c});
        area.push_back(4.0);
        group.push_back(0);
        centre.push_back({0.0, 0.0, 0.0});
    }
    const auto h = RegularTrimesh::conservativeRemap(src, fmap, cmap, 1, area, group, centre, 0);
    double moved = 0.0;
    for (size_t j = 0; j < h.size(); ++j) {
        moved += h[j] * area[j];
    }
    require(std::abs(moved - total) <= 1e-12 * total, "area * thickness sums to the fine volume");
    std::cout << "conservative remap: " << moved << " of " << total << "\n";
}
} // namespace

int main()
{
    try {
        parentNeedsTwoChildren();
        remapConservesVolume();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
