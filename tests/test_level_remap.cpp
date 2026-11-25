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
#include <cmath>
#include <iostream>
#include <map>
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

bool same(double a, double b)
{
    return std::abs(a - b) <= 1e-12 * std::max(std::abs(a), std::abs(b));
}

// Point-wise transfer as Fracture::redistribute_values does it.
std::vector<double>
transfer(const std::vector<double>& v, const Map& map1, const Map& map2, int level, int map1_level)
{
    std::vector<double> sum(map2.size(), 0.0);
    std::vector<double> weight(map2.size(), 0.0);
    for (const auto& [i, j, w] : RegularTrimesh::createGridToGridMap(map1, map2, level, map1_level)) {
        sum[j] += v[i] * w;
        weight[j] += w;
    }
    for (size_t j = 0; j < sum.size(); ++j) {
        sum[j] = weight[j] > 0.0 ? sum[j] / weight[j] : 0.0;
    }
    return sum;
}

Map cellMap(const RegularTrimesh& mesh)
{
    Map m;
    for (const auto& cell : mesh.cellIndices()) {
        m.push_back({cell});
    }
    return m;
}

RegularTrimesh fineMesh()
{
    std::vector<CellRef> cells;
    for (int i = -6; i < 6; ++i) {
        for (int j = -6; j < 6; ++j) {
            cells.push_back({i, j, 0});
            cells.push_back({i, j, 1});
        }
    }
    return RegularTrimesh(cells.begin(), cells.end());
}

// Values tagged by cell index, so a wrong donor shows as a wrong value.
std::vector<double> tagged(size_t n)
{
    std::vector<double> v(n);
    for (size_t i = 0; i < n; ++i) {
        v[i] = 1.0 + static_cast<double>(i);
    }
    return v;
}

// A second round on the same level-1 grid must give back what the first round put there.
void sameCoarseLevelIsIdentity()
{
    const Map coarse = cellMap(fineMesh().coarsen(true));
    const auto v = tagged(coarse.size());
    const auto out = transfer(v, coarse, coarse, 1, 1);
    for (size_t i = 0; i < v.size(); ++i) {
        require(same(out[i], v[i]), "level 1 -> level 1 keeps every value");
    }
    // Negative control: reading the level-1 refs as level 0 coarsens them twice.
    const auto bad = transfer(v, coarse, coarse, 1, 0);
    size_t wrong = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        wrong += !same(bad[i], v[i]);
    }
    require(wrong > 0, "negative control: level-unaware transfer changes values");
    std::cout << "same level: " << coarse.size() << " cells kept, " << wrong
              << " changed without the level\n";
}

// Refining back to level 0 gives each fine cell its coarse parent's value.
void coarseToFineTakesParent()
{
    const RegularTrimesh fine = fineMesh();
    const Map coarse = cellMap(fine.coarsen(true));
    const Map fmap = cellMap(fine);
    std::map<CellRef, double> parent_value;
    const auto v = tagged(coarse.size());
    for (size_t i = 0; i < coarse.size(); ++i) {
        parent_value[coarse[i].front()] = v[i];
    }
    const auto out = transfer(v, coarse, fmap, 0, 1);
    size_t covered = 0;
    for (size_t j = 0; j < fmap.size(); ++j) {
        const auto it = parent_value.find(RegularTrimesh::fine_to_coarse(fmap[j].front(), 1));
        if (it == parent_value.end()) {
            require(out[j] == 0.0, "a fine cell outside the coarse mesh gets no donor");
            continue;
        }
        ++covered;
        require(same(out[j], it->second), "a fine cell takes its parent's value");
    }
    require(covered > 0, "the coarse mesh covers some fine cells");
    // Negative control: without the level the coarse refs are matched as fine refs.
    const auto bad = transfer(v, coarse, fmap, 0, 0);
    size_t wrong = 0;
    for (size_t j = 0; j < fmap.size(); ++j) {
        wrong += !same(bad[j], out[j]);
    }
    require(wrong > 0, "negative control: level-unaware refinement misplaces values");
    std::cout << "coarse to fine: " << covered << " covered fine cells, " << wrong
              << " wrong without the level\n";
}

// A piecewise-constant field survives 0 -> 1 -> 1 -> 0, the sequence of one propagation.
void roundTripWithinPropagation()
{
    const RegularTrimesh fine = fineMesh();
    const Map fmap = cellMap(fine);
    const Map coarse = cellMap(fine.coarsen(true));
    std::map<CellRef, double> value_of_parent;
    for (size_t i = 0; i < coarse.size(); ++i) {
        value_of_parent[coarse[i].front()] = 10.0 + static_cast<double>(i);
    }
    std::vector<double> v0(fmap.size(), 0.0);
    for (size_t j = 0; j < fmap.size(); ++j) {
        const auto it = value_of_parent.find(RegularTrimesh::fine_to_coarse(fmap[j].front(), 1));
        v0[j] = (it == value_of_parent.end()) ? 0.0 : it->second;
    }
    const auto v1 = transfer(v0, fmap, coarse, 1, 0);
    const auto v1b = transfer(v1, coarse, coarse, 1, 1);
    const auto back = transfer(v1b, coarse, fmap, 0, 1);
    for (size_t j = 0; j < fmap.size(); ++j) {
        require(same(back[j], v0[j]), "round trip keeps a field that is constant per coarse cell");
    }
    std::cout << "round trip: " << fmap.size() << " fine cells unchanged\n";
}
} // namespace

int main()
{
    try {
        sameCoarseLevelIsIdentity();
        coarseToFineTakesParent();
        roundTripWithinPropagation();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
