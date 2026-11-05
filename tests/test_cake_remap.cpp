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
#include <numeric>
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

bool close(double a, double b)
{
    return std::abs(a - b) <= 1e-12 * std::max({1.0, std::abs(a), std::abs(b)});
}

double total(const std::vector<double>& h, const std::vector<double>& area)
{
    return std::inner_product(h.begin(), h.end(), area.begin(), 0.0);
}

// The non-conservative overlap average, as the reference to improve on.
std::vector<double> overlapAverage(const std::vector<double>& h, const Map& map1, const Map& map2, int level)
{
    std::vector<double> sum(map2.size(), 0.0);
    std::vector<double> weight(map2.size(), 0.0);
    for (const auto& [i, j, w] : RegularTrimesh::createGridToGridMap(map1, map2, level)) {
        sum[j] += h[i] * w;
        weight[j] += w;
    }
    for (size_t j = 0; j < sum.size(); ++j) {
        sum[j] = weight[j] > 0.0 ? sum[j] / weight[j] : 0.0;
    }
    return sum;
}

RegularTrimesh::CellVolumes
source(const std::vector<double>& h, double area, const std::vector<int>& group, const std::vector<Coord3D>& centre)
{
    RegularTrimesh::CellVolumes s;
    for (double x : h) {
        s.volume.push_back(x * area);
    }
    s.group = group;
    s.centre = centre;
    return s;
}

const CellRef c0 {0, 0, 0}, c1 {1, 0, 0}, c2 {0, 1, 0}, c3 {0, 0, 1};

void fullCoarsening()
{
    const std::vector<double> h {1, 2, 3, 6};
    const auto s = source(h, 1.0, {7, 7, 7, 7}, std::vector<Coord3D>(4));
    const auto out = RegularTrimesh::conservativeRemap(s, {{c0}, {c1}, {c2}, {c3}}, {{c0}}, 1, {4.0}, {7}, {Coord3D {}});
    require(close(out[0], 3.0), "full coarsening gives the mean");
}

void rimInSameGroup()
{
    // The fifth cell's parent is not in the target mesh.
    const std::vector<double> h {1, 1, 1, 1, 8};
    const Map map1 {{c0}, {c1}, {c2}, {c3}, {CellRef {2, 0, 0}}};
    const Map map2 {{c0}};
    const auto s = source(h, 1.0, {7, 7, 7, 7, 7}, std::vector<Coord3D>(5));
    const auto out = RegularTrimesh::conservativeRemap(s, map1, map2, 1, {4.0}, {7}, {Coord3D {}});
    require(close(total(out, {4.0}), 12.0), "rim cell's volume is kept");

    const auto old = overlapAverage(h, map1, map2, 1);
    require(total(old, {4.0}) < 12.0 - 7.0, "negative control: the overlap average drops the rim");
}

void rimToNearestWhenGroupMissing()
{
    const auto s = source({5.0}, 1.0, {3}, {Coord3D {10, 0, 0}});
    const Map map1 {{CellRef {8, 0, 0}}};
    const Map map2 {{c0}, {c1}};
    const auto out = RegularTrimesh::conservativeRemap(
        s, map1, map2, 0, {2.0, 2.0}, {1, 2}, {Coord3D {0, 0, 0}, Coord3D {9, 0, 0}});
    require(close(out[0], 0.0) && close(out[1], 2.5), "unmapped volume goes to the nearest cell");
}

void newCellsStartEmpty()
{
    const auto s = source({2.0}, 1.0, {1}, {Coord3D {}});
    const auto out = RegularTrimesh::conservativeRemap(
        s, {{c0}}, {{c0}, {c1}}, 0, {1.0, 1.0}, {1, 1}, {Coord3D {}, Coord3D {}});
    require(close(out[0], 2.0) && out[1] == 0.0, "a cell with no overlap stays at zero");
}

void partlyCoveredCellIsDiluted()
{
    const std::vector<double> h {4, 4, 4};
    const Map map1 {{c0}, {c1}, {c2}};
    const Map map2 {{c0}};
    const auto s = source(h, 1.0, {1, 1, 1}, std::vector<Coord3D>(3));
    const auto out = RegularTrimesh::conservativeRemap(s, map1, map2, 1, {4.0}, {1}, {Coord3D {}});
    require(close(out[0], 3.0), "three of four children: volume over the whole area");
    require(close(overlapAverage(h, map1, map2, 1)[0], 4.0), "negative control: the overlap average gains volume");
}

void coarseSourceSplit()
{
    // A cell made of four fine references maps onto those four as separate cells.
    const auto s = source({2.0}, 4.0, {1}, {Coord3D {}});
    const Map map1 {{c0, c1, c2, c3}};
    const Map map2 {{c0}, {c1}, {c2}, {c3}};
    const std::vector<double> area(4, 1.0);
    const auto out = RegularTrimesh::conservativeRemap(
        s, map1, map2, 0, area, {1, 1, 1, 1}, std::vector<Coord3D>(4));
    for (double x : out) {
        require(close(x, 2.0), "uniform thickness stays uniform");
    }
}

void strictCoarseningOfARealMesh()
{
    std::vector<CellRef> cells;
    for (int i = -5; i < 6; ++i) {
        for (int j = -5; j < 6; ++j) {
            cells.push_back({i, j, 0});
            cells.push_back({i, j, 1});
        }
    }
    const RegularTrimesh fine(cells.begin(), cells.end());
    const RegularTrimesh coarse = fine.coarsen(true);

    Map map1;
    for (const auto& cell : fine.cellIndices()) {
        map1.push_back({cell});
    }
    Map map2;
    for (const auto& cell : coarse.cellIndices()) {
        map2.push_back({cell});
    }
    const size_t n1 = map1.size();
    const size_t n2 = map2.size();
    const std::vector<double> h(n1, 1.0);
    const auto s = source(h, 1.0, std::vector<int>(n1, 0), std::vector<Coord3D>(n1));
    const std::vector<double> area(n2, 4.0);
    const auto out = RegularTrimesh::conservativeRemap(
        s, map1, map2, 1, area, std::vector<int>(n2, 0), std::vector<Coord3D>(n2));

    const double volume = static_cast<double>(n1);
    require(close(total(out, area), volume), "real strict coarsening keeps the volume");
    for (double x : out) {
        require(x >= 0.0, "no negative thickness");
    }
    const auto g2g = RegularTrimesh::createGridToGridMap(map1, map2, 1);
    double mapped = 0.0;
    for (const auto& e : g2g) {
        mapped += std::get<2>(e);
    }
    require(mapped < volume, "negative control: strict coarsening leaves cells unmapped");
    std::cout << "strict coarsening: " << n1 << " -> " << n2 << " cells, mapped " << mapped << "\n";
}
} // namespace

int main()
{
    try {
        fullCoarsening();
        rimInSameGroup();
        rimToNearestWhenGroupMissing();
        newCellsStartEmpty();
        partlyCoveredCellIsDiluted();
        coarseSourceSplit();
        strictCoarseningOfARealMesh();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
