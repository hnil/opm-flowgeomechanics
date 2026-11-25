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
#include <opm/geomech/MechSector.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

using namespace Opm::MechSector;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// One cell whose I = 0 edge pillars point towards each other in J, so shifting
// each along its own direction makes them cross.
Grdecl convergingBox()
{
    Grdecl g;
    g.dims = {1, 1, 1};
    const double xy[4][2] = {{0.0, 0.0}, {100.0, -50.0}, {0.0, 100.0}, {100.0, 150.0}};
    for (const auto& p : xy) {
        g.coord.insert(g.coord.end(), {p[0], p[1], 1000.0, p[0], p[1], 1100.0});
    }
    g.zcorn.resize(8);
    for (int c = 0; c < 8; ++c) {
        g.zcorn[zIndex(g.dims, 0, 0, 0, c & 1, (c >> 1) & 1, c >> 2)] = (c >> 2) ? 1100.0 : 1000.0;
    }
    return g;
}

double topArea(const Grdecl& g, int i)
{
    const int px = g.dims[0] + 1;
    const int ring[4] = {i, i + 1, i + 1 + px, i + px};
    double a = 0.0;
    for (int k = 0; k < 4; ++k) {
        const double* p = &g.coord[6 * ring[k]];
        const double* q = &g.coord[6 * ring[(k + 1) % 4]];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

int invertedColumns(const Grdecl& g)
{
    const double ref = topArea(g, g.dims[0] / 2);
    int n = 0;
    for (int i = 0; i < g.dims[0]; ++i) {
        n += topArea(g, i) * ref <= 0.0;
    }
    return n;
}

Padding lateral(bool side)
{
    Padding p;
    p.cells = 3;
    p.size = 600.0;
    p.sideDirection = side;
    return p;
}

// Every padding pillar on a side is its edge pillar moved by distance times one unit vector.
void sideDirectionTranslatesEdge()
{
    const Grdecl in = convergingBox();
    const Grdecl out = padLateral(in, 0, lateral(true), 1000.0, 1100.0);
    const int px = out.dims[0] + 1;
    const double dist[4] = {0.0, 200.0, 400.0, 600.0};
    for (int j = 0; j <= 1; ++j) {
        for (int m = 1; m <= 3; ++m) {
            for (int side = 0; side < 2; ++side) {
                const int q = (side == 0 ? 3 - m : 4 + m) + px * j;
                const double* e = &in.coord[6 * (side + 2 * j)];
                const double* o = &out.coord[6 * q];
                const double ux = side == 0 ? -1.0 : 1.0;
                for (int t = 0; t < 2; ++t) {
                    require(std::abs(o[3 * t] - (e[3 * t] + dist[m] * ux)) < 1e-9
                                && std::abs(o[3 * t + 1] - e[3 * t + 1]) < 1e-9
                                && o[3 * t + 2] == e[3 * t + 2],
                            "side_direction shifts by the mean outward direction (here +-x)");
                }
            }
        }
    }
    require(invertedColumns(out) == 0, "side_direction keeps every padding column oriented");
    std::cout << "side direction: " << out.dims[0] << " columns, none inverted\n";
}

// Negative control: the default per-pillar direction crosses the converging pillars.
void perPillarDirectionCrosses()
{
    const Grdecl out = padLateral(convergingBox(), 0, lateral(false), 1000.0, 1100.0);
    const int bad = invertedColumns(out);
    require(bad > 0, "negative control: per-pillar shift inverts columns on converging pillars");
    std::cout << "per-pillar direction: " << bad << " inverted columns\n";
}

void specReadsSideDirection()
{
    const auto file = std::filesystem::temp_directory_path()
        / ("test_mech_side_direction_" + std::to_string(::getpid()) + ".json");
    {
        std::ofstream f(file);
        f << R"({"lateral": {"cells": 3, "width": 600, "side_direction": true}})";
    }
    const Spec on = readSpec(file.string(), {1, 1, 1});
    {
        std::ofstream f(file);
        f << R"({"lateral": {"cells": 3, "width": 600}})";
    }
    const Spec off = readSpec(file.string(), {1, 1, 1});
    std::filesystem::remove(file);
    require(on.lateral.sideDirection && !off.lateral.sideDirection,
            "lateral.side_direction is read, default off");
    std::cout << "spec: side_direction read, default off\n";
}
} // namespace

int main()
{
    try {
        sideDirectionTranslatesEdge();
        perPillarDirectionCrosses();
        specReadsSideDirection();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
