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
#include <opm/geomech/MechanicalConstraintValidation.hpp>
#include <opm/grid/CpGrid.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool value, const char* message)
{
    if (!value) {
        throw std::runtime_error(message);
    }
}

void localChecks()
{
    using namespace Opm::Elasticity;
    const std::vector<std::pair<std::int64_t, std::int64_t>> edges{{0, 1}, {1, 2}};
    require(rigidBodyConstraintRanks(edges, {})[0].rank == 0, "Free body must have rank zero");
    require(rigidBodyConstraintRanks(edges, {{0, {0, 0, 0}, 7}})[0].rank == 3,
            "A fixed point leaves three rotations");
    require(rigidBodyConstraintRanks(edges, {{0, {0, 0, 0}, 7}, {1, {1, 0, 0}, 7}})[0].rank == 5,
            "Two fixed points leave rotation around their line");
    const std::vector<ConstrainedMechanicalVertex> supports{
        {0, {0, 0, 0}, 7}, {1, {1, 0, 0}, 6}, {2, {0, 1, 0}, 4}};
    require(rigidBodyConstraintRanks(edges, supports)[0].rank == 6,
            "Six independent scalar supports must constrain a body");
    auto disconnected = edges;
    disconnected.emplace_back(10, 11);
    const auto result = rigidBodyConstraintRanks(disconnected, supports);
    require(result.size() == 2 && result[0].rank == 6 && result[1].rank == 0,
            "A constrained body cannot hide a disconnected free body");
    auto duplicates = supports;
    duplicates[0].mask = 1;
    duplicates.push_back({0, {0, 0, 0}, 6});
    duplicates.push_back(supports[1]);
    require(rigidBodyConstraintRanks(edges, duplicates)[0].rank == 6,
            "Overlap copies must merge constraint masks by global vertex ID");
    for (double scale : {1e-6, 1.0, 1e6}) {
        auto shifted = supports;
        for (auto& node : shifted) {
            for (double& x : node.position) {
                x = scale*(x + 1e6);
            }
        }
        require(rigidBodyConstraintRanks(edges, shifted)[0].rank == 6,
                "Rank must be independent of origin and physical units");
    }
    require(rigidBodyConstraintRanks({}, {}).empty(), "Empty topology has no components");
}
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    Dune::CpGrid grid;
    std::string error;
    try {
        localChecks();
    } catch (const std::exception& e) {
        error = e.what();
    }
    if (grid.comm().max(static_cast<int>(!error.empty()))) {
        std::cerr << error << '\n';
        return 1;
    }
    try {
        grid.createCartesian({8, 2, 2}, {1, 1, 1});
        if (grid.comm().size() > 1) {
            std::vector<int> parts(32);
            for (int k = 0; k < 2; ++k) {
                for (int j = 0; j < 2; ++j) {
                    for (int i = 0; i < 8; ++i) {
                        parts[i + 8*j + 16*k] = i*grid.comm().size()/8;
                    }
                }
            }
            grid.loadBalance(parts, false, true, 1);
        }
        const auto gv = grid.leafGridView();
        using Nodes = std::vector<std::tuple<std::size_t, Opm::MechBCValue>>;
        Nodes supports, line, rollers;
        for (const auto& vertex : vertices(gv)) {
            const auto index = gv.indexSet().index(vertex);
            const auto p = vertex.geometry().center();
            Opm::MechBCValue fixed;
            fixed.fixeddir = {true, true, true};
            if (grid.comm().rank() == 0 && p[0] == 0) {
                if ((p[1] == 0 && p[2] == 0) || (p[1] == 1 && p[2] == 0)
                    || (p[1] == 0 && p[2] == 1)) {
                    supports.emplace_back(index, fixed);
                }
                if (p[1] == 0) {
                    line.emplace_back(index, fixed);
                }
            }
            Opm::MechBCValue roller;
            roller.fixeddir = {p[0] == 0, p[1] == 0, p[2] == 0};
            if (roller.fixeddir[0] || roller.fixeddir[1] || roller.fixeddir[2]) {
                rollers.emplace_back(index, roller);
            }
        }
        auto check = [&](const Nodes& nodes, int rank, const char* message) {
            const auto result = Opm::Elasticity::mechanicalConstraintRanks(grid, nodes);
            require(result.size() == 1 && result[0].rank == rank, message);
        };
        check(supports, 6, "Globally valid supports must pass even with no local supports on other ranks");
        check(line, 5, "Collinear supports must be rejected consistently on all ranks");
        check({}, 0, "Unconstrained distributed body must be rejected on every rank");
        check(rollers, 6, "Distributed rollers and overlap copies must constrain all modes");
        Nodes invalid;
        if (grid.comm().rank() == 0) {
            invalid.emplace_back(gv.size(3), Opm::MechBCValue{});
        }
        bool threw = false;
        try {
            Opm::Elasticity::mechanicalConstraintRanks(grid, invalid);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        require(threw, "Invalid local BC indices must fail collectively without hanging");
    } catch (const std::exception& e) {
        error = e.what();
    }
    const int failed = grid.comm().max(static_cast<int>(!error.empty()));
    if (!error.empty()) {
        std::cerr << error << '\n';
    }
    return failed;
}
