#include <config.h>
#include <opm/geomech/FractureConnectionId.hpp>
#include <opm/geomech/FractureWell.hpp>
#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    try {
        Opm::WellConnections connections;
        for (int grid : {0, 1, 2}) {
            connections.addConnection(0, 0, 0, 6562, Opm::Connection::State::OPEN,
                0.0, Opm::Connection::CTFProperties{}, 1,
                Opm::Connection::Direction::Z, Opm::Connection::CTFKind::DeckValue,
                grid, grid);
        }
        for (int grid : {0, 1, 2}) {
            const auto* connection = Opm::FractureConnectionId{6562, grid}.find(connections);
            require(connection && connection->get_lgr_level() == grid,
                    "Equal numeric indices in different grids must remain distinct");
        }
        require(!Opm::FractureConnectionId{6562, 3}.find(connections),
                "Missing grid must not fall back to another grid");
        require(!Opm::FractureConnectionId{1875, 0}.find(connections),
                "WSEED parent index must not resolve as an LGR completion");

        // Only two completions are local, stored in a different order from
        // the schedule. Schedule ordinals cannot index this list.
        Opm::FractureWell::Connection first, second;
        first.cell = 42;
        first.identity = {6562, 2};
        second.cell = 17;
        second.identity = {6562, 1};
        const Opm::FractureWell well(".", "identity_test", "TEST", {first, second},
            {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}}, {{0, 1}, {1, 2}});
        require(well.reservoirCellForConnection({6562, 1}) == 17,
                "Local mapping must use identity rather than schedule order");
        require(well.reservoirCellForConnection({6562, 2}) == 42,
                "Local mapping must preserve the grid identity");
        require(well.reservoirCellForConnection({6562, 0}) == -1,
                "A nonlocal completion must not index the local list");
        require(well.connectionIdentity(1) == Opm::FractureConnectionId{6562, 1},
                "Fracture creation must retain the selected completion identity");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
