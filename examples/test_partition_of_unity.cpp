//==============================================================================
//!
//! \file test_partition_of_unity.cpp
//!
//! \brief Standalone test that the node (vertex) degrees of freedom form a
//!        partition of unity over the MPI ranks after load balancing a deck.
//!
//! For a partition of unity every shared node must be owned by exactly one
//! rank. Opm::makeEntityEntityCommunication<3>() builds the owner/copy parallel
//! index set for the vertices (codim 3) by selecting a unique owner per shared
//! node and then verifies, through an All-All communication, that every node is
//! owned exactly once. If the ownership is not a partition of unity it throws
//! "Owner is not a partition of unity". Opm::entityToDofIndexSet() expands the
//! node index set to the per-node block of mechanical degrees of freedom, so the
//! same partition-of-unity property is exercised on the DOF level as well.
//!
//! The number of MPI processes is the input to the test: run it under mpirun and
//! the communicator size determines how the deck is partitioned, e.g.
//!
//!   mpirun -np 4 ./test_partition_of_unity gridfilename=deck.grdecl \
//!       --add-corners=true --edge-conformal=true --edge-weights=uniform
//!
//==============================================================================
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define DISABLE_ALUGRID_SFC_ORDERING 1

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/version.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/GridEnums.hpp>

#include <opm/geomech/DuneCommunicationHelpers.hpp>

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

//! \brief Command line configurable parameters
struct Params
{
    std::string file;
    bool add_corners {true};
    bool edge_conformal {true};
    //! \brief Edge weight method for the partitioner: uniform|trans|logtrans
    Dune::EdgeWeightMethod edge_weights {Dune::EdgeWeightMethod::uniformEdgeWgt};
    //! \brief Partition method: simple|zoltan|metis|zoltanwell (zoltanwell ==
    //!        zoltanGoG, the well-aware method the simulator uses by default).
    Dune::PartitionMethod partition_method {Dune::PartitionMethod::zoltan};
    int overlap_layers {1};
    bool verbose {false};
};

void
syntax(char** argv)
{
    std::cerr
        << "Usage: mpirun -np <N> " << argv[0] << " gridfilename=deck.grdecl \\\n"
        << "           [--add-corners=true] [--edge-conformal=true] \\\n"
        << "           [--edge-weights=uniform|trans|logtrans] \\\n"
        << "           [--overlap-layers=1] [--verbose=false]\n\n"
        << "\t The number of MPI processes (mpirun -np <N>) is the input that\n"
        << "\t controls the partitioning. The test verifies that the node (vertex)\n"
        << "\t degrees of freedom form a partition of unity across the ranks.\n";
}

//! \brief Parse "true"/"false"/"1"/"0"/"yes"/"no" into a bool.
bool
parseBool(const std::string& value)
{
    return value == "true" || value == "1" || value == "yes" || value == "on";
}

Dune::EdgeWeightMethod
parseEdgeWeights(const std::string& value)
{
    if (value == "uniform") {
        return Dune::EdgeWeightMethod::uniformEdgeWgt;
    }
    if (value == "trans" || value == "default" || value == "defaultTrans") {
        return Dune::EdgeWeightMethod::defaultTransEdgeWgt;
    }
    if (value == "logtrans" || value == "logTrans") {
        return Dune::EdgeWeightMethod::logTransEdgeWgt;
    }
    throw std::runtime_error("Unknown --edge-weights value: " + value
                             + " (expected uniform|trans|logtrans)");
}

//! \brief Mirror FlowGenericVanguard's string -> PartitionMethod mapping.
Dune::PartitionMethod
parsePartitionMethod(const std::string& value)
{
    if (value == "simple") {
        return Dune::PartitionMethod::simple;
    }
    if (value == "zoltan") {
        return Dune::PartitionMethod::zoltan;
    }
    if (value == "metis") {
        return Dune::PartitionMethod::metis;
    }
    if (value == "zoltanwell" || value == "zoltanGoG") {
        return Dune::PartitionMethod::zoltanGoG;
    }
    throw std::runtime_error("Unknown --partition-method value: " + value
                             + " (expected simple|zoltan|metis|zoltanwell)");
}

//! \brief Accept both "key=value" and "--key=value" forms on the command line.
Params
parseCommandLine(int argc, char** argv)
{
    Params p;
    bool got_file = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        // strip a leading "--"
        if (arg.rfind("--", 0) == 0) {
            arg = arg.substr(2);
        }
        const auto eq = arg.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        const std::string key = arg.substr(0, eq);
        const std::string value = arg.substr(eq + 1);

        if (key == "gridfilename" || key == "grid") {
            p.file = value;
            got_file = true;
        } else if (key == "add-corners" || key == "add_corners") {
            p.add_corners = parseBool(value);
        } else if (key == "edge-conformal" || key == "edge_conformal") {
            p.edge_conformal = parseBool(value);
        } else if (key == "edge-weights" || key == "edge_weights") {
            p.edge_weights = parseEdgeWeights(value);
        } else if (key == "partition-method" || key == "partition_method") {
            p.partition_method = parsePartitionMethod(value);
        } else if (key == "overlap-layers" || key == "overlap_layers") {
            p.overlap_layers = std::stoi(value);
        } else if (key == "verbose") {
            p.verbose = parseBool(value);
        }
    }
    if (!got_file) {
        throw std::runtime_error("Missing required parameter gridfilename=<deck>");
    }
    return p;
}

const char*
edgeWeightName(Dune::EdgeWeightMethod m)
{
    switch (m) {
    case Dune::EdgeWeightMethod::uniformEdgeWgt:
        return "uniform";
    case Dune::EdgeWeightMethod::defaultTransEdgeWgt:
        return "trans";
    case Dune::EdgeWeightMethod::logTransEdgeWgt:
        return "logtrans";
    default:
        return "unknown";
    }
}

const char*
partitionMethodName(Dune::PartitionMethod m)
{
    switch (m) {
    case Dune::PartitionMethod::simple:
        return "simple";
    case Dune::PartitionMethod::zoltan:
        return "zoltan";
    case Dune::PartitionMethod::metis:
        return "metis";
    case Dune::PartitionMethod::zoltanGoG:
        return "zoltanwell";
    default:
        return "unknown";
    }
}

} // anonymous namespace

int
main(int argc, char** argv)
{
    Dune::MPIHelper& mpi = Dune::MPIHelper::instance(argc, argv);
    auto world_comm = mpi.getCommunication();

    if (argc < 2 || std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        if (world_comm.rank() == 0) {
            syntax(argv);
        }
        return 1;
    }

    try {
        const Params p = parseCommandLine(argc, argv);

        if (world_comm.rank() == 0) {
            std::cout << "Partition-of-unity test on " << world_comm.size() << " process(es)\n"
                      << "  deck           : " << p.file << '\n'
                      << "  add-corners    : " << std::boolalpha << p.add_corners << '\n'
                      << "  edge-conformal : " << std::boolalpha << p.edge_conformal << '\n'
                      << "  edge-weights   : " << edgeWeightName(p.edge_weights) << '\n'
                      << "  partition      : " << partitionMethodName(p.partition_method) << '\n'
                      << "  overlap-layers : " << p.overlap_layers << std::endl;
        }

        // --- build the grid from the deck (honouring --edge-conformal) ---------
        // Use a lenient ParseContext (WARN instead of THROW) so that minor deck
        // quirks irrelevant to the grid topology - e.g. an extra record
        // terminator - do not abort the test, matching how the simulator runs.
        Opm::Parser parser;
        Opm::ParseContext parseContext(Opm::InputErrorAction::WARN);
        Opm::ErrorGuard errorGuard;
        const auto deck = parser.parseFile(p.file, parseContext, errorGuard);
        Opm::EclipseState eclState(deck);

        // Pass &eclState (not nullptr) so the grid undergoes the same PORV /
        // pinch based cell processing as the simulator's vanguard
        // (GenericCpGridVanguard::doCreateGrids_). Passing nullptr skips that
        // and yields a different active-cell set / node topology.
        Dune::CpGrid grid;
        grid.processEclipseFormat(&eclState.getInputGrid(),
                                  /*eclipseState*/ &eclState,
                                  /*periodic*/ false,
                                  /*clip_z*/ false,
                                  /*pinchActive*/ false,
                                  /*edge_conformal*/ p.edge_conformal);

        if (world_comm.rank() == 0) {
            std::cout << "Global grid size (cells): " << grid.size(0) << std::endl;
        }

        // The well-aware partitioner (zoltanwell == zoltanGoG) needs the wells:
        // it collapses all cells perforated by a well into a single graph
        // vertex. This is exactly what the simulator passes (schedule.getWellsatEnd()).
        std::vector<Opm::Well> wells;
        if (p.partition_method == Dune::PartitionMethod::zoltanGoG) {
            auto python = std::make_shared<Opm::Python>();
            Opm::Schedule schedule(deck, eclState, parseContext, errorGuard, python);
            wells = schedule.getWellsatEnd();
            if (world_comm.rank() == 0) {
                std::cout << "Wells passed to partitioner: " << wells.size() << std::endl;
            }
        }

        // --- load balance with the requested edge weights and corner cells -----
        if (world_comm.size() > 1) {
            grid.loadBalance(p.edge_weights,
                             /*wells*/ wells.empty() ? nullptr : &wells,
                             /*possibleFutureConnections*/ {},
                             /*transmissibilities*/ nullptr,
                             /*ownersFirst*/ false,
                             /*addCornerCells*/ p.add_corners,
                             /*overlapLayers*/ p.overlap_layers,
                             /*partitionMethod*/ p.partition_method,
                             /*imbalanceTol*/ 1.1,
                             /*level*/ -1);
            grid.switchToDistributedView();
        }

        for (int rank = 0; rank < world_comm.size(); ++rank) {
            if (rank == world_comm.rank()) {
                std::cout << "  rank " << rank << " local cells " << grid.size(0) << " local nodes "
                          << grid.size(3) << std::endl;
            }
            world_comm.barrier();
        }

        // --- the actual partition-of-unity check on the node DOFs --------------
        // makeEntityEntityCommunication<3> selects a unique owner for every
        // shared vertex and throws "Owner is not a partition of unity" if any
        // node ends up with zero or more than one owner.
        constexpr int nodeCodim = 3;
        auto node_indexset
            = Opm::makeEntityEntityCommunication<nodeCodim>(grid, p.verbose);

        // Expand the node partition of unity to the mechanical DOFs (3 per node)
        // so the per-DOF index set is exercised as well.
        constexpr int dofsPerNode = 3;
        auto dof_indexset = Opm::entityToDofIndexSet(node_indexset, dofsPerNode, p.verbose);

        // Sanity: the number of owned nodes summed over ranks must equal the
        // global number of nodes (every node owned exactly once).
        std::size_t owned_nodes = 0;
        for (const auto& ind : node_indexset) {
            if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                ++owned_nodes;
            }
        }
        std::size_t owned_dofs = 0;
        for (const auto& ind : dof_indexset) {
            if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                ++owned_dofs;
            }
        }
        const std::size_t total_owned_nodes = world_comm.sum(owned_nodes);
        const std::size_t total_owned_dofs = world_comm.sum(owned_dofs);

        if (world_comm.rank() == 0) {
            std::cout << "Owned nodes (sum over ranks): " << total_owned_nodes << std::endl;
            std::cout << "Owned node DOFs (sum over ranks): " << total_owned_dofs << std::endl;
            if (total_owned_dofs != dofsPerNode * total_owned_nodes) {
                std::cout << "FAILED: owned DOF count is inconsistent with owned node count"
                          << std::endl;
            } else {
                std::cout << "PASSED: node degrees of freedom form a partition of unity"
                          << std::endl;
            }
        }

        if (total_owned_dofs != dofsPerNode * total_owned_nodes) {
            return 2;
        }
        return 0;
    } catch (const Dune::Exception& e) {
        std::cerr << "Rank " << world_comm.rank() << " Dune error: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Rank " << world_comm.rank() << " error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Rank " << world_comm.rank() << " unknown exception" << std::endl;
        return 1;
    }
}
