#include <config.h>
#include <opm/geomech/MechanicalVertexStencilValidation.hpp>
#if HAVE_MPI
#include <opm/geomech/DuneCommunicationHelpers.hpp>
#endif
#include <opm/geomech/vem/vemutils.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace {
using Entries = std::map<std::pair<int, int>, double>;
struct System { Entries matrix; std::map<int, double> load; };
std::vector<int> physicalDofs(const Dune::CpGrid& grid)
{
    const auto gv = grid.leafGridView();
    std::vector<int> dofs(3*gv.size(3));
    for (const auto& vertex : vertices(gv)) {
        const auto p = vertex.geometry().center();
        const int v = static_cast<int>(std::lround(p[0]))
            + 5*static_cast<int>(std::lround(p[1])) + 25*static_cast<int>(std::lround(p[2]));
        for (int d = 0; d < 3; ++d) {
            dofs[3*gv.indexSet().index(vertex)+d] = 3*v+d;
        }
    }
    return dofs;
}
System assemble(const Dune::CpGrid& grid)
{
    std::vector<double> coords;
    std::vector<int> nf, nv, nodes;
    std::vector<char> mask;
    vem::getGridVectors(grid, coords, nf, nv, nodes, &mask);
    const auto gv = grid.leafGridView();
    std::vector<double> young(gv.size(0), 2.3e9), poisson(gv.size(0), 0.23);
    std::vector<double> body(3*gv.size(0));
    for (const auto& cell : elements(gv)) {
        const auto p = cell.geometry().center();
        const int c = gv.indexSet().index(cell);
        body[3*c] = 1+p[0]; body[3*c+1] = 2+p[1]; body[3*c+2] = 3+p[2];
    }
    std::vector<std::tuple<int, int, double>> entries;
    std::vector<double> load;
    vem::assemble_mech_system_3D_dune(grid, coords.data(), gv.size(0), nf.data(),
        nv.data(), nodes.data(), young.data(), poisson.data(), body.data(),
        0, nullptr, nullptr, 0, nullptr, nullptr, entries, load, vem::SIMPLE, false, mask.data());
    const auto dofs = physicalDofs(grid);
    System result;
    for (const auto& [i, j, value] : entries) {
        result.matrix[{dofs[i], dofs[j]}] += value;
    }
    for (std::size_t i = 0; i < load.size(); ++i) {
        result.load[dofs[i]] = load[i];
    }
    return result;
}
int mismatchedRows(const System& reference, const System& local, const std::set<int>& owners)
{
    int bad = 0;
    for (const int row : owners) {
        double error = 0, scale = 1;
        std::set<int> columns;
        for (const auto* a : {&reference.matrix, &local.matrix}) {
            for (auto it = a->lower_bound({row, -1}); it != a->end() && it->first.first == row; ++it) {
                columns.insert(it->first.second);
            }
        }
        for (const int col : columns) {
            const auto value = [&](const Entries& a) {
                const auto it = a.find({row, col});
                return it == a.end() ? 0.0 : it->second;
            };
            const double ref = value(reference.matrix);
            scale = std::max(scale, std::abs(ref));
            error = std::max(error, std::abs(ref-value(local.matrix)));
        }
        const double ref = reference.load.at(row);
        const double loadError = std::abs(ref-local.load.at(row));
        bad += error > 1e-10*scale || loadError > 1e-10*std::max(1.0, std::abs(ref));
    }
    return bad;
}
void check(const System& reference, bool corners)
{
    Dune::CpGrid grid;
    grid.createCartesian({4, 4, 4}, {1, 1, 1});
    if (grid.comm().size() > 1) {
        if (grid.comm().size() != 4) {
            throw std::runtime_error("This partition fixture requires one or four ranks");
        }
        std::vector<int> parts(64);
        for (int k = 0; k < 4; ++k)
            for (int j = 0; j < 4; ++j)
                for (int i = 0; i < 4; ++i)
                    parts[i+4*j+16*k] = i >= 2 ? 1 : j >= 2 ? 2 : k >= 2 ? 3 : 0;
        grid.setCornerCellsByVertex(corners);
        grid.loadBalance(parts, false, corners, 2);
    }
    std::set<int> owners;
    bool rejected = false;
#if HAVE_MPI
    if (grid.comm().size() > 1) {
        using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
        Comm comm(grid.comm(), Dune::SolverCategory::overlapping, false);
        comm.indexSet() = Opm::entityToDofIndexSet(Opm::makeEntityEntityCommunication<3>(grid), 3);
        comm.remoteIndices().rebuild<false>();
        const auto dofs = physicalDofs(grid);
        for (const auto& index : comm.indexSet()) {
            if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                owners.insert(dofs[index.local().local()]);
            }
        }
        try {
            Opm::Elasticity::validateMechanicalVertexStencils(grid, comm);
        } catch (const std::runtime_error& e) {
            rejected = std::string(e.what()).find("Incomplete VEM vertex stencils") != std::string::npos;
        }
        bool preflightRejected = false;
        try {
            Opm::Elasticity::validateMechanicalVertexStencils(grid);
        } catch (const std::runtime_error& e) {
            preflightRejected = std::string(e.what()).find("Incomplete VEM vertex stencils") != std::string::npos;
        }
        if (preflightRejected != rejected) {
            throw std::runtime_error("Initialization and solver stencil checks must agree");
        }
    } else
#endif
    {
        for (const auto& [row, value] : reference.load) {
            owners.insert(row);
        }
    }
    const int bad = grid.comm().sum(mismatchedRows(reference, assemble(grid), owners));
    const bool shouldReject = !corners && grid.comm().size() > 1;
    if (grid.comm().min(static_cast<int>(rejected == shouldReject)) != 1
        || (shouldReject ? bad == 0 : bad != 0)) {
        throw std::runtime_error("Stencil validation must track the actual assembled owner-row mismatch");
    }
    if (grid.comm().rank() == 0) {
        std::cout << "corners=" << corners << " mismatched owner rows=" << bad
                  << " rejected=" << rejected << '\n';
    }
}
}
int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    try {
#if HAVE_MPI
        Dune::CpGrid serial(MPI_COMM_SELF);
#else
        Dune::CpGrid serial;
#endif
        serial.createCartesian({4, 4, 4}, {1, 1, 1});
        const auto reference = assemble(serial);
        check(reference, false);
        check(reference, true);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
