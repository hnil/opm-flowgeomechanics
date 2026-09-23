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
#ifndef OPM_MECH_GRID_PARTITION_HPP
#define OPM_MECH_GRID_PARTITION_HPP

#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgrid/RetainedCornerPointInput.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{

/// Partitioning for a run where the mechanics grid is a coarsening of the flow
/// grid: the coarse cells are partitioned, and every flow cell follows its
/// coarse cell. Both grids then have the same rank boundaries, so a mechanics
/// cell and all of its flow cells are always on one rank.
namespace MechPartition
{

/// COARSEN records, one per line: I1 I2 J1 J2 K1 K2 NX NY NZ, 1-based and
/// inclusive. '--' and '#' start a comment.
inline std::vector<Coarsening::CoarsenRequest> readRecords(const std::string& filename)
{
    std::ifstream is(filename);
    if (!is) {
        OPM_THROW(std::runtime_error, "Cannot open mechanics coarsening file " + filename);
    }

    std::vector<Coarsening::CoarsenRequest> requests;
    std::string line;
    while (std::getline(is, line)) {
        const auto comment = line.find_first_of("-#");
        if (comment != std::string::npos
            && (line[comment] == '#' || line.compare(comment, 2, "--") == 0)) {
            line.erase(comment);
        }
        std::istringstream ls(line);
        std::vector<int> v;
        for (int x; ls >> x;) {
            v.push_back(x);
        }
        if (v.empty()) {
            continue;
        }
        if (v.size() != 9) {
            OPM_THROW(std::runtime_error,
                      "A coarsening record needs 9 numbers (I1 I2 J1 J2 K1 K2 NX NY NZ): " + line);
        }
        Coarsening::CoarsenRequest r;
        for (int d = 0; d < 3; ++d) {
            r.startIJK[d] = v[2*d] - 1;
            r.endIJK[d] = v[2*d + 1];
            r.cellsPerDim[d] = v[6 + d];
        }
        requests.push_back(r);
    }
    return requests;
}

namespace detail
{

/// Recursive coordinate bisection over the coarse index box.
inline void bisect(std::vector<int>& cells, const std::array<int,3>& dims,
                   const std::vector<double>& weight, int firstRank, int numRanks,
                   std::vector<int>& part)
{
    if (numRanks <= 1) {
        for (const int c : cells) {
            part[c] = firstRank;
        }
        return;
    }

    const auto coord = [&dims](int c, int axis) {
        switch (axis) {
        case 0:  return c % dims[0];
        case 1:  return (c/dims[0]) % dims[1];
        default: return c/(dims[0]*dims[1]);
        }
    };

    int axis = 0;
    int widest = -1;
    for (int a = 0; a < 3; ++a) {
        int lo = std::numeric_limits<int>::max(), hi = std::numeric_limits<int>::min();
        for (const int c : cells) {
            lo = std::min(lo, coord(c, a));
            hi = std::max(hi, coord(c, a));
        }
        if (hi - lo > widest) {
            widest = hi - lo;
            axis = a;
        }
    }

    std::sort(cells.begin(), cells.end(), [&](int a, int b) {
        const int ca = coord(a, axis), cb = coord(b, axis);
        return (ca != cb) ? (ca < cb) : (a < b);
    });

    const int leftRanks = numRanks/2;
    double total = 0.0;
    for (const int c : cells) {
        total += weight[c];
    }

    // Split at the weighted median, but leave every rank at least one cell.
    const double target = total*leftRanks/numRanks;
    std::size_t split = 0;
    double acc = 0.0;
    while (split < cells.size() && acc + weight[cells[split]] <= target) {
        acc += weight[cells[split]];
        ++split;
    }
    const std::size_t minLeft = static_cast<std::size_t>(leftRanks);
    const std::size_t maxLeft = cells.size() - static_cast<std::size_t>(numRanks - leftRanks);
    split = std::clamp(split, minLeft, maxLeft);

    std::vector<int> left(cells.begin(), cells.begin() + split);
    std::vector<int> right(cells.begin() + split, cells.end());
    bisect(left, dims, weight, firstRank, leftRanks, part);
    bisect(right, dims, weight, firstRank + leftRanks, numRanks - leftRanks, part);
}

} // namespace detail

/// Rank of each coarse Cartesian cell. Cells with zero weight carry no flow
/// cells and are left out of the balance; they follow their neighbours.
inline std::vector<int> partitionCoarseCells(const std::array<int,3>& coarseDims,
                                             const std::vector<double>& weight,
                                             const int numRanks)
{
    std::vector<int> part(weight.size(), 0);
    if (numRanks <= 1) {
        return part;
    }

    std::vector<int> cells;
    cells.reserve(weight.size());
    for (std::size_t c = 0; c < weight.size(); ++c) {
        if (weight[c] > 0.0) {
            cells.push_back(static_cast<int>(c));
        }
    }
    if (static_cast<int>(cells.size()) < numRanks) {
        throw std::runtime_error("Fewer coarse mechanics cells than MPI ranks");
    }
    detail::bisect(cells, coarseDims, weight, 0, numRanks, part);

    // An empty cell takes the rank of the nearest one along the Cartesian
    // ordering, so the mechanics grid has no rank without cells either.
    int last = part[cells.front()];
    for (std::size_t c = 0; c < part.size(); ++c) {
        if (weight[c] > 0.0) {
            last = part[c];
        } else {
            part[c] = last;
        }
    }
    return part;
}

/// The rank of each cell of the flow grid's Cartesian space, computed once by
/// the load balancer and read back when the mechanics grid is distributed, so
/// the two cannot disagree.
class Registry
{
public:
    static void set(std::vector<int> partOfCoarseCell)
    { instance() = std::move(partOfCoarseCell); }

    static const std::vector<int>& get() { return instance(); }

    static bool empty() { return instance().empty(); }

private:
    static std::vector<int>& instance()
    {
        static std::vector<int> part;
        return part;
    }
};

/// The rank of every flow cell: the rank of the coarse cell it belongs to.
/// Install with GenericCpGridVanguard::setExternalLoadBalancer so that flow is
/// partitioned this way from the start.
template <class Grid>
std::vector<int> flowPartition(const Grid& grid,
                               const std::vector<Coarsening::CoarsenRequest>& requests)
{
    const auto size = grid.logicalCartesianSize();
    const std::array<int,3> dims{size[0], size[1], size[2]};
    // The blocks, not the coarse grid: this has to work for a coarsening the
    // mechanics builds by merging cells as well as for one written as
    // COORD/ZCORN.
    const auto layout = Coarsening::blockLayout(dims, requests);

    const auto anchor = [&layout, &dims](int cartesian) {
        const auto& box = layout.boxes[layout.blockOfCartesian[cartesian]];
        return box[0] + dims[0]*(box[1] + static_cast<std::size_t>(dims[1])*box[2]);
    };

    // Each block is weighted where its first corner sits, so the bisection
    // runs over the Cartesian box as usual.
    const auto& globalCell = grid.globalCell();
    std::vector<double> weight(layout.blockOfCartesian.size(), 0.0);
    for (const int cartesian : globalCell) {
        weight[anchor(cartesian)] += 1.0;
    }

    const auto partOfAnchor = partitionCoarseCells(dims, weight, grid.comm().size());

    // Rank of every Cartesian cell, which the mechanics grid reads back.
    std::vector<int> partOfCartesian(layout.blockOfCartesian.size(), 0);
    for (std::size_t c = 0; c < partOfCartesian.size(); ++c) {
        partOfCartesian[c] = partOfAnchor[anchor(static_cast<int>(c))];
    }
    Registry::set(partOfCartesian);

    std::vector<int> flowPart(globalCell.size(), 0);
    for (std::size_t c = 0; c < globalCell.size(); ++c) {
        flowPart[c] = partOfCartesian[globalCell[c]];
    }
    return flowPart;
}

/// The value of --mech-coarsen-file on the command line, or empty. Read
/// straight from argv because the partition has to be chosen before the
/// parameter system is up and the grid is built.
inline std::string coarsenFileFromArgs(int argc, char** argv)
{
    const std::string flag = "--mech-coarsen-file";
    for (int a = 1; a < argc; ++a) {
        const std::string arg(argv[a]);
        if (arg.rfind(flag + "=", 0) == 0) {
            return arg.substr(flag.size() + 1);
        }
        if (arg == flag && a + 1 < argc) {
            return argv[a + 1];
        }
    }
    return {};
}

/// Partition the flow grid by the mechanics coarsening, when one is asked for.
/// Call from main(), before the simulator is created.
template <class TypeTag>
void installFlowPartition(int argc, char** argv)
{
    const auto file = coarsenFileFromArgs(argc, argv);
    if (file.empty()) {
        return;
    }
    const auto requests = readRecords(file);
    // The mechanics grid is a coarsening of what flow was built from, so the
    // processing has to keep that description.
    RetainCornerPointInput::enable();
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    Vanguard::setExternalLoadBalancer(
        [requests](const Dune::CpGrid& grid) { return flowPartition(grid, requests); });
}

} // namespace MechPartition
} // namespace Opm

#endif // OPM_MECH_GRID_PARTITION_HPP
