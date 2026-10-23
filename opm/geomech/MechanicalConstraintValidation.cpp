#include <config.h>
#include <opm/geomech/MechanicalConstraintValidation.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>
#include <dune/grid/common/rangegenerators.hh>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

namespace Opm::Elasticity {
namespace {
class Components
{
public:
    std::int64_t find(std::int64_t id)
    {
        parent_.emplace(id, id);
        auto root = id;
        while (parent_.at(root) != root) {
            root = parent_.at(root);
        }
        while (id != root) {
            const auto next = parent_.at(id);
            parent_[id] = root;
            id = next;
        }
        return root;
    }
    void join(std::int64_t a, std::int64_t b)
    {
        a = find(a);
        b = find(b);
        parent_[std::max(a, b)] = std::min(a, b);
    }
    const auto& entries() const { return parent_; }
private:
    std::map<std::int64_t, std::int64_t> parent_;
};

int constraintRank(const std::vector<ConstrainedMechanicalVertex>& nodes)
{
    if (nodes.empty()) {
        return 0;
    }
    auto low = nodes.front().position;
    auto high = low;
    for (const auto& node : nodes) {
        for (int d = 0; d < 3; ++d) {
            low[d] = std::min(low[d], node.position[d]);
            high[d] = std::max(high[d], node.position[d]);
        }
    }
    // Center and nondimensionalize before rank estimation. Translation of the
    // coordinate origin and a uniform change of units must not change the rank.
    double scale = 0.0;
    for (int d = 0; d < 3; ++d) {
        scale = std::max(scale, high[d] - low[d]);
    }
    if (scale == 0.0) {
        scale = 1.0;
    }
    std::vector<std::array<double, 6>> basis;
    for (const auto& node : nodes) {
        std::array<double, 3> p;
        for (int d = 0; d < 3; ++d) {
            p[d] = ((node.position[d] - low[d]) - 0.5*(high[d] - low[d])) / scale;
        }
        const std::array<std::array<double, 6>, 3> rows{{
            {{1, 0, 0, 0, p[2], -p[1]}},
            {{0, 1, 0, -p[2], 0, p[0]}},
            {{0, 0, 1, p[1], -p[0], 0}}
        }};
        for (int d = 0; d < 3; ++d) {
            if (!(node.mask & (1u << d))) {
                continue;
            }
            auto row = rows[d];
            // Reorthogonalized Gram-Schmidt, avoiding normal equations and
            // their squared condition number. Six columns bound the work.
            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& q : basis) {
                    const double dot = std::inner_product(row.begin(), row.end(), q.begin(), 0.0);
                    for (int j = 0; j < 6; ++j) {
                        row[j] -= dot*q[j];
                    }
                }
            }
            const double norm = std::sqrt(std::inner_product(row.begin(), row.end(), row.begin(), 0.0));
            if (norm > 1e-10) {
                for (double& value : row) {
                    value /= norm;
                }
                basis.push_back(row);
                if (basis.size() == 6) {
                    return 6;
                }
            }
        }
    }
    return static_cast<int>(basis.size());
}
} // namespace

std::vector<MechanicalConstraintComponent> rigidBodyConstraintRanks(
    const std::vector<std::pair<std::int64_t, std::int64_t>>& edges,
    const std::vector<ConstrainedMechanicalVertex>& constraints)
{
    Components components;
    for (const auto& [a, b] : edges) {
        components.join(a, b);
    }
    std::map<std::int64_t, ConstrainedMechanicalVertex> unique;
    for (const auto& node : constraints) {
        if (node.mask == 0) {
            continue;
        }
        for (double coordinate : node.position) {
            if (!std::isfinite(coordinate)) {
                throw std::runtime_error("Mechanical constraint has a nonfinite vertex coordinate");
            }
        }
        auto [it, inserted] = unique.emplace(node.id, node);
        if (!inserted) {
            it->second.mask |= node.mask;
        }
    }
    std::map<std::int64_t, std::vector<ConstrainedMechanicalVertex>> groups;
    for (const auto& [id, parent] : components.entries()) {
        groups[components.find(id)];
    }
    for (const auto& [id, node] : unique) {
        if (components.entries().count(id) == 0) {
            throw std::runtime_error("Mechanical constraint references a vertex outside the mechanical grid");
        }
        groups[components.find(id)].push_back(node);
    }
    std::vector<MechanicalConstraintComponent> result;
    for (const auto& [id, nodes] : groups) {
        result.push_back({id, constraintRank(nodes)});
    }
    return result;
}

std::vector<MechanicalConstraintComponent> mechanicalConstraintRanks(
    const Dune::CpGrid& grid,
    const std::vector<std::tuple<std::size_t, MechBCValue>>& boundaryNodes)
{
    const auto gv = grid.leafGridView();
    const auto& comm = grid.comm();
    const std::size_t numVertices = gv.size(3);
    std::vector<std::int64_t> ids(numVertices);
    std::vector<std::array<double, 3>> positions(numVertices);
    for (const auto& vertex : vertices(gv)) {
        const auto index = gv.indexSet().index(vertex);
        ids[index] = grid.globalIdSet().id(vertex);
        const auto p = vertex.geometry().center();
        positions[index] = {p[0], p[1], p[2]};
    }
    // Collapse rank-local cell connectivity before communication: one pair
    // per local vertex, rather than transmitting every cell's corner list.
    Components local;
    for (const auto& cell : elements(gv)) {
        const auto first = ids[gv.indexSet().subIndex(cell, 0, 3)];
        for (unsigned corner = 0; corner < cell.subEntities(3); ++corner) {
            local.join(first, ids[gv.indexSet().subIndex(cell, corner, 3)]);
        }
    }
    std::vector<std::int64_t> localEdges;
    for (const auto& [id, parent] : local.entries()) {
        localEdges.push_back(id);
        localEdges.push_back(local.find(id));
    }
    int invalid = 0;
    std::vector<std::int64_t> localMasks;
    std::vector<double> localPositions;
    for (const auto& [index, bc] : boundaryNodes) {
        if (index >= numVertices) {
            invalid = 1;
            continue;
        }
        unsigned mask = 0;
        for (int d = 0; d < 3; ++d) {
            if (bc.fixeddir[d]) {
                mask |= 1u << d;
            }
        }
        if (mask == 0) {
            continue;
        }
        localMasks.insert(localMasks.end(), {ids[index], static_cast<std::int64_t>(mask)});
        localPositions.insert(localPositions.end(), positions[index].begin(), positions[index].end());
    }
    if (comm.max(invalid)) {
        throw std::runtime_error("Mechanical boundary condition has an invalid local vertex index");
    }
    const auto gatheredEdges = Opm::gatherv(localEdges, comm, 0).first;
    const auto gatheredMasks = Opm::gatherv(localMasks, comm, 0).first;
    const auto gatheredPositions = Opm::gatherv(localPositions, comm, 0).first;
    std::vector<MechanicalConstraintComponent> result;
    std::string error;
    if (comm.rank() == 0) {
        try {
            std::vector<std::pair<std::int64_t, std::int64_t>> edges;
            for (std::size_t i = 0; i < gatheredEdges.size(); i += 2) {
                edges.emplace_back(gatheredEdges[i], gatheredEdges[i+1]);
            }
            std::vector<ConstrainedMechanicalVertex> constraints;
            for (std::size_t i = 0; i < gatheredMasks.size()/2; ++i) {
                constraints.push_back({gatheredMasks[2*i],
                    {gatheredPositions[3*i], gatheredPositions[3*i+1], gatheredPositions[3*i+2]},
                    static_cast<unsigned>(gatheredMasks[2*i+1])});
            }
            result = rigidBodyConstraintRanks(edges, constraints);
        } catch (const std::exception& e) {
            error = e.what();
        }
    }
    int errorSize = static_cast<int>(error.size());
    comm.broadcast(&errorSize, 1, 0);
    error.resize(errorSize);
    if (errorSize) {
        comm.broadcast(error.data(), errorSize, 0);
        throw std::runtime_error(error);
    }
    int count = static_cast<int>(result.size());
    comm.broadcast(&count, 1, 0);
    result.resize(count);
    for (auto& component : result) {
        comm.broadcast(&component.id, 1, 0);
        comm.broadcast(&component.rank, 1, 0);
    }
    return result;
}
} // namespace Opm::Elasticity
