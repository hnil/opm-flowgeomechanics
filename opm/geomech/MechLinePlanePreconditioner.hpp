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

#ifndef OPM_GEOMECH_MECH_LINE_PLANE_PRECONDITIONER_HPP
#define OPM_GEOMECH_MECH_LINE_PLANE_PRECONDITIONER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvercategory.hh>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm
{

// Registers type "mech_lines" for the mechanics operators; coords must outlive the solver.
void addMechLinePlanePreconditioner(const std::vector<double>& coords);

// Symmetric block Gauss-Seidel over vertical lines and lateral plane blocks, around an inner
// preconditioner. Dofs are three per node, coords three per node. Blocks hold owned nodes only.
template <class Matrix, class Vector, class Comm>
class MechLinePlanePreconditioner : public Dune::PreconditionerWithUpdate<Vector, Vector>
{
public:
    using Inner = Dune::PreconditionerWithUpdate<Vector, Vector>;

    MechLinePlanePreconditioner(const Matrix& M,
                                const std::vector<double>& coords,
                                const PropertyTree& prm,
                                std::shared_ptr<Inner> inner,
                                const Comm& comm)
        : M_(M)
        , coords_(coords)
        , inner_(std::move(inner))
        , comm_(comm)
        , plane_theta_(prm.get<double>("plane_theta", 0.3))
        , plane_comps_(prm.get<int>("plane_comps", 3))
        , plane_cap_(prm.get<int>("plane_cap", 200))
    {
        setup();
    }

    void pre(Vector& v, Vector& d) override
    {
        inner_->pre(v, d);
    }

    void apply(Vector& v, const Vector& d) override
    {
        v = 0.0;
        // Ghost values stay frozen within a sweep, so the backward sweep is the forward sweep's transpose.
        for (std::size_t i = 0; i < blocks_.size(); ++i)
            sweep(blocks_[i], v, d);
        comm_.copyOwnerToAll(v, v);
        r_ = d;
        M_.mmv(v, r_);
        comm_.project(r_);
        c_ = 0.0; // the inner preconditioner takes this as its initial guess
        inner_->apply(c_, r_);
        v += c_;
        for (std::size_t i = blocks_.size(); i-- > 0;)
            sweep(blocks_[i], v, d);
        comm_.copyOwnerToAll(v, v);
    }

    void post(Vector& v) override
    {
        inner_->post(v);
    }

    void update() override
    {
        inner_->update();
        setup();
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

    Dune::SolverCategory::Category category() const override
    {
        return std::is_same_v<Comm, Dune::Amg::SequentialInformation> ? Dune::SolverCategory::sequential
                                                                      : Dune::SolverCategory::overlapping;
    }

private:
    struct Block {
        std::vector<int> dofs;
        int band = 0;
        std::vector<double> chol; // row i holds columns i-band..i
    };

    void setup()
    {
        Dune::Timer timer;
        blocks_.clear();
        const int nn = int(M_.N()) / 3;
        owned_.assign(nn, 1);
        if constexpr (!std::is_same_v<Comm, Dune::Amg::SequentialInformation>) {
            for (const auto& index : comm_.indexSet())
                if (index.local().attribute() != Dune::OwnerOverlapCopyAttributeSet::owner)
                    owned_[int(index.local().local()) / 3] = 0;
        }
        where_.assign(3 * nn, -1);
        addLines(nn);
        const std::size_t lines = blocks_.size();
        addPlanes(nn);
        r_.resize(M_.N());
        c_.resize(M_.N());

        const auto& cc = comm_.communicator();
        double counts[3] = {double(lines), double(blocks_.size() - lines), double(plane_nodes_)};
        cc.sum(counts, 3);
        const double seconds = cc.max(timer.elapsed());
        OpmLog::info("mech_lines: " + std::to_string(long(counts[0])) + " lines, " + std::to_string(long(counts[1]))
                     + " plane blocks over " + std::to_string(long(counts[2])) + " nodes, largest band "
                     + std::to_string(cc.max(max_band_)) + ", setup " + std::to_string(seconds) + " s");
    }

    void addLines(const int nn)
    {
        std::vector<int> up(nn, -1), down(nn, -1);
        std::vector<double> sup(nn, 0.0), sdown(nn, 0.0);
        std::map<int, double> coupling;
        for (int n = 0; n < nn; ++n) {
            if (!owned_[n])
                continue;
            coupling.clear();
            for (int d = 0; d < 3; ++d)
                for (auto c = M_[3 * n + d].begin(); c != M_[3 * n + d].end(); ++c)
                    coupling[int(c.index()) / 3] += double(*c) * double(*c);
            for (const auto& [o, s] : coupling) {
                if (o == n || !owned_[o] || lateral(n, o))
                    continue;
                const double dz = coords_[3 * o + 2] - coords_[3 * n + 2];
                if (dz > 0 && s > sup[n]) {
                    sup[n] = s;
                    up[n] = o;
                }
                if (dz < 0 && s > sdown[n]) {
                    sdown[n] = s;
                    down[n] = o;
                }
            }
        }
        for (int n = 0; n < nn; ++n) {
            if (!owned_[n] || (down[n] >= 0 && up[down[n]] == n))
                continue; // not the lowest node of its line
            Block b;
            for (int k = n;; k = up[k]) {
                for (int d = 0; d < 3; ++d)
                    b.dofs.push_back(3 * k + d);
                if (up[k] < 0 || down[up[k]] != k)
                    break;
            }
            finish(b);
        }
    }

    // Lateral uz couplings that are strong and positive mark modes a diagonal smoother cannot reach.
    void addPlanes(const int nn)
    {
        plane_nodes_ = 0;
        if (plane_theta_ <= 0)
            return;
        std::vector<std::vector<int>> adj(nn);
        for (int n = 0; n < nn; ++n) {
            if (!owned_[n])
                continue;
            for (auto c = M_[3 * n + 2].begin(); c != M_[3 * n + 2].end(); ++c) {
                const int o = int(c.index()) / 3;
                if (int(c.index()) % 3 != 2 || o == n || !owned_[o])
                    continue;
                const double dg = double(M_[3 * n + 2][3 * n + 2]) * double(M_[3 * o + 2][3 * o + 2]);
                if (double(*c) > plane_theta_ * std::sqrt(dg) && lateral(n, o))
                    adj[n].push_back(o);
            }
        }
        std::vector<char> seen(nn, 0);
        for (int n = 0; n < nn; ++n) {
            if (seen[n] || adj[n].empty())
                continue;
            std::vector<int> q {n};
            seen[n] = 1;
            for (std::size_t h = 0; h < q.size() && int(q.size()) < plane_cap_; ++h)
                for (int o : adj[q[h]])
                    if (!seen[o] && int(q.size()) < plane_cap_) {
                        seen[o] = 1;
                        q.push_back(o);
                    }
            if (q.size() < 4)
                continue;
            Block b;
            for (auto it = q.rbegin(); it != q.rend(); ++it)
                for (int d = 3 - plane_comps_; d < 3; ++d)
                    b.dofs.push_back(3 * *it + d);
            finish(b);
            plane_nodes_ += q.size();
        }
    }

    bool lateral(const int n, const int o) const
    {
        const double dx = coords_[3 * o] - coords_[3 * n];
        const double dy = coords_[3 * o + 1] - coords_[3 * n + 1];
        const double dz = coords_[3 * o + 2] - coords_[3 * n + 2];
        return std::abs(dz) <= std::hypot(dx, dy);
    }

    // Band Cholesky of the block's dofs in the order given.
    void finish(Block& b)
    {
        const int m = int(b.dofs.size());
        for (int i = 0; i < m; ++i)
            where_[b.dofs[i]] = i;
        for (int i = 0; i < m; ++i)
            for (auto c = M_[b.dofs[i]].begin(); c != M_[b.dofs[i]].end(); ++c)
                if (where_[c.index()] >= 0)
                    b.band = std::max(b.band, i - where_[c.index()]);
        const int w = b.band + 1;
        b.chol.assign(m * w, 0.0);
        for (int i = 0; i < m; ++i)
            for (auto c = M_[b.dofs[i]].begin(); c != M_[b.dofs[i]].end(); ++c) {
                const int j = where_[c.index()];
                if (j >= 0 && j <= i)
                    b.chol[i * w + j - i + b.band] = double(*c);
            }
        for (int i = 0; i < m; ++i)
            for (int j = std::max(0, i - b.band); j <= i; ++j) {
                double s = b.chol[i * w + j - i + b.band];
                for (int k = std::max(0, i - b.band); k < j; ++k)
                    if (j - k <= b.band)
                        s -= b.chol[i * w + k - i + b.band] * b.chol[j * w + k - j + b.band];
                if (j == i) {
                    if (s <= 0)
                        OPM_THROW(std::runtime_error, "mech_lines: block not positive definite");
                    b.chol[i * w + b.band] = std::sqrt(s);
                } else {
                    b.chol[i * w + j - i + b.band] = s / b.chol[j * w + b.band];
                }
            }
        for (int i : b.dofs)
            where_[i] = -1;
        max_band_ = std::max(max_band_, b.band);
        blocks_.push_back(std::move(b));
    }

    void sweep(const Block& b, Vector& x, const Vector& d)
    {
        const int m = int(b.dofs.size()), w = b.band + 1;
        y_.resize(m);
        for (int i = 0; i < m; ++i) {
            double s = d[b.dofs[i]];
            for (auto c = M_[b.dofs[i]].begin(); c != M_[b.dofs[i]].end(); ++c)
                s -= double(*c) * double(x[c.index()]);
            for (int k = std::max(0, i - b.band); k < i; ++k)
                s -= b.chol[i * w + k - i + b.band] * y_[k];
            y_[i] = s / b.chol[i * w + b.band];
        }
        for (int i = m - 1; i >= 0; --i) {
            for (int k = i + 1; k < std::min(m, i + w); ++k)
                y_[i] -= b.chol[k * w + i - k + b.band] * y_[k];
            y_[i] /= b.chol[i * w + b.band];
            x[b.dofs[i]] += y_[i];
        }
    }

    const Matrix& M_;
    const std::vector<double>& coords_;
    std::shared_ptr<Inner> inner_;
    const Comm& comm_;
    double plane_theta_;
    int plane_comps_;
    int plane_cap_;
    std::vector<Block> blocks_;
    std::vector<char> owned_;
    std::vector<int> where_;
    std::vector<double> y_;
    Vector r_, c_;
    std::size_t plane_nodes_ = 0;
    int max_band_ = 0;
};

} // namespace Opm

#endif // OPM_GEOMECH_MECH_LINE_PLANE_PRECONDITIONER_HPP
