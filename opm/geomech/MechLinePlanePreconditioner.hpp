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
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvercategory.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm
{

// Registers "mech_lines" and "mech_hypre" for the mechanics operators; all arguments must
// outlive the solver; axis is per node (pillar direction, layer normal) or empty.
void addMechPreconditioners(const std::vector<double>& coords,
                            const std::vector<Dune::BlockVector<Dune::FieldVector<double, 1>>>& modes,
                            const std::vector<double>& axis);

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
                                const Comm& comm,
                                const std::vector<double>* axis = nullptr)
        : M_(M)
        , coords_(coords)
        , inner_(std::move(inner))
        , comm_(comm)
        , plane_theta_(prm.get<double>("plane_theta", 0.3))
        , plane_comps_(prm.get<int>("plane_comps", 3))
        , plane_cap_(prm.get<int>("plane_cap", 200))
        , axis_(prm.get<bool>("structure_axis", false) && axis && !axis->empty() ? axis : nullptr)
        , emulate_ranks_(prm.get<int>("emulate_ranks", 1))
        , emulate_cut_(prm.get<std::string>("emulate_cut", "xy"))
        , additive_(prm.get<std::string>("mode", "multiplicative") == "additive")
        , cheb_degree_(prm.get<int>("cheb_degree", 2))
        , cheb_ratio_(prm.get<double>("cheb_ratio", 30.0))
        , line_ratio_(prm.get<double>("line_ratio", 0.0))
    {
        setup();
    }

    void pre(Vector& v, Vector& d) override
    {
        inner_->pre(v, d);
    }

    void apply(Vector& v, const Vector& d) override
    {
        if (additive_) {
            v = 0.0;
            chebyshev(v, d);
            residual(v, d, r_);
            c_ = 0.0;
            inner_->apply(c_, r_);
            v += c_;
            chebyshev(v, d);
            return;
        }
        v = 0.0;
        if (!part_.empty())
            x0_ = v;
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
        if (!part_.empty())
            x0_ = v;
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
        int part = 0;
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
        makeParts(nn);
        addLines(nn);
        std::size_t lines = blocks_.size();
        addPlanes(nn);
        if (blockParts(nn)) {
            blocks_.clear();
            max_band_ = 0;
            addLines(nn);
            lines = blocks_.size();
            addPlanes(nn);
        }
        r_.resize(M_.N());
        c_.resize(M_.N());
        z_.resize(M_.N());
        if (additive_)
            estimateLambdaMax();

        const auto& cc = comm_.communicator();
        double counts[4] = {double(lines), double(blocks_.size() - lines), double(plane_nodes_), double(line_nodes_)};
        cc.sum(counts, 4);
        const double seconds = cc.max(timer.elapsed());
        OpmLog::info("mech_lines: " + std::to_string(long(counts[0])) + " lines, " + std::to_string(long(counts[1]))
                     + " plane blocks over " + std::to_string(long(counts[2])) + " nodes, largest band "
                     + std::to_string(cc.max(max_band_)) + ", " + std::to_string(long(counts[3]))
                     + " nodes in multi-node lines, setup " + std::to_string(seconds) + " s"
                     + (part_.empty() ? std::string() : ", emulating " + std::to_string(emulate_ranks_)
                                                             + " ranks, cut " + emulate_cut_)
                     + (additive_ ? ", additive Chebyshev degree " + std::to_string(cheb_degree_)
                                        + ", lambda_max " + std::to_string(lambda_max_)
                                  : std::string()));
    }

    // r = d - M x, with x made consistent first.
    void residual(Vector& x, const Vector& d, Vector& r)
    {
        comm_.copyOwnerToAll(x, x);
        r = d;
        M_.mmv(x, r);
        comm_.project(r);
    }

    // Additive Schwarz over all blocks; lines and planes overlap, which the Chebyshev bounds absorb.
    void schwarz(const Vector& r, Vector& z)
    {
        z = 0.0;
        for (const auto& b : blocks_)
            solveAdd(b, r, z);
        comm_.copyOwnerToAll(z, z);
    }

    double dot(const Vector& a, const Vector& b) const
    {
        double s = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (owned_[i / 3])
                s += double(a[i]) * double(b[i]);
        return comm_.communicator().sum(s);
    }

    // Power iteration on S^{-1} M; Chebyshev targets [lambda_max / cheb_ratio, lambda_max].
    void estimateLambdaMax()
    {
        Vector q(M_.N()), w(M_.N());
        for (std::size_t i = 0; i < q.size(); ++i)
            q[i] = 1.0 + 0.5 * std::sin(0.7 * double(i));
        double lambda = 1.0;
        for (int it = 0; it < 15; ++it) {
            q *= 1.0 / std::sqrt(dot(q, q));
            comm_.copyOwnerToAll(q, q);
            M_.mv(q, w);
            comm_.project(w);
            schwarz(w, z_);
            lambda = std::sqrt(dot(z_, z_));
            q = z_;
        }
        lambda_max_ = 1.1 * lambda;
    }

    void chebyshev(Vector& x, const Vector& d)
    {
        const double lmin = lambda_max_ / cheb_ratio_;
        const double theta = 0.5 * (lambda_max_ + lmin), delta = 0.5 * (lambda_max_ - lmin);
        const double sigma = theta / delta;
        double rho = 1.0 / sigma;
        p_.resize(M_.N());
        res_.resize(M_.N());
        residual(x, d, res_);
        schwarz(res_, z_);
        p_ = z_;
        p_ *= 1.0 / theta;
        x += p_;
        for (int k = 1; k < cheb_degree_; ++k) {
            const double rho_new = 1.0 / (2.0 * sigma - rho);
            residual(x, d, res_);
            schwarz(res_, z_);
            p_ *= rho_new * rho;
            p_.axpy(2.0 * rho_new / delta, z_);
            x += p_;
            rho = rho_new;
        }
        comm_.copyOwnerToAll(x, x);
    }

    void addLines(const int nn)
    {
        std::vector<int> up(nn, -1), down(nn, -1);
        std::vector<double> sup(nn, 0.0), sdown(nn, 0.0);
        std::vector<double> dlat(nn, std::numeric_limits<double>::max());
        std::map<int, double> coupling;
        for (int n = 0; n < nn; ++n) {
            if (!owned_[n])
                continue;
            coupling.clear();
            for (int d = 0; d < 3; ++d)
                for (auto c = M_[3 * n + d].begin(); c != M_[3 * n + d].end(); ++c)
                    coupling[int(c.index()) / 3] += double(*c) * double(*c);
            for (const auto& [o, s] : coupling) {
                if (o == n || !owned_[o] || !samePart(n, o))
                    continue;
                if (lateral(n, o)) {
                    dlat[n] = std::min(dlat[n], dist(n, o));
                    continue;
                }
                const double dz = along(n, o);
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
        // line_ratio > 0: chain only nodes whose nearest lateral neighbour is line_ratio times
        // farther than the nearest vertical one (flat cells); the others become single-node blocks.
        // Geometric, since the stabilization couples all nodes of a flat cell strongly.
        if (line_ratio_ > 0) {
            std::vector<char> flat(nn);
            for (int n = 0; n < nn; ++n) {
                double dv = std::numeric_limits<double>::max();
                if (up[n] >= 0)
                    dv = dist(n, up[n]);
                if (down[n] >= 0)
                    dv = std::min(dv, dist(n, down[n]));
                flat[n] = dv < std::numeric_limits<double>::max() && dlat[n] > line_ratio_ * dv;
            }
            for (int n = 0; n < nn; ++n) {
                if (up[n] >= 0 && !(flat[n] && flat[up[n]]))
                    up[n] = -1;
                if (down[n] >= 0 && !(flat[n] && flat[down[n]]))
                    down[n] = -1;
            }
        }
        line_nodes_ = 0;
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
            if (b.dofs.size() > 3)
                line_nodes_ += b.dofs.size() / 3;
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
                if (int(c.index()) % 3 != 2 || o == n || !owned_[o] || !samePart(n, o))
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

    // Diagnostic, sequential only: emulate_ranks parts by coordinate bisection (emulate_cut "xy":
    // vertical cuts, "z": horizontal). Blocks stay inside a part and read other parts' values as
    // frozen within a sweep, as the parallel blocks read ghosts.
    void makeParts(const int nn)
    {
        part_.clear();
        if (emulate_ranks_ <= 1 || !std::is_same_v<Comm, Dune::Amg::SequentialInformation>
            || emulate_cut_ == "blocks")
            return;
        part_.assign(nn, 0);
        std::vector<int> all(nn);
        std::iota(all.begin(), all.end(), 0);
        bisect(all, 0, emulate_ranks_, part_);
        x0_.resize(M_.N());
    }

    // emulate_cut "blocks": parts made of whole connected groups of blocks; a group larger than a
    // part's share is split by coordinates. Returns true when the blocks must be rebuilt.
    bool blockParts(const int nn)
    {
        if (emulate_ranks_ <= 1 || emulate_cut_ != "blocks"
            || !std::is_same_v<Comm, Dune::Amg::SequentialInformation>)
            return false;
        std::vector<int> root(nn);
        std::iota(root.begin(), root.end(), 0);
        auto find = [&root](int a) {
            while (root[a] != a)
                a = root[a] = root[root[a]];
            return a;
        };
        for (const auto& b : blocks_)
            for (int i : b.dofs)
                root[find(i / 3)] = find(b.dofs[0] / 3);
        std::map<int, std::vector<int>> groups;
        for (int n = 0; n < nn; ++n)
            groups[find(n)].push_back(n);
        const double share = double(nn) / emulate_ranks_;
        std::vector<std::vector<int>> items;
        std::size_t largest = 0;
        int split = 0;
        std::vector<int> label(nn, 0);
        for (auto& [r, g] : groups) {
            largest = std::max(largest, g.size());
            if (double(g.size()) <= share) {
                items.push_back(std::move(g));
                continue;
            }
            const int k = int(std::ceil(g.size() / share));
            bisect(g, 0, k, label);
            std::vector<std::vector<int>> pieces(k);
            for (int n : g)
                pieces[label[n]].push_back(n);
            for (auto& q : pieces)
                items.push_back(std::move(q));
            ++split;
        }
        part_.assign(nn, 0);
        std::vector<int> idx(items.size());
        std::iota(idx.begin(), idx.end(), 0);
        assignItems(items, idx, 0, emulate_ranks_);
        x0_.resize(M_.N());
        OpmLog::info("mech_lines emulation: " + std::to_string(groups.size()) + " block groups, largest "
                     + std::to_string(largest) + " of " + std::to_string(nn) + " nodes, " + std::to_string(split)
                     + " split by coordinates");
        return true;
    }

    // Weighted bisection of the items' centroids along the longer horizontal axis.
    void assignItems(const std::vector<std::vector<int>>& items, std::vector<int> idx, const int first,
                     const int count)
    {
        if (count <= 1 || idx.empty()) {
            for (int i : idx)
                for (int n : items[i])
                    part_[n] = first;
            return;
        }
        std::vector<std::array<double, 2>> c(items.size());
        double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300}, total = 0;
        for (int i : idx) {
            c[i] = {0, 0};
            for (int n : items[i])
                for (int a = 0; a < 2; ++a)
                    c[i][a] += coords_[3 * n + a] / items[i].size();
            for (int a = 0; a < 2; ++a) {
                lo[a] = std::min(lo[a], c[i][a]);
                hi[a] = std::max(hi[a], c[i][a]);
            }
            total += items[i].size();
        }
        const int axis = hi[0] - lo[0] >= hi[1] - lo[1] ? 0 : 1;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return c[a][axis] < c[b][axis]; });
        const int left = count / 2;
        double acc = 0;
        std::size_t m = 0;
        while (m < idx.size() && acc + 0.5 * items[idx[m]].size() < total * left / count)
            acc += items[idx[m++]].size();
        assignItems(items, std::vector<int>(idx.begin(), idx.begin() + m), first, left);
        assignItems(items, std::vector<int>(idx.begin() + m, idx.end()), first + left, count - left);
    }

    void bisect(std::vector<int> nodes, const int first, const int count, std::vector<int>& label)
    {
        if (count <= 1 || nodes.empty()) {
            for (int n : nodes)
                label[n] = first;
            return;
        }
        int axis = 2;
        if (emulate_cut_ != "z") {
            double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
            for (int n : nodes)
                for (int a = 0; a < 2; ++a) {
                    lo[a] = std::min(lo[a], coords_[3 * n + a]);
                    hi[a] = std::max(hi[a], coords_[3 * n + a]);
                }
            axis = hi[0] - lo[0] >= hi[1] - lo[1] ? 0 : 1;
        }
        const int left = count / 2;
        const auto mid = nodes.begin() + std::ptrdiff_t(nodes.size() * left / count);
        std::nth_element(nodes.begin(), mid, nodes.end(),
                         [&](int a, int b) { return coords_[3 * a + axis] < coords_[3 * b + axis]; });
        bisect(std::vector<int>(nodes.begin(), mid), first, left, label);
        bisect(std::vector<int>(mid, nodes.end()), first + left, count - left, label);
    }

    bool samePart(const int n, const int o) const
    {
        return part_.empty() || part_[n] == part_[o];
    }

    // With axis_ (per node: pillar direction t, layer normal n) an offset d splits as d = a t + (in-layer
    // part), a = d.n / t.n, so lines follow the pillars and planes the layers; otherwise t = n = z.
    void split(const int n, const int o, double& a, double& inlayer) const
    {
        double d[3], t[3] = {0, 0, 1}, m[3] = {0, 0, 1};
        for (int i = 0; i < 3; ++i)
            d[i] = coords_[3 * o + i] - coords_[3 * n + i];
        if (axis_) {
            for (int i = 0; i < 3; ++i) {
                t[i] = (*axis_)[6 * n + i] + (*axis_)[6 * o + i];
                m[i] = (*axis_)[6 * n + 3 + i] + (*axis_)[6 * o + 3 + i];
            }
            const double tl = std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
            const double ml = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
            if (tl > 0 && t[0] * m[0] + t[1] * m[1] + t[2] * m[2] > 1e-3 * tl * ml) {
                for (double& x : t)
                    x /= tl;
            } else {
                t[0] = t[1] = m[0] = m[1] = 0;
                t[2] = m[2] = 1;
            }
        }
        a = (d[0] * m[0] + d[1] * m[1] + d[2] * m[2]) / (t[0] * m[0] + t[1] * m[1] + t[2] * m[2]);
        double r2 = 0;
        for (int i = 0; i < 3; ++i)
            r2 += (d[i] - a * t[i]) * (d[i] - a * t[i]);
        inlayer = std::sqrt(r2);
    }

    double along(const int n, const int o) const
    {
        if (!axis_)
            return coords_[3 * o + 2] - coords_[3 * n + 2];
        double a, r;
        split(n, o, a, r);
        return a;
    }

    double dist(const int n, const int o) const
    {
        return std::hypot(coords_[3 * o] - coords_[3 * n], coords_[3 * o + 1] - coords_[3 * n + 1],
                          coords_[3 * o + 2] - coords_[3 * n + 2]);
    }

    bool lateral(const int n, const int o) const
    {
        if (!axis_) {
            const double dx = coords_[3 * o] - coords_[3 * n];
            const double dy = coords_[3 * o + 1] - coords_[3 * n + 1];
            return std::abs(coords_[3 * o + 2] - coords_[3 * n + 2]) <= std::hypot(dx, dy);
        }
        double a, r;
        split(n, o, a, r);
        return std::abs(a) <= r;
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
        if (!part_.empty())
            b.part = part_[b.dofs[0] / 3];
        max_band_ = std::max(max_band_, b.band);
        blocks_.push_back(std::move(b));
    }

    // z[b] += B^{-1} r[b]
    void solveAdd(const Block& b, const Vector& r, Vector& z)
    {
        const int m = int(b.dofs.size()), w = b.band + 1;
        y_.resize(m);
        for (int i = 0; i < m; ++i) {
            double s = r[b.dofs[i]];
            for (int k = std::max(0, i - b.band); k < i; ++k)
                s -= b.chol[i * w + k - i + b.band] * y_[k];
            y_[i] = s / b.chol[i * w + b.band];
        }
        for (int i = m - 1; i >= 0; --i) {
            for (int k = i + 1; k < std::min(m, i + w); ++k)
                y_[i] -= b.chol[k * w + i - k + b.band] * y_[k];
            y_[i] /= b.chol[i * w + b.band];
            z[b.dofs[i]] += y_[i];
        }
    }

    void sweep(const Block& b, Vector& x, const Vector& d)
    {
        const int m = int(b.dofs.size()), w = b.band + 1;
        y_.resize(m);
        for (int i = 0; i < m; ++i) {
            double s = d[b.dofs[i]];
            for (auto c = M_[b.dofs[i]].begin(); c != M_[b.dofs[i]].end(); ++c)
                s -= double(*c) * double(part_.empty() || part_[c.index() / 3] == b.part ? x[c.index()] : x0_[c.index()]);
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
    const std::vector<double>* axis_;
    std::vector<Block> blocks_;
    std::vector<char> owned_;
    std::vector<int> where_;
    std::vector<double> y_;
    Vector r_, c_;
    std::size_t plane_nodes_ = 0;
    int emulate_ranks_;
    std::string emulate_cut_;
    std::vector<int> part_;
    Vector x0_;
    int max_band_ = 0;
    bool additive_;
    int cheb_degree_;
    double cheb_ratio_;
    double lambda_max_ = 1.0;
    double line_ratio_;
    std::size_t line_nodes_ = 0;
    Vector res_, z_, p_;
};

} // namespace Opm

#endif // OPM_GEOMECH_MECH_LINE_PLANE_PRECONDITIONER_HPP
