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
#ifndef OPM_MECH_FLOW_MAP_HPP
#define OPM_MECH_FLOW_MAP_HPP

#include <cassert>
#include <cstddef>
#include <vector>

namespace Opm
{

/// Which mechanics cell each flow cell belongs to, for a mechanics grid that
/// is a coarsening of the flow grid: a mechanics cell is the union of its flow
/// cells, so restriction is a volume-weighted mean and prolongation is a copy.
///
/// The default is the identity — one mechanics grid cell per flow cell — and
/// then every caller keeps its flow indexing and nothing is transferred.
class MechFlowMap
{
public:
    MechFlowMap() = default;

    /// Mechanics and flow share the grid.
    static MechFlowMap identity(std::size_t numFlowCells)
    {
        MechFlowMap map;
        map.num_mech_cells_ = numFlowCells;
        return map;
    }

    /// flowToMech[flowCell] is the mechanics cell it belongs to; weights are
    /// the flow cells' bulk volumes.
    MechFlowMap(std::vector<int> flowToMech, std::vector<double> weights,
                std::size_t numMechCells)
        : identity_(false)
        , num_mech_cells_(numMechCells)
        , flow_to_mech_(std::move(flowToMech))
        , weight_(std::move(weights))
        , mech_weight_(numMechCells, 0.0)
    {
        assert(flow_to_mech_.size() == weight_.size());
        for (std::size_t f = 0; f < flow_to_mech_.size(); ++f) {
            if (flow_to_mech_[f] >= 0) {
                mech_weight_[flow_to_mech_[f]] += weight_[f];
            }
        }
    }

    bool isIdentity() const { return identity_; }
    std::size_t numMechCells() const { return num_mech_cells_; }
    std::size_t numFlowCells() const
    { return identity_ ? num_mech_cells_ : flow_to_mech_.size(); }

    /// The mechanics cell of a flow cell, or -1 if the flow cell has no
    /// mechanics counterpart (a smaller mechanics domain).
    int mechCell(std::size_t flowCell) const
    { return identity_ ? static_cast<int>(flowCell) : flow_to_mech_[flowCell]; }

    /// Volume-weighted mean over each mechanics cell's flow cells, which
    /// preserves the integral of a force density.
    template <class Vector>
    void restrict(const Vector& flow, Vector& mech) const
    {
        if (identity_) {
            mech = flow;
            return;
        }
        mech.resize(num_mech_cells_);
        mech = 0.0;
        for (std::size_t f = 0; f < flow_to_mech_.size(); ++f) {
            const int m = flow_to_mech_[f];
            if (m >= 0) {
                mech[m] += flow[f]*weight_[f];
            }
        }
        for (std::size_t m = 0; m < num_mech_cells_; ++m) {
            if (mech_weight_[m] > 0.0) {
                mech[m] /= mech_weight_[m];
            }
        }
    }

    /// Piecewise constant: every flow cell takes its mechanics cell's value.
    template <class Vector>
    void prolong(const Vector& mech, Vector& flow) const
    {
        if (identity_) {
            flow = mech;
            return;
        }
        flow.resize(flow_to_mech_.size());
        for (std::size_t f = 0; f < flow_to_mech_.size(); ++f) {
            const int m = flow_to_mech_[f];
            if (m >= 0) {
                flow[f] = mech[m];
            }
        }
    }

private:
    bool identity_{true};
    std::size_t num_mech_cells_{0};
    std::vector<int> flow_to_mech_;
    std::vector<double> weight_;
    std::vector<double> mech_weight_;
};

} // namespace Opm

#endif // OPM_MECH_FLOW_MAP_HPP
