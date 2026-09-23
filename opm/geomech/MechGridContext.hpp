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
#ifndef OPM_MECH_GRID_CONTEXT_HPP
#define OPM_MECH_GRID_CONTEXT_HPP

#include <opm/geomech/MechFlowMap.hpp>
#include <opm/geomech/MechGridPartition.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <fstream>
#include <type_traits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{

/// The mechanics grid when it is not the flow grid: a coarsening of the deck's
/// grid, built with the same corner points, plus the map to the flow cells.
///
/// In parallel the flow grid must have been partitioned by
/// MechPartition::flowPartition, which puts every flow cell on the rank of its
/// coarse cell; the mechanics grid is then distributed the same way, so a
/// mechanics cell and all of its flow cells share a rank.
class MechGridContext
{
public:
    static std::vector<Coarsening::CoarsenRequest> readRecords(const std::string& filename)
    { return MechPartition::readRecords(filename); }

    /// Build the mechanics grid from the deck's grid and the records, and the
    /// map from the flow grid's cells to its own.
    template <class FlowGrid, class CartesianMapper>
    MechGridContext(const EclipseGrid* inputGrid,
                    const std::vector<Coarsening::CoarsenRequest>& requests,
                    const FlowGrid& flowGrid,
                    const CartesianMapper& flowCartesianMapper)
    {
        if (flowGrid.comm().size() > 1) {
            // The balancer runs on rank 0 only, so ask it there and let every
            // rank throw together: a throw on one rank alone would leave the
            // others in a collective call.
            int haveIt = (flowGrid.comm().rank() != 0) || !MechPartition::Registry::empty();
            haveIt = flowGrid.comm().min(haveIt);
            if (!haveIt) {
                OPM_THROW_NOLOG(std::runtime_error,
                                "In parallel the flow grid must be partitioned by the "
                                "mechanics coarsening: install "
                                "MechPartition::installFlowPartition in main() before the "
                                "simulator is created.");
            }
        }
        // Only rank 0 holds the deck's grid, so it decides and tells the others.
        int pinchOrMinpv = 0;
        if (inputGrid != nullptr) {
            pinchOrMinpv = (inputGrid->isPinchActive()
                            || inputGrid->getMinpvMode() != MinpvMode::Inactive) ? 1 : 0;
        }
        pinchOrMinpv = flowGrid.comm().max(pinchOrMinpv);
        if (pinchOrMinpv != 0) {
            OPM_THROW_NOLOG(std::runtime_error,
                            "A separate mechanics grid needs the deck's geometry as flow sees "
                            "it; MINPV/PINCH decks need the edge-conformal processed grid to "
                            "be retained first, which is not implemented.");
        }

        const auto size = flowGrid.logicalCartesianSize();
        const std::array<int,3> fineDims{size[0], size[1], size[2]};
        // The index side needs no geometry, so every rank can have it.
        const auto cartesian = Coarsening::cartesianMap(fineDims, requests);

        // Corner-point processing is a rank-0 job with collective steps, so
        // every rank calls it and only rank 0 passes the description.
        std::unique_ptr<EclipseGrid> coarseGrid;
        if (inputGrid != nullptr) {
            Coarsening::Grdecl fine;
            fine.dims = fineDims;
            fine.coord = inputGrid->getCOORD();
            fine.zcorn = inputGrid->getZCORN();
            fine.actnum = inputGrid->getACTNUM();

            Coarsening::Options options;
            // Mechanics wants rock everywhere the block has volume, so inactive
            // cells with real volume are absorbed instead of left as holes.
            options.activity = Coarsening::Activity::FillHoles;
            const auto result = Coarsening::coarsenCornerPoint(fine, requests, options);
            for (const auto& note : result.report.notes) {
                OpmLog::info("Mechanics grid: " + note);
            }
            coarseGrid = std::make_unique<EclipseGrid>(result.grid.dims, result.grid.coord,
                                                       result.grid.zcorn,
                                                       result.grid.actnum.data());
        }

        grid_ = std::make_unique<Dune::CpGrid>(flowGrid.comm());
        grid_->processEclipseFormat(coarseGrid.get(), /*ecl_state*/ nullptr,
                                    /*periodic_extension*/ false, /*turn_normals*/ false,
                                    /*clip_z*/ false, /*pinchActive*/ false,
                                    /*edge_conformal*/ true);

        if (flowGrid.comm().size() > 1) {
            distribute();
        }
        buildMap(cartesian, flowGrid, flowCartesianMapper);
    }

    const Dune::CpGrid& grid() const { return *grid_; }
    Dune::CpGrid& grid() { return *grid_; }
    const MechFlowMap& map() const { return map_; }

    /// Fill the overlap cells of a mechanics-cell vector from their owners.
    /// Restriction only fills cells whose flow children are local.
    template <class Vector>
    void communicate(Vector& mechCellData) const
    {
        if (grid_->comm().size() == 1) {
            return;
        }
        const auto& indexSet = grid_->leafGridView().indexSet();
        VectorHandle<Vector, std::decay_t<decltype(indexSet)>> handle(mechCellData, indexSet);
        grid_->communicate(handle, Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication);
    }

private:
    /// Copies one value per cell, owner to copy.
    template <class Vector, class IndexSet>
    class VectorHandle
    {
    public:
        VectorHandle(Vector& data, const IndexSet& indexSet)
            : data_(data), index_set_(indexSet) {}

        using DataType = double;
        bool contains(int /*dim*/, int codim) const { return codim == 0; }
        bool fixedSize(int /*dim*/, int /*codim*/) const { return true; }
        template <class Entity> std::size_t size(const Entity&) const { return 1; }

        template <class Buffer, class Entity>
        void gather(Buffer& buffer, const Entity& e) const
        { buffer.write(value(data_[index_set_.index(e)])); }

        template <class Buffer, class Entity>
        void scatter(Buffer& buffer, const Entity& e, std::size_t)
        {
            double v;
            buffer.read(v);
            data_[index_set_.index(e)] = v;
        }

    private:
        template <class T>
        static double value(const T& v)
        {
            if constexpr (std::is_arithmetic_v<T>) {
                return v;
            } else {
                return v[0];
            }
        }

        Vector& data_;
        const IndexSet& index_set_;
    };

    void distribute()
    {
        std::vector<int> parts;
        if (grid_->comm().rank() == 0) {
            const Dune::CartesianIndexMapper<Dune::CpGrid> mapper(*grid_);
            const auto& partOfCoarseCell = MechPartition::Registry::get();
            parts.resize(grid_->leafGridView().size(0));
            for (std::size_t c = 0; c < parts.size(); ++c) {
                parts[c] = partOfCoarseCell[mapper.cartesianIndex(static_cast<int>(c))];
            }
        }
        grid_->loadBalance(parts, /*ownersFirst*/ false, /*addCornerCells*/ true,
                           /*overlapLayers*/ 1);
    }

    template <class FlowGrid, class CartesianMapper>
    void buildMap(const Coarsening::CartesianMap& cartesian, const FlowGrid& flowGrid,
                  const CartesianMapper& flowCartesianMapper)
    {
        const Dune::CartesianIndexMapper<Dune::CpGrid> mechMapper(*grid_);
        const auto& mechView = grid_->leafGridView();
        std::vector<int> cartesianToMech(static_cast<std::size_t>(cartesian.coarseDims[0])
                                         * cartesian.coarseDims[1]*cartesian.coarseDims[2], -1);
        for (const auto& cell : elements(mechView)) {
            const int idx = mechView.indexSet().index(cell);
            cartesianToMech[mechMapper.cartesianIndex(idx)] = idx;
        }

        const auto& flowView = flowGrid.leafGridView();
        const std::size_t numFlowCells = flowView.size(0);
        std::vector<int> flowToMech(numFlowCells, -1);
        std::vector<double> weight(numFlowCells, 0.0);

        std::size_t unmappedInterior = 0, unmappedOverlap = 0;
        for (const auto& cell : elements(flowView)) {
            const int f = flowView.indexSet().index(cell);
            weight[f] = cell.geometry().volume();
            const int coarseCartesian = cartesian.fineToCoarse[flowCartesianMapper.cartesianIndex(f)];
            const int m = (coarseCartesian < 0) ? -1 : cartesianToMech[coarseCartesian];
            flowToMech[f] = m;
            if (m < 0) {
                if (cell.partitionType() == Dune::InteriorEntity) {
                    ++unmappedInterior;
                } else {
                    ++unmappedOverlap;
                }
            }
        }
        if (unmappedInterior > 0) {
            // In parallel this means the partitions disagree: an interior flow
            // cell must always be on the rank that owns its mechanics cell.
            OPM_THROW(std::runtime_error,
                      std::to_string(unmappedInterior) + " interior flow cells have no "
                      "mechanics cell on this rank; the mechanics grid must cover the flow "
                      "grid, and the two partitions must correspond.");
        }
        if (unmappedOverlap > 0) {
            OpmLog::info("Mechanics grid: " + std::to_string(unmappedOverlap)
                         + " overlap flow cells have no local mechanics cell; they take no "
                         "part in the mechanics.");
        }

        map_ = MechFlowMap(std::move(flowToMech), std::move(weight), mechView.size(0));

        std::ostringstream os;
        os << "Mechanics grid: " << mechView.size(0) << " cells for " << numFlowCells
           << " flow cells on this rank";
        OpmLog::info(os.str());
    }

    std::unique_ptr<Dune::CpGrid> grid_;
    MechFlowMap map_;
};

} // namespace Opm

#endif // OPM_MECH_GRID_CONTEXT_HPP
