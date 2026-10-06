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
#include <opm/grid/cpgrid/RetainedCornerPointInput.hpp>
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
    /// `processed` is what the flow grid was built from (MINPV and PINCH
    /// included), and exists on rank 0 only.
    enum class Method { Auto, Grdecl, Merge, Collapse };

    static Method method(const std::string& name)
    {
        if (name == "auto")   { return Method::Auto; }
        if (name == "grdecl") { return Method::Grdecl; }
        if (name == "merge")  { return Method::Merge; }
        if (name == "collapse") { return Method::Collapse; }
        OPM_THROW(std::runtime_error,
                  "Unknown mechanics coarsening method '" + name
                  + "'; use grdecl, merge, collapse or auto");
    }

    template <class FlowGrid, class CartesianMapper>
    MechGridContext(const RetainedCornerPointInput* processed,
                    const std::vector<Coarsening::CoarsenRequest>& requests,
                    const FlowGrid& flowGrid,
                    const CartesianMapper& flowCartesianMapper,
                    const Method wanted = Method::Auto)
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
        // Only rank 0 holds the processed description, so it decides for all.
        const bool isRoot = flowGrid.comm().rank() == 0;
        int problem = 0;
        if (isRoot) {
            problem = (processed == nullptr) ? 1 : 0;
        }
        problem = flowGrid.comm().max(problem);
        if (problem != 0) {
            OPM_THROW_NOLOG(std::runtime_error,
                            "The mechanics grid needs the description the flow grid was built "
                            "from: call RetainCornerPointInput::enable() before the grid is "
                            "created (MechPartition::installFlowPartition does).");
        }

        const auto size = flowGrid.logicalCartesianSize();
        const std::array<int,3> fineDims{size[0], size[1], size[2]};

        // Corner-point processing is a rank-0 job with collective steps, so
        // every rank calls it and only rank 0 passes the description.
        // The grdecl route gives clean six-faced cells, but pillars run
        // through every layer, so it cannot coarsen part of a column. The
        // merge can, at the price of a coarse cell keeping all the faces it
        // has towards finer neighbours.
        int useMerge = (wanted == Method::Merge || wanted == Method::Collapse) ? 1 : 0;
        if (isRoot && wanted == Method::Auto) {
            try {
                Coarsening::cartesianMap(fineDims, requests);
            } catch (const std::invalid_argument& e) {
                OpmLog::info(std::string("Mechanics grid: cannot be written as a corner-point "
                                         "description (") + e.what()
                             + "); merging cells of the flow grid instead.");
                useMerge = 1;
            }
        }
        useMerge = flowGrid.comm().max(useMerge);
        merged_ = useMerge != 0;

        // The index side needs no geometry, so every rank can have it.
        Coarsening::CartesianMap cartesian;
        if (!merged_) {
            cartesian = Coarsening::cartesianMap(fineDims, requests);
        }

        std::unique_ptr<EclipseGrid> coarseGrid;
        if (processed != nullptr && !merged_) {
            Coarsening::Grdecl fine;
            fine.dims = fineDims;
            fine.coord = processed->coord;
            fine.zcorn = processed->zcorn;
            fine.actnum = processed->actnum;

            Coarsening::Options options;
            // Mechanics wants rock everywhere the block has volume, so cells
            // that are inactive or that MINPV removed are absorbed instead of
            // left as holes.
            options.activity = Coarsening::Activity::FillHoles;
            // Without edge-conformal processing, MINPV leaves the removed
            // cells' volume as a gap between their neighbours; the mechanics
            // body takes it back as rock.
            options.allowVerticalGaps = !processed->edgeConformal;
            const auto result = Coarsening::coarsenCornerPoint(fine, requests, options);
            for (const auto& note : result.report.notes) {
                OpmLog::info("Mechanics grid: " + note);
            }
            if (result.report.absorbedGapVolume > 0.0) {
                OpmLog::warning("Mechanics grid: took "
                                + std::to_string(result.report.absorbedGapVolume)
                                + " m3 of gaps left by cell removal back as rock. Run with "
                                "--edge-conformal=true so the removal merges the cells "
                                "geometrically; gaps under cells the mechanics does not "
                                "coarsen stay holes in the body.");
            }
            coarseGrid = std::make_unique<EclipseGrid>(result.grid.dims, result.grid.coord,
                                                       result.grid.zcorn,
                                                       result.grid.actnum.data());
        }

        grid_ = std::make_unique<Dune::CpGrid>(flowGrid.comm());
        if (merged_) {
            grdecl input{};
            Coarsening::BlockLayout layout;
            if (processed != nullptr) {
                input.dims[0] = fineDims[0];
                input.dims[1] = fineDims[1];
                input.dims[2] = fineDims[2];
                input.coord = processed->coord.data();
                input.zcorn = processed->zcorn.data();
                input.actnum = processed->actnum.data();
                layout = Coarsening::blockLayout(fineDims, requests);
                std::ostringstream os;
                os << "Mechanics grid: merging " << fineDims[0]*fineDims[1]*fineDims[2]
                   << " cells into " << layout.boxes.size() << " blocks";
                OpmLog::info(os.str());
            }
            grid_->processEclipseFormatCoarsened(input, layout.blockOfCartesian, layout.boxes,
                                                 /*edge_conformal*/ true,
                                                 /*collapse_coarse_faces*/ wanted == Method::Collapse);
        } else {
            grid_->processEclipseFormat(coarseGrid.get(), /*ecl_state*/ nullptr,
                                        /*periodic_extension*/ false, /*turn_normals*/ false,
                                        /*clip_z*/ false, /*pinchActive*/ false,
                                        /*edge_conformal*/ true);
        }

        if (flowGrid.comm().size() > 1) {
            // The balancer gave every Cartesian cell a rank; the merge keeps
            // that index space, the grdecl route has its own.
            const auto& partOfFine = MechPartition::Registry::get();
            std::vector<int> partOfMech;
            if (merged_) {
                partOfMech = partOfFine;
            } else if (flowGrid.comm().rank() == 0) {
                partOfMech.assign(static_cast<std::size_t>(cartesian.coarseDims[0])
                                  * cartesian.coarseDims[1]*cartesian.coarseDims[2], 0);
                for (std::size_t f = 0; f < cartesian.fineToCoarse.size(); ++f) {
                    if (cartesian.fineToCoarse[f] >= 0) {
                        partOfMech[cartesian.fineToCoarse[f]] = partOfFine[f];
                    }
                }
            }
            distribute(partOfMech);
        }
        if (merged_) {
            // Anchor index of each cell's block, in the fine index space.
            const auto layout = Coarsening::blockLayout(fineDims, requests);
            std::vector<int> anchorOfCartesian(layout.blockOfCartesian.size(), -1);
            for (std::size_t c = 0; c < layout.blockOfCartesian.size(); ++c) {
                const auto& box = layout.boxes[layout.blockOfCartesian[c]];
                anchorOfCartesian[c] = box[0] + fineDims[0]*(box[1] + fineDims[1]*box[2]);
            }
            buildMap(anchorOfCartesian, fineDims, flowGrid, flowCartesianMapper);
        } else {
            buildMap(cartesian.fineToCoarse, cartesian.coarseDims, flowGrid, flowCartesianMapper);
        }
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

    /// `partOfMechCartesian` gives the rank of each cell of the mechanics
    /// grid's own Cartesian space.
    void distribute(const std::vector<int>& partOfMechCartesian)
    {
        std::vector<int> parts;
        if (grid_->comm().rank() == 0) {
            const Dune::CartesianIndexMapper<Dune::CpGrid> mapper(*grid_);
            parts.resize(grid_->leafGridView().size(0));
            for (std::size_t c = 0; c < parts.size(); ++c) {
                parts[c] = partOfMechCartesian[mapper.cartesianIndex(static_cast<int>(c))];
            }
        }
        grid_->loadBalance(parts, /*ownersFirst*/ false, /*addCornerCells*/ true,
                           /*overlapLayers*/ 1);
    }

    /// `toMech` takes a fine Cartesian index to the mechanics grid's own
    /// Cartesian index: the coarse index for the grdecl route, the block's
    /// anchor for the merge, which keeps the fine index space.
    template <class FlowGrid, class CartesianMapper, class ToMech>
    void buildMap(const ToMech& toMech, const std::array<int,3>& mechDims,
                  const FlowGrid& flowGrid, const CartesianMapper& flowCartesianMapper)
    {
        const Dune::CartesianIndexMapper<Dune::CpGrid> mechMapper(*grid_);
        const auto& mechView = grid_->leafGridView();
        std::vector<int> cartesianToMech(static_cast<std::size_t>(mechDims[0])
                                         * mechDims[1]*mechDims[2], -1);
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
            const int coarseCartesian = toMech[flowCartesianMapper.cartesianIndex(f)];
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
    bool merged_{false};
};

} // namespace Opm

#endif // OPM_MECH_GRID_CONTEXT_HPP
