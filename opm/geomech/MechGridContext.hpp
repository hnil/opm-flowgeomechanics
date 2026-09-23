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

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{

/// The mechanics grid when it is not the flow grid: a coarsening of the deck's
/// grid, built with the same corner points, plus the map to the flow cells.
///
/// Serial only so far. In parallel the two grids would need partitions that
/// correspond, which is not implemented.
class MechGridContext
{
public:
    /// COARSEN records, one per line: I1 I2 J1 J2 K1 K2 NX NY NZ, 1-based and
    /// inclusive. '--' and '#' start a comment.
    static std::vector<Coarsening::CoarsenRequest> readRecords(const std::string& filename)
    {
        std::ifstream is(filename);
        if (!is) {
            OPM_THROW(std::runtime_error, "Cannot open mechanics coarsening file " + filename);
        }

        std::vector<Coarsening::CoarsenRequest> requests;
        std::string line;
        while (std::getline(is, line)) {
            const auto comment = line.find_first_of("-#");
            if (comment != std::string::npos && line.compare(comment, 2, "--") == 0) {
                line.erase(comment);
            } else if (comment != std::string::npos && line[comment] == '#') {
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

    /// Build the mechanics grid from the deck's grid and the records, and the
    /// map from the flow grid's cells to its own.
    template <class FlowGrid, class CartesianMapper>
    MechGridContext(const EclipseGrid& inputGrid,
                    const std::vector<Coarsening::CoarsenRequest>& requests,
                    const FlowGrid& flowGrid,
                    const CartesianMapper& flowCartesianMapper)
    {
        if (flowGrid.comm().size() > 1) {
            OPM_THROW(std::runtime_error,
                      "A separate mechanics grid is not supported in parallel yet: the two "
                      "grids would need corresponding partitions.");
        }
        if (inputGrid.isPinchActive() || inputGrid.getMinpvMode() != MinpvMode::Inactive) {
            OPM_THROW(std::runtime_error,
                      "A separate mechanics grid needs the deck's geometry as flow sees it; "
                      "MINPV/PINCH decks need the edge-conformal processed grid to be "
                      "retained first, which is not implemented.");
        }

        Coarsening::Grdecl fine;
        fine.dims = {static_cast<int>(inputGrid.getNX()),
                     static_cast<int>(inputGrid.getNY()),
                     static_cast<int>(inputGrid.getNZ())};
        fine.coord = inputGrid.getCOORD();
        fine.zcorn = inputGrid.getZCORN();
        fine.actnum = inputGrid.getACTNUM();

        Coarsening::Options options;
        // Mechanics wants rock everywhere the block has volume, so inactive
        // cells with real volume are absorbed instead of left as holes.
        options.activity = Coarsening::Activity::FillHoles;
        const auto result = Coarsening::coarsenCornerPoint(fine, requests, options);

        for (const auto& note : result.report.notes) {
            OpmLog::info("Mechanics grid: " + note);
        }

        grdecl input;
        input.dims[0] = result.grid.dims[0];
        input.dims[1] = result.grid.dims[1];
        input.dims[2] = result.grid.dims[2];
        input.coord = result.grid.coord.data();
        input.zcorn = result.grid.zcorn.data();
        input.actnum = result.grid.actnum.data();
        grid_ = std::make_unique<Dune::CpGrid>();
        grid_->processEclipseFormat(input, /*remove_ij_boundary*/ false,
                                    /*turn_normals*/ false, /*edge_conformal*/ true);

        buildMap(result, flowGrid, flowCartesianMapper);
    }

    const Dune::CpGrid& grid() const { return *grid_; }
    Dune::CpGrid& grid() { return *grid_; }
    const MechFlowMap& map() const { return map_; }

private:
    template <class FlowGrid, class CartesianMapper>
    void buildMap(const Coarsening::Result& result, const FlowGrid& flowGrid,
                  const CartesianMapper& flowCartesianMapper)
    {
        const Dune::CartesianIndexMapper<Dune::CpGrid> mechMapper(*grid_);
        const auto& mechView = grid_->leafGridView();
        std::vector<int> cartesianToMech(static_cast<std::size_t>(result.grid.dims[0])
                                         * result.grid.dims[1] * result.grid.dims[2], -1);
        for (const auto& cell : elements(mechView)) {
            const int idx = mechView.indexSet().index(cell);
            cartesianToMech[mechMapper.cartesianIndex(idx)] = idx;
        }

        const auto& flowView = flowGrid.leafGridView();
        const std::size_t numFlowCells = flowView.size(0);
        std::vector<int> flowToMech(numFlowCells, -1);
        std::vector<double> weight(numFlowCells, 0.0);

        std::size_t unmapped = 0;
        for (const auto& cell : elements(flowView)) {
            const int f = flowView.indexSet().index(cell);
            weight[f] = cell.geometry().volume();
            const int coarseCartesian = result.fineToCoarse[flowCartesianMapper.cartesianIndex(f)];
            const int m = (coarseCartesian < 0) ? -1 : cartesianToMech[coarseCartesian];
            flowToMech[f] = m;
            unmapped += (m < 0) ? 1 : 0;
        }
        if (unmapped > 0) {
            OPM_THROW(std::runtime_error,
                      std::to_string(unmapped) + " flow cells have no mechanics cell; the "
                      "mechanics grid must cover the flow grid.");
        }

        map_ = MechFlowMap(std::move(flowToMech), std::move(weight), mechView.size(0));

        std::ostringstream os;
        os << "Mechanics grid: " << mechView.size(0) << " cells for " << numFlowCells
           << " flow cells";
        OpmLog::info(os.str());
    }

    std::unique_ptr<Dune::CpGrid> grid_;
    MechFlowMap map_;
};

} // namespace Opm

#endif // OPM_MECH_GRID_CONTEXT_HPP
