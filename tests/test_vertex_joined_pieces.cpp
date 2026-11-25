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

#include <config.h>
#include <opm/geomech/RegularTrimesh.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}
} // namespace

// A piece joined to the rest at one node must be its own face-connected piece in the Dune grid,
// with every cell an outer cell: the case Fracture::restoreLeakoffOfIsolatedPieces_ handles.
int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    try {
        std::vector<CellRef> cells;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cells.push_back({i, j, 0});
                cells.push_back({i, j, 1});
            }
        }
        cells.push_back({4, 4, 0}); // shares only node (4,4) with cell {3,3,1}
        cells.push_back({4, 4, 1});
        const RegularTrimesh mesh(cells.begin(), cells.end());
        for (const CellRef n : {CellRef {3, 4, 1}, CellRef {4, 3, 1}, CellRef {4, 4, 0}}) {
            require(n == CellRef {4, 4, 0} || !mesh.isActive(n), "the piece has no edge neighbour in the block");
        }

        auto [grid, cellmap, bmap] = mesh.createDuneGrid(0, {}, false);
        const auto gv = grid->leafGridView();
        Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView> mapper(gv, Dune::mcmgElementLayout());
        const std::size_t nc = gv.size(0);
        require(nc == cells.size(), "no smoothing or coarsening: one Dune cell per trimesh cell");

        std::vector<std::vector<int>> nbs(nc);
        std::vector<bool> on_boundary(nc, false);
        for (const auto& e : elements(gv)) {
            const int c = mapper.index(e);
            for (const auto& is : intersections(gv, e)) {
                if (is.neighbor()) {
                    nbs[c].push_back(mapper.index(is.outside()));
                } else {
                    on_boundary[c] = true;
                }
            }
        }
        // Face-connected pieces, the same search as Fracture::restoreLeakoffOfIsolatedPieces_.
        std::vector<int> piece(nc, -1);
        std::vector<std::size_t> sizes;
        for (std::size_t s = 0; s < nc; ++s) {
            if (piece[s] >= 0) {
                continue;
            }
            std::vector<int> todo {static_cast<int>(s)};
            piece[s] = static_cast<int>(sizes.size());
            for (std::size_t k = 0; k < todo.size(); ++k) {
                for (const int n : nbs[todo[k]]) {
                    if (piece[n] < 0) {
                        piece[n] = piece[s];
                        todo.push_back(n);
                    }
                }
            }
            sizes.push_back(todo.size());
        }
        require(sizes.size() == 2, "a vertex-only join gives two face-connected pieces");
        const std::size_t small = sizes[0] < sizes[1] ? 0 : 1;
        require(sizes[small] == 2 && sizes[1 - small] == 32, "pieces of 2 and 32 cells");
        for (std::size_t c = 0; c < nc; ++c) {
            if (piece[c] == static_cast<int>(small)) {
                require(on_boundary[c], "every cell of the small piece is an outer cell");
            }
        }
        std::cout << "pieces: " << sizes[0] << " + " << sizes[1] << " cells, small piece all on the boundary\n";
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
