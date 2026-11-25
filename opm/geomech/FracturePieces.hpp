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
#ifndef OPM_FRACTURE_PIECES_HPP
#define OPM_FRACTURE_PIECES_HPP

#include <cstddef>
#include <vector>

namespace Opm
{
//! Cells of the face-connected pieces that contain no fed cell and no cell with positive leak-off.
inline std::vector<int>
cellsInPiecesWithoutSink(const std::vector<std::vector<int>>& neighbours,
                         const std::vector<bool>& fed,
                         const std::vector<double>& leakoff)
{
    const std::size_t nc = neighbours.size();
    std::vector<bool> seen(nc, false);
    std::vector<int> result;
    for (std::size_t start = 0; start < nc; ++start) {
        if (seen[start]) {
            continue;
        }
        std::vector<int> cells {static_cast<int>(start)};
        seen[start] = true;
        bool has_sink = false;
        for (std::size_t k = 0; k < cells.size(); ++k) {
            const int c = cells[k];
            has_sink = has_sink || fed[c] || leakoff[c] > 0.0;
            for (const int n : neighbours[c]) {
                if (!seen[n]) {
                    seen[n] = true;
                    cells.push_back(n);
                }
            }
        }
        if (!has_sink) {
            result.insert(result.end(), cells.begin(), cells.end());
        }
    }
    return result;
}
} // namespace Opm

#endif // OPM_FRACTURE_PIECES_HPP
