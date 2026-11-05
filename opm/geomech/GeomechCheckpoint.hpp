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
#ifndef OPM_GEOMECH_CHECKPOINT_HPP
#define OPM_GEOMECH_CHECKPOINT_HPP

#include <cstddef>
#include <vector>

namespace Opm::GeomechCheckpoint
{

//! Serialize a container of Dune::FieldVector, which the Serializer has no packer for.
template <class Serializer, class Container>
void fieldVectors(Serializer& serializer, Container& values)
{
    using Block = typename Container::value_type;
    constexpr std::size_t N = Block::dimension;

    std::vector<double> flat;
    if (serializer.isSerializing()) {
        flat.reserve(values.size() * N);
        for (const auto& block : values) {
            for (std::size_t i = 0; i < N; ++i) {
                flat.push_back(block[i]);
            }
        }
    }
    serializer(flat);
    if (!serializer.isSerializing()) {
        values.resize(flat.size() / N);
        std::size_t pos = 0;
        for (auto& block : values) {
            for (std::size_t i = 0; i < N; ++i) {
                block[i] = flat[pos++];
            }
        }
    }
}

} // namespace Opm::GeomechCheckpoint

#endif
