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

#ifndef OPM_GEOMECH_MECH_OPTIONS_HPP
#define OPM_GEOMECH_MECH_OPTIONS_HPP

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <stdexcept>
#include <string>

namespace Opm {

//! Reservoir mechanics options, read once from the fracture parameter tree:
//! a key in its `mech` block wins over the same key at the top level.
struct MechOptions
{
    int vem_stability_choice = 1;
    bool vem_source = true;
    bool vem_stress = true;
    bool stab_on_stress = false;
    bool diagonal_scaling = false;
    std::string face_quadrature = "kite";
    bool smooth_force = false;
    double smooth_force_length = 0.0;
    bool patch_recovery_stress = false;

    static MechOptions read(const PropertyTree& param)
    {
        const auto mech = param.get_child_optional("mech");
        auto get = [&](const std::string& key, auto fallback) {
            using T = decltype(fallback);
            const T top = param.get<T>(key, fallback);
            return mech ? mech->template get<T>(key, top) : top;
        };
        MechOptions o;
        if (get("reduce_boundary", false)) {
            throw std::invalid_argument("reduce_boundary=true was removed; "
                                        "Dirichlet nodes are eliminated symmetrically");
        }
        o.vem_stability_choice = get("vem_stability_choice", o.vem_stability_choice);
        o.vem_source = get("vem_source", o.vem_source);
        o.vem_stress = get("vem_stress", o.vem_stress);
        o.stab_on_stress = get("stab_on_stress", o.stab_on_stress);
        o.diagonal_scaling = get("mech_diagonal_scaling", o.diagonal_scaling);
        o.face_quadrature = get("vem_face_quadrature", o.face_quadrature);
        o.smooth_force = get("smooth_force", o.smooth_force);
        o.smooth_force_length = get("smooth_force_length", o.smooth_force_length);
        o.patch_recovery_stress = get("patch_recovery_stress", o.patch_recovery_stress);
        return o;
    }
};

} // namespace Opm

#endif
