#!/usr/bin/env python3
# Copyright 2026 Equinor ASA.
#
# This file is part of the Open Porous Media project (OPM).  OPM is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.  See
# <http://www.gnu.org/licenses/>.
"""The LGR-inner layout: one uniform grid, refined only where it matters.

A base grid of equal, well-shaped cells covers the whole domain and the
reservoir is a refined box inside it — the arrangement a CARFIN would give,
without the cross that a padded tensor grid leaves behind.

Writes T6_LGR_INNER.DATA (flow on the refined resolution everywhere) and
mech_coarsen_T6_level0.txt, which merges everything outside the box back to
base cells: the mechanics grid is then the base grid outside and the refined
grid inside, i.e. what mechanics on level zero plus a refined region would be.

The mechanics grid is built by merging cells, since a box refined in part of a
column cannot be written as COORD/ZCORN.
"""

import os
import sys

import make_coarsen_tests as base

HERE = os.path.dirname(os.path.abspath(__file__))

# Base cells are REFINE times the refined ones in each direction.
REFINE = 5
FINE_DX = 200.0
FINE_DZ = 100.0
NX = NY = 30                     # 6 base columns each way
NZ = 30                          # 6 base layers
TOP = 500.0

# The refined box, in fine indices (1-based, inclusive).
BOX_I = (11, 20)
BOX_J = (11, 20)
BOX_K = (11, 20)
RES_K = (13, 18)                 # reservoir layers inside the box
WELL = (15, 15, 15)


def layout():
    xs = base.cumulative(0.0, [FINE_DX]*NX)
    ys = base.cumulative(0.0, [FINE_DX]*NY)
    zs = base.cumulative(TOP, [FINE_DZ]*NZ)
    return base.Layout(xs, ys, zs,
                       core_i=(1, NX), core_j=(1, NY), core_k=(1, NZ),
                       res_k=RES_K, well=WELL)


def records():
    """Coarsen everything outside the box by REFINE, leave the box alone."""
    def split(n, lo, hi):
        """Ranges before, inside and after the box, as (i1, i2) 1-based."""
        out = []
        if lo > 1:
            out.append((1, lo - 1))
        out.append((lo, hi))
        if hi < n:
            out.append((hi + 1, n))
        return out

    groups = []
    for i1, i2 in split(NX, *BOX_I):
        for j1, j2 in split(NY, *BOX_J):
            for k1, k2 in split(NZ, *BOX_K):
                inside = (i1, i2) == BOX_I and (j1, j2) == BOX_J and (k1, k2) == BOX_K
                if inside:
                    continue
                n = [i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1]
                # as close to base cells as the range allows
                counts = [max(1, c//REFINE) for c in n]
                groups.append((i1, i2, j1, j2, k1, k2) + tuple(counts))
    return groups


def main():
    with open(base.BASE_DECK) as f:
        deck = f.read()

    lay = layout()
    grid = "T6_LGR_INNER_GRID.INC"
    base.write_grid_include(os.path.join(HERE, grid), lay)
    with open(os.path.join(HERE, "T6_LGR_INNER.DATA"), "w") as f:
        f.write(base.make_deck(deck, lay, grid, 3))

    with open(os.path.join(HERE, "mech_coarsen_T6_level0.txt"), "w") as f:
        f.write("-- Mechanics on the base grid outside the refined box, refined\n"
                "-- resolution inside it. Needs the merge route: a box refined in\n"
                "-- part of a column cannot be written as COORD/ZCORN.\n"
                "-- I1 I2 J1 J2 K1 K2 NX NY NZ\n")
        for r in records():
            f.write(" ".join("%3d" % x for x in r) + "\n")

    cells = NX*NY*NZ
    print("T6_LGR_INNER %dx%dx%d = %d cells, refined box i%s j%s k%s, well %s"
          % (NX, NY, NZ, cells, BOX_I, BOX_J, WELL[2:], WELL))
    print("mech_coarsen_T6_level0.txt: %d records" % len(records()))


if __name__ == "__main__":
    sys.exit(main())
