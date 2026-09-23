#!/usr/bin/env python3
# Copyright 2026 Equinor ASA.
#
# This file is part of the Open Porous Media project (OPM).  OPM is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.  See
# <http://www.gnu.org/licenses/>.
"""Generate the padded SIMPLE decks used by the corner-point coarsening work.

Writes, next to this script:

  T1_PAD_FINE.DATA      padded SIMPLE, fine everywhere
  T1_PAD_COARSE.DATA    the same model with the padding merged (the exact
                        reference the coarsening tool must reproduce)
  T2_CAPROCK.DATA       T1_PAD_FINE with thin low-perm layers on top of the
                        reservoir (flow keeps them, mechanics coarsens them away)
  coarsen_spec_T1.json  fine -> coarse description, as axis groups and as
                        COARSEN records
  coarsen_spec_T2.json  the mechanics-only coarsening of T2

Physics (PROPS/SOLUTION/SUMMARY/SCHEDULE) is taken verbatim from the base deck,
so only geometry and the index-dependent keywords are generated here.

Plan: opm-gridrefined/docs/COARSENING-FOR-MECHANICS.md
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# Works with both the BCPROP and the BCMECH form of the base deck: only BCCON
# is rewritten here, and that keyword is the same in both.
BASE_DECK = os.path.join(HERE, "..", "SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA")

# --- the model ---------------------------------------------------------------

CORE_NX = CORE_NY = 11
CORE_DX = CORE_DY = 181.0
CORE_NZ = 25
CORE_DZ = 20.0
CORE_TOP = 1900.0

RES_K1, RES_K2 = 11, 20          # reservoir layers, core-local 1-based
WELL_I = WELL_J = 6              # core-local
WELL_K = 16

# Side padding: fine cell widths outward from the core, merged in pairs for the
# coarse reference.  Growth factor 2 keeps the outer boundary far (~5.4 km)
# without many cells.
PAD_LATERAL = [362.0, 724.0, 1448.0, 2896.0]
# Over- and underburden: fine layer thicknesses away from the core, likewise
# merged in pairs.  The overburden reaches the surface.
PAD_ABOVE = [20.0, 50.0, 100.0, 200.0, 400.0, 1130.0]   # sums to CORE_TOP
PAD_BELOW = [20.0, 50.0, 100.0, 200.0, 400.0, 1130.0]

# How fine cells are merged for the coarse reference: pairs in the padding,
# untouched in the core.
PAD_GROUPING = 2

# T2: the bottom of the seal, immediately above the reservoir, is resolved as
# thin low-permeability layers.  3 x 0.5 m at 181 m lateral size is aspect ~360.
CAPROCK_N = 3
CAPROCK_DZ = 0.5

PROPS_ROCK = {"PERMX": 10.0, "PORO": 0.10}
PROPS_RES = {"PERMX": 1000.0, "PORO": 0.28}
PROPS_PAD = {"PERMX": 0.001, "PORO": 0.01}
PROPS_CAPROCK = {"PERMX": 0.0001, "PORO": 0.05}


def cumulative(start, deltas):
    """Breakpoints from a start value and a list of increments."""
    out = [start]
    for d in deltas:
        out.append(out[-1] + d)
    return out


class Layout:
    """Tensor grid: coordinate breakpoints plus where the core sits in them.

    All index ranges are 1-based and inclusive, as in the deck.
    """

    def __init__(self, xs, ys, zs, core_i, core_j, core_k, res_k, well,
                 caprock_k=None):
        self.xs, self.ys, self.zs = xs, ys, zs
        self.nx, self.ny, self.nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
        self.core_i, self.core_j, self.core_k = core_i, core_j, core_k
        self.res_k = res_k
        self.well = well
        self.caprock_k = caprock_k


def build_fine(with_caprock=False):
    npad = len(PAD_LATERAL)
    xs = cumulative(-sum(PAD_LATERAL),
                    PAD_LATERAL[::-1] + [CORE_DX] * CORE_NX + PAD_LATERAL)
    ys = cumulative(-sum(PAD_LATERAL),
                    PAD_LATERAL[::-1] + [CORE_DY] * CORE_NY + PAD_LATERAL)

    core_dz = [CORE_DZ] * CORE_NZ
    caprock_k = None
    if with_caprock:
        # Split the last seal layer (core-local RES_K1 - 1) into a thick
        # remainder plus CAPROCK_N thin layers sitting on the reservoir.
        host = RES_K1 - 2                      # 0-based index of that layer
        core_dz = (core_dz[:host]
                   + [CORE_DZ - CAPROCK_N * CAPROCK_DZ]
                   + [CAPROCK_DZ] * CAPROCK_N
                   + core_dz[host + 1:])
        first_thin = len(PAD_ABOVE) + host + 2
        caprock_k = (first_thin, first_thin + CAPROCK_N - 1)

    zs = cumulative(CORE_TOP - sum(PAD_ABOVE), PAD_ABOVE + core_dz + PAD_BELOW)
    extra = CAPROCK_N if with_caprock else 0
    kcore0 = len(PAD_ABOVE)

    return Layout(xs, ys, zs,
                  core_i=(npad + 1, npad + CORE_NX),
                  core_j=(npad + 1, npad + CORE_NY),
                  core_k=(kcore0 + 1, kcore0 + CORE_NZ + extra),
                  res_k=(kcore0 + RES_K1 + extra, kcore0 + RES_K2 + extra),
                  well=(npad + WELL_I, npad + WELL_J, kcore0 + WELL_K + extra),
                  caprock_k=caprock_k)


def group_sizes(fine_count, core_range, grouping):
    """Merge `grouping` fine cells at a time outside the core, 1:1 inside."""
    sizes = []
    lo, hi = core_range
    i = 1
    while i <= fine_count:
        if lo <= i <= hi:
            sizes.append(1)
            i += 1
        else:
            n = min(grouping, (lo - i) if i < lo else (fine_count - i + 1))
            sizes.append(max(n, 1))
            i += max(n, 1)
    return sizes


def coarsen_layout(fine, grouping=PAD_GROUPING):
    """Coarse layout: pad cells merged in groups, core untouched."""
    gi = group_sizes(fine.nx, fine.core_i, grouping)
    gj = group_sizes(fine.ny, fine.core_j, grouping)
    gk = group_sizes(fine.nz, fine.core_k, grouping)

    def breaks(fine_breaks, sizes):
        out, idx = [fine_breaks[0]], 0
        for s in sizes:
            idx += s
            out.append(fine_breaks[idx])
        return out

    xs, ys, zs = (breaks(fine.xs, gi), breaks(fine.ys, gj), breaks(fine.zs, gk))
    # nesting: every coarse breakpoint must be a fine one
    for c, f in ((xs, fine.xs), (ys, fine.ys), (zs, fine.zs)):
        assert all(v in f for v in c), "coarse grid is not nested in the fine one"

    # A fine index maps to the coarse group containing it.
    def to_coarse(sizes, fine_index):
        return sum_to_index(sizes, fine_index) + 1

    core_i = (to_coarse(gi, fine.core_i[0]), to_coarse(gi, fine.core_i[1]))
    core_j = (to_coarse(gj, fine.core_j[0]), to_coarse(gj, fine.core_j[1]))
    core_k = (to_coarse(gk, fine.core_k[0]), to_coarse(gk, fine.core_k[1]))
    res_k = (to_coarse(gk, fine.res_k[0]), to_coarse(gk, fine.res_k[1]))
    well = tuple(to_coarse(g, idx)
                 for g, idx in zip((gi, gj, gk), fine.well))

    return Layout(xs, ys, zs, core_i, core_j, core_k, res_k, well), (gi, gj, gk)


# --- grdecl ------------------------------------------------------------------


def fmt(values, per_line=8):
    out, line = [], []
    for v in values:
        line.append("%.4f" % v if isinstance(v, float) else str(v))
        if len(line) == per_line:
            out.append(" " + " ".join(line))
            line = []
    if line:
        out.append(" " + " ".join(line))
    return "\n".join(out)


def write_grid_include(path, lay):
    coord = []
    for y in lay.ys:
        for x in lay.xs:
            coord += [x, y, lay.zs[0], x, y, lay.zs[-1]]

    zcorn = []
    for k in range(lay.nz):
        for surface in (lay.zs[k], lay.zs[k + 1]):
            for _j in range(lay.ny):
                for _row in range(2):
                    for _i in range(lay.nx):
                        zcorn += [surface, surface]

    with open(path, "w") as f:
        f.write("-- generated by make_coarsen_tests.py -- do not edit\n")
        f.write("SPECGRID\n %d %d %d 1 F /\n\n" % (lay.nx, lay.ny, lay.nz))
        f.write("COORD\n%s\n/\n\n" % fmt(coord, 6))
        f.write("ZCORN\n%s\n/\n\n" % fmt(zcorn))
        f.write("ACTNUM\n %d*1 /\n" % (lay.nx * lay.ny * lay.nz))


# --- deck keywords that depend on the index space ---------------------------


def box(i, j, k):
    return "%d %d %d %d %d %d" % (i[0], i[1], j[0], j[1], k[0], k[1])


def props_block(lay):
    whole = ((1, lay.nx), (1, lay.ny), (1, lay.nz))
    lines = ["EQUALS"]
    for key, val in PROPS_PAD.items():
        lines.append("  %-8s %-8g %s /" % (key, val, box(*whole)))
    lines.append("/\n")

    lines.append("EQUALS")
    for key, val in PROPS_ROCK.items():
        lines.append("  %-8s %-8g %s /"
                     % (key, val, box(lay.core_i, lay.core_j, lay.core_k)))
    lines.append("/\n")

    lines.append("EQUALS")
    for key, val in PROPS_RES.items():
        lines.append("  %-8s %-8g %s /"
                     % (key, val, box(lay.core_i, lay.core_j, lay.res_k)))
    lines.append("/\n")

    if lay.caprock_k:
        lines.append("-- thin seal layers on top of the reservoir (T2)")
        lines.append("EQUALS")
        for key, val in PROPS_CAPROCK.items():
            lines.append("  %-8s %-8g %s /"
                         % (key, val, box(lay.core_i, lay.core_j, lay.caprock_k)))
        lines.append("/\n")

    return "\n".join(lines)


def bccon_block(lay):
    full_i, full_j, full_k = (1, lay.nx), (1, lay.ny), (1, lay.nz)
    faces = [
        (1, (1, 1), full_j, full_k, "X-"),
        (2, (lay.nx, lay.nx), full_j, full_k, "X+"),
        (3, full_i, (1, 1), full_k, "Y-"),
        (4, full_i, (lay.ny, lay.ny), full_k, "Y+"),
        (5, full_i, full_j, (1, 1), "Z-"),
        (6, full_i, full_j, (lay.nz, lay.nz), "Z+"),
    ]
    lines = ["BCCON"]
    for idx, i, j, k, direction in faces:
        lines.append("  %d %s %s /" % (idx, box(i, j, k), direction))
    lines.append("/")
    return "\n".join(lines)


# --- deck assembly -----------------------------------------------------------


def substitute(text, pattern, replacement, what, flags=re.S):
    new, n = re.subn(pattern, lambda _m: replacement, text, flags=flags)
    if n != 1:
        raise RuntimeError("expected exactly one %s in the base deck, found %d"
                           % (what, n))
    return new


def make_deck(base, lay, grid_include, steps):
    text = base

    text = substitute(text, r"DIMENS\n(?:--[^\n]*\n)*\s*\d+\s+\d+\s+\d+\s*/[^\n]*",
                      "DIMENS\n  %d %d %d /" % (lay.nx, lay.ny, lay.nz),
                      "DIMENS")

    # geometry: the block-centred description becomes an INCLUDE of the grdecl
    text = substitute(text, r"TOPS\n.*?DZV\n\s*25\*20 /",
                      "INCLUDE\n  '%s' /" % grid_include,
                      "TOPS/DXV/DYV/DZV block")

    # the two property EQUALS blocks (rock default, reservoir box)
    text = substitute(text,
                      r"EQUALS\n  PERMX 10 6\* /.*?PORO  0\.28 4\* 11 20 / \n/",
                      props_block(lay), "property EQUALS blocks")

    text = substitute(text, r"BCCON\n.*?\n/", bccon_block(lay), "BCCON")

    i, j, k = lay.well
    for kw in ("BPR", "BTEMP", "BSTRSSXX"):
        text = substitute(text, kw + r"\n  6 6 16 /",
                          "%s\n  %d %d %d /" % (kw, i, j, k), kw)
    text = substitute(text, r"'B-3H' 6 6 16 16 'OPEN'",
                      "'B-3H' %d %d %d %d 'OPEN'" % (i, j, k, k), "COMPDAT")
    text = substitute(text, r"'B-3H' 6 6 16 0 1 0 5 0 1e-4/",
                      "'B-3H' %d %d %d 0 1 0 5 0 3e-4/" % (i, j, k), "WSEED")

    # The base deck stops after 15 days, which is too short for the fracture to
    # grow past its seed.
    text = substitute(text, r"TSTEP\n -- 18\*5\n 3\*5\n/",
                      "TSTEP\n %d*5\n/" % steps, "TSTEP")

    header = ("-- Generated by data/coarsen_tests/make_coarsen_tests.py.\n"
              "-- Padded SIMPLE for the corner-point coarsening work; see\n"
              "-- opm-gridrefined/docs/COARSENING-FOR-MECHANICS.md\n")
    # The base deck is full of trailing whitespace; don't carry it into a new file.
    return header + "\n".join(line.rstrip() for line in text.splitlines()) + "\n"


def spec_records(fine, groups):
    """The axis grouping, also expanded as COARSEN records."""
    gi, gj, gk = groups

    def spans(sizes):
        out, start = [], 1
        for s in sizes:
            out.append((start, start + s - 1))
            start += s
        return out

    # One record per (i,j,k) run of merged cells; runs that stay 1:1 are skipped.
    def runs(sizes):
        out, start, i = [], 1, 0
        while i < len(sizes):
            j, n = i, 0
            while j < len(sizes) and sizes[j] == sizes[i]:
                n += sizes[j]
                j += 1
            out.append((start, start + n - 1, j - i, sizes[i]))
            start += n
            i = j
        return out

    records = []
    for i1, i2, ni, si in runs(gi):
        for j1, j2, nj, sj in runs(gj):
            for k1, k2, nk, sk in runs(gk):
                if si == 1 and sj == 1 and sk == 1:
                    continue
                records.append([i1, i2, j1, j2, k1, k2, ni, nj, nk])

    return {
        "axis_groups": {"i": gi, "j": gj, "k": gk},
        "coarsen_records": records,
        "fine_dims": [fine.nx, fine.ny, fine.nz],
        "spans": {"i": spans(gi), "j": spans(gj), "k": spans(gk)},
    }


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--steps', type=int, default=18,
                    help='number of 5-day steps (default 18 = 90 days)')
    args = ap.parse_args()

    with open(BASE_DECK) as f:
        base = f.read()

    fine = build_fine()
    coarse, groups = coarsen_layout(fine)
    caprock = build_fine(with_caprock=True)

    cases = [
        ("T1_PAD_FINE", fine),
        ("T1_PAD_COARSE", coarse),
        ("T2_CAPROCK", caprock),
    ]
    for name, lay in cases:
        grid = "%s_GRID.INC" % name
        write_grid_include(os.path.join(HERE, grid), lay)
        with open(os.path.join(HERE, "%s.DATA" % name), "w") as f:
            f.write(make_deck(base, lay, grid, args.steps))
        print("%-14s %3d x %3d x %3d = %6d cells   core i%s j%s k%s   well %s"
              % (name, lay.nx, lay.ny, lay.nz, lay.nx * lay.ny * lay.nz,
                 lay.core_i, lay.core_j, lay.core_k, lay.well))

    spec = spec_records(fine, groups)
    spec["coarse_dims"] = [coarse.nx, coarse.ny, coarse.nz]
    with open(os.path.join(HERE, "coarsen_spec_T1.json"), "w") as f:
        json.dump(spec, f, indent=2)

    # T2: the mechanics-only coarsening is T1's padding grouping plus merging
    # the thin layers back into their host layer.
    cap_groups = coarsen_layout(caprock)[1]
    gk = list(cap_groups[2])
    k1 = caprock.caprock_k[0] - 1        # host layer, 1-based
    at = sum_to_index(gk, k1)
    gk[at:at + 1 + CAPROCK_N] = [1 + CAPROCK_N]
    cap_spec = spec_records(caprock, (cap_groups[0], cap_groups[1], gk))
    cap_spec["note"] = ("mechanics-only: pad merged as in T1, plus the %d thin "
                        "seal layers merged back into their host layer"
                        % CAPROCK_N)
    with open(os.path.join(HERE, "coarsen_spec_T2.json"), "w") as f:
        json.dump(cap_spec, f, indent=2)
    print("wrote coarsen_spec_T1.json, coarsen_spec_T2.json")


def sum_to_index(sizes, target_cell):
    """Index in `sizes` of the group holding 1-based fine cell `target_cell`."""
    acc = 0
    for idx, s in enumerate(sizes):
        acc += s
        if acc >= target_cell:
            return idx
    raise IndexError(target_cell)


if __name__ == "__main__":
    sys.exit(main())
