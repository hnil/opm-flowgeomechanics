#!/usr/bin/env python3
# Copyright 2026 Equinor ASA.
#
# This file is part of the Open Porous Media project (OPM).  OPM is free
# software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.  See
# <http://www.gnu.org/licenses/>.
"""Very thin layers at the base of the overburden, on SIMPLE with or without a CARFIN.

The overburden layer directly above the reservoir (k=10) is split into
(20 m - N*eps) + N layers of eps, all with the overburden's rock, so the model is
the same as SIMPLE and only the grid differs. Flow keeps the thin layers; the
records file merges them back into their host for the mechanics.

    make_thin_layers.py OUTDIR [--eps 0.01] [--count 3] [--carfin]

writes OUTDIR/THIN_<eps>[_LGR].DATA and OUTDIR/mech_merge_thin.txt, and with
--eps 0 the unsplit reference deck.
"""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..", "SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA")
HOST_K = 10                      # overburden layer right above the reservoir
CARFIN = " 'LGRW'   5  7  5  7 %d %d   9  9  3 /"


def sub(text, pattern, replacement, what):
    out, n = re.subn(pattern, replacement, text, flags=re.M)
    if n == 0:
        raise SystemExit("no match for " + what)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--count", type=int, default=3)
    ap.add_argument("--carfin", action="store_true")
    ap.add_argument("--steps", default="2*5")
    args = ap.parse_args()

    n = args.count if args.eps > 0 else 0
    kk = lambda k: k + n if k > HOST_K else k     # noqa: E731

    with open(BASE) as f:
        deck = f.read()

    dz = [20.0]*25
    dz[HOST_K - 1:HOST_K] = [20.0 - n*args.eps] + [args.eps]*n
    deck = sub(deck, r"^(\s*)11 11 25 /", r"\g<1>11 11 %d /" % len(dz), "DIMENS")
    deck = sub(deck, r"^DZV\s*\n\s*25\*20 /",
               "DZV\n  " + " ".join("%.6g" % d for d in dz) + " /", "DZV")
    deck = sub(deck, r"4\* 11 20 /", "4* %d %d /" % (kk(11), kk(20)), "reservoir EQUALS")
    deck = sub(deck, r"^(\s*6\s+2\*\s+2\*\s+)25(\s+)25(\s+Z\+)",
               r"\g<1>%d\g<2>%d\g<3>" % (len(dz), len(dz)), "BCCON Z+")
    well_k = kk(16)
    deck = sub(deck, r"6 6 16 16 'OPEN'", "6 6 %d %d 'OPEN'" % (well_k, well_k), "COMPDAT")
    deck = sub(deck, r"^(\s*)6 6 16 /", r"\g<1>6 6 %d /" % well_k, "BPR/BTEMP/BSTRSS")
    deck = sub(deck, r"'B-3H' 6 6 16 0", "'B-3H' 6 6 %d 0" % well_k, "WSEED")
    deck = sub(deck, r"^TSTEP\s*\n -- 18\*5\n 3\*5", "TSTEP\n -- 18*5\n " + args.steps,
               "TSTEP")
    if args.carfin:
        deck = sub(deck, r"^RPTGRID", "CARFIN\n-- name  I1 I2 J1 J2 K1 K2  NX NY NZ\n"
                   + CARFIN % (kk(15), kk(17)) + "\nENDFIN\n\nRPTGRID", "CARFIN")

    os.makedirs(args.outdir, exist_ok=True)
    name = "THIN_%g" % args.eps if n else "THIN_NONE"
    name = name.replace(".", "p") + ("_LGR" if args.carfin else "")
    with open(os.path.join(args.outdir, name + ".DATA"), "w") as f:
        f.write(deck)
    with open(os.path.join(args.outdir, "mech_merge_thin.txt"), "w") as f:
        f.write("-- Merge the thin layers back into the overburden layer they were cut from.\n"
                "-- I1 I2 J1 J2 K1 K2 NX NY NZ\n"
                " 1 11 1 11 %d %d 11 11 1\n" % (HOST_K, HOST_K + n))
    print(name, "layers", len(dz), "thin", n, "x", args.eps, "m, well k", well_k)


if __name__ == "__main__":
    sys.exit(main())
