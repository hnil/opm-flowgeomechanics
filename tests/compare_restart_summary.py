#!/usr/bin/env python3
# Copyright 2026 Equinor ASA.
#
# This file is part of the Open Porous Media project (OPM).
#
# OPM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPM.  If not, see <http://www.gnu.org/licenses/>.
"""Compare a restarted run's summary with the full run's.

The restart must have the same connection keys with the given prefixes, and equal values
at every time both runs report.
"""
import argparse, re, subprocess, sys


# `summary -l` prints fixed-width columns that run long keys together, so match whole connection keys.
def keys(summary, smspec, prefixes):
    out = subprocess.run([summary, "-l", smspec], capture_output=True, text=True, check=True).stdout
    pattern = "(?:" + "|".join(map(re.escape, prefixes)) + r")[A-Z0-9_]*:[A-Za-z0-9_\-]+:\d+,\d+,\d+"
    return set(re.findall(pattern, out))


def series(summary, smspec, key):
    out = subprocess.run([summary, smspec, "TIME", key], capture_output=True, text=True, check=True).stdout
    vals = {}
    for line in out.splitlines():
        p = line.split()
        try:
            vals[round(float(p[0]), 6)] = float(p[1])
        except (IndexError, ValueError):
            pass
    return vals


def main():
    p = argparse.ArgumentParser()
    p.add_argument("summary")
    p.add_argument("full")
    p.add_argument("restart")
    p.add_argument("--prefix", nargs="+", required=True)
    p.add_argument("--rel-tol", type=float, default=1e-5)
    p.add_argument("--abs-tol", type=float, default=1e-5)
    a = p.parse_args()

    full, rst = a.full + ".SMSPEC", a.restart + ".SMSPEC"
    kf, kr = keys(a.summary, full, a.prefix), keys(a.summary, rst, a.prefix)
    if not kf:
        print("no keys with prefixes", a.prefix, "in", full)
        return 1
    if kf != kr:
        print("keys missing in restart:", sorted(kf - kr)[:10], "extra:", sorted(kr - kf)[:10])
        return 1
    bad = 0
    points = 0
    for k in sorted(kf):
        f, r = series(a.summary, full, k), series(a.summary, rst, k)
        for t in r.keys() & f.keys():
            points += 1
            d = abs(f[t] - r[t])
            if d > a.abs_tol and d > a.rel_tol * max(abs(f[t]), abs(r[t])):
                bad += 1
                if bad <= 10:
                    print(f"{k} t={t}: full {f[t]} restart {r[t]}")
    print(f"{len(kf)} keys, {points} points compared, {bad} differ")
    return 1 if bad or points == 0 else 0


if __name__ == "__main__":
    sys.exit(main())
