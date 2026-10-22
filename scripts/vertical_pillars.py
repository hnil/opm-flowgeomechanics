#!/usr/bin/env python3
"""Make every pillar of a corner-point grid vertical (no dependencies).

usage: vertical_pillars.py GRID.GRDECL [ACTNUM.prop] OUT.GRDECL

Each pillar becomes the vertical line through its own point at a reference depth: the mean
depth of the active cell corners on it, or of all its corners if none is active. ZCORN is kept,
so layer depths are unchanged and only the horizontal corner positions move. Slanted pillars
cross when they are extended past the reservoir; vertical ones never do.
"""
import math
import re
import sys


def read_keywords(path, wanted):
    text = open(path).read()
    text = re.sub(r"--[^\n]*", "", text)
    out, pos = {}, 0
    for m in re.finditer(r"^\s*([A-Z][A-Z0-9_]{1,7})\s*$", text, re.M):
        if m.group(1) in wanted:
            end = text.index("/", m.end())
            vals = []
            for tok in text[m.end():end].split():
                if "*" in tok:
                    n, v = tok.split("*")
                    vals += [float(v)] * int(n)
                else:
                    try:
                        vals.append(float(tok))
                    except ValueError:  # SPECGRID's coordinate flag
                        pass
            out[m.group(1)] = vals
    return out, text


def main():
    args = sys.argv[1:]
    grid_path, out_path = args[0], args[-1]
    kw, text = read_keywords(grid_path, {"SPECGRID", "COORD", "ZCORN"})
    nx, ny, nz = (int(v) for v in kw["SPECGRID"][:3])
    coord, zcorn = kw["COORD"], kw["ZCORN"]
    actnum = [1] * (nx * ny * nz)
    if len(args) == 3:
        actnum = [int(v) for v in read_keywords(args[1], {"ACTNUM"})[0]["ACTNUM"]]

    npil = (nx + 1) * (ny + 1)
    zsum, zn, zall, zalln = [0.0] * npil, [0] * npil, [0.0] * npil, [0] * npil
    zof = [[] for _ in range(npil)]  # active corner depths, for the move report
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                active = actnum[i + nx * (j + ny * k)] != 0
                for c in range(8):
                    ci, cj, ck = c & 1, (c >> 1) & 1, c >> 2
                    z = zcorn[(2 * i + ci) + 2 * nx * (2 * j + cj) + 4 * nx * ny * (2 * k + ck)]
                    p = (i + ci) + (nx + 1) * (j + cj)
                    zall[p] += z
                    zalln[p] += 1
                    if active:
                        zsum[p] += z
                        zn[p] += 1
                        zof[p].append(z)

    moved, degenerate = 0.0, 0
    new = list(coord)
    for p in range(npil):
        x1, y1, z1, x2, y2, z2 = coord[6 * p: 6 * p + 6]
        if z2 == z1:
            degenerate += 1
            continue
        zref = zsum[p] / zn[p] if zn[p] else zall[p] / max(zalln[p], 1)
        at = lambda z, a, b: a + (z - z1) / (z2 - z1) * (b - a)
        x, y = at(zref, x1, x2), at(zref, y1, y2)
        for z in (min(zof[p]), max(zof[p])) if zof[p] else ():
            moved = max(moved, math.hypot(at(z, x1, x2) - x, at(z, y1, y2) - y))
        new[6 * p: 6 * p + 6] = [x, y, z1, x, y, z2]

    m = re.search(r"^\s*COORD\s*$", text, re.M)
    end = text.index("/", m.end()) + 1
    body = "\n".join(" ".join(f"{v:.3f}" for v in new[q: q + 6]) for q in range(0, len(new), 6))
    open(out_path, "w").write(text[: m.start()] + "COORD\n" + body + "\n/" + text[end:])
    print(f"{npil} pillars made vertical ({degenerate} degenerate left alone); "
          f"active corners moved up to {moved:.1f} m")


if __name__ == "__main__":
    main()
