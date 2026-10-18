#!/usr/bin/env python3
"""Last-step restart arrays of two runs, compared over their common Cartesian cells
(no dependencies). The runs may have different active cells, e.g. a USEFLUX sector
against its parent.

usage: cmp_cells.py [--cell I,J,K] REF TEST [KW ...]   (case paths without extension)
"""
import struct, sys

def records(path):
    with open(path, "rb") as f:
        data = f.read()
    pos = 0
    def rec():
        nonlocal pos
        n, = struct.unpack(">i", data[pos:pos+4]); body = data[pos+4:pos+4+n]; pos += n + 8
        return body
    while pos < len(data):
        h = rec(); name = h[:8].decode().strip(); cnt, = struct.unpack(">i", h[8:12]); typ = h[12:16].decode()
        size = {"INTE": 4, "REAL": 4, "DOUB": 8, "LOGI": 4, "CHAR": 8, "MESS": 0}.get(typ, 8)
        fmt = {"INTE": "i", "REAL": "f", "DOUB": "d", "LOGI": "i"}.get(typ)
        vals = []
        while len(vals) < cnt and size:
            body = rec()
            if fmt:
                vals.extend(struct.unpack(">%d%s" % (len(body)//size, fmt), body))
            else:
                vals.extend(body[i:i+8].decode() for i in range(0, len(body), 8))
        yield name, vals

def grid(base):
    dims = None
    for name, v in records(base + ".EGRID"):
        if name == "GRIDHEAD":
            dims = v[1:4]
        if name == "ACTNUM":
            return dims, [g for g, a in enumerate(v) if a != 0]
    raise RuntimeError("no ACTNUM in " + base)

def last_step(base, wanted):
    out, cur = {}, {}
    for name, v in records(base + ".UNRST"):
        if name == "SEQNUM":
            cur = {}; out = cur
        elif name in wanted:
            cur[name] = v
    return out

args = sys.argv[1:]
cell = None
if args and args[0] == "--cell":
    cell = [int(x) - 1 for x in args[1].split(",")]
    args = args[2:]
ref, test = args[0], args[1]
kws = args[2:] or ["PRESSURE", "TEMP", "STRESSXX", "STRESSYY", "STRESSZZ", "STRESSXY", "DELSTRXX", "DELSTRZZ"]
dims, aref = grid(ref)
_, atest = grid(test)
iref = {g: i for i, g in enumerate(aref)}
common = [(iref[g], j, g) for j, g in enumerate(atest) if g in iref]
r, t = last_step(ref, set(kws)), last_step(test, set(kws))
well = None if cell is None else cell[0] + dims[0]*(cell[1] + dims[1]*cell[2])
print(f"common cells {len(common)} (ref {len(aref)}, test {len(atest)})")
print(f"{'kw':10} {'max|d|':>11} {'max|ref|':>11} {'rel':>9} {'at cell ref':>12} {'test':>12}")
for k in kws:
    if k not in r or k not in t:
        print(f"{k:10} missing"); continue
    d = max(abs(r[k][i] - t[k][j]) for i, j, _ in common)
    s = max(abs(r[k][i]) for i, j, _ in common) or 1.0
    w = [(r[k][i], t[k][j]) for i, j, g in common if g == well]
    ws = f"{w[0][0]:12.5g} {w[0][1]:12.5g}" if w else ""
    print(f"{k:10} {d:11.4g} {s:11.4g} {d/s:9.2e} {ws}")
