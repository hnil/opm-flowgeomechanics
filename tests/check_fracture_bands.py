#!/usr/bin/env python3
"""Band check for best-practice fracture runs (self-contained, no deps).

Guards the qualitative result - the fracture grows, conducts, and the well
pressure lands in a physical band - without pinning bits. Bands are per-config
(see CMakeLists); tighten only with a validated re-baseline.
"""
import argparse, glob, math, re, sys

def last_vtu_area(vtu_glob):
    files = sorted(glob.glob(vtu_glob))
    if not files:
        return None, None
    txt = open(files[-1]).read()
    m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt, re.S)
    pts = [float(x) for x in m.group(1).split()]
    conn = [int(x) for x in re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S).group(1).split()]
    offs = [int(x) for x in re.search(r'Name="offsets"[^>]*>(.*?)</DataArray>', txt, re.S).group(1).split()]
    fw = [float(x) for x in re.search(r'Name="FractureWidth"[^>]*>(.*?)</DataArray>', txt, re.S).group(1).split()]
    tot = 0.0; open_a = 0.0; prev = 0
    for ci, o in enumerate(offs):
        idx = conn[prev:o]; prev = o
        a = 0.0
        ax, ay, az = pts[3*idx[0]], pts[3*idx[0]+1], pts[3*idx[0]+2]
        for k in range(1, len(idx)-1):
            bx, by, bz = pts[3*idx[k]], pts[3*idx[k]+1], pts[3*idx[k]+2]
            cx, cy, cz = pts[3*idx[k+1]], pts[3*idx[k+1]+1], pts[3*idx[k+1]+2]
            ux, uy, uz = bx-ax, by-ay, bz-az
            vx, vy, vz = cx-ax, cy-ay, cz-az
            wx, wy, wz = uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx
            a += 0.5*math.sqrt(wx*wx+wy*wy+wz*wz)
        tot += a
        if ci < len(fw) and fw[ci] > 1e-6:
            open_a += a
    return tot, open_a

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--log", required=True)
    p.add_argument("--vtu-glob", required=True)
    p.add_argument("--bhp-min", type=float, required=True, help="bar")
    p.add_argument("--bhp-max", type=float, required=True)
    p.add_argument("--area-min", type=float, required=True, help="m2")
    p.add_argument("--area-max", type=float, required=True)
    p.add_argument("--open-frac-min", type=float, default=0.0)
    a = p.parse_args()

    fails = []
    bhps = re.findall(r'Inj bhp: ([0-9.eE+-]+)', open(a.log).read())
    if not bhps:
        fails.append("no 'Inj bhp:' lines in the log (run failed before the fracture solve?)")
    else:
        bhp = float(bhps[-1]) / 1e5
        print(f"final Inj bhp: {bhp:.1f} bar (band {a.bhp_min}-{a.bhp_max})")
        if not (a.bhp_min <= bhp <= a.bhp_max):
            fails.append(f"BHP {bhp:.1f} outside [{a.bhp_min}, {a.bhp_max}]")
    tot, open_a = last_vtu_area(a.vtu_glob)
    if tot is None:
        fails.append(f"no vtu files match {a.vtu_glob}")
    else:
        of = open_a / tot if tot > 0 else 0.0
        print(f"final fracture area: {tot:.0f} m2 (band {a.area_min}-{a.area_max}), open fraction {of:.2f}")
        if not (a.area_min <= tot <= a.area_max):
            fails.append(f"area {tot:.0f} outside [{a.area_min}, {a.area_max}] - "
                         "under-growth means fracturing was suppressed, over-growth means runaway")
        if of < a.open_frac_min:
            fails.append(f"open fraction {of:.2f} < {a.open_frac_min}")
    if fails:
        print("BAND TEST FAILED:\n  " + "\n  ".join(fails))
        return 1
    print("BAND TEST PASSED")
    return 0

if __name__ == "__main__":
    sys.exit(main())
