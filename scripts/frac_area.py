#!/usr/bin/env python3
"""Area, open area and volume of fracture VTU files (no dependencies).

usage: frac_area.py RUN_DIR/*nr-*.vtu
The fracture files are written with --enable-write-all-solutions=true.
"""
import math
import sys
import xml.etree.ElementTree as ET


def data_array(piece, name):
    for da in piece.iter('DataArray'):
        if da.get('Name') == name:
            return [float(x) for x in da.text.split()]
    return None


def area_of(fn, wtol=1e-9):
    piece = ET.parse(fn).getroot().find('.//Piece')
    flat = data_array(piece, 'Coordinates')
    pts = [flat[i:i + 3] for i in range(0, len(flat), 3)]
    conn = [int(x) for x in data_array(piece, 'connectivity')]
    offsets = [int(x) for x in data_array(piece, 'offsets')]
    width = data_array(piece, 'FractureWidth')
    total = open_area = volume = 0.0
    start = 0
    for i, end in enumerate(offsets):
        idx = conn[start:end]
        start = end
        if len(idx) < 3:
            continue
        p0 = pts[idx[0]]
        a = 0.0
        for k in range(1, len(idx) - 1):
            u = [pts[idx[k]][d] - p0[d] for d in range(3)]
            v = [pts[idx[k + 1]][d] - p0[d] for d in range(3)]
            c = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
            a += 0.5 * math.sqrt(sum(x * x for x in c))
        total += a
        if width and width[i] > wtol:
            open_area += a
            volume += a * width[i]
    return total, open_area, volume


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        t, o, v = area_of(fn)
        print(f"{fn.split('/')[-1]}: total_area={t:.1f} m^2  open_area={o:.1f} m^2  volume={v:.4f} m^3")
