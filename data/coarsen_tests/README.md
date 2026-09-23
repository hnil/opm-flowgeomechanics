# Padded SIMPLE test cases for corner-point coarsening

Test decks for the coarse-mechanics-grid work. Plan:
`~/Documents/OPM/opm_vscode_clean/opm-gridrefined/docs/COARSENING-FOR-MECHANICS.md`.

Everything here is generated — edit `make_coarsen_tests.py`, not the decks:

```bash
python3 make_coarsen_tests.py            # 18 steps of 5 days (90 days)
python3 make_coarsen_tests.py --steps 3  # 15 days, quick smoke run
```

Physics (PROPS/SOLUTION/SUMMARY/SCHEDULE) is taken verbatim from
`../SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA`; only the geometry and the
index-dependent keywords (DIMENS, property EQUALS boxes, BCCON, BPR/BTEMP/BSTRSSXX,
COMPDAT, WSEED, TSTEP) are generated. The base deck's BCPROP and BCMECH forms both
work, since only BCCON is rewritten.

## Cases

| Deck | Grid | What it is |
|---|---|---|
| `T1_PAD_FINE` | 19 × 19 × 37 = 13 357 | SIMPLE's 11×11×25 core (181 m, 20 m layers) with 4 lateral padding cells per side (362 … 2896 m) and 6 layers of over- and underburden (20 … 1130 m), reaching the surface and 4300 m |
| `T1_PAD_COARSE` | 15 × 15 × 31 = 6 975 | the same model with every padding pair merged; the **exact geometric reference** the coarsening tool must reproduce from T1_PAD_FINE |
| `T2_CAPROCK` | 19 × 19 × 40 = 14 440 | T1_PAD_FINE with the seal layer above the reservoir split into 18.5 m + 3 × 0.5 m low-perm layers (aspect ≈ 360), which flow keeps and mechanics is meant to coarsen away |

Both coarse grids are strictly nested in their fine grid: every coarse cell is a union
of fine cells, and the generator asserts that every coarse breakpoint is a fine one.

`coarsen_spec_T1.json` / `coarsen_spec_T2.json` describe the fine → coarse mapping in
two equivalent forms: `axis_groups` (how many fine cells each coarse cell merges, per
axis) and `coarsen_records` (the same thing as COARSEN records
`I1 I2 J1 J2 K1 K2 NX NY NZ`). T2's spec is the mechanics-only one: T1's padding
grouping plus the thin layers merged back into their host.

The padding is rock for the mechanics and nearly inert for flow: PERMX 0.001 mD,
PORO 0.01. Mechanical properties are the base deck's, uniform over the whole grid.
BCCON sits on the outer faces of the padded grid.

## Running

```bash
BIN=~/Documents/OPM/opm_geomech/bcmech/build/opm-flowgeomechanics/bin/flow_energy_geomech
$BIN --output-dir=/tmp/coarsen/T1_PAD_FINE \
     --fracture-param-file=sequential_implicit \
     --enable-write-all-solutions=true T1_PAD_FINE.DATA
```

Check with `summary <case>.SMSPEC TIME WBHP:B-3H WWIRFRAC:B-3H 'BSTRSSXX:*'` and
`python3 ../test_oss/refine/frac_area.py <dir>/*_nr-*.vtu`.

## Coarsening the fine grid

The coarsening library lives in opm-gridrefined
(`opm/grid/cpgrid/coarsening/CornerPointCoarsening.{hpp,cpp}`), with the CLI
`examples/coarsen_grdecl.cpp`. Feeding it T1's records reproduces
`T1_PAD_COARSE_GRID.INC` exactly, and feeding it T2's mechanics-only records gives that
same grid, so the thin layers are merged away exactly:

```bash
ARGS=$(python3 -c "import json;print(' '.join('--coarsen \"%s\"' % ' '.join(map(str,r))
       for r in json.load(open('coarsen_spec_T1.json'))['coarsen_records']))")
eval coarsen_grdecl T1_PAD_FINE_GRID.INC /tmp/T1_COARSENED.INC $ARGS --lgr
```

`--lgr` prints the CARFIN-style requests that refine the coarse grid back to the fine
one (26 boxes for T1; the 1:1 core is not listed).

## Running the mechanics on its own grid

`mech_coarsen_T1.txt` holds the same records as `coarsen_spec_T1.json`, one per line,
for `--mech-coarsen-file`: flow stays on the fine grid, the mechanics runs on the
coarse one.

```bash
$BIN --output-dir=/tmp/coarsen/twogrid \
     --fracture-param-file=sequential_implicit \
     --mech-coarsen-file=mech_coarsen_T1.txt T1_PAD_FINE.DATA
```

The boundary of that grid is not the one BCCON describes, so it gets its own
constraint: `"mech_grid_bc"` in the parameter JSON, one of `fixed`, `roller` or
`roller_free_top` (the default, which matches these decks' free top).

In parallel the flow grid is partitioned by the coarsening too, so each mechanics
cell and all of its flow cells share a rank:

```bash
mpirun -np 4 $BIN --threads-per-process=1 --output-dir=/tmp/coarsen/twogrid_np4 \
     --fracture-param-file=sequential_implicit \
     --mech-coarsen-file=mech_coarsen_T1.txt T1_PAD_FINE.DATA
```

Not for decks with MINPV/PINCH yet.

## Baseline (2026-09-23, bcmech build, serial, 90 days, `sequential_implicit`)

Values at the last step. BHP is on its 290 bar limit in all three.

| Run | BSTRSSXX | WWIRFRAC | fracture area | fracture volume | run time |
|---|---|---|---|---|---|
| T1_PAD_FINE (flow fine, mech fine) | 175.86 | 719 | 3879.8 m² | 92.60 m³ | 114 s |
| T1_PAD_FINE + `--mech-coarsen-file` | 176.88 | 719 | 3879.8 m² | 91.92 m³ | 70 s |
| T1_PAD_COARSE (both coarse) | 176.96 | 691 | 3879.8 m² | 91.66 m³ | 62 s |
| T2_CAPROCK (flow fine, mech fine) | 190.38 | 527 | 3879.8 m² | 79.74 m³ | 149 s |

On 2 ranks, the two-grid run gives 175.92 / 718.98 / 3879.8 m² / 92.71 m³ in 64 s,
against 175.87 / 718.69 / 3879.8 m² / 92.78 m³ in 110 s for the single-grid run on the
same ranks. np=4 gives 175.04 / 718.73 / 3879.8 m² / 93.43 m³ in 62 s.

Coarsening the padding — for flow and mechanics together — moves the stress at the
well by 1.1 bar (0.6 %) and the fracture volume by 1.0 %, with the same fracture area,
and halves the run time.

The two-grid run sits between the two, at 0.6 % in stress and 0.7 % in fracture volume
from the all-fine run, with the same fracture area and the same fracture injection
rate, in 61 % of the time. T2 is a different model, not a coarsening of T1, so its
numbers are its own reference.


With `--steps 3` (15 days) the fracture is still at its 74.7 m² seed in every case, so
that schedule is a smoke test only.
