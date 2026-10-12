# WELLREF: refinement around wells, with fractures

`WELLREF` (GRID section) refines nested boxes ("rings") around the completions of the
wells it names. It is expanded **at parse time** into ordinary `CARFIN ... ENDFIN`
blocks, so everything downstream is the normal LGR path: grid builder, EGRID/INIT/UNRST
per LGR, ResInsight, and the fracture-on-LGR code.

## Where it works

Branch **`geomech-gridrefined-2026-09`** in all five modules of the
`opm_geomech/gridrefined/` superbuild (pushed to hnil/opm-common, hnil/opm-simulators,
hnil/opm-upscaling, hnil/opm-flowgeomechanics; grid module on hnil/opm-gridrefined).
Binary: `gridrefined/build/opm-flowgeomechanics/bin/flow_energy_geomech`.

- opm-common `a18f30435` (WELLREF keyword and expansion), `f9ed4651a` (COMPDAT
  connections inside a box moved into the innermost LGR).
- opm-simulators `0a5a17260` (host transmissibility of nested cells).
- opm-flowgeomechanics: fracture seeding in a **nested** ring needs `LgrLookup::father`
  to walk up to level 0 (2026-10-12); a single ring works without it.

Not on the plain `geomech` / bcmech branches (no WELLREF, no fracture-on-LGR).

## Syntax

```
GRID
...
WELLREF
-- well    NX NY NZ  RING1 RING2 RING3  LAYERS  KLAYERS
  'B-3H'    3  3  1    1     2     0     PERF     3 /
/
```

- `NX NY NZ`: refinement of every ring relative to the ring outside it.
- `RING1..3`: ring widths in coarse cells, outermost first (0 = absent). The outermost
  box is the completion dilated by the sum of the widths, each inner ring by one width
  less.
- `LAYERS PERF` (default): the innermost ring covers the completed layers plus
  `KLAYERS` above and below, each outer ring one layer more. `ALL`: the whole column.
- Seeds are the `COMPDAT` cells of the matching wells (any report step). Write the well
  with plain `WELSPECS`/`COMPDAT`; each connection is moved into the innermost ring,
  split over the refined cells along the well with its factor divided between them.
- The log prints each box: `WELLREF: CARFIN 'WR1R1B1' on GLOBAL: cells [4,8]x[4,8]x[13,19] -> 25x25x21`.

## Rules for geomechanics

- Run with `--edge-conformal=true`; the mechanics needs a hanging-node-free leaf.
- `WELLDIMS` item 2 (connections per well) must hold the refined and fracture
  connections; the examples use 800.
- `WSEED` keeps the **level-zero** cell; the seed is found through the refined
  connections' ancestor.
- **The innermost ring must contain the whole fracture.** A fracture connection outside
  its well's LGR stops the run, naming the cell and the box to enlarge.
- Block summary vectors (`BPR`, `BSTRSSXX`, ...) on a refined-away cell report zero
  (warning); use the `LB*` forms.
- `EDIT` modifiers on `TRANX/TRANY/TRANZ` with a refined grid: serial only.

## Examples

All with `best_practice/seq_implicit_fb_robust_v2.json`, `--edge-conformal=true`.

| deck | grid | how to run |
|---|---|---|
| `SIMPLE_MECH_FRAC_WELLREF.DATA` | one ring, 5x5x3, = the hand-written `..._FRAC_LGR_CONTAINED.DATA` box | 4 ranks |
| `SIMPLE_MECH_FRAC_WELLREF_NESTED.DATA` | two rings 3x3x1 each: inner cells ~20 m cubes | 4 ranks |
| `make_model2_wellref.py IN.DATA OUTDIR` | model2 with one 5x5x1 ring around INJ1 | serial, see below |

```bash
BIN=gridrefined/build/opm-flowgeomechanics/bin/flow_energy_geomech
D=gridrefined/opm-flowgeomechanics/data
mpirun -np 4 $BIN $D/wellref/SIMPLE_MECH_FRAC_WELLREF.DATA \
  --fracture-param-file=$D/best_practice/seq_implicit_fb_robust_v2.json \
  --linear-solver-mech=$D/amgmechsolverhypre_geosxcpu_new_nv.json \
  --edge-conformal=true --parsing-strictness=low --threads-per-process=1 \
  --enable-write-all-solutions=true --output-dir=/tmp/wellref_simple
```

Results (2026-10-12, 106 days, 4 ranks):

| case | leaf cells | BHP day 106 | WWIRFRAC/WWIR | fracture area (open) | wall |
|---|---|---|---|---|---|
| hand-written CARFIN + COMPDATL | 15975 | 272.0 | 65.0 % | 12416 (12416) m² | 166 s |
| WELLREF one ring | 15975 | 240.9 | 64.9 % | 13136 (12998) m² | 210 s |
| WELLREF two rings | 19153 | 247.3 | 61.0 % | 3014 (2903) m² | 141 s |

- WELLREF builds exactly the hand-written box. The BHP differs because the
  hand-written deck completes one refined cell (`COMPDATL ... 11 11`), while WELLREF
  spreads the coarse connection over the three refined layers.
- The nested run's fracture is a quarter of the single ring's on ~5x finer cells near
  the well: the fracture-reservoir coupling depends on cell size (see the
  grid-controlled-initiation notes). Not a truncation: it sits well inside the ring.

## model2: not geometric yet

`make_model2_wellref.py` writes the model2 deck with `'INJ1' 5 5 1 1 0 0 PERF 3`
(I 1-3, J 12-14, K 10-19; fault F2 runs along the box's east face). The leaf is **not
edge-conformal**: the mechanics check finds 61 cells with hanging nodes and unpaired
boundary edges (with `KLAYERS 2`: 45 cells, a crack). The refined box against faulted
corner-point neighbours is the cause; SIMPLE (Cartesian) is clean. It runs only with
`--mech-grid-check=warn`, which puts VEM on that grid, so its numbers are a
demonstration, not results. Also serial only (model2 has `EDIT` TRANX modifiers).

```bash
python3 $D/wellref/make_model2_wellref.py <model2 dir>/OPM_TEST_MODEL2_NOPARTICLES_MECH_FRAC.DATA /tmp/m2wr
cd /tmp/m2wr && $BIN --parameter-file=<tests_2808 params.ini> OPM_TEST_MODEL2_NOPARTICLES_MECH_FRAC_WELLREF.DATA \
  --mech-grid-check=warn --matrix-add-well-contributions=true --shut-unsolvable-wells=false \
  --solver-max-time-step-in-days=5 --solver-min-time-step=0.005 \
  --solver-continue-on-convergence-failure=true --output-dir=/tmp/m2wr/out
```

2026-10-12, serial, 334 days, 56 min (unrefined np4: 6 min), 7594 leaf cells: BHP
384.5 at day 10 and 411.2 at the end, WWIRFRAC/WWIR 48 % -> 65 %, area 2030 (1605
open) m². Unrefined: 468.8 at day 10, 414.4, 64 %, 1322 m²; reference: 429.7, 411.1,
1571 m².

## Known gaps

- LGR boxes on faulted corner-point grids are not edge-conformal (model2 above).
- In parallel, the mechanics grid check sees level 0 only (the global view of an LGR
  grid is unrefined); the refined leaf is checked only in serial.
- ResInsight does not draw a nested inner ring correctly; the outer ring is fine.
