# Sector tests

Flow on a FLUXNUM sector (USEFLUX), mechanics on a larger grid.

## T0: SIMPLE_MECH, sector inside the box

The decks are copies of `../SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA` with three lines
added to GRID:

- **`FLUXNUM`** is 1 in I 3–9, J 3–9 and all layers (7×7×25 cells, injector B-3H at
  6,6,16), and 0 elsewhere.
- **`T0_PARENT`** adds DUMPFLUX with `FLUXTYPE 'PRESSURE'`. It runs the whole box and
  writes `T0_PARENT.FLUX`.
- **`T0_SECTOR`** adds `USEFLUX 'T0_PARENT'` and `FLUXREG 1`. Flow covers the sector only.
- **`T0_IDENTITY.coarsen`** gives both runs a mechanics grid of their own, one cell per
  cell of the whole box. In the sector run, the 1800 mechanics cells outside the flow
  sector take YMODULE and PRATIO from the deck and carry no pressure or temperature load.

The FLUX file is looked up next to the sector deck, so run the sector from a copy:

```bash
D=opm-flowgeomechanics/data
run() {
  build/opm-flowgeomechanics/bin/flow_energy_geomech --parameter-file=$D/myparamsfrac_seq.ini \
    --fracture-param-file=$D/fracture_simple_bhp_seq_final.json \
    --linear-solver-mech=$D/amgmechsolverhypre_geosxcpu_new_nv.json \
    --mech-coarsen-file=$D/sector_tests/T0_IDENTITY.coarsen "$@"
}
run $D/sector_tests/T0_PARENT.DATA --output-dir=/tmp/t0/parent
cp $D/sector_tests/T0_SECTOR.DATA /tmp/t0/parent/
run /tmp/t0/parent/T0_SECTOR.DATA --output-dir=/tmp/t0/sector
python3 opm-flowgeomechanics/scripts/cmp_cells.py --cell 6,6,16 /tmp/t0/parent/T0_PARENT /tmp/t0/sector/T0_SECTOR
```

### Results (2026-10-18, `sector-geomech` stack)

The parent (T0_PARENT) is identical to the plain deck. That is, the mechanics grid with
`mech_grid_bc = roller_free_top` reproduces the deck's BCCON/BCMECH exactly, and
DUMPFLUX does not change the run.

Sector against parent, over the 1225 sector cells:

| | max diff | at the injector cell |
|---|---|---|
| PRESSURE | 4.6e-5 bar | 220.93 / 220.93 |
| STRESSXX | 0.16 bar | 238.16 / 238.12 |
| STRESSZZ | 0.10 bar | 414.52 / 414.51 |
| WBHP:B-3H, day 15 | 1.5e-4 bar | 257.0247 / 257.0246 |
| fracture area / volume | – | 93.3 / 93.3 m², 0.1047 / 0.1049 m³ |

The stress difference comes from the cells outside the sector. They are loaded in the
parent and unloaded in the sector run.

For contrast, the same sector run with mechanics on the flow grid only (no
`--mech-coarsen-file`):

- the lateral BCCON boxes fall on inactive cells, and the run reports an improper
  mechanical boundary;
- stress at the injector is off by 2.0 (XX) and 4.6 (ZZ) bar;
- fracture volume is 0.0940 m³, 10 % low.

## T1: mechanics box and padding built at run time

`--mech-sector-file=<json>` replaces the coarsening file. The mechanics grid is cut out
of the deck's grid (`"box"`, I1 I2 J1 J2 K1 K2, one-based and inclusive; the whole grid
if left out) and padded around it:

- `"lateral"`: `cells` over `width` metres on every side;
- `"top"`: up to `depth`, or over `thickness`;
- `"bottom"`: over `thickness`, or down to `depth`;
- cell sizes grow by `growth` away from the box.

With `"vertical_pillars": true` in `"lateral"`, the padding pillars are vertical instead of
parallel to the box's boundary pillars. The first one clears the slanted boundary pillar by
the first cell's width at every depth of the grid, so padding pillars never cross. Pillars
above and below the box are the box's own and keep their slant. A coarsening that merges
the first padding column with the next one is refused, because the dropped pillar is not on
the line between the two it lies between.

The run warns when a column's pillars cross or fold within its depths.

With `"bend_pillars": true` (top level), each box pillar continues vertically above the
shallowest and below the deepest box corner on it, where the fault-healing padding layer
ends. The top and bottom padding then hangs on vertical pillars too. The dumped GRDECL
keeps straight pillars, since the format cannot hold the bend; look at the grid through
`--vtk-write-mech-grid` instead.

Padding cells, and box cells outside the flow grid, take YMODULE and PRATIO from the
nearest box cell and carry no load. The outer boundary takes `mech_grid_bc` (default
`roller_free_top`).

| spec | mechanics grid |
|---|---|
| `{}` | the deck's 11×11×25 box, no padding |
| `T1_PAD.json` | 17×17×31: 3 cells / 2000 m laterally, 3 cells up to the surface, 3 cells / 1000 m below |
| `T1_PAD_WIDE.json` | 21×21×35: 5 cells / 6000 m, 5 cells to the surface, 5 cells / 3000 m |
| `T1_SECTOR_PAD.json` | the 9×9×25 box at (2,2,1), one cell beyond the T0 flow sector, padded to 15×15×31 |

```bash
D=opm-flowgeomechanics/data
S=$D/sector_tests
runm() {
  build/opm-flowgeomechanics/bin/flow_energy_geomech --parameter-file=$D/myparamsfrac_seq.ini \
    --fracture-param-file=$D/fracture_simple_bhp_seq_final.json \
    --linear-solver-mech=$D/amgmechsolverhypre_geosxcpu_new_nv.json "$@"
}
runm $S/T0_PARENT.DATA --mech-sector-file=$S/T1_PAD.json --output-dir=/tmp/t1/pad
runm $S/T0_PARENT.DATA --mech-sector-file=$S/T1_PAD_WIDE.json --output-dir=/tmp/t1/wide
mkdir -p /tmp/t1/sector && cp /tmp/t1/pad/T0_PARENT.FLUX $S/T0_SECTOR.DATA /tmp/t1/sector/
runm /tmp/t1/sector/T0_SECTOR.DATA --mech-sector-file=$S/T1_SECTOR_PAD.json --output-dir=/tmp/t1/sector/out
python3 opm-flowgeomechanics/scripts/cmp_cells.py --cell 6,6,16 /tmp/t1/pad/T0_PARENT /tmp/t1/sector/out/T0_SECTOR
```

The sector's FLUX file comes from the padded parent here. Mechanics does not feed back
into flow, so it is the same file as T0's.

### Results (2026-10-19)

Stress at the injector cell (6,6,16) after 15 days, and the fracture:

| run | STRESSXX | STRESSZZ | area m² | volume m³ |
|---|---|---|---|---|
| `{}` (same as the plain deck, bit-identical) | 238.16 | 414.52 | 93.3 | 0.1047 |
| `T1_PAD` | 237.69 | 413.58 | 95.1 | 0.1098 |
| `T1_PAD_WIDE` | 237.64 | 413.48 | 95.5 | 0.1106 |
| `T1_SECTOR_PAD` (flow on the T0 sector) | 237.86 | 414.04 | 94.2 | 0.1075 |

What the numbers show:

- **Padding size.** Padding moves the injector stress by about 1 bar; the rollers at the
  deck's edge were holding it. Widening the padding further changes it by 0.1 bar, so
  `T1_PAD` is enough here.
- **Sector against full flow.** With the same padding, the sector run is within 0.5 bar
  of `T1_PAD`; pressure agrees to 5e-5 bar.
