# Sector tests

Flow on a FLUXNUM sector (USEFLUX), mechanics on a larger grid. The design is in
`opm_geomech/fixed_tests/SECTOR_GEOMECH_DESIGN_2026-10-16.md` (not in this repository).

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
