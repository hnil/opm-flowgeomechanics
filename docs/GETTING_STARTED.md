# Getting started with opm-flowgeomechanics (October 2026)

`flow_energy_geomech` couples OPM Flow (thermal, black-oil) with a VEM mechanics solve and
fracture growth from injection wells. This note gets you from nothing to a run you can look
at: build the right branches, run the small SIMPLE model, check and plot the result.

Current line: branch **`lgr-status-geomech`** in all five repositories on
`github.com/hnil`. It is upstream OPM master (13 Oct 2026) plus the LGR work (both LGR
implementations behind `--lgr-backend`) and the geomechanics work. Geomechanics always runs
on **edge-conformal** grids with the **conforming** LGR backend; `flow_energy_geomech` sets
both by default and stops if you override them.

## 1. Build

### Prerequisites

Tested on macOS (Apple Silicon) with Apple Clang and Homebrew:

```bash
brew install cmake ninja open-mpi boost fmt suite-sparse hypre cgal gmp libomp python
```

DUNE **2.11** is needed: dune-common, dune-geometry, dune-grid, dune-istl,
dune-localfunctions and **dune-foamgrid** (the fracture grid), installed into one prefix
(`$DUNE` below). They build with plain CMake, in that order, each with
`-DCMAKE_INSTALL_PREFIX=$DUNE -DCMAKE_PREFIX_PATH=$DUNE`. Use the same compiler for DUNE
and OPM (on macOS: Apple Clang, not Homebrew's llvm/gcc).

On Linux the same packages exist (`libboost-all-dev libfmt-dev libsuitesparse-dev
libhypre-dev libcgal-dev libopenmpi-dev`); the OpenMP flag below is then not needed. This
note has only been tested on macOS.

### Get the code

```bash
mkdir -p ~/opm-geomech && cd ~/opm-geomech
for m in opm-common opm-grid opm-simulators opm-upscaling opm-flowgeomechanics; do
  git clone -b lgr-status-geomech https://github.com/hnil/$m.git
done
cp opm-flowgeomechanics/superbuild/CMakeLists.txt .
```

The copied `CMakeLists.txt` builds all five as one project, so nothing has to be installed.

### Configure and build

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$DUNE;/opt/homebrew" \
  -DUSE_MPI=ON -DUSE_HYPRE=ON -DUSE_OPENMP=ON \
  -DOpenMP_ROOT=/opt/homebrew/opt/libomp
cmake --build build --target flow_energy_geomech summary compareECL -j 6
```

- Check that OpenMP was found: `grep '^OpenMP_CXX_FLAGS:' build/CMakeCache.txt` must not say
  `NOTFOUND`.
- If the build is killed with `[code=137]`, it ran out of memory: use fewer jobs (`-j 3`).
- The first build takes a while; after that only what you change is rebuilt. Build the
  target you need, not everything: `flow_energy_geomech` is the simulator, `summary` a small
  tool for reading results, `compareECL` what the regression tests compare with.

Quick check (under a minute):

```bash
cd build/opm-flowgeomechanics && ctest -R REGTEST_geomech -j 4; cd ../..
```

The test driver rewrites a few files under `opm-flowgeomechanics/result/`; undo that with
`git -C opm-flowgeomechanics checkout -- result/`.

## 2. Run the SIMPLE model

An 11 x 11 x 25 box with one water injector (`B-3H`, 20 000 m³/day of 15 °C water into a
90 °C reservoir) and a fracture seeded at the well (`WSEED`). 106 days.

```bash
BIN=build/opm-flowgeomechanics/bin/flow_energy_geomech
DATA=opm-flowgeomechanics/data
OUT=~/opm-geomech/runs/simple
$BIN $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.DATA \
     --fracture-param-file=$DATA/best_practice/seq_implicit_fb_robust_v2.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --enable-write-all-solutions=true \
     --output-dir=$OUT
```

About 45 seconds on one core. What the options do:

| option | |
|---|---|
| `--fracture-param-file` | fracture solver settings; the files in `data/best_practice/` are the curated ones (see its README) |
| `--linear-solver-mech` | linear solver for the mechanics (AMG from HYPRE) |
| `--enable-write-all-solutions=true` | write the fracture geometry (`*.vtu`) at every step |
| `--output-dir` | always give one: otherwise the output lands next to the deck |

In `$OUT` you get the usual OPM output (`CASE.EGRID`, `.INIT`, `.UNRST`, `.SMSPEC`/`.UNSMRY`,
`.PRT` log, `.DBG`), the reservoir solution as VTK (`CASE-NNNNN.vtu`, `CASE.pvd`), and the
fracture per step (`CASEFracture_on_B-3H_perf_0_nr-NNNNN.vtu`, collected in a `.pvd`).

## 3. Check and plot

**Never judge a run by "it finished".** A run that finishes without fracturing is a failed
run. Check that the fracture grew and that it carries the injection.

Numbers, with the `summary` tool you built:

```bash
S=build/opm-common/bin/summary
$S -l $OUT/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.SMSPEC        # list vectors
$S $OUT/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.SMSPEC TIME WBHP:B-3H WWIR:B-3H WWIRFRAC:B-3H
python3 opm-flowgeomechanics/scripts/frac_area.py $(ls $OUT/*nr-*.vtu | tail -1)
```

`WWIRFRAC/WWIR` is the share of the injection that goes through the fracture. Expected at
day 106: BHP about 221 bar, `WWIRFRAC` about 19 300 of 20 000 m³/day (96 %), fracture area
about 11 700 m² (the seed is about 300 m²). The fracture starts to grow between day 14 and
21: the BHP drops from 261 to 226 bar while the fracture share rises from 7 % to 88 %.

Pictures:

- **ResInsight** (https://resinsight.org): open `CASE.EGRID` for pressure, temperature and
  stresses in the grid; summary plots from `CASE.SMSPEC`.
- **ParaView**: open the fracture `...Fracture_on_B-3H_perf_0_nr.pvd` and colour by
  `FractureWidth` or pressure; `CASE.pvd` holds the reservoir fields.

### Plotting through Claude (optional)

The `resinsight-live` plugin lets Claude Code drive ResInsight: ask for a picture or a plot
in plain words and it renders it. Its own guide is at
https://github.com/hnil/resinsight-live (GETTING_STARTED.md). In short:

```bash
claude plugin marketplace add https://github.com/hnil/resinsight-live
claude plugin install resinsight-live@hnil-reservoir-tools
mkdir -p ~/.config/resinsight-live
echo "RESINSIGHT_LIVE_ALLOWED_DIRS=$HOME/opm-geomech/runs" >> ~/.config/resinsight-live/env
```

Restart Claude Code. It only opens files inside `RESINSIGHT_LIVE_ALLOWED_DIRS`, so write
your runs there. Then, for example:

- *"Show PRESSURE and TEMP on the slice j=6 of ~/opm-geomech/runs/simple at the first and
  last step."*
- *"Plot WBHP, WWIR and WWIRFRAC for B-3H."*
- With your own ResInsight open (start the program itself on macOS:
  `/Applications/ResInsight.app/Contents/MacOS/ResInsight &`): *"What are the values in the
  cell I just clicked?"*

## 4. Next steps

- **Local grid refinement around the well:** `data/wellref/SIMPLE_MECH_FRAC_WELLREF.DATA`
  (one refinement ring from the `WELLREF` keyword). Run it on 4 processes with
  `mpirun -np 4 $BIN ... --threads-per-process=1 --parsing-strictness=low` and the same
  options as above; about 5 minutes. Expected at day 106: BHP about 220 bar, 99 % through
  the fracture, about 10 900 m². See `data/wellref/README.md`.
- Fracture settings: `data/best_practice/README.md` (which file for which case) and
  `data/fracture_settings_REFERENCE.md` (every option).
- Things that bite:
  - `WSEED` item 10 (seed width) overrides `min_width` in the JSON; too small a value
    starves the fracture.
  - Rate-controlled decks: add `--solver-max-time-step-in-days=2` and use a
    `seq_implicit_rate*` file.
  - JSON settings files cannot hold comments; the `"__..."` keys are dividers.
