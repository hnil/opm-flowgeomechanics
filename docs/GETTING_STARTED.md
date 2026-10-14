# Getting started with opm-flowgeomechanics (October 2026)

`flow_energy_geomech` couples OPM Flow (thermal, black-oil) with a VEM mechanics solve and
fracture growth from injection wells. This note gets you from nothing to a run you can look
at: build the right branches, run the small SIMPLE model, check and plot the result.
For macOS see [GETTING_STARTED_MACOS.md](GETTING_STARTED_MACOS.md).

Current line: branch **`lgr-status-geomech`** in all five repositories on
`github.com/hnil`. It is upstream OPM master (13 Oct 2026) plus the LGR work (both LGR
implementations behind `--lgr-backend`) and the geomechanics work. Geomechanics always runs
on **edge-conformal** grids with the **conforming** LGR backend; `flow_energy_geomech` sets
both by default and stops if you override them.

## 1. Build (Ubuntu 24.04)

Use Ubuntu 24.04 or newer: the code is C++20 and needs GCC 13 or newer. The DUNE and OPM
steps below were tried end to end on macOS; the Ubuntu package list has not been tried on
Ubuntu itself, so if configure reports something missing, install it and tell Halvor.

### System packages

```bash
sudo apt update
sudo apt install build-essential gfortran git cmake ninja-build pkg-config python3 \
  libopenmpi-dev openmpi-bin libboost-all-dev libfmt-dev libsuitesparse-dev \
  libblas-dev liblapack-dev libhypre-dev libtrilinos-zoltan-dev \
  libcgal-dev libgmp-dev paraview
```

HYPRE is the mechanics linear solver, Zoltan the parallel partitioner, CGAL the fracture
geometry.

### DUNE 2.11 from source

Ubuntu packages DUNE 2.9 without dune-foamgrid (the fracture grid), so build these six
modules into a folder of your own, in this order:

```bash
DUNE=$HOME/opm-geomech/dune-2.11
mkdir -p ~/opm-geomech/dune-src && cd ~/opm-geomech/dune-src
for m in core/dune-common core/dune-geometry core/dune-grid core/dune-istl \
         core/dune-localfunctions extensions/dune-foamgrid; do
  git clone -b v2.11.0 --depth 1 https://gitlab.dune-project.org/$m.git
done
for m in dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-foamgrid; do
  cmake -S $m -B build-$m -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$DUNE -DCMAKE_PREFIX_PATH=$DUNE \
    -DDUNE_ENABLE_PYTHONBINDINGS=OFF -DDUNE_PYTHON_USE_VENV=OFF
  cmake --build build-$m -j 8 && cmake --install build-$m
done
```

A few minutes. Put `export DUNE=$HOME/opm-geomech/dune-2.11` in your `~/.bashrc`.

### Get the code

```bash
cd ~/opm-geomech
for m in opm-common opm-grid opm-simulators opm-upscaling opm-flowgeomechanics; do
  git clone -b lgr-status-geomech https://github.com/hnil/$m.git
done
cp opm-flowgeomechanics/superbuild/CMakeLists.txt .
```

The copied `CMakeLists.txt` builds all five as one project, so nothing has to be installed.

### Configure and build

```bash
cd ~/opm-geomech
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$DUNE \
  -DUSE_MPI=ON -DUSE_HYPRE=ON -DUSE_OPENMP=ON 2>&1 | tee configure.log
grep -i "requested but not found" configure.log    # must print nothing
cmake --build build --target flow_energy_geomech summary compareECL -j 8
```

- The `grep` catches a missing HYPRE: configure only warns, and the mechanics solver
  settings used below then fail at run time. Install `libhypre-dev` and reconfigure.
- If the build is killed (`[code=137]` or "Killed"), it ran out of memory: the heaviest
  files need several GB each, so use fewer jobs (`-j 4`).
- The first build takes a while (it compiles all of OPM); after that only what you change
  is rebuilt. Build the targets you need, not everything: `flow_energy_geomech` is the
  simulator, `summary` a small tool for reading results, `compareECL` what the regression
  tests compare with.

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

About 45 seconds on one core of a recent laptop. What the options do:

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

- **ResInsight**: download the Linux release from
  https://github.com/OPM/ResInsight/releases and unpack it. Open `CASE.EGRID` for pressure,
  temperature and stresses in the grid; summary plots from `CASE.SMSPEC`.
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
echo "RESINSIGHT_EXECUTABLE=/path/to/ResInsight" >> ~/.config/resinsight-live/env
```

`RESINSIGHT_EXECUTABLE` points at the `ResInsight` program you unpacked (not needed if
`ResInsight` is on your `PATH`). The plugin needs Python 3.11 or newer (24.04 has 3.12).

Restart Claude Code. It only opens files inside `RESINSIGHT_LIVE_ALLOWED_DIRS`, so write
your runs there. Then, for example:

- *"Show PRESSURE and TEMP on the slice j=6 of ~/opm-geomech/runs/simple at the first and
  last step."*
- *"Plot WBHP, WWIR and WWIRFRAC for B-3H."*
- With your own ResInsight open: *"What are the values in the cell I just clicked?"*

## 4. Next steps

- **An offline AI agent** on a local model can run this example and report the numbers:
  [LOCAL_AI_AGENT.md](LOCAL_AI_AGENT.md).
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
