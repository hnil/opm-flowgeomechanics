# Geomechanics physics baselines

A small, dependency-free Python harness for the SIMPLE showcase and one-ring
WELLREF case. It checks that a run reaches the intended date, grows a fracture,
and routes the expected injection through it. It also records fracture open area,
volume, cumulative well injection, and reported solver failures/rescues. No
simulator equations or existing regression reference files are changed.

These are **regression acceptance bands**, based on `docs/GETTING_STARTED.md`,
`data/wellref/README.md` and the October 2026 Linux setup. They are not evidence of
mesh convergence, mass/energy conservation, or agreement with field measurements.
The bands intentionally allow small numerical variation; passing them does not
establish that every configuration or MPI rank count is equivalent.

[`reference_results.json`](reference_results.json) records the fresh 15 October
2026 observations, input and binary fingerprints, solver-baseline revisions and
warning counts. It is a compact review artifact, not an automatically updated
golden file. Full logs and per-step histories remain in the run directories.

## Run fresh baselines

Requires Python 3.9+, a configured five-repository superbuild and the targets
`flow_energy_geomech` and `summary`. The two decks are self-contained; their
fracture settings and mechanics-solver JSON files are specified in `cases.json`.
The harness does not build anything or change the decks/configuration.

From the directory containing the five repositories:

```bash
python3 opm-flowgeomechanics/scripts/physics_baselines.py run \
  --build-dir build \
  --output-root ../runs/physics-baselines-001
```

Both cases run sequentially: SIMPLE on one rank, WELLREF on four ranks, each with
one thread per rank. Expect several minutes overall, depending on hardware.
`--output-root` must not already exist. Simulator output and logs stay outside
the source tree. A failing case returns a nonzero suite exit status; its log and
manifest remain available. The default timeout is 900 seconds per case and
terminates the MPI process group if exceeded.

Useful options:

```bash
# Run only the short case.
python3 opm-flowgeomechanics/scripts/physics_baselines.py run \
  --build-dir build --output-root ../runs/simple-002 --cases simple

# Override the MPI launcher, or add a launcher argument (repeatable).
python3 opm-flowgeomechanics/scripts/physics_baselines.py run \
  --build-dir build --output-root ../runs/physics-baselines-003 \
  --mpi-launcher /usr/bin/mpiexec --launcher-arg=--oversubscribe
```

For MPI builds, the default launcher, process-count flag and pre/postflags come
from `CMakeCache.txt`; the simulator is not launched as an MPI singleton. A serial
build can run SIMPLE directly; the four-rank WELLREF profile requires MPI.
MPI needs working local sockets, including when launched from a restricted tool
environment. Raw stderr is captured in `simulator.log`.

## Inspect existing outputs

Analysis needs the matching summary tool, the original output directory, and a
complete run log (captured stdout or the simulator PRT file containing the solver
statistics and end-of-simulation marker). Reports go into a **new** directory:

```bash
python3 opm-flowgeomechanics/scripts/physics_baselines.py analyze \
  --case simple --run-dir ../runs/simple \
  --summary-tool build/opm-common/bin/summary \
  --log ../runs/simple/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.PRT \
  --report-dir ../runs/simple-analysis-001
```

This mode hashes the analyzed artifacts but does not claim to know which source,
build or configuration originally produced them. It never assigns current Git
revisions to historical results. Only `run` creates a fresh-run provenance manifest.

## Artifacts and interpretation

Every fresh case writes:

- `physics-manifest.json`: all five Git revisions, tracked-diff fingerprints,
  untracked source hashes, deck/configuration hashes, executable hashes, launcher
  command, selected CMake settings, tool/DUNE versions, process/thread counts,
  platform information, return code and elapsed wall time. Commit dirty inputs
  before publishing results; a fingerprint identifies changes but cannot restore
  their contents. Git revisions describe the source tree at launch; binary hashes
  identify what actually executed. The harness does not prove the build is current.
- `physics-summary.csv`: the full summary time series with BHP in bar, injection
  rates in m³/day, cumulative injection in m³ and fracture-flow fraction. Fraction
  is empty when total injection is zero, rather than reporting a misleading zero.
- `physics-geometry.csv`: every snapshot referenced by the fracture PVD, with its
  timestamp, total area, open area, volume and filename. Geometry is read from
  ASCII VTU files using the existing `scripts/frac_area.py`; open means width
  greater than `1e-9` m. These are geometric aperture metrics, not an independent
  hydraulic-aperture or conservation measurement.
- `physics-report.json`: observed final metrics, acceptance bands, pass/fail
  checks, log diagnostics, and hashes of summary, log and geometry inputs.
- `simulator.log`: unmodified stdout/stderr from the simulator/launcher.

`physics-suite.json` at the output root collects case results and errors. The
harness never updates acceptance bands automatically.

**The final summary and geometry timestamps are distinct.** Current reference
summaries end at day 106, but the last PVD geometry is labelled day 101. The report
records this five-day gap and checks that it does not increase. It does not relabel
the last geometry as day 106. Duplicate initial PVD timestamps are preserved.

Solver counts are the maximum cumulative counts emitted in the log, not a sum
over ranks. A nonconvergence message, reported linear failure, rescue, absent
solver diagnostics or missing completion marker fails the baseline. This deliberately
strict diagnostic policy can reject a run that ultimately recovered; inspect the
log rather than changing physics to obtain a pass. These checks do not yet certify
that every rejected timestep or relaxed CNV acceptance has been detected.

Known warnings are visible as `pass_with_warnings`: unsupported deck keywords in
both cases, plus WELLREF's local mechanical-BC diagnostics,
refined-away summary vectors and unavailable LGR well-average pressures. Their
occurrences and other warning samples are retained. This is an explicit historical
allowlist, not a declaration that those warnings are harmless. A new warning type
fails; `--strict-warnings` also fails the known types. Occurrences count emitted
messages and must not be interpreted as the number of affected MPI ranks. Warning
classification currently recognizes `Warning:` messages and the known BC error;
other launcher diagnostics remain available in the raw log.
Missing-connection warnings now fail both cases: LGR completions are resolved by
(grid number, completion index), separately from the WSEED parent-cell identity.

## Fast harness checks

```bash
python3 -m unittest discover \
  -s opm-flowgeomechanics/tests/physics_baselines -p 'test_*.py' -v

# After CMake configuration, the same checks are registered with CTest.
ctest --test-dir build/opm-flowgeomechanics -R '^Unit.PhysicsBaselineHarness$' \
  --output-on-failure
```

The tests cover an analytic triangle's area/volume, PVD timestamps, absent files,
truncated output, non-finite/nonmonotone summaries, no growth, low fracture flow,
solver failure and warning handling. They do not run the simulator. Full physics
runs remain explicit because they take minutes and need MPI.

When adding a new physical regime, define a new profile and independently justify
its bands. For this first harness the two supported cases each contain exactly
one fracture collection. Multiple fractures, different units, binary VTU, nested
WELLREF, restart and independent refinement studies require further work. Do not
silently aggregate asynchronous fracture snapshots or tune physical parameters to
make a regression pass.
