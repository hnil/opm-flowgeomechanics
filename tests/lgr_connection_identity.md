# LGR fracture completion identity

Branch: `fix/lgr-fracture-connection-identity` (based on `0ab3be4`).

WELLREF's WSEED cell 1875 is a level-zero parent cell. Runtime seed-selection
instrumentation confirmed that its chosen schedule completion is index 6562 in
LGR grid 1. Looking up completion 1875 therefore emitted repeated warnings and
skipped reservoir/filter-cake property work. The temporary instrumentation has
been removed.

The fracture now retains two identities: the parent WSEED index for activation,
and `(completion index, LGR grid number)` for schedule properties. `Fracture::init` requires callers to supply the completion identity explicitly;
there is no fallback from the WSEED parent index. A grid number
is essential because different grids can contain the same numeric index.
Seed selection also maps the chosen completion into the rank-local connection
list by identity; a full-schedule ordinal is invalid after nonlocal completions
have been removed. A direct WSEED match is restricted to grid zero.

The existing policy of choosing the middle completion under a refined WSEED
parent remains in place. This does not introduce multiple fractures per WSEED,
prove behavior for nested refinement, or certify MPI decomposition invariance.

## Checks

From the workspace root:

```sh
cmake --build src/build --target flow_energy_geomech test_fracture_connection_identity -j6
ctest --test-dir src/build/opm-flowgeomechanics \
  -R 'Unit\.(FractureConnectionIdentity|PhysicsBaselineHarness)' --output-on-failure
python3 src/opm-flowgeomechanics/scripts/physics_baselines.py run \
  --build-dir src/build --output-root runs/lgr-connection-identity-NEW \
  --cases simple wellref
```

The output directory must be new. MPI requires working local socket support.
The C++ regression checks grid/index collisions, a missing grid, parent versus
completion identity, reordered local connections, and a nonlocal completion.
The Python harness checks that missing-connection warnings fail both profiles.
All nine existing `REGTEST_geomech` tests passed; tracked result files were
restored after the run.

Full baseline results are recorded below after validation. The original
`physics_baselines/reference_results.json` remains the historical observation
before this fix, including its 440 WELLREF missing-connection warnings.

## Filter-cake follow-up

An additional four-rank smoke run completed using a temporary copy of WELLREF:
activate its commented WINJDAM (`'B-3H' 'LINRAD' 10.0 0.0`) and WINJFCNC
(`'B-3H' 10.0000`) records, and replace both final TSTEP blocks with `2*1 /`.
Use the baseline fracture/mechanics JSON files and the usual
`--threads-per-process=1 --parsing-strictness=low` flags.

On this machine the deck is `/tmp/opm-wellref-nonzero-cake.DATA`, output is
`/tmp/opm-wellref-nonzero-cake-output`, and log is
`/tmp/opm-wellref-nonzero-cake.log`. There were no missing-connection warnings,
but all saved FilterCakeThickness values were zero. This run establishes a
completed execution, not positive filter-cake accumulation. It also reported
one relaxed CNV acceptance. Investigate concentration propagation, completion
property preservation during schedule updates, rate scaling, and output timing
before adding a quantitative filter-cake acceptance test. These are candidate
investigation points, not diagnosed causes.

## Full baseline validation (2026-10-15)

Outputs: `runs/lgr-connection-identity-2026-10-15/` under the workspace root.

Both cases completed day 106 and passed all acceptance checks, with only the
remaining allowlisted warnings. Geometry is from day 101.

| Case | BHP (bar) | Fracture flow fraction | Area (m²) | Volume (m³) | Missing connections |
| --- | ---: | ---: | ---: | ---: | ---: |
| simple | 220.976288 | 0.9638071289 | 11694.944 | 293.231735 | 0 |
| wellref | 219.740097 | 0.9939622070 | 11085.292 | 221.388269 | 0 |

All reported physics metrics match the historical observations exactly. WELLREF
missing-connection messages dropped from 440 to zero. Its remaining warning counts
are mechanical BC 3, unsupported keywords 1, refined summary 3, and LGR well-average
pressure 1. SIMPLE retains one unsupported-keyword warning. Both cases report
zero nonconvergence messages and zero linear-solve failures.

## Review follow-up

Removed the optional completion-identity argument and its implicit grid-zero
fallback after review. Omitting the identity now fails at compile time rather
than silently treating a WSEED parent as a completion. All three production
creation paths already supply the identity explicitly. The obsolete six-argument
call in `examples/single_fracture_sim.cpp` predates this branch and remains outside
this change.
