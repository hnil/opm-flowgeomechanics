# Global mechanical constraint validation

Branch: `fix/mpi-mechanics-constraint-validation`, based on `73e75ab`.

## Problem and behavior

The old initialization check counted fixed displacement directions on each MPI
rank and guessed rotational constraints from node counts, without coordinates.
A valid partition with no local Z support could warn despite the globally
constrained body. Conversely, fixing two points constrains all translations but
leaves a rotation about the line through them; the old heuristic could accept it.
It also skipped validation when a separate mechanics grid was configured.

Initialization now builds connected components of the actual mechanics grid's
vertex DOFs and checks the constraint rank on the six rigid-body modes for every
component. Local indices never cross ranks: communication uses CpGrid global
vertex IDs. Overlap copies merge masks by ID. No rank needs local supports as long
as its mechanical body is globally supported.

For normalized coordinates `(x,y,z)`, a fixed displacement component contributes
one of these rows, with columns `(tx,ty,tz,rx,ry,rz)`:

```text
ux: 1 0 0  0  z -y
uy: 0 1 0 -z  0  x
uz: 0 0 1  y -x  0
```

Coordinates are centered and scaled per component using the constrained points.
Reorthogonalized Gram-Schmidt computes the rank without normal equations, with a
normalized residual threshold of `1e-10`. Each body must have rank six. Rank 0
logs deficient component IDs/ranks, and all processes reject the configuration
before the mechanics solve. Valid configurations produce one informational message.
The check also runs on a separate mechanics grid rather than skipping it.

## Tests and reproduction

From the workspace root:

```sh
cmake --build src/build --target flow_energy_geomech test_mechanical_constraints -j6
ctest --test-dir src/build/opm-flowgeomechanics \
  -R 'Unit.(MechanicalConstraints|PhysicsBaselineHarness)' --output-on-failure
python3 src/opm-flowgeomechanics/scripts/physics_baselines.py run \
  --build-dir src/build --output-root runs/mpi-constraints-NEW --cases simple wellref
```

The output directory must be new and MPI must have local socket support.

All three native tests passed with 1, 2 and 4 MPI processes. They cover:

- Valid supports located entirely on rank 0, with no local supports elsewhere.
- Two/collinear fixed points (rank five), a fixed point (rank three), and a free body.
- Six independent scalar supports and distributed roller constraints.
- Overlap copies with duplicate or complementary masks.
- A supported component alongside a separate free component.
- Coordinate translations and changes of physical units.
- Invalid local BC indices: all ranks fail without a collective hang.

The Python harness tests also pass. Mechanical-BC warnings are no longer allowed
in either physics profile; a regression verifies their rejection.

An end-to-end negative experiment used a temporary copy of the SIMPLE showcase:
within its BCMECH block, replace each `FIXED` direction triple with `0 0 0`,
leaving all other inputs unchanged. Run with the baseline fracture and mechanics
JSON files, `--threads-per-process=1 --parsing-strictness=low`, and MPI np1/2/4.
All three runs exited with code 1, reported `0/6 rigid-body modes`, and rejected
the configuration before completing a simulation. Local reproduction files are
`/tmp/opm-unconstrained-mechanics.DATA` and
`/tmp/opm-unconstrained-mechanics-np{1,2,4}.log`.

## Scope and limits

This checks the elimination of rigid-body modes, not the full stiffness matrix's
positive definiteness, conditioning, material validity or arbitrary hinge modes.
Connectivity follows shared vertex DOFs; bodies attached only at a point/line may
have additional mechanisms beyond the six whole-component rigid motions.
Near-dependent supports can be rejected by the numerical rank tolerance.

The check gathers compressed local connectivity and constrained coordinates on
rank 0 at initialization, with memory proportional to the gathered vertex records.
It does not assume the retained serial grid matches distributed refinement, but
its communication/memory cost has not been benchmarked on very large meshes.
It validates the masks supplied to the mechanics solver; it does not introduce
new boundary-condition semantics or validate contradictory prescribed values.
Scheduled changes to mechanical constraints are not revalidated by this
initialization-only hook. TPSA is outside this VEM-specific change.

## Full baseline results (2026-10-15)

Both full runs completed day 106 and passed. Every reported physics metric
matched the historical baseline exactly (geometry saved at day 101).

| Case | BHP (bar) | Fracture flow fraction | Area (m²) | Volume (m³) | Mechanical BC warnings |
| --- | ---: | ---: | ---: | ---: | ---: |
| simple | 220.976288 | 0.9638071289 | 11694.944 | 293.231735 | 0 |
| wellref | 219.740097 | 0.9939622070 | 11085.292 | 221.388269 | 0 |

WELLREF mechanical-BC messages decreased from three emitted messages to zero.
Both cases also retain zero missing-connection warnings and zero reported linear
solve failures/nonconvergence messages. Remaining allowlisted warnings are the
unsupported-keyword warning in both cases, and WELLREF refined-summary (3) and
LGR well-average-pressure (1) warnings.

Run evidence: `runs/mpi-constraint-validation-2026-10-15/{simple,wellref}/`.
Because the warning policy was tightened during these runs, their outputs were
also reanalyzed with the final profile into
`runs/mpi-constraint-final-policy-2026-10-15/{simple,wellref}/`; both passed.
The reanalysis preserves the original run manifests rather than rewriting them.

All nine existing geomechanics regression tests passed (51.66 seconds); tracked
result files were restored afterward. Final native tests passed on np1/2/4, and
the Python harness passed all 11 checks.
