# WELLREF partitioning and rank-specific validation

The VEM geomechanics executable defaults to `--refine-before-redistribute=false`.
This invokes the existing conforming LGR partition-cell-group path: each outer
CARFIN region, its overlap halo, and touching groups stay on one rank, using
Zoltan graph partitioning and at least two overlap layers. The mesh and physical
settings are not enlarged or retuned.

The generic conforming auto policy refines before distributing a flat leaf and
deliberately splits refinement boxes. The fracture completion-to-schedule mapper
uses per-level grids. An explicit parallel fracture/LGR request for the flat path
is rejected during initialization with an actionable message; serial use is
unaffected. Supporting distributed fractures on a split flat leaf is separate work.

## Validation

Run SIMPLE, one-rank WELLREF and four-rank WELLREF with the same dependency matrix,
fracture JSON, mechanics solver and input deck. Require completion, growth,
conduction, zero missing-completion/mechanical-BC diagnostics, and numerical checks.
An explicit unsupported flat-path request must fail before a timestep.

## Refinement-quality follow-up

The colleague request for nicer refinement remains an explicit separate review:

- Measure cell aspect ratios, skew/degeneracy and the mechanics geometric-grid check.
- Check refinement transitions and clearance from the grown fracture.
- Measure cell count, partition imbalance, halo overhead and peak memory.
- Compare pressure, conduction and fracture geometry with preserved baselines.

Propose a new quality option only after measurements define its behavior and
tradeoffs. The partitioning switch controls ownership/order, not a smoother
transition or a new mesh-quality algorithm.
