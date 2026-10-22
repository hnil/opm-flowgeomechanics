# Decomposition plan for Fracture / FractureModel (design note, 2026-08-25)

Current state: `Fracture` is 6,470 lines over four files with 114 member
variables and ~78 public methods; `FractureModel` is ~1,250 lines mixing
registry, orchestration and well plumbing. The 114 members fall into seven
natural clusters — the class is seven objects wearing one hat. The proposal is
**extraction into member structs/components, not a rewrite**: every step is a
mechanical move gated by the 9 REGTEST + 2 BANDTEST suite.

## Fracture → seven components

1. **`FractureGeometry`** — the mesh and its history:
   `grid_, grid_prev_, trimesh_, trimesh_prev_, grid_stretcher_*,
   grid_mesh_map_*, axis_, naxis_, origo_, cell_normals_, layers_, nlinear_,
   out_indices_`. Owns expansion/remap (`redistribute_values`, the
   expand_to_criterion glue) and the prev/current swap discipline that today is
   scattered ad-hoc — the swap becoming one method removes a whole bug class
   (the OOB remaps we fixed were exactly missed prev/current pairs).

2. **`ReservoirSample`** — everything sampled from the flow/mech state per cell:
   `reservoir_{cells,pressure,mobility,perm,dist,stress,cstress,density,cell_z}_,
   all_reservoir_*_, map_reservoir_*_, fracture_water_property_evaluator_`.
   One `resize(ncf)` + one `sampleFrom(simulator)`; kills the family of
   "vector sized to the old grid" bugs by construction (single resize site).

3. **`WellCoupling`** — the well side of the fracture:
   `wellinfo_, perfinj_, well_source_*, well_perf_cells_, well_perf_rate_,
   well_rate_(+mix), total_wellindex_(+mix), total_WI_well_, wi_dz_,
   wi_respress_, well_ref_depth_, perf_ref_depth_, perf_pressure_(+mix,
   +step_start), dp_perf_, density_perf_, mobility_water_perf_,
   well_control_is_rate_, effective_control_latched_, well_indices_(+accel)`.
   This is the piece the **generalized well-constraint row** will rebuild —
   extract it FIRST, so the row lands in a component instead of in the god
   class. `effectiveControlType()`, `setWellProps`, the WI upscaling and the
   CTF/`RuntimePerforation` production all live here.

4. **`FilterCake`** — `filtercake_{thikness,thikness_prev,perm,poro}_,
   has_filtercake_` + the update/rollback logic (already explicit-per-step;
   small, clean extraction).

5. **`FractureSolver`** — the width/pressure solve (most of
   `Fracture_fullSystemIteration.cpp` is already this in TU form):
   `A_, I_, S_, S_linop_, pressure_matrix_, coupling_matrix_,
   pressure_operator_, pressure_solver_, psolver_, frac_flow_precond_,
   rhs_pressure_, rhs_width_, closed_cells_, cell_flip_counts_,
   fracture_matrix_ (BEM), E_, nu_, last_solve_stats_, prmpressure_`.
   Interface: `solve(state...) -> stats`. The contact policies (sticky/FB) are
   strategy details inside it. Formalizing the TU boundary into a class is
   mostly moving the members it already monopolizes.

6. **`CouplingState`** — the outer-loop/per-step state (was CLEANUP_PLAN B5):
   `area_step_start_, volume_step_start_, perf_pressure_step_start_,
   width_round_prev_, width_round_valid_, coupling_dt_cap_,
   coupling_quiet_steps_, max_flow_time_step_` + the checkpoint capture
   (`moveForwardInTime`) and the dt-controller/onset-hold logic
   (`maxFlowTimeStep`, `onsetHoldActive`). Small, high-clarity win.

7. **`FractureOutput`** — `vtkwriter_, vtkmultiwriter_` + `writeIterationSnapshot`
   and the debug dumps. Pure sink, trivial to extract, immediately shrinks
   Fracture.cpp.

What remains in `Fracture`: identity (`prm_, active_, gravity_`), the solve
orchestration (the `if_propagate_trimesh` loop, which after extracting the
propagation-round hooks — old CLEANUP_PLAN B6 — reads as: sample → solve →
score front → expand → remap), and the physics glue (`solvePressure`,
`normalFractureTraction`, K1). Target ≈ 1,500 lines and ~15 members.

## FractureModel → registry + three helpers

`FractureModel` should keep only the registry (wells_, well_fractures_,
search tree) and per-step orchestration (solve loop with the failure policy,
moveForwardInTime fan-out, maxFlowTimeStep min-reduction). Extract:

- **`WellStateSampler`** (template) — today's `updateWellProperties` body:
  reads the well state (rates, WI, cmode, depths) and feeds each fracture's
  `WellCoupling`. This is also where the generalized row gets its inputs.
- **`ConnectionsUpdater`** — `addConnectionsToSchedual`/`assignGeoMechWellState`
  and the CTF-push gating (topology hysteresis): the one place that mutates the
  schedule/well state.
- **`SeedFactory`** — WSEED parsing + fracture construction (the seed-size /
  min_width wiring now documented at the creation site).

## Sequencing and constraints

1. Do NOTHING before the generalized well-constraint row is designed — then
   extract **WellCoupling first** and build the row inside it (one refactor,
   not two).
2. Then, in order of value/risk: FractureOutput (trivial) → CouplingState (B5)
   → ReservoirSample → FractureGeometry (+ propagation hooks, B6) →
   FractureSolver formalization. One extraction per commit, full gate each.
3. LGR: `FractureModel*`/`FractureMechHost` are on the `geomech_lgr` conflict
   surface — perform each extraction on `geomech-fracture-auxcells`, merge to
   `geomech-lgr-auxcells` immediately (small merges stay trivial; a batched
   refactor merge would not).
4. The BANDTESTs are the behavioural guard (the coupled runs are not
   run-to-run reproducible, so REGTEST bit-identity only covers the frozen
   legacy paths); a pure mechanical move should still keep REGTEST green.
