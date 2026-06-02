# Fracture parameter-file reference (`SeqMechFrac`)

Companion to `fracture_simple_bhp_seq_final.json`. JSON (boost property_tree) has **no comments**, so
the JSON uses harmless `"__…"` divider keys (ignored by the reader) and this file documents every option.

## 0. Recommended values (validated June 2026)

Validated on the SEQ regression deck (`SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ`, converges, full growth)
and the longer `model2` BHP/RATE decks (complete, 0 non-convergence, 0 linear-solver failures, fracture
grows to grid ~1556). The four committed regression tests still pass unchanged.

**The settings that actually matter:**
1. `closed_cell` toggle guard must stay effectively unlimited — `max_closed_cell_toggle_count=1e9`,
   `max_closed_cell_toggle_fraction=1.0`. Too-tight values freeze the contact set and the solve never
   converges. (Also fixed in code: the guard now throttles only re-openings, never contact growth.)
2. `reopen_width_tolerance=0.0` (force-driven reopen). Any value >0 locks closed cells (their width is
   pinned at 0) and suppresses fracture growth. (Also fixed in code.)
3. **Robust linear solver beats the fancy one.** `use_ad_pressure_assembly=false`,
   `drop_fluid_mech_linearization=true`, preconditioner `mode_policy=manual` / `mech_first=false`
   (mech-last), `linsolver.tol=1e-10`, `linsolver.max_iter=100`. The AD + `coupling_auto` + `mech_first`
   path converges fast on short cases but **fails the coupled linear solve on longer runs** (it aborts).
   The mech-last / decoupled path completes long runs with zero failures.
4. Physics on: `config.gravity_off=false`, `guard_state=true`.

Growth fidelity (your top priority) comes from (1)+(2); robustness on long runs comes from (3).


Two property trees are parsed from the same file:
- `getGeomechParam()` = the **whole** file → the outer SeqMechFrac driver reads `solver.*` here.
- `getFractureParam()` = the **`fractureparam`** subtree → the `Fracture` object reads `solver.*` (i.e. `fractureparam.solver.*`).

⚠️ Consequence: the CTF / perf-pressure change thresholds are read from **`fractureparam.solver.*`**, *not*
the top-level `solver` block. A `ctf_change_threshold` placed in the top-level `solver` block (as in the old
configs) is **dead**. It now lives under `fractureparam.solver` in the final config.

---

## 1. Top-level `solver` (outer SeqMechFrac driver — `BlackoilModelGeomech.hpp`)

| key | default | meaning |
|---|---|---|
| `method` | — | `SeqMechFrac` (sequential implicit: flow → mech → fracture, iterated) · `PostSolve` (solve fracture after the converged step) · `SeqMech` |
| `implicit_flow` | — | `true`: only do mech+fracture once the reservoir flow step has converged; gates outer iteration on `report.converged` |
| `max_mech_it` | 10 | max outer iterations in which the mechanics solve is still done |
| `max_frac_it` | 10 | max outer iterations in which the fracture solve is still done |

## 2. `fractureparam` (outer coupling / well-index update gating)

| key | default | meaning |
|---|---|---|
| `addconnections` | — | allow the fracture to add/modify well connections (CTF) in the schedule |
| `require_converged_fracture_for_wi_update` | true | only push new well indices to the reservoir when the fracture solve converged (avoids feeding garbage CTF) |
| `topology_hysteresis.enable` | false | enable hysteresis so connection-set changes don't chatter between outer iterations |
| `topology_hysteresis.on_threshold` | 1.0 | composite coupling-change norm above which a topology update is *armed* |
| `topology_hysteresis.off_threshold` | 0.3 | norm below which the pending update is disarmed |
| `topology_hysteresis.confirm_iterations` | 1 | consecutive armed iterations required before actually updating |
| `topology_hysteresis.cooldown_iterations` | 0 | iterations to suppress further updates after one is applied |
| `topology_hysteresis.emergency_threshold` | 10.0 | norm that forces an immediate update regardless of hysteresis |

## 3. `fractureparam.solver` — propagation / growth  *(group `__1`)*
**This is the qualitative behaviour to protect.** Controls trimesh expansion per step/iteration.

| key | meaning |
|---|---|
| `area_change_fac` | max fracture-area growth factor allowed per expansion step (growth limiter ↔ your timestep coupling) |
| `dt_limit` | growth/step limiter |
| `efac`, `rfac` | expansion heuristics (target-expansion gain `efac`; max grow as fraction of box length `rfac`) |
| `max_expand_iter` | max grid-expansion iterations per fracture solve |
| `max_iter_on_same_level` | cap on solves at one refinement level |
| `target_cellcount`, `cellcount_threshold`, `numcell_threshold`, `max_num_coarsening` | trimesh coarsening / size control |
| `smooth_boundary`, `full_intersections`, `divide_wellidx`, `no_leakof_outercells`, `remap_solution` | geometry / remap options (`remap_solution=true` carries width/pressure across grid changes instead of re-initialising) |

## 4. `fractureparam.solver` — inner nonlinear (width/pressure) solve  *(group `__2`)*

| key | meaning |
|---|---|
| `max_iter` | max Newton iterations per fracture solve (non-convergence → warning, or throw if `failure_on_nonconvergence`) |
| `tolerance` | convergence tol on the (mech, flow) residual infinity-norms |
| `damping` | Newton step scaling (1 = full step) |
| `min_width`,`max_width` | aperture bounds used in flow transmissibility |
| `max_dwidth`,`max_dp` | per-iteration clamps on width / pressure update (stability) |
| `max_change` | misc step clamp |
| `force_limit` | traction floor in closed-cell force test |
| `failure_on_nonconvergence` | throw (true) vs warn-and-continue (false) when `max_iter` is hit |

## 5. `fractureparam.solver` — closed-cell (contact) handling  *(group `__3`)*
**Drives both convergence and growth.** A closed cell has its mechanics row replaced by `w=0`.

| key | default | meaning |
|---|---|---|
| `closing_type` | `org` | `org` (force = rhs − A·w − I·p) · `simple` (per-face force vs pressure) · `open` (never close) |
| `closed_cell_policy` | `sticky` | `legacy` (stateless test each iter) · `sticky` (hysteresis: stay closed unless reopen test passes) |
| `close_force_tolerance` | 0.0 | compressive force above which an open cell closes |
| `reopen_force_tolerance` | =close_tol | (sticky) opening force below which a closed cell may reopen |
| `reopen_width_tolerance` | 0.0 | **(sticky) keep at 0.** A closed cell's width is pinned at 0, so any value >0 makes the reopen test `width≥tol` unsatisfiable → cells lock closed → fracture under-grows. 0 ⇒ force-driven reopen. |
| `max_closed_cell_toggle_count` | 1e9 | **keep large.** Max cells allowed to flip open/closed per iteration; too small (e.g. 32) freezes the closed set at empty and the solve never converges |
| `max_closed_cell_toggle_fraction` | 1.0 | **keep at 1.0.** Same guard as a fraction of N; effective limit is `min(count, fraction·N)` |

## 6. `fractureparam.solver` — fracture flow–mech coupling block `C`  *(group `__4`)*

| key | default | meaning |
|---|---|---|
| `drop_fluid_mech_linearization` | true | `true` ⇒ always zero the coupling block `C` (decoupled / Picard). `false` ⇒ keep `C` (true Newton) subject to the drop tolerances below |
| `drop_tol_h` | 1e1 | drop `C` when the **mech traction** infinity-norm exceeds this. The traction is ~confining stress (O(1e6–1e7 Pa)), so values like 100 mean `C` is *always* dropped even with `drop_fluid_mech_linearization=false`. To actually keep coupling, set ≳1e8 — **but** the coupled linear solver currently fails on these systems (see notes). |
| `drop_tol_p` | 1.0 | same, on the pressure residual |
| `use_ad_pressure_assembly` | false | use the consistent AD assembler for the pressure+coupling matrices |

## 7. `fractureparam.solver` — state guards  *(group `__5`)*

| key | default | meaning |
|---|---|---|
| `guard_state` | true | validate the nonlinear state each iteration (throw on NaN / blow-up) |
| `guard_max_abs_width` | 1e30 | max plausible aperture before declaring blow-up |
| `guard_max_abs_pressure` | 1e30 | max plausible fracture pressure |

## 8. `fractureparam.solver` — well-index / perf-pressure coupling stabilization  *(group `__6`)*
**These are the accelerators for the outer (well-index-from-fracture) convergence.**
Change-detection thresholds decide *when* the reservoir well CTF/perf-pressure get updated from the fracture;
the mixing options decide *how* (damped / accelerated) the perf-pressure and (optionally) well-index are blended
across outer iterations. Implemented in `Fracture::applyCouplingUpdate`.

| key | default | meaning |
|---|---|---|
| `ctf_change_threshold` | 1e-15 | abs CTF change above which connections are considered changed (triggers schedule/WI update) |
| `ctf_change_threshold_rel` | 0.0 | optional relative CTF change gate (AND-ed with abs when >0) |
| `ctf_relative_scale` | 1.0 | floor scale for the relative CTF change |
| `perf_pressure_change_threshold` | 1e30 | abs perf-pressure change gate (1e30 ⇒ effectively never gates on perf pressure) |
| `perf_pressure_change_threshold_rel` | 0.0 | relative perf-pressure gate |
| `perf_pressure_relative_scale` | 1.0 | floor scale for relative perf-pressure change |
| `coupling_update_mode` | `legacy` | how perf-pressure (and WI, if enabled) are updated each outer iter: `legacy` (fixed damped average) · `aitken` (adaptive relaxation ω) · `anderson` (1-history Anderson mixing, Aitken/legacy fallback) |
| `damping_factor_perf` | 2.0 | legacy damping weight for perf-pressure: `new=(old·d+target)/(1+d)` |
| `damping_factor_wi` | 2.0 | legacy damping weight for well-index / CTF averaging |
| `enable_wi_coupling_update` | false | if true, also relax `well_rate`/`total_wellindex` through `applyCouplingUpdate`; if false they are taken directly |
| `coupling_relax_init` | 1/(1+damping) | initial relaxation ω for aitken/anderson |
| `coupling_relax_min`,`coupling_relax_max` | 0.05, 1.25 | clamp on ω |
| `coupling_acceptance_factor` | 1.0 | accept an accelerated step only if it reduces the residual by at least this factor, else fall back |
| `coupling_fallback_legacy` | true | on rejection, fall back to the legacy damped update |
| `coupling_anderson_alpha_min`,`coupling_anderson_alpha_max` | -1.0, 2.0 | clamp on the Anderson mixing coefficient |
| `coupling_legacy_first_damp` | false | emulate pre-`0e6bb48` behaviour (damp even on the first call) |
| `legacy_coupling_change_logic` | false | revert to old "always update connections" behaviour (bypasses hysteresis + convergence gate) |
| `legacy_parent_setup_iteration` | false | **physically incorrect** legacy reset of the parent first-iteration after a connection update; `false` preserves reservoir state (`runParentFirstIterationPreservingState`) |

### Practical guidance for the WI stabilization
- Default (`legacy`, `damping_factor_perf=2`) is the proven, stable baseline.
- To **accelerate** the outer well-index/perf-pressure convergence, try `coupling_update_mode=anderson`
  (keep `coupling_fallback_legacy=true`, `coupling_acceptance_factor≈1`). It only changes *rate*, not the
  fixed point, and falls back to the damped update when it would overshoot.
- `enable_wi_coupling_update=true` additionally relaxes the well-index itself — use if the CTF (not just
  perf pressure) oscillates between outer iterations.

## 9. `fractureparam.solver` — diagnostics  *(group `__7`)*

| key | default | meaning |
|---|---|---|
| `verbosity` | 0 | 2 ⇒ per-iteration `Nonlinear iteration: K rhs: <mech> <flow> closed_cells: n`; 3 ⇒ x/dx norms; >2 ⇒ coupling-update logging |
| `dump_linear_system_on_failure` | true | **off in the final config**; dumps `.mtx`/`.txt` snapshots on a linear-solver failure (large, slow) |
| `write_coupled_linear_system` | false | dump every coupled system (debug only) |

## 10. `fractureparam.solver.linsolver` (coupled fracture linear solve)
`solver` (`bicgstab`/`gmres`/`fgmres`), `tol`, `atol`, `max_iter`, `verbosity`, `restart`, and
`preconditioner` (`FractureMechanicsPreconditioner`): `mode_policy` (`manual`/`coupling_auto`),
`mech_first`, `fixed_stress`, `mode_switch_coupling_threshold`, `diag_mech`, `diag_flow`,
`mech_press_coupling`, `update_mech_on_reuse`, and a nested `flow_solver` (Dune FlexibleSolver) block.
`coupling_auto` estimates a coupling indicator and switches to `fixed_stress` when it exceeds the threshold.
