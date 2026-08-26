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
- `getGeoMechParam()` = the **whole** file → the outer SeqMechFrac driver reads `solver.*` here.
- `getFractureParam()` = the **`fractureparam`** subtree → the `Fracture` object reads `solver.*` (i.e. `fractureparam.solver.*`).

⚠️ Consequence: the CTF / perf-pressure change thresholds are read from **`fractureparam.solver.*`**, *not*
the top-level `solver` block. A `ctf_change_threshold` placed in the top-level `solver` block (as in the old
configs) is **dead**. It now lives under `fractureparam.solver` in the final config.

---

## 1. Top-level `solver` (outer SeqMechFrac driver — `NonlinearSystemBlackOilReservoirGeoMech.hpp`)

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
| `min_width`,`max_width` | numerical width clamps inside the width solve. NOT the flow floor — that is `config.min_width`, and for deck runs it is SET BY THE DECK: `WSEED` item 10 (WIDTH, default 1e-4 m) overrides the JSON (`FractureModel.cpp`, logged at seed creation). Recommended deck value 3e-4 (≈ Reveal's conductivity floor); `config.min_width` in JSON matters only for standalone/test setups without a deck `WSEED`. |
| `max_dwidth`,`max_dp` | per-iteration clamps on width / pressure update (stability) |
| `max_change` | misc step clamp |
| `force_limit` | traction floor in closed-cell force test |
| `failure_on_nonconvergence` | throw (true) vs warn-and-continue (false) when `max_iter` is hit |

## 5. `fractureparam.solver` — closed-cell (contact) handling  *(group `__3`)*
**Drives both convergence and growth.** A closed cell has its mechanics row replaced by `w=0`.

| key | default | meaning |
|---|---|---|
| `closing_type` | `org` | `org` (force = rhs − A·w − I·p) · `simple` (per-face force vs pressure) · `open` (never close) |
| `closed_cell_policy` | `sticky` | `legacy` (stateless test each iter) · `sticky` (hysteresis: stay closed unless reopen test passes) · `pdas` (semismooth-Newton / primal-dual active set: one consistently-scaled criterion `closed ⟺ tmp_i − \|A_ii\|·w_i > tol`). `pdas` cuts contact chatter ~20% on the decoupled path with identical growth; it does NOT make keep-C true Newton converge (that non-convergence is the contact complementarity itself, not the active-set rule). |
| `close_force_tolerance` | 0.0 | compressive force above which an open cell closes |
| `reopen_force_tolerance` | =close_tol | (sticky) opening force below which a closed cell may reopen |
| `reopen_width_tolerance` | 0.0 | **(sticky) keep at 0.** A closed cell's width is pinned at 0, so any value >0 makes the reopen test `width≥tol` unsatisfiable → cells lock closed → fracture under-grows. 0 ⇒ force-driven reopen. |
| `max_closed_cell_toggle_count` | 1e9 | **keep large.** Max cells allowed to flip open/closed per iteration; too small (e.g. 32) freezes the closed set at empty and the solve never converges |
| `max_closed_cell_toggle_fraction` | 1.0 | **keep at 1.0.** Same guard as a fraction of N; effective limit is `min(count, fraction·N)` |

## 6. `fractureparam.solver` — fracture flow–mech coupling block `C`  *(group `__4`)*

| key | default | meaning |
|---|---|---|
| `drop_fluid_mech_linearization` | true | `true` ⇒ always zero the coupling block `C` (decoupled / Picard). `false` ⇒ keep `C` (true Newton) subject to the drop tolerances below |
| `drop_tol_h` | 1e1 | drop `C` when the **mech traction** infinity-norm exceeds this. The traction is ~confining stress (O(1e6–1e7 Pa)), so values like 100 mean `C` is *always* dropped even with `drop_fluid_mech_linearization=false` — setting that flag alone reproduces drop-C bit for bit. To actually keep coupling, lift this (≳1e8) **and** `drop_tol_p`. ⚠️ The old note here said the coupled solver fails on these systems; that was measured before the 2026-08-26 out-of-bounds fix and no longer holds — keep-C now runs with zero linear-solve failures and fewer nonlinear iterations than drop-C on both the SIMPLE deck and 92-day model2 (see `fixed_tests/NEXT_WORK.md` item 1c). Defaults unchanged pending a full matrix re-run. |
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
| `wi_upscaling` | `legacy` | KEEP BOTH options; the default stays `legacy` (the safe one). `legacy` zeroes CTF contributions where p_well < p_cell — a discontinuity, but it implicitly limits over-injection through cells the fracture cannot actually feed. `conductivity` is continuous but MISSES the in-fracture pressure drop (it upscales as if the whole fracture sat at the perforation pressure), over-stating conduction on long fractures. The proper fix arrives with the well-in-the-loop integration, where the fracture pressure profile enters the well equations directly. |
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
| `coupling_convergence_mode` | `threshold` | how the OUTER SeqMechFrac loop declares fracture-coupling convergence: `threshold` (legacy booleans from the change thresholds above) · `residual` (opt-in: converged iff inner fracture solve converged AND no structure change AND max relative change of any CTF / perf pressure ≤ `coupling_tolerance`; logged each outer iteration) |
| `coupling_tolerance` | 1e-3 | dimensionless outer residual tolerance used by `coupling_convergence_mode=residual` |
| `wi_vector_acceleration` | `none` | opt-in acceleration of the per-cell **CTF vector** fixed point (the quantity `wellIndicesAvrg` damps): `none` (legacy damped average) · `aitken` (vector Aitken relaxation ω from the residual history) · `anderson` (depth-1 Anderson mixing of the last two targets). Reuses `coupling_relax_*`/`coupling_anderson_alpha_*`/`coupling_acceptance_factor`/`coupling_fallback_legacy`; CTFs kept ≥ 0; rejected steps fall back to the legacy damped update; history reset at timestep boundaries and on fracture reset |
| `value_only_wi_update` | false | opt-in cheap connection-update path in the outer loop: when the coupling change has **no structure change** (same cells, only CTF values moved), write the CTFs in place on the existing well perforations and skip the schedule rebuild + well-model `beginTimeStep` + `prepareTimeStep` + parent setup iteration. Structure changes always take the full legacy path |
| `k1_extraction` | `aperture_cell` | how the propagation criterion's K1 is evaluated at front cells: `aperture_cell` (legacy: single boundary-cell width at its centroid distance — mesh sensitive) · `displacement_correlation` (opt-in: least-squares fit of the asymptotic opening w(s)=C·√s over a corridor behind the front edge, converted with the same calibrated constant; reduces to legacy for a single point, falls back to legacy when <2 corridor cells) |
| `toggle_guard_mode` | `count` | anti-chatter strategy for the contact set: `count` (legacy: global re-opening count vs `max_closed_cell_toggle_*`; **WARNING:** a finite count limit wholesale-rejects legitimate en-masse re-opening — observed ~100× fracture-area suppression on the refine decks with 32/0.1) · `chatter` (opt-in: per-cell flip counting within a solve; first state changes always pass so monotone opening is never blocked; a cell is pinned only after `chatter_flip_limit` changes — true chatter) |
| `chatter_flip_limit` | 3 | per-cell contact state changes allowed per solve before the cell is pinned (`toggle_guard_mode=chatter`) |
| `k1_correlation_length` | 0.0 | corridor length L in meters for `displacement_correlation`; 0 ⇒ auto = `k1_correlation_cells` × local cell size |
| `k1_correlation_cells` | 3.0 | corridor length in local cell sizes when `k1_correlation_length` is 0 |

### Practical guidance for the WI stabilization
- Default (`legacy`, `damping_factor_perf=2`) is the proven, stable baseline.
- To **accelerate** the outer well-index/perf-pressure convergence, try `coupling_update_mode=anderson`
  (keep `coupling_fallback_legacy=true`, `coupling_acceptance_factor≈1`). It only changes *rate*, not the
  fixed point, and falls back to the damped update when it would overshoot.
- `enable_wi_coupling_update=true` additionally relaxes the well-index itself — use if the CTF (not just
  perf pressure) oscillates between outer iterations.
- `wi_vector_acceleration`: use **`anderson`**; do **not** use `aitken`. Validated on
  vegard model2_refpad (native JSON, 64-cell CTF coupling, 334 days): anderson cut the inner
  fracture NL iterations 19228 → 6813 (−65%) and the inner non-convergence warnings 208 → **0**
  at equal solve time, with physics within a few % of the (partially unconverged) baseline and
  conduction ~94%. Aitken's single relaxation factor over the sign-mixed CTF vector *destabilizes*
  the same case (53099 NL iters, 774 non-convergence warnings). On decks whose CTF fixed point is
  already trivial (e.g. model2 rate_short, SIMPLE SEQ: residual ~1e-15) the acceleration correctly
  stays dormant and results are bit-identical.

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

## 11. New robustness / performance levers (all opt-in, defaults reproduce old behaviour)

| key | location | default | meaning |
|---|---|---|---|
| `linsolver.scaling` | `fractureparam.solver.linsolver` | `none` | `blockmax`: symmetric block scaling `D S D` of the coupled system (s_m=1/√max\|A\|, s_p=1/√max\|M\|) before the Krylov solve. Strong, cheap conditioning win; dumps stay unscaled (physical). |
| `linsolver.failure_policy` | `fractureparam.solver.linsolver` | `throw` | `ladder`: on linear non-convergence, try a rescue ladder (configured→fgmres→fixed-stress preconditioner→**drop-C / decoupled**) before giving up. The drop-C rung is essentially always solvable, so the run progresses instead of aborting. Since 2026-08-24 a Krylov *abort* (bicgstab breakdown / `defect=nan` from overflow on the ill-scaled system) also falls into the ladder instead of failing the nonlinear iteration — this plus `scaling=blockmax` + `preconditioner.lu_pivoting` is what makes the iterative rate_well solve survive establishment (109 s vs direct 458 s on the caked SHORT10D, same answer). With `throw` (default) aborts still rethrow unchanged. |
| `linsolver.preconditioner.lu_pivoting` | `…linsolver.preconditioner` | `false` | Partially pivoted dense LU for the mechanics block of `FractureMechanicsPreconditioner` (+ a named recoverable error on a zero/non-finite pivot). The unpivoted legacy LU can overflow when O(1) pinned closed-cell rows mix with O(1e9) BEM rows. Off by default (pivoting changes round-off → REGTEST bit-identity); on in `seq_implicit_ratewell.json`. |
| `linsolver.debug_checksum` | `fractureparam.solver.linsolver` | `false` | Per coupled solve, print a bit-exact hash of the system matrix, the residual and the preconditioner's action on a fixed probe. Two runs of the same case must print identical streams; the first differing line localizes a reproducibility break to assembly, residual or preconditioner. (Found the 2026-08-26 out-of-bounds `dR/dw` read this way; see also `-DGEOMECH_ASAN=ON`.) |
| `linsolver.solver` | `fractureparam.solver.linsolver` | `bicgstab` | also accepts `direct`: assemble the dense 2N×2N coupled system and LU-solve it. For SMALL fractures this is a robust *linear* solve (the dense-BEM coupled system has cond ~1e11 and defeats Krylov; the only good iterative preconditioner — exact Schur — is dense, i.e. same cost as direct). Gated by `direct_max_cells` (default 4000); above it, falls back to `direct_fallback_solver` (default `fgmres`). NOTE: a robust linear solve does **not** by itself make keep-C *true Newton* converge — the nonlinear coupled step is erratic for this contact problem; decoupled Picard (drop-C) remains the robust path. |
| `keep_zero_ctf_topology` | `fractureparam.solver` | `false` | Keep fracture-generated dynamic completions in the schedule even when their current fracture WI/CTF is non-positive. Default `false` skips these topology-only connections because an open zero-CF schedule completion is **not** inert: it still enters `initializeWellPerfData()` / `initializeWellState()`, changing well topology, perf ordering and first-perforation pressure seeding before the fracture WI is added later. Enable only for topology studies or debugging. |
| `conservative_propagation` | `fractureparam.solver` | `false` | Don't grow the fracture (trimesh `if_propagate_trimesh`) off a non-converged inner solve — K1 from an unconverged state is unreliable, and opening on it can irreversibly alter the result. |
| `retract_established_compressed` | `fractureparam.solver` | `false` | In `removeNewZeroWithCells`, also retract **established** cells (not just newly-grown) that end zero-width/compressed, instead of only closing them. ⚠️ lets the fracture shrink → can oscillate against growth; experimental. |
| `dump_case_min_nonlin_iter` | `fractureparam.solver` | `0` | If >0, dump a physical (unscaled) coupled-system snapshot when a solve reaches this nonlinear-iteration count — captures a *relevant* hard case for offline study (`replay_fracture_linear_system`, incl. `--contact-report`). |
| `mech_skip_coupling_threshold` | **top-level `solver`** (NOT `fractureparam.solver` — read by the SeqMechFrac driver, like `max_mech_it`) | `0` | If >0, skip the outer mechanics solve when the previous outer iteration's fracture→well coupling `composite_norm` is below this, reusing the prior stress. Mechanics depends on the fracture only via well flow, so this is exact when the wells are stable. Self-correcting. **Caveat:** a purely thermal stress transient that doesn't move the wells is not detected — keep the threshold conservative on strongly thermal/rate-controlled decks. Validated on the BHP SEQ deck: skipped 9/14 mech solves, faster, bit-identical fracture area/volume. **Calibration:** `composite_norm` scale ~0..1e3. Use a CONSERVATIVE threshold (~1.0 = skip only near-zero coupling change): on model2 RATE that cut mech solves 184→76 (−59%) with the fracture solution preserved to ~1%. An AGGRESSIVE threshold (e.g. 1e30) diverges on rate-controlled decks (area +50%) — never use it there. MPI-consistent (decision uses the global comm.max'd norm). |
| `max_growth_iterations` | **top-level `solver`** (SeqMechFrac driver) | `0` | If >0, when the outer iteration's fracture solve leaves the coupling unconverged (typically: the fracture is still propagating during establishment), iterate [flow-to-convergence → mech → fracture + connection update] **inside the same outer Newton iteration**, up to this many rounds, with the Newton iteration counter frozen (`SetupIterationContextGuard`). This makes the trajectory insensitive to `area_change_fac` (a smaller extension cap only costs more rounds) and stops establishment from exhausting `NewtonMaxIterations` — whose dt chops cannot help, since the required growth is a ratio per solve, not per unit time. Respects `max_mech_solves_per_step`. Validated on the model_1208 model2 RATE deck: the first timestep completes at full size instead of dt-death-spiralling at day 1.7. |
| `growth_flow_iterations` | **top-level `solver`** | `12` | Flow Newton iterations allowed per growth round (`max_growth_iterations`) to reconverge the flow after the round's connection update; if flow does not reconverge within this, the growth loop returns control to the outer Newton loop. |
| `coupling_state_relaxation` | `fractureparam.solver` | `1.0` | Damped fixed point on the fracture **width state** exposed to the coupling (ω = weight of the new solve; 1 = off). WI, volume and K1 all derive from the width, so this is a trust region on the whole coupling round — it tames the contact open/close bang-bang between rounds without the conduction delay that CTF damping causes. On the no-cake model_1208 short deck ω=0.5 reproduced the undamped answer with ~65 % fewer growth rounds; on caked decks it lowers the BHP overshoot (546→460 with accept) but does not by itself break the contact limit cycle. History is per-timestep (reset in `moveForwardInTime` and on rollback); blending is skipped on rounds where the mesh changed size. |
| `growth_driven_switch` | **top-level `solver`** | `false` | With `max_growth_iterations`: keep iterating only while a fracture is at its per-solve growth cap (finite `maxFlowTimeStep`); otherwise accept the lagged coupling (PostSolve-style) — the seq-implicit rounds are spent only where growth is happening. Essential partner of the contact hysteresis (`close_force_tolerance`/`reopen_force_tolerance`): hysteresis alone ratchets (476 rounds, 41 000 m²), switch alone still flaps once the fracture is large; together they gave the robust caked-model2 result (`best_practice/seq_implicit_robust.json`). |
| `disable_propagation` | `fractureparam.solver` | `false` | Diagnostic — freezes the trimesh at the seed (no expansion ever). For isolating the width/contact/leak-off machinery from the growth criterion, e.g. against a fixed-geometry reference. Note the WSEED seed (~0.3 m²) conducts essentially nothing, so on a rate-controlled caked well this runs the BHP away by design. |
| `propagate_open_front_only` | `fractureparam.solver` | `false` | Only allow expansion at front cells that are OPEN in the current contact set — a closed front cell carries no mechanical opening load, so growth from it would ride on the width-floor conductivity rather than on K1. Pair with `conservative_propagation` so new cells are only added off a converged open/closed set. Establishment behaviour under both gates should be validated per case (during establishment everything is open, so the open-front gate is usually inert there). |
| `onset_dt_limit` | `fractureparam.solver` | `0` | If >0 (days): via the `maxFlowTimeStep` hook, hold dt at this value while the fracture is still an unopened seed (area < `onset_seed_area`, default 2 m²) and the perforation pressure rose more than `onset_pressure_rise` (default 1 bar) since the step checkpoint. Makes the seed be evaluated near its opening pressure instead of after the well has run away (caked model2: opens at 430–460 bar instead of 977) — without prescribing a small dt globally. |
| `max_area_growth_per_step` | `fractureparam.solver` | `0` | If >0: per-step growth guard — after the seq-implicit outer loop converges, if any active fracture's area exceeds `max(factor·A0, A0 + min_area_growth_guard)` (A0 at the step checkpoint; absolute floor default 200 m² so a seed can open) a `NumericalProblem` is raised inside the retry scope → dt chopped, fracture rolled back from checkpoint. In `endTimeStep` (PostSolve) it can only WARN (post-acceptance). Gives a Reveal-paced onset (0.4→173→225→282→435→791 m² at 440–490 bar) but at ~30× cost and it also blocks legitimate growth later (stall, BHP 842) — the fracture solve grows by a fixed ratio per solve regardless of dt, so a per-step guard can only work by making dt tiny. Use as a diagnostic of over-growth; prefer `volume_paced_propagation`. |
| `volume_paced_propagation` | `fractureparam.solver` | `false` | Mass-balance pacing of the expansion rounds inside one fracture solve: stop expanding when the fracture volume created in this step exceeds what the well can deliver minus leak-off, `V(k)−V0 > (q − Q_leak(k))·dt·f`. See `opm/geomech/VOLUME_PACED_PROPAGATION.md` for the formulation, the first test and its caveats. Partner keys: `volume_pacing_factor` (`f`, default 1), `volume_pacing_pressure` (also draw the BC pressure down by the over-draw over the secant storage compliance, eq. 3–4), `volume_pacing_source` (`perf` = this perforation's CTF share — circular for an establishing fracture, stalls at 252 m²; `well` = whole well rate — the natural single-fracture statement), `volume_pacing_verbosity`. Storage is never the limiter (V ≤ 3 m³ vs q·dt of 10²–10³ m³); the bound is carried by leak-off, so its magnitude inherits the fracture solve's leak-off at high dp — which is what stalls the caked model2 at ~1000 m² for days 2–7. |
| `write_iteration_snapshots` | `fractureparam.solver` | `false` | Diagnostic: after every fracture solve in the seq-implicit outer loop, append one line to `<iteration_snapshot_dir>/<fracture>_iter_summary.csv` (step, round, tag outer/growth, perf pressure, well & perf rate, ncells, area, open area, volume, Q_leak, WI sum, mean p_frac, max width, n_closed) and, if `iteration_snapshot_cells` (default true), per-cell rows to `<fracture>_iter.csv` (x,y,z, area, p_frac, p_res, traction, width, closed, K1, leak-off, WI). Lets the open/closed branches the coupling visits be plotted round by round (see `model_1208/.../README.md`: the (p_frac, Q_leak/q) branch diagram shows there is no open state with Q=q under a fixed-pressure BC). |
| `fractureWI_mode` | `fractureparam.solver` | `constant` | How the well→fracture connection strength `fWI` per well-source cell is set: `constant` (the config `fractureWI`, legacy) · `estimate` (radial flow in the fracture plane towards the wellbore with cubic-law permeability, `WI = 2π(w³/12)/ln(r_e/r_w)`, `w = max(current width, min flow width)`, `r_e` from the cell area, `r_w = fractureWI_rw` default 0.1 m — the same estimate the aux-cell path uses as `embedded_perf_wi_mode=estimate`). On the model_1208 decks the seed connection is not the limiter (100–2000× the fracture's leak-off conductance), so this mostly matters as physics/consistency, not stability. |
| `well_source_all_perfs` | `fractureparam.solver` | `false` | Feed the fracture from **every** reservoir cell the well perforates (each fracture cell lying in a perforated reservoir cell becomes a well source), not only the seed ring. Needed when the fracture grows across a multi-perforation interval: with the seed-only sources the other perforations stay pure matrix perforations and (under `rate_well`) still take most of the rate — on model2 INJ1 (4 perfs, seed on k=13) ~60 % went through caked matrix perfs while Reveal's 31-m-tall fracture takes everything. Also sizes the implicit-BCRS overflow for the many-entry well row. |
| **`WSEED` `WIDTH` item (deck!)** | deck | `1e-4` | Not a JSON key: `FractureModel.cpp:574` puts the WSEED width into `config.min_width`, the cubic-law flow-width floor, overriding the JSON. The decks that give only two sizes get the default 0.1 mm; Reveal's cubic law carries a 2.39 D·m conductivity floor ≡ 0.31 mm. With 0.1 mm the fracture beyond a few metres is hydraulically starved and collapses → the onset bistability; **1 mm gives a flat, Reveal-like onset** (model2: 418–431 bar, no snap). Set the WIDTH item explicitly. |
| `control.type = rate_well` | `fractureparam.control` | `perf_pressure` | **Recommended for single-fracture, rate-controlled injectors** (`best_practice/seq_implicit_ratewell.json`), *not* the default. Adds one DOF `p_w` with `Σ_matrix WI·λ·(p_w−p_res) + Σ_sources fWI·(p_w−p_frac) = q_well`, so the fracture pressure is *solved for* under the rate constraint: it lands on the near-closure, leak-off-equilibrium state (closure+8 bar, ~500 m² on model2) by itself, is insensitive to the flow-width floor (0.3 mm ≡ 1 mm), and gives Reveal's decaying BHP shape after opening — the fixed-pressure BC cannot reach that state and needs a 3× floor to compensate. Use the **iterative** coupled solver (bicgstab, mech-last, umfpack flow block; 4.6× faster than `direct`): the earlier NaN was the establishment death spiral, cured by `max_growth_iterations`+`growth_driven_switch`, not a solver defect. **Known problems, be explicit:** (1) the well row is a *duplicated, static, single-well, single-phase* well model — matrix WI/λ/p_res frozen from the last flow solve, λ=1, one fracture per well, not in sync with the well model's BHP/THP/group logic; the clean version is the well perforating fracture DOFs (embedded aux-cell path). (2) The fracture is fed only where it exists: on multi-perforation wells the perforations outside the fracture remain matrix perfs and take their share (model2 4-perf: 30–40 %), holding the BHP 40–60 bar above Reveal until the fracture is tall enough. (3) The well DOF sits high (~650 bar) until the fracture is big because the caked matrix perfs carry the rate — physics, not a bug, but it looks like late opening. (4) One day of testing vs months for perf_pressure; not in the regression suite. |

Diagnostics added to the per-solve stats log (`eclgeomechmodel.hh`): `closed_cell_toggles`
(contact chatter), `linear_solve_failures`, `ladder_rescues`. The contact/opening state of a
dumped case can be inspected offline with `replay_fracture_linear_system <prefix> --contact-report`.

| `control.type = well` | `fractureparam.control` | — | GENERALIZED WELL ROW (v2): the well DOF p_w is kept for every constraint; the single well equation expresses the currently ACTIVE one — rate-like (RATE/RESV/GRUP) → the rate balance (`rate_well`), pressure-constrained (BHP) → the constraint row scale·(p_w − BHP@perf-datum) = 0 (`bhp_well`, same scale as the rate row's WI diagonal; cell↔well fWI couplings identical in both). A mid-run control switch changes ONE ROW's content — no structural change (validated: 29 switches in one run, zero events). THP still throws (needs a VFP linearization). `bhp_well` can also be set explicitly. Physics note: a fixed BHP far above closure has no finite fracture equilibrium — growth continues until a rate cap or pressure support limits it; on the CASE_PRESS dev deck the fracture's injectivity at 420 bar immediately exceeds the deck's 200,000 m³/d cap, so the well correctly ends up rate-capped. |
