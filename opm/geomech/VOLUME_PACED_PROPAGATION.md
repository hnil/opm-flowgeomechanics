# Volume-paced fracture propagation (`solver.volume_paced_propagation`)

Status: opt-in, off by default (2026-08-15). Implemented in `Fracture_impl.hpp`
(the `if_propagate_trimesh` branch) as an additional stop criterion for
`expand_to_criterion`, plus an optional draw-down of the fracture-pressure BC.

## Why

The fracture solve propagates by repeated *expansion rounds* inside one call:
solve the coupled width/pressure system on the current mesh, evaluate
K1/K1c on the front, add every front cell above criterion, repeat. Every round
is evaluated against the **same** perforation pressure `p_perf` (set once from
the flow state before the solve). Nothing inside the loop knows how much fluid
the well can actually deliver in the timestep, so the fracture grows by a fixed
*ratio per solve* (`area_change_fac`) regardless of `dt`:

    435 → 1127 m²  in a 1.34-day step,   435 → 791 m²  in the 0.44-day retry

(caked model2, `model_1208`). Reveal, on the same deck, grows 0.3 → 100 → 389 m²
over its first three steps, at the same opening pressure (~500 bar). The
missing physics is that a fracture can only open the volume it is fed: a round
that would need more fluid than the well supplies in `dt` must draw the
fracture pressure down, which stops the front. This note states that balance
at the timestep scale.

## Quantities (all per fracture, SI)

| symbol | meaning | where it comes from |
|---|---|---|
| `dt` | length of the current timestep [s] | `simulator.timeStepSize()` |
| `q_perf` | reservoir-volume rate the well delivers into this fracture's perforation [m³/s] | `well_perf_rate_` (`setPerfProps`) |
| `V(k)` | fracture pore volume after expansion round `k`, `Σ_cells area·width` [m³] | `calculateFractureProperties().volume` |
| `V0` | fracture volume at the start of the timestep (checkpoint) [m³] | stored in `moveForwardInTime()` |
| `Q_leak(k)` | total leak-off out of the fracture after round `k`, `Σ_cells leakof·(p_frac−p_res)` [m³/s] | `calculateFractureProperties().flux` |
| `V_avail` | fluid volume available for opening new fracture volume in this step [m³] | see below |

## Balance

Fluid entering the fracture during the step either leaks off or is stored as
new fracture volume:

    q_perf · dt  =  Q_leak · dt  +  (V(k) − V0)

so the volume that the expansion may create in this step is bounded by

    V(k) − V0  ≤  V_avail := max( 0,  ( q_perf − Q_leak(k) ) · dt ) · f_vol      (1)

with a dimensionless safety/tuning factor `f_vol` (`solver.volume_pacing_factor`,
default 1.0). `Q_leak(k)` is the leak-off *of the current round's solution*, so a
larger, more leaky fracture is credited with less storage — a fracture that
already leaks off everything the well gives it cannot grow further this step.

## Stop criterion (implemented)

After every expansion round `k` (i.e. every call of the score function),
compute `V(k)` and `Q_leak(k)` and stop expanding when (1) is violated:

    stop_volume(k)  :=  ( V(k) − V0 )  >  V_avail(k)                             (2)

The criterion is OR-ed with the existing area-ratio criterion
(`expantionMax`, `area_change_fac`), so the old cap still bounds any single
solve; the new one is the one that makes growth per step scale with `dt`.
When (2) triggers, `max_flow_time_step_` is left untouched (no dt limit is
implied by the volume balance itself — a smaller `dt` simply gives a smaller
`V_avail`, which is the intended scaling).

## Optional pressure draw-down (implemented, off unless `volume_pacing_pressure=true`)

Rather than only *stopping* the loop when the volume budget is spent, the
fracture-pressure BC can be relaxed so the *next* round's K1 is evaluated at
the pressure the well would actually be at. With the well modelled as a
constant-rate source feeding a compressible fracture + leak-off system, the
first-order response of the perforation pressure to an over-drawn volume is

    p_perf(k+1)  =  p_perf(0)  −  ( V(k) − V0 − V_avail(k) )⁺ / C_eff              (3)

    C_eff  :=  Σ_cells area · ( ∂width/∂p )  ≈  V(k) / ( p_perf(0) − p_res )       (4)

i.e. `C_eff` is the fracture's storage compliance estimated secant-wise from
the current state (width is set by the elastic response to net pressure). (3)
is only applied when the over-draw is positive; `p_perf(k+1)` is floored at the
reservoir pressure at the perforation. This is the "pressure falls as the
fracture takes volume" feedback that the plain stop criterion only enforces
implicitly.

## What this is and is not

- It is a **mass-balance pacing at the timestep scale**, not a new propagation
  criterion: K1/K1c still decides *where* the front moves; (1)–(2) decide *how
  much* can open in `dt`.
- It makes growth per solve **dt-consistent** by construction: halving `dt`
  halves `V_avail`, so the guard/chop machinery is no longer needed to pace the
  onset.
- It does **not** change the converged, steady state (once the fracture leaks
  off what it receives, `V_avail → 0` and the criterion is inactive), so
  long-time results are unaffected.
- Uncertainty (user's): whether `q_perf·dt` is a *tight enough* bound. On the
  caked model2 onset, `q_perf ≈ 8000 Sm³/d ≈ 0.09 m³/s`, `dt = 0.1 d` gives
  `V_avail ≈ 800 m³` — vs the ~50 m³ a 1000 m² fracture at 5 cm... At the
  observed widths (mm), 800 m³ would allow ~10⁵ m² — far too loose. This is
  why the leak-off term matters: with `Q_leak(k)` credited round by round the
  bound tightens as the fracture opens; whether it tightens *enough* is the
  thing to measure. If it does not, the pressure draw-down (3) is the term that
  supplies the missing feedback, and `f_vol < 1` is the honest tuning knob.

## First test (caked model2 SHORT10D, 2026-08-15) — what the numbers say

- The storage term is irrelevant, as suspected: `V(k)` is 0.0006–3 m³ while
  `q_perf·dt` is 200–5000 m³. The bound (1) is carried entirely by the leak-off
  term: it fires exactly when `Q_leak(k) > q_perf` (e.g. 0.034 > 0.031 m³/s), i.e.
  when the fracture already leaks off more than its perforation delivers, so
  `V_avail = 0` and the front halts. The draw-down (3) then does what it says
  (477 → 372 bar). So the *physics* works; the *magnitudes* hinge on `q_perf`.
- With `q_perf = well_perf_rate_` (this perforation's current CTF share of the
  well) the fracture stalls at 252 m² for the whole 10 days (BHP 465 → 517):
  the share is **circular** for an establishing fracture — a small fracture gets
  a small share, the small share halts it. Reveal's open fracture takes
  essentially the whole well. Hence `volume_pacing_source = well` uses the
  whole well's reservoir rate as the budget (single-fracture wells; the
  natural first-order statement).
- Consequence for the formulation: `q_perf` in (1) should be read as *the rate
  the well would deliver into the fracture at the drawn-down pressure*, which
  for one fracture per well is the well rate; a proper multi-fracture split
  would need the well model in the loop (option (b) in the discussion).

## Options

| key (`fractureparam.solver.*`) | default | meaning |
|---|---|---|
| `volume_paced_propagation` | `false` | enable the stop criterion (2) |
| `volume_pacing_factor` | `1.0` | `f_vol` in (1) |
| `volume_pacing_pressure` | `false` | also apply the draw-down (3)–(4) between rounds |
| `volume_pacing_verbosity` | `0` | >0 logs `V(k)`, `V0`, `Q_leak`, `V_avail`, and any draw-down per round |

## Review 2026-09-04 (HMN + Claude): measured on the Reveal-comparison boxes

**Verdict: the volume budget never binds on the thermal boxes, and the physics it
encodes is not what paces Reveal.** Keep it as an opt-in guard; do not rely on it
for establishment.

Measured, dz20 caked (20 m cells, 3000 Sm3/d), candX with and without
`volume_paced_propagation=true, volume_pacing_source=well, volume_pacing_pressure=true`:
**bit-identical** — BHP, area (0.3 → 6327 m2 in the single step day 70.5→75.5),
`WWIRFRAC`, everything. The budget is `q·dt ≈ 3062 m3/d × 5 d ≈ 15 000 m3` against a
fracture volume of ~800 m3 (8400 m2 × 0.1 m): a factor ~20 too loose, and
`Q_leak(k)` never exceeds `q_perf` during the burst, so (1) is inactive. This is the
"far too loose" case the uncertainty paragraph above anticipated.

**What actually paces growth in the reference.** Reveal's per-node output on the same
deck (`..._10ppm_nodes.csv`) shows its closure stress at the seed relaxing with
cooling — 312 → 265 → 225 → 213 bar over days 0–90 as T falls 78 → 33 °C — the
crack opening where the cake-driven fracture pressure (223 → 232) crosses it, and
from then on **every perimeter node sits inside the cooled zone** (node T 32–45 °C at
day 90, 17–20 °C by day 300, far field 90 °C) with closure 183–215 bar against a
fracture pressure of ~220. The front tracks the isotherm at which σ(T) drops below p.
Growth is limited by how far the thermal front has reached, not by fluid volume. On
the isothermal-looking fine box (dxyz2, 2 m cells, injection 30 °C) the same holds
with a very small front: the seed node is already at 30 °C on day 0, `p_frac ≈
p_res` to 0.1 bar, and the fracture creeps 0.3 → 974 m2 over 180 days.

**Why OPM candX bursts instead.** Inside one step the sub-iterations run the
fracture to its K1 = K1c equilibrium against a reservoir frozen at the start of the
step; the feedback that exists (the duplicated well row's BHP falls 227.5 → 214.9
bar while the area walks 0.7 → 6327 m2) is pressure-only and local to the fracture
solve. The stress state it equilibrates against is the cell-averaged one:
`THERMEXR=3e-5, YMODULE=14 GPa` give ≈5.6 bar/°C, and OPM relieves −125 to −189 bar
at the seed cell by day 60–90 (Reveal's node: −87 bar by day 60, at a colder node),
applied to whole 20 m cells. Reveal's closure responds to the cooled *volume*
(sub-grid). So the reach of the cooled region — the thing that limits Reveal — is
grid-resolved in OPM and a full cell face opens at once. That is a discretisation
issue, not missing physics; the criterion basis is the same (the fracture samples the
total stress, `GeoMechModel::stress()` = effective + `initStress` + pore term, and
OPM's initial closure, 0.167 bar/m ≈ 369 bar, is *higher* than Reveal's 312).

**Consequence for this note.** A per-step pacing that would matter here is one tied
to the *stress state change per step* at the front, not to fluid volume: e.g. do not
propagate into a front cell whose sampled closure has not changed since the last
mech solve of this step (equivalently: one propagation round per mech solve, which
is what PostSolve does and which reproduces Reveal's BHP to ~1 bar on both boxes).
The volume balance (1)–(4) stays correct and harmless; it is simply not the active
constraint on injection-fracturing cases with leak-off.
