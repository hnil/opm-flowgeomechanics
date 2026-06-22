# Best-practice fracture parameter files

Curated `--fracture-param-file` configurations for `flow_energy_geomech`, distilled from
the 2026-06 validation sweeps (see `../../STATUS_2026-06-12.md` and
`data/mac_tests_result/*/FINDINGS.md`). Every file here shares the **validated robust
inner block** (decoupled Picard, `closed_cell_policy=sticky`, mech-last preconditioner,
`reopen_width_tolerance=0`, full anti-chatter guard) and differs only in the coupling
*method*, the well *control*, and a few case-specific levers.

Full key-by-key documentation of every option is in `../fracture_settings_REFERENCE.md`.

> **Prime rule:** none of these trades fracturing away for stability. Always check
> conduction (`WWIRFRAC/WWIR %`) and fracture area, not just completion — use the harness
> `data/mac_tests_result/analysis/compare_runs.py --bands`.

---

## The files

| file | method | control | use it for |
|---|---|---|---|
| `seq_implicit_bhp.json` | SeqMechFrac (implicit_flow) | perf_pressure (BHP) | **general default**: tight flow↔mech↔fracture coupling, highest conduction |
| `seq_implicit_rate.json` | SeqMechFrac | rate_well | rate-controlled wells (the hard case) — adds Anderson CTF acceleration |
| `postsolve.json` | PostSolve | perf_pressure | cheaper per step, fracture solved once after the converged flow+mech step; coupling is *lagged* |
| `seq_implicit_chatter.json` | SeqMechFrac | perf_pressure | compressed / vemstab decks with contact chatter (large closed-cell populations) |
| `seq_implicit_performance.json` | SeqMechFrac | perf_pressure | large / grid-refined runs where the VEM mech solve dominates |

### `seq_implicit_bhp.json` — general sequential-implicit default
The robust base. `method=SeqMechFrac`, `implicit_flow=true`: flow, mechanics and fracture
are iterated each timestep until the coupling converges. Decoupled Picard
(`drop_fluid_mech_linearization=true`), `closed_cell_policy=sticky`, mech-last
preconditioner, `use_ad_pressure_assembly=true`, gravity on, full toggle guard
(no finite re-opening cap — it would suppress fracture growth, see STATUS §3.5).
**Highest fracture conduction on every validated deck (~85–99 %).**

### `seq_implicit_rate.json` — rate-controlled
Same base, `control.type=rate_well`. Rate control is the documented stressor (the well must
push the per-step growth), so this turns on `wi_vector_acceleration=anderson` — vector
Anderson on the per-cell CTF fixed point, validated to cut inner fracture NL iterations
~65 % and remove non-convergence on the hard coupling case, at equal time. For very stiff
rate growth also pass `--solver-max-time-step-in-days=2..5`.

### `postsolve.json` — post-solve (lagged) coupling
`method=PostSolve`, `implicit_flow=false`, `max_mech_it=2`: the fracture is solved **once
after** the converged flow+mech step (the original scheme). Cheaper per step, no outer
coupling iteration; physically correct when fracture growth per step is small. Coupling is
lagged, so conduction develops more slowly on short timescales and catches up as the
fracture matures.

### `seq_implicit_chatter.json` — chatter-prone / compressed cases
General base + `toggle_guard_mode=chatter`, `chatter_flip_limit=3`. The per-cell
anti-chatter guard pins only cells that flip contact state more than `chatter_flip_limit`
times within a solve; first-time openings/closings always pass, so it breaks the
nonsmooth-contact chatter loop (e.g. vegard_0806: vemstab + large compressed zone) **without
blocking fracture establishment**. Use when you see runaway fracture-NL-iteration counts /
huge `closed_cell_toggles`.
> Requires a build that includes the chatter guard (`toggle_guard_mode` option). On older
> binaries the key is simply ignored (legacy `count` guard, unlimited cap) — still safe.

### `seq_implicit_performance.json` — large / refined runs
General base + top-level `solver.mech_skip_coupling_threshold=1.0` (conservative: skip the
mechanics solve only when the fracture→well coupling change is near zero) + Anderson CTF
acceleration. On grid-refined decks the VEM mechanics solve dominates runtime; the
conservative skip cut mech solves up to ~60 % with the solution preserved to ~1 %.
**Do not raise the threshold on rate-controlled decks** (an aggressive value diverges under
sharp rate-driven pressure drops).

---

## Correspondence to the built-in aliases

`--fracture-param-file` accepts two **built-in alias names** (no file needed); they are
compiled into the binary in `opm/geomech/BuiltinFractureParams.cpp`:

| alias | method | equivalent file here |
|---|---|---|
| `sequential_implicit` | SeqMechFrac (implicit_flow) | ≈ `seq_implicit_bhp.json` |
| `standard` | PostSolve | ≈ `postsolve.json` (the alias uses the legacy `closed_cell_policy=legacy` + `vem_stability_choice=3`; the file here uses the newer sticky/robust block) |

So for the two main schemes you can skip the file entirely and pass the alias. The files in
this directory are the place to start when you need to change control (rate), turn on a
case-specific lever (chatter, performance), or pin settings for reproducibility.

The committed regression decks use the frozen reference configs in the parent directory
(`fracture_simple_bhp_seq_final.json`, `fracture_simple_bhp_seq_legacy_full.json`,
`fracture_simple_rate_seq_full.json`, alias `sequential_implicit`) — do **not** edit those;
edit copies here instead.

---

## Example commands

Paths below assume the repo root `/Users/hnil/Documents/OPM/opm_geomech`. The mech linear
solver JSON and `myparamsfrac*.ini` live in `opm-flowgeomechanics/data/`.

```bash
BIN=builds/release/opm-flowgeomechanics/bin/flow_energy_geomech
DATA=opm-flowgeomechanics/data
BP=$DATA/best_practice

# --- General sequential-implicit (BHP), explicit file ---
$BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA \
     --fracture-param-file=$BP/seq_implicit_bhp.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --enable-write-all-solutions=true --output-dir=/tmp/run_seq_bhp

# --- Same thing via the built-in alias (no file) ---
$BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA \
     --fracture-param-file=sequential_implicit \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --output-dir=/tmp/run_seq_alias

# --- Rate-controlled well (the hard case) + small max step ---
$BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA \
     --fracture-param-file=$BP/seq_implicit_rate.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --solver-max-time-step-in-days=2 --output-dir=/tmp/run_seq_rate

# --- Post-solve (cheaper, lagged) ---
$BIN --parameter-file=$DATA/myparamsfrac.ini \
     $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC.DATA \
     --fracture-param-file=$BP/postsolve.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --output-dir=/tmp/run_post
#   or:  --fracture-param-file=standard

# --- Chatter-prone / compressed deck ---
$BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     <your_compressed_deck>.DATA \
     --fracture-param-file=$BP/seq_implicit_chatter.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --output-dir=/tmp/run_chatter

# --- Large / grid-refined run (mech-skip + Anderson) ---
$BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     data/test_oss/refine/CASE_REFINE_nx_31_ny_31_nz_25.DATA \
     --fracture-param-file=$BP/seq_implicit_performance.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --output-dir=/tmp/run_perf

# --- MPI: same options, just launch under mpirun ---
mpirun -np 2 $BIN --parameter-file=$DATA/myparamsfrac_seq.ini \
     $DATA/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA \
     --fracture-param-file=$BP/seq_implicit_bhp.json \
     --linear-solver-mech=$DATA/amgmechsolverhypre_geosxcpu_new_nv.json \
     --threads-per-process=1 --output-dir=/tmp/run_mpi
```

### Always verify conduction afterwards
```bash
# per-run table incl. WWIRFRAC% and fracture area/volume
python3 data/mac_tests_result/analysis/compare_runs.py <output-dir-parent> \
        --bands data/mac_tests_result/analysis/acceptance_bands.json
# or, for a single fracture VTU:
python3 data/test_oss/refine/frac_area.py /tmp/run_seq_bhp/*nr-*.vtu
```

## Verified validation cases (fast, SIMPLE-based, strong well–fracture coupling)
- **Pressure-driven:** `data/test_oss/refine/CASE_PRESS_FRAC_nx_21_ny_21_nz_25.DATA`
  (isothermal, BHP 420 > closure; ~56 s; ~98 % conduction).
- **Temperature-driven:** `data/test_oss/refine/CASE_REFINE_nx_21_ny_21_nz_25.DATA`
  (cold 15 °C injection into 90 °C reservoir, BHP below closure → thermal-stress opening).

---

## How these were verified (reproducible)

Both decks live in `data/test_oss/refine/` and use that directory's
`myparamsfrac_seq.ini` (adaptive time-stepping on). The fracture config used for
verification is the general default, `seq_implicit_bhp.json` (identical to the parent
`fracture_simple_bhp_seq_final.json`). Run **from `data/test_oss/refine/`**:

```bash
ROOT=/Users/hnil/Documents/OPM/opm_geomech
BIN=$ROOT/builds/release/opm-flowgeomechanics/bin/flow_energy_geomech
BP=$ROOT/opm-flowgeomechanics/data/best_practice
MECH=$ROOT/opm-flowgeomechanics/data/amgmechsolverhypre_geosxcpu_new_nv.json
cd $ROOT/data/test_oss/refine

# (A) PRESSURE-driven case
$BIN --parameter-file=myparamsfrac_seq.ini CASE_PRESS_FRAC_nx_21_ny_21_nz_25.DATA \
     --fracture-param-file=$BP/seq_implicit_bhp.json \
     --linear-solver-mech=$MECH \
     --enable-write-all-solutions=true --threads-per-process=2 \
     --output-dir=/tmp/press_frac

# (B) TEMPERATURE-driven case (same command, other deck)
$BIN --parameter-file=myparamsfrac_seq.ini CASE_REFINE_nx_21_ny_21_nz_25.DATA \
     --fracture-param-file=$BP/seq_implicit_bhp.json \
     --linear-solver-mech=$MECH \
     --enable-write-all-solutions=true --threads-per-process=2 \
     --output-dir=/tmp/temp_frac
```

`--enable-write-all-solutions=true` is **required** for the geometry check — it writes the
per-step fracture `*.vtu` that `frac_area.py` reads.

Then three checks per run (the same ones the harness automates):

```bash
SUM=$ROOT/builds/release/opm-common/bin/summary

# 1) fracture growth: total_area (grid extent) must NOT stay at the ~299 m^2 seed
python3 frac_area.py /tmp/press_frac/*nr-*.vtu | tail -1
#   -> total_area=567559 m^2  open_area=0.0 m^2  volume=0.0
#   NB: on the PRESSURE case the final snapshot reports open_area~0 even though
#   97.8% of injection routes through the fracture -- the "all cells closed but
#   conducting" contact state (open/close vs conduction inconsistency). So judge
#   GROWTH by total_area, and CONDUCTION by WWIRFRAC (check 2), not by open_area.
#   On the TEMPERATURE case open_area is large (~2.68e4 m^2).

# 2) conduction  frac% = WWIRFRAC/WWIR  (must stay high on a fracture-dominated case)
P=$(ls /tmp/press_frac/*.SMSPEC); P=${P%.SMSPEC}
for v in WBHP:B-3H WWIR:B-3H WWIRFRAC:B-3H; do
  echo "$v = $($SUM $P $v 2>/dev/null | awk 'NF' | tail -1)"
done
#   WBHP=420 (BHP-limited), WWIR=11483, WWIRFRAC=11233  -> frac% = 97.8 %

# 3) no non-convergence
grep -c "did not converge" /tmp/press_frac/../press_frac.log   # or the tee'd run log; expect 0
```

What "verified" means for these cases, and the **option each check exercises**:

| check | tool | option it pins down |
|---|---|---|
| completes, 0 timestep cuts / `did not converge` | run log `grep` | the inner solver block (Picard, sticky, mech-last) is stable |
| fracture **grows** (`total_area` ≫ seed ~299 m²) | `frac_area.py` last VTU | the toggle guard does **not** suppress opening (the whole reason `seq_implicit_bhp.json` keeps the *unlimited* guard, not 32/0.1) |
| **conduction** `WWIRFRAC/WWIR` high | `summary` CLI | `addconnections=true` + control type actually route flow through the fracture |
| pass/fail vs known-good | `compare_runs.py --bands` | the whole config, scored against `acceptance_bands.json` |

**Measured baselines** with `seq_implicit_bhp.json` (= what the bands encode):
- pressure case: area ≈ 5.68e5 m², frac% ≈ 97.8 %, 0 non-convergence, ~56 s.
- temperature case (CASE_REFINE L21): area ≈ 2.93e4 m², open ≈ 2.68e4 m², frac% ≈ 99.8 %, 0 non-convergence.

**Option sensitivity worth knowing (pressure case):** the drive is the injector BHP limit
relative to the closure stress (~334 bar + poroelastic back-stress). It is a sharp
threshold on this uniform-stress deck — `BHP 380` produces *no* opening (area stays at the
299 m² seed), `BHP 420` opens the whole domain. So when adapting this case, the BHP limit
(in `WCONINJE`) is the knob that turns fracturing on/off, not the fracture JSON. The
temperature case instead keeps `WCONINJE` BHP *below* closure and drives opening purely by
the injected/reservoir temperature contrast (`WTEMP` 15 °C vs 90 °C) — so there the
`WTEMP` value is the drive.

Single-config pass/fail against the recorded baselines:

```bash
# point compare_runs.py at a directory of <case>/output/*.SMSPEC runs, e.g. a sweep dir,
# or use frac_area.py + summary directly as above for one-off runs.
python3 $ROOT/data/mac_tests_result/analysis/compare_runs.py <run-dir> \
        --bands $ROOT/data/mac_tests_result/analysis/acceptance_bands.json
# FAILS if any case aborts, frac% drops below its floor, or (refine/pressure cases)
# fracture area collapses -- i.e. it catches a config that gained stability by not fracturing.
```
