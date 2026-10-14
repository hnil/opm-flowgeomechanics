# Prompts for the offline agent

Tested on 14 Oct 2026 with Qwen3-8B and the SIMPLE showcase run (setup:
[LOCAL_AI_AGENT.md](LOCAL_AI_AGENT.md)). Give each one with `opm-agent-offline fast --prompt "..."`.
`$D` is `~/opm-geomech/opm-flowgeomechanics/data`; RUN is the run folder the first answer names
(the one holding the `.SMSPEC` file).

## Works

**Run and report** (80 s): 221.0 bar, 20 000 and 19 276 m³/day, 96.4 %.

```text
Run the deck $D/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.DATA with binary geomech and these extra flow flags: --fracture-param-file=$D/best_practice/seq_implicit_fb_robust_v2.json --linear-solver-mech=$D/amgmechsolverhypre_geosxcpu_new_nv.json --enable-write-all-solutions=true . When it has finished, report WBHP, WWIR and WWIRFRAC for well B-3H at the last report step, and WWIRFRAC/WWIR in percent.
```

**Report from a finished run:** the same values.

```text
From the run RUN, report WBHP, WWIR and WWIRFRAC for well B-3H at the last report step, and WWIRFRAC/WWIR in percent.
```

**Convergence:** 35 time steps, 3.23 Newton per step, none wasted (flow only, not the
mechanics or fracture iterations).

```text
How did the run RUN converge? Give Newton iterations per time step, wasted iterations and time-step chops.
```

**A table:** the values are copied exactly, but it may add a row outside the range. Needs
opm-agentic from 14 Oct 2026 or later (`opm_summary` rows, `9b085ff7a`); before that it
shifted a column by one row.

```text
From the run RUN, give a table of time (days), WBHP:B-3H and WWIRFRAC:B-3H at every report step from day 10 to day 25. Copy the numbers exactly.
```

## Does not work

**Reasoning over a series:** the fracture takes over half the injection at day 17, with BHP
going from 244.8 bar (day 16) to 237.5 bar (day 17). The model answered day 12 (day 17 with
`--think`), each time with the wrong BHP values. Ask for the table instead and read it.

```text
In the run RUN, at which time does WWIRFRAC for well B-3H first exceed half of WWIR, and what is WBHP just before and just after that time?
```

Rules of thumb: name the deck, `binary geomech` and every flag; ask for values, not
conclusions; check the numbers against the tool results it prints.
