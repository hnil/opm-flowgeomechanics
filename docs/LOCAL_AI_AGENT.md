# Running the geomechanics example with an offline AI agent

The OPM agent from opm-agentic runs a small open model (Qwen3) on your own machine, with no
internet access, and drives simulations through its `opm_*` tools. This note sets it up for
`flow_energy_geomech` and shows what to ask it.

Before you start, build `flow_energy_geomech` as in [GETTING_STARTED.md](GETTING_STARTED.md)
(section 1). The paths below assume that layout (`~/opm-geomech`).

Tested on 14 Oct 2026 with Qwen3-8B (the smallest model the installer offers), on macOS; the
Ubuntu installer itself is described in opm-agentic.

## What it can and cannot do

| works with the small model | not reliable, or not possible |
|---|---|
| run a deck with given options and say whether it succeeded | conclusions from a time series, e.g. when a value crosses a threshold (wrong with and without `--think`) |
| report a vector at the last report step, and simple arithmetic on it | pictures or plots: the agent has no plotting tool |
| summarise how the run converged | runs on several processes (`mpirun`) |
| copy a table of vectors over time (needs the opm-agentic `summary-rows` change) | |

The prompts that were tried, with their answers: [LOCAL_AI_PROMPTS.md](LOCAL_AI_PROMPTS.md).

Check every number against the tool results it prints (the `->` lines). For series and
pictures, use the `summary` tool and ResInsight as in GETTING_STARTED.md section 3.

## 1. Install the agent

opm-agentic is a private repository: ask Halvor to add your machine's deploy key. Then follow
its guide `deploy/linux/README.md`, sections 2 to 4:

1. get the code with the deploy key (section 2);
2. `bash deploy/linux/install-box.sh` (section 3): it builds the model server and downloads
   the models, choosing what fits your hardware. `--cpu-only` works without a graphics card;
3. `sudo bash deploy/linux/setup-isolation.sh`, once (section 3), and
   `opm-agent-offline check` to see that nothing outside the machine is reachable (section 4).

The guide's `setup-opm-build.sh` step builds plain OPM Flow; you do not need it for this.

## 2. Tell it about flow_energy_geomech

Open `~/.config/opm-agent/config.toml` and add the binary under `[binaries]`, with your own
home folder:

```toml
scratch = "/home/yourname/opm-agent/scratch"   # where runs are written

[binaries]
geomech = "/home/yourname/opm-geomech/build/opm-flowgeomechanics/bin/flow_energy_geomech"
```

Keep any lines the installer already wrote.

## 3. Run the SIMPLE model

One question, then exit (the prompt is built from your own paths):

```bash
D=$HOME/opm-geomech/opm-flowgeomechanics/data
opm-agent-offline fast --prompt "Run the deck $D/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SHOWCASE.DATA with binary geomech and these extra flow flags: --fracture-param-file=$D/best_practice/seq_implicit_fb_robust_v2.json --linear-solver-mech=$D/amgmechsolverhypre_geosxcpu_new_nv.json --enable-write-all-solutions=true . When it has finished, report WBHP, WWIR and WWIRFRAC for well B-3H at the last report step, and WWIRFRAC/WWIR in percent."
```

Expected (80 s on a laptop, of which 45 s is the simulation): it calls `opm_run`, then
`opm_summary`, and answers WBHP about 221 bar, WWIR 20 000 m³/day, WWIRFRAC about
19 300 m³/day, 96 %.

Name everything in the question: the deck path, `binary geomech`, the flags. A small model
does not guess them well. The answer also gives the run folder (under `scratch`); use that
folder, the one with the `.SMSPEC` file in it, in later questions.

Then, for example (replace RUN with that folder):

```bash
opm-agent-offline fast --prompt "How did the run RUN converge? Give Newton iterations per time step, wasted iterations and time-step chops."
```

Expected: 35 time steps, 3.2 Newton iterations per step, none wasted. These are the flow
iterations only; the mechanics and fracture iterations are in the run's `.DBG` log
("Mechanical solve stats", "Fracture solve stats").

`fast` is the model the installer chose for chat: the largest that fits on your graphics card,
which on a small card is Qwen3-8B. `big` is a larger, slower one.

For a conversation instead of single questions, start `opm-agent-offline` and type at the
`you>` prompt; `/help` lists the commands.

## 4. Look at the result

The run folder holds the same files as in GETTING_STARTED.md section 2. Look at them as in
its section 3: numbers with `build/opm-common/bin/summary`, fracture area with
`scripts/frac_area.py`, pictures in ResInsight and ParaView.
