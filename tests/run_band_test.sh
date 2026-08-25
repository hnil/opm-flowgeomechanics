#!/usr/bin/env bash
# run_band_test.sh <simulator> <deck> <fracture-config> <mech-solver-json> <ini> <outdir> <checker> [checker args...]
set -u
SIM=$1; DECK=$2; CFG=$3; MECH=$4; INI=$5; OUT=$6; CHECKER=$7; shift 7
rm -rf "$OUT"; mkdir -p "$OUT"
"$SIM" --parameter-file="$INI" "$DECK" \
  --fracture-param-file="$CFG" --linear-solver-mech="$MECH" \
  --solver-max-time-step-in-days=2 --solver-min-time-step=0.005 \
  --solver-continue-on-convergence-failure=true \
  --enable-write-all-solutions=true --output-dir="$OUT" > "$OUT/run.log" 2>&1
status=$?
if [ $status -ne 0 ]; then
  echo "simulator exited with $status"; tail -20 "$OUT/run.log"; exit 1
fi
python3 "$CHECKER" --log "$OUT/run.log" "$@"
