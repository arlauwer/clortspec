#!/usr/bin/env bash
#
# grid.sh
# -------
# Build a model grid and its XSPEC table from start to finish.
#
# Usage:
#   ./grid.sh torus  [options]     shared stage-1 grid       -> torus.fits
#   ./grid.sh TARGET [options]     torus + shell, e.g. 3C234 -> runs/TARGET/combined.fits
#
# Options:
#   -j JOBS          concurrent SKIRT simulations (default: 1)
#   -t THREADS       threads per SKIRT simulation (default: nproc / JOBS)
#   -c JOBS          concurrent CLOUDY runs       (default: nproc)
#   --from STEP      start at STEP
#   --to STEP        stop after STEP
#   --only STEP      run just STEP
#   --dry-run        skirt-run only lists what is left to simulate
#   -h, --help       this text
#
# Steps, in order (scripts in steps/):
#   cloudy-gen   cloudy_generator.py   CLOUDY inputs            (not for torus)
#   cloudy-run   cloudy_run.sh         run CLOUDY               (not for torus)
#   skirt-gen    skirt_generator.py    SKIRT models; validates the CLOUDY runs
#   skirt-run    skirt_run.sh          run SKIRT
#   table        build_xspec.py        XSPEC table
#
# The two run steps skip what is already finished, so an interrupted grid is
# continued by calling the script again.
#
# Which grid is what:
#   torus    TORUS_GRID in config.py. One table for ALL targets (rest-frame
#            energies, redshift is a fit parameter); fit.py uses it for stage 1.
#            Build it once; rebuild only after changing config.py.
#   TARGET   COMBINED_GRID in config/TARGET.py, centred on that target's
#            stage-1 result.
#
# After editing a grid, just run this script again: runs whose parameters are
# unchanged keep their results, the others are reset and redone. (Inserting a
# value into an axis renumbers the runs, so most of them are redone.) If the
# grid became smaller, the generator tells you which run_* directories to delete.
#
# Examples:
#   ./grid.sh 3C234 -j 4 -t 5                 everything, 4 SKIRT runs x 5 threads
#   ./grid.sh 3C234 --only skirt-run --dry-run     what is left?
#   ./grid.sh 3C234 --to skirt-gen            CLOUDY here ...
#   nohup ./grid.sh 3C234 --from skirt-run -j 3 -t 8 > 3C234.out 2>&1 &    ... SKIRT on a server
#
# Memory: with 30 equal-column shell cells a combined run peaks at ~17 GB at
# logNHs 24 (was up to ~60 GB with 60 log-spaced cells). Before skirt-run,
#   python3 scripts/predict_memory.py RESOURCE_DIR runs/TARGET/combined_grid/run_* --ram GB
# prints the largest expected peak and how many JOBS fit; THREADS ~ cores / JOBS.
#
# The whole procedure for a target:
#   python3 fit.py TARGET        stage 1, needs only torus.fits
#   edit config/TARGET.py        centre COMBINED_GRID on the stage-1 result
#   ./grid.sh TARGET             this script
#   python3 fit.py TARGET        stages 1-3

set -euo pipefail
cd "$(dirname "$0")"

usage() {
    sed -n '7,56p' "$0" | sed 's/^# \{0,1\}//'
    echo; echo "Targets: $(ls config/*.py | xargs -n1 basename -s .py | grep -vx -e __init__ -e model | xargs)"
    exit 1
}

GRID=""
SKIRT_ARGS=()
CLOUDY_JOBS=$(nproc)
FROM=""
TO=""
ONLY=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|-t) SKIRT_ARGS+=("$1" "$2"); shift 2 ;;
        -c) CLOUDY_JOBS="$2"; shift 2 ;;
        --from) FROM="$2"; shift 2 ;;
        --to) TO="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        --dry-run) SKIRT_ARGS+=(--dry-run); shift ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1" >&2; usage ;;
        *) GRID="$1"; shift ;;
    esac
done
[[ -z "$GRID" ]] && usage

if [[ "$GRID" == torus ]]; then
    STEPS=(skirt-gen skirt-run table)
else
    [[ -f "config/$GRID.py" && "$GRID" != model && "$GRID" != __init__ ]] || { echo "error: unknown target '$GRID' (no config/$GRID.py)" >&2; exit 1; }
    STEPS=(cloudy-gen cloudy-run skirt-gen skirt-run table)
fi
for s in "$FROM" "$TO" "$ONLY"; do
    [[ -z "$s" || " ${STEPS[*]} " == *" $s "* ]] || { echo "error: no step '$s' for $GRID (steps: ${STEPS[*]})" >&2; exit 1; }
done

run_step() {
    case "$1" in
        cloudy-gen) python3 steps/cloudy_generator.py "$GRID" ;;
        cloudy-run) steps/cloudy_run.sh "$GRID" -j "$CLOUDY_JOBS" ;;
        skirt-gen)  python3 steps/skirt_generator.py "$GRID" ;;
        skirt-run)  steps/skirt_run.sh "$GRID" "${SKIRT_ARGS[@]}" ;;
        table)      python3 steps/build_xspec.py "$GRID" ;;
    esac
}

started=0
[[ -z "$FROM" ]] && started=1
for step in "${STEPS[@]}"; do
    [[ "$step" == "$FROM" ]] && started=1
    [[ "$started" -eq 1 ]] || continue
    if [[ -z "$ONLY" || "$step" == "$ONLY" ]]; then
        echo "[$(date '+%T')] $GRID: $step"
        run_step "$step" || { echo "error: step '$step' failed -- fix it and continue with --from $step" >&2; exit 1; }
    fi
    [[ "$step" == "$TO" ]] && break
done
echo "[$(date '+%T')] $GRID: done"
