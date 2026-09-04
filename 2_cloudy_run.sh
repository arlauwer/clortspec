#!/usr/bin/env bash
#
# 2_cloudy_run.sh
# ----------------
# Runs `cloudy < sim.in > sim.out` for every run_*/cloudy/ dir produced by
# 1_cloudy_generator.py, keeping N instances running at once.
#
# Usage:
#   ./2_cloudy_run.sh [N_JOBS] [ROOT_DIR] [--force]
#
#   N_JOBS    number of concurrent Cloudy instances (default: nproc)
#   ROOT_DIR  directory containing run_* subdirs (default: cloudy_grid)
#
# Skips a run if sim.out already exists and is non-empty (resume-friendly).
# Use --force to re-run everything anyway.

set -u

N_JOBS="${1:-}"
ROOT="${2:-cloudy_grid}"
FORCE=0

# Handle --force wherever it appears
for arg in "$@"; do
    if [[ "$arg" == "--force" ]]; then
        FORCE=1
    fi
done

if [[ -z "$N_JOBS" || "$N_JOBS" == "--force" ]]; then
    N_JOBS=$(nproc)
fi

if ! command -v cloudy >/dev/null 2>&1; then
    echo "error: 'cloudy' executable not found on PATH" >&2
    exit 1
fi

if [[ ! -d "$ROOT" ]]; then
    echo "error: root directory '$ROOT' not found" >&2
    exit 1
fi

# Resolve to an absolute path so relative paths still work correctly after
# each background job cd's into its own run directory.
ROOT="$(realpath "$ROOT")"

# Collect the list of run directories that actually need running
pending=()

for dir in "$ROOT"/run_*; do
    [[ -d "$dir" ]] || continue

    rundir="$dir/cloudy"

    [[ -f "$rundir/sim.in" ]] || continue

    # Resume-friendly: don't rerun successful/existing simulations unless forced.
    if [[ "$FORCE" -eq 0 && -s "$rundir/sim.out" ]]; then
        continue
    fi

    pending+=("$rundir")
done

total=${#pending[@]}

if [[ "$total" -eq 0 ]]; then
    if [[ "$FORCE" -eq 1 ]]; then
        echo "Nothing to do -- no sim.in files found."
    else
        echo "Nothing to do -- all runs already have a non-empty sim.out (use --force to redo)."
    fi
    exit 0
fi

echo "Running $total Cloudy simulation(s) under $ROOT, $N_JOBS at a time..."

log="$ROOT/run_log.txt"
rm -f "$log"

# --- job pool -------------------------------------------------------------
for rundir in "${pending[@]}"; do

    # Wait for a free slot.
    while (( $(jobs -pr | wc -l) >= N_JOBS )); do
        sleep 0.2
    done

    (
        cd "$rundir" || exit 1

        cloudy < sim.in > sim.out
        st=$?

        if [[ "$st" -eq 0 ]]; then
            echo "OK    $rundir" >> "$log"
        else
            echo "FAIL  $rundir" >> "$log"
        fi

        exit "$st"
    ) &
done

# --- wait for the rest ----------------------------------------------------
wait

# --- tally results --------------------------------------------------------
ok_count=$(grep -c '^OK' "$log" 2>/dev/null || true)
fail_count=$(grep -c '^FAIL' "$log" 2>/dev/null || true)

echo "Done: $ok_count OK, $fail_count failed (see $log)."
