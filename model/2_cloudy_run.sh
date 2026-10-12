#!/usr/bin/env bash
#
# 2_cloudy_run.sh
# ---------------
# Run `cloudy < sim.in > sim.out` for every run_*/ of a target's shell CLOUDY
# grid written by 1_cloudy_generator.py, JOBS at a time (each CLOUDY run is
# single-threaded, ~100-400 s). The neutral torus needs no CLOUDY.
#
# Usage:
#   ./2_cloudy_run.sh TARGET [-j JOBS] [--force]
#
#   TARGET     e.g. 3C234 (runs/3C234/shell_cloudy), or a path to a grid directory
#   -j JOBS    concurrent CLOUDY runs (default: nproc)
#   --force    also re-run finished runs
#
# A run counts as finished when its sim.out ends with "Cloudy exited OK", so
# the script can be interrupted and restarted. Summary: <grid>/cloudy_run_log.txt

set -u
cd "$(dirname "$0")" || exit 1

JOBS=$(nproc)
FORCE=0
ROOT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) JOBS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) sed -n '9,15p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
        *) ROOT="$1"; shift ;;
    esac
done
[[ -z "$ROOT" ]] && { sed -n '9,15p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
[[ -d "runs/$ROOT/shell_cloudy" ]] && ROOT="runs/$ROOT/shell_cloudy"
[[ -d "$ROOT" ]] || { echo "error: '$ROOT' not found (run 1_cloudy_generator.py TARGET first)" >&2; exit 1; }
ROOT="$(realpath "$ROOT")"
command -v cloudy >/dev/null 2>&1 || { echo "error: 'cloudy' not on PATH" >&2; exit 1; }

finished() { grep -q "Cloudy exited OK" "$1/sim.out" 2>/dev/null; }

pending=()
for dir in "$ROOT"/run_*; do
    [[ -f "$dir/sim.in" ]] || continue
    if [[ "$FORCE" -eq 0 ]] && finished "$dir"; then continue; fi
    pending+=("$dir")
done
total=${#pending[@]}
echo "$ROOT: $total CLOUDY run(s) to do, $JOBS at a time"
(( total == 0 )) && exit 0

log="$ROOT/cloudy_run_log.txt"
echo "# $(date '+%F %T') start: $total runs, $JOBS jobs" >> "$log"
base=$(wc -l < "$log")
start=$SECONDS

for dir in "${pending[@]}"; do
    while (( $(jobs -pr | wc -l) >= JOBS )); do
        sleep 2
        n=$(( $(wc -l < "$log") - base ))
        printf '\r  %d/%d done, %d min elapsed   ' "$n" "$total" $(( (SECONDS - start) / 60 ))
    done
    (
        cd "$dir" || exit 1
        t0=$SECONDS
        cloudy < sim.in > sim.out 2>&1
        if finished "$dir"; then
            echo "OK   ${dir##*/} $((SECONDS - t0))s" >> "$log"
        else
            echo "FAIL ${dir##*/}" >> "$log"
        fi
    ) &
done
wait
echo

n_fail=$(tail -n +"$((base + 1))" "$log" | grep -c '^FAIL ')
echo "Done in $(( (SECONDS - start) / 60 )) min: $((total - n_fail)) OK, $n_fail failed (log: $log)"
(( n_fail == 0 ))
