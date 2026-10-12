#!/usr/bin/env bash
#
# 4_skirt_run.sh
# --------------
# Run SKIRT for every run_*/model.ski of a grid, JOBS simulations at a time.
#
# Usage:
#   ./4_skirt_run.sh GRID [-j JOBS] [-t THREADS] [--force] [--dry-run]
#
#   GRID        torus (shared stage-1 grid) | TARGET, e.g. 3C234 (its combined
#               grid, runs/3C234/combined_grid) | <path to a grid directory>
#   -j JOBS     concurrent SKIRT simulations          (default: 1)
#   -t THREADS  threads per simulation                 (default: nproc / JOBS)
#   --force     also re-run simulations that already finished
#   --dry-run   only list what would be run
#
# A run counts as finished when it has an SED file for every instrument in its
# model.ski, so the script can be interrupted and restarted at any time.
# SKIRT output goes to run_*/skirt_stdout.txt, a summary to
# <grid>/skirt_run_log.txt. On a server, e.g.:
#
#   nohup ./4_skirt_run.sh 3C234 -j 8 -t 8 > 3C234.out 2>&1 &

set -u
cd "$(dirname "$0")" || exit 1

usage() { sed -n '7,15p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

GRID=""
JOBS=1
THREADS=""
FORCE=0
DRY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) JOBS="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) usage ;;
        -*) echo "unknown option $1" >&2; usage ;;
        *) GRID="$1"; shift ;;
    esac
done
[[ -z "$GRID" ]] && usage

if [[ "$GRID" == torus ]]; then
    ROOT=torus_grid
elif [[ -d "runs/$GRID/combined_grid" ]]; then
    ROOT="runs/$GRID/combined_grid"
else
    ROOT="$GRID"
fi
[[ -d "$ROOT" ]] || { echo "error: '$ROOT' not found (run 3_skirt_generator.py first)" >&2; exit 1; }
ROOT="$(realpath "$ROOT")"

[[ -z "$THREADS" ]] && THREADS=$(( $(nproc) / JOBS ))
(( THREADS < 1 )) && THREADS=1
command -v skirt >/dev/null 2>&1 || { echo "error: 'skirt' not on PATH" >&2; exit 1; }

finished() {
    local name
    for name in $(grep -o 'instrumentName="[^"]*"' "$1/model.ski" | cut -d'"' -f2); do
        [[ -s "$1/model_${name}_sed.dat" ]] || return 1
    done
}

pending=()
n_all=0
for dir in "$ROOT"/run_*; do
    [[ -f "$dir/model.ski" ]] || continue
    n_all=$((n_all + 1))
    if [[ "$FORCE" -eq 0 ]] && finished "$dir"; then continue; fi
    pending+=("$dir")
done
total=${#pending[@]}

echo "$ROOT: $total of $n_all runs to do, $JOBS job(s) x $THREADS thread(s)"
if [[ "$DRY" -eq 1 || "$total" -eq 0 ]]; then
    (( total > 0 )) && printf '  %s\n' "${pending[@]##*/}"
    exit 0
fi

log="$ROOT/skirt_run_log.txt"
echo "# $(date '+%F %T') start: $total runs, $JOBS jobs x $THREADS threads" >> "$log"
base=$(wc -l < "$log")
start=$SECONDS

progress() {
    local n elapsed
    n=$(( $(wc -l < "$log") - base ))
    elapsed=$(( SECONDS - start ))
    (( n > 0 )) && printf '\r  %d/%d done, %d min elapsed, ~%d min left   ' \
        "$n" "$total" $((elapsed / 60)) $((elapsed * (total - n) / n / 60))
}

wait_for_slot() {
    while (( $(jobs -pr | wc -l) >= $1 )); do
        sleep 2
        progress
    done
}

for dir in "${pending[@]}"; do
    wait_for_slot "$JOBS"
    (
        cd "$dir" || exit 1
        t0=$SECONDS
        if skirt -t "$THREADS" model.ski > skirt_stdout.txt 2>&1 && finished "$dir"; then
            echo "OK   ${dir##*/} $((SECONDS - t0))s" >> "$log"
        else
            echo "FAIL ${dir##*/} (see ${dir##*/}/skirt_stdout.txt)" >> "$log"
        fi
    ) &
done
wait_for_slot 1
wait
echo

n_fail=$(tail -n +"$((base + 1))" "$log" | grep -c '^FAIL ')
echo "Done in $(( (SECONDS - start) / 60 )) min: $((total - n_fail)) OK, $n_fail failed (log: $log)"
(( n_fail == 0 ))
