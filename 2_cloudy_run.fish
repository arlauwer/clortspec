#!/usr/bin/fish
#
# 2_cloudy_run.fish
# ------------------
# Runs `cloudy < sim.in > sim.out` for every run_*/cloudy/ dir produced by
# 1_cloudy_generator.py, keeping N instances running at once (a new one
# starts as soon as a slot frees up).
#
# Usage:
#   ./2_cloudy_run.fish [N_JOBS] [ROOT_DIR]
#
#   N_JOBS    number of concurrent Cloudy instances (default: nproc)
#   ROOT_DIR  directory containing run_* subdirs   (default: cloudy_grid)
#
# Skips a run if sim.out already exists and is non-empty (resume-friendly).
# Use --force as an extra argument to re-run everything anyway.

set -l N_JOBS $argv[1]
if test -z "$N_JOBS"
    set N_JOBS (nproc)
end

set -l ROOT $argv[2]
if test -z "$ROOT"
    set ROOT cloudy_grid
end

set -l FORCE 0
if contains -- --force $argv
    set FORCE 1
end

if not command -q cloudy
    echo "error: 'cloudy' executable not found on PATH" >&2
    exit 1
end

if not test -d $ROOT
    echo "error: root directory '$ROOT' not found" >&2
    exit 1
end

# Resolve to an absolute path so relative paths still work correctly after
# each background job cd's into its own run directory.
set ROOT (realpath $ROOT)

# Collect the list of run directories that actually need running
set -l pending
for dir in $ROOT/run_*
    set -l rundir $dir/cloudy
    if not test -f $rundir/sim.in
        continue
    end
    set pending $pending $rundir
end

set -l total (count $pending)
if test $total -eq 0
    echo "Nothing to do -- all runs already have a sim.out (use --force to redo)."
    exit 0
end

echo "Running $total Cloudy simulation(s) under $ROOT, $N_JOBS at a time..."

set -l log $ROOT/run_log.txt
rm -f $log

# --- job pool -------------------------------------------------------------
for rundir in $pending
    # wait for a free slot
    while test (jobs -p | count) -ge $N_JOBS
        sleep 0.2
    end

    fish -c "
        cd $rundir
        and cloudy < sim.in > sim.out
        set -l st \$status
        if test \$st -eq 0
            echo OK  $rundir >> $log
        else
            echo FAIL $rundir >> $log
        end
    " &
end

# --- wait for the rest, then tally results ---------------------------------
for job in (jobs -p)
    wait $job
end

set -l ok_count (grep -c '^OK' $log)
set -l fail_count (grep -c '^FAIL' $log)
echo "Done: $ok_count OK, $fail_count failed (see $log)."