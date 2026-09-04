#!/usr/bin/env bash

cd cloudy_grid || exit 1

for run in run_*; do
    echo "$run"
    (
        cd "$run/skirt" || exit 1
        skirt model.ski >/dev/null
    )
done
