# Model: SKIRT grids + XSPEC fit

Neutral torus + photoionised polar shell (CLOUDY), in one SKIRT model, fitted
to the XMM spectra in three stages. The model is described in the docstring of
`config/model.py`.

```
config/model.py      model settings, shared by all targets
config/<TARGET>.py   redshift, Galactic N_H, COMBINED_GRID of one target
grid.sh              build a grid + XSPEC table
fit.py               fit a target
torus.fits           shared stage-1 table (torus_grid/ on the server)
runs/<TARGET>/       shell_cloudy/, combined_grid/, combined.fits, fit_stage*
```

## Workflow

```bash
python3 fit.py 3C445           # stage 1 (torus only); stops if combined.fits is missing
# edit config/3C445.py: centre COMBINED_GRID on runs/3C445/fit_stage1.txt
./grid.sh 3C445 -j 2 -t 16     # CLOUDY + SKIRT + table -> runs/3C445/combined.fits
python3 fit.py 3C445           # stages 1-3
```

`fit.py` output in `runs/TARGET/`:
- `fit_stage1.txt`: torus values to centre the combined grid on;
- `fit_stage2.txt`, `fit_stage3.txt`: every parameter with its 90% interval;
  `AT LIMIT` means it hit the edge of the grid;
- `fit_stage*.png`: the plots; `fit_stage3.xcm`: load it in XSPEC.

The stages:
1. torus only, 3-10 keV;
2. combined table, 0.3-10 keV, torus/Γ/norm frozen at stage 1;
3. everything free.

## grid.sh

```bash
./grid.sh TARGET [-j JOBS] [-t THREADS] [-c CLOUDY_JOBS] [--from|--to|--only STEP] [--dry-run]
./grid.sh torus  ...                    # shared torus grid -> torus.fits (only after changing it)
```

| option | meaning |
|---|---|
| `-j` | SKIRT runs at the same time (default 1) |
| `-t` | threads per SKIRT run (default cores / JOBS) |
| `-c` | CLOUDY runs at the same time (default all cores) |
| `--from/--to/--only` | `cloudy-gen`, `cloudy-run`, `skirt-gen`, `skirt-run`, `table` |
| `--dry-run` | list the SKIRT runs still to do |

- Finished runs are skipped, so after an interruption just run it again.
- After editing a grid, rerun it: runs with unchanged inputs keep their
  results. Adding a value to an axis renumbers the runs, so most are redone.

Memory: a combined SKIRT run takes up to ~60 GB, so use JOBS ≈ RAM / 60 GB and
THREADS ≈ cores / JOBS. About 2/3 of the run time is building the cells.

Single steps:

```bash
steps/cloudy_run.sh TARGET|DIR [-j JOBS] [--force]
steps/skirt_run.sh TARGET|torus|DIR [-j JOBS] [-t THREADS] [--force] [--dry-run]
```

## On a server

The server needs `skirt`, `cloudy` and Python with numpy and astropy (no HEASoft).

```bash
rsync -av --exclude 'runs/*/*.tar.gz' --exclude __pycache__ model/ server:model/
ssh server
cd model && nohup ./grid.sh 3C445 -j 2 -t 16 > 3C445.out 2>&1 &
tail -f 3C445.out
# afterwards, copy back runs/3C445/combined.fits (and torus.fits if rebuilt)
```

To run CLOUDY locally and SKIRT on the server:
1. Locally: `./grid.sh 3C445 --to skirt-gen`.
2. Copy `runs/3C445/combined_grid/` to the server.
3. On the server: `./grid.sh 3C445 --from skirt-run`.

## Adding a target

1. Reduce it (root `README.md`).
2. Copy `config/3C234.py` to `config/<TARGET>.py`; set `TITLE`, `REDSHIFT`,
   `NH_GALACTIC_1E22` (https://www.swift.ac.uk/analysis/nhtot/).
3. Run `python3 fit.py <TARGET>` (stage 1), then centre `COMBINED_GRID`, run
   `./grid.sh <TARGET>`, and fit again.

Grid rules:
- A single-valued axis is held fixed and left out of the table.
- The cone must fit inside the torus opening:
  `shell_opening_angle + SHELL_WALL_DEG/2 <= 90 - torus_opening_angle`.
  This is checked when the target is loaded.

Optional entries in a target file: `INCLINATION` (default 90), `SPECTRA_DIR`.

## Helper scripts

```bash
python3 scripts/predict_memory.py RESOURCE_DIR runs/3C445/combined_grid/run_00012 [--min-fraction 1e-4]
# memory estimate per run; RESOURCE_DIR = .../SKIRT9_Resources_AtomsMolecules/XRayIonicGas

python3 scripts/radial_test.py cloudy|skirt|compare 3C445   # SHELL_N_RADIAL 30 vs 60, see its --help
```
