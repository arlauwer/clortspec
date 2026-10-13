# SKIRT X-ray pipeline: neutral torus + photoionised polar shell

Neutral torus + photoionised hollow polar bi-cone ("conical shell", n ∝ r⁻²),
simulated together in one SKIRT model and fitted to XMM spectra in stages
(Bert's suggestion). The model (geometry, energy range, torus grid) is in
`config.py`; its docstring describes it. Everything that belongs to one AGN is
in `targets/<name>.py`.

```
model/
├── grid.sh              build a grid + XSPEC table      } the two commands
├── fit.py               fit a target with the tables    } you run
├── steps/               the individual steps grid.sh calls
├── config.py            the model, shared by all targets
├── targets/3C234.py     redshift, Galactic N_H, combined grid of one target
├── torus_grid/          stage-1 SKIRT grid   } shared: the tables are in the
├── torus.fits           stage-1 XSPEC table  } rest frame, z is a parameter
└── runs/3C234/          everything computed for one target (not in git)
    ├── shell_cloudy/    CLOUDY grid
    ├── combined_grid/   SKIRT grid
    ├── combined.fits    XSPEC table
    └── fit_stage*.png, fit_stage3.xcm
```

Spectra come from `../reduction/xmm_reduce.sh <target>` (→ `../<target>/spectra/`).

## Running it

```bash
./grid.sh torus -j JOBS -t THREADS     # once for all targets (already done) -> torus.fits
python3 fit.py 3C445                   # stage 1 (3–10 keV); stops there until combined.fits exists
# edit targets/3C445.py: centre the torus axes (and Γ) of COMBINED_GRID on the stage-1 result
./grid.sh 3C445 -j JOBS -t THREADS     # CLOUDY + SKIRT + table -> runs/3C445/combined.fits
python3 fit.py 3C445                   # stages 1–3 -> runs/3C445/fit_stage*.png, fit_stage3.xcm
```

Both have `--help`. `grid.sh` runs these steps in order; `--from`, `--to` and `--only` select a
part, and the two run steps skip what is already finished, so an interrupted
grid is continued by calling it again.

| step | script | what |
|---|---|---|
| `cloudy-gen` | `steps/cloudy_generator.py TARGET` | CLOUDY inputs for every shell configuration → `runs/TARGET/shell_cloudy/` (not for `torus`) |
| `cloudy-run` | `steps/cloudy_run.sh TARGET -j JOBS` | run CLOUDY, ≈100–400 s per run (not for `torus`) |
| `skirt-gen` | `steps/skirt_generator.py torus\|TARGET` | SKIRT models → `torus_grid/` or `runs/TARGET/combined_grid/`; validates CLOUDY (`shell_cloudy/cloudy_check.csv`) |
| `skirt-run` | `steps/skirt_run.sh torus\|TARGET -j JOBS -t THREADS` | run SKIRT (`--dry-run` lists what is left, `--force` reruns finished models) |
| `table` | `steps/build_xspec.py torus\|TARGET` | → `torus.fits` or `runs/TARGET/combined.fits` |

## The two grids

- **`torus`**: `TORUS_GRID` in `config.py`, built with `./grid.sh torus` into
  `torus_grid/` and `torus.fits`. It is shared by all targets, because the
  tables are in rest-frame energies and the redshift is a fit parameter.
  `fit.py` uses it for stage 1. Build it once (or copy `torus.fits` from where
  it was built; `.fits` files are not in git) and rebuild only after changing
  `config.py`, e.g. when a target's stage-1 fit runs into the edge of the grid.
- **a target**: `COMBINED_GRID` in `targets/<name>.py`, built with
  `./grid.sh <name>` into `runs/<name>/`.

After editing a grid, run `grid.sh` again. Runs whose parameters are unchanged
keep their results; the others are reset and redone. Inserting a value into an
axis renumbers the runs, so then most of them are redone. If the grid became
smaller, the generator lists the `run_*` directories to delete.

## Adding a target

1. Reduce it (`../README.md`, "Adding a target") so that `../<name>/spectra/`
   holds `{pn,mos1,mos2}_grp.fits`.
2. Copy `targets/3C234.py` to `targets/<name>.py`; set `TITLE`, `REDSHIFT` and
   `NH_GALACTIC_1E22` (https://www.swift.ac.uk/analysis/nhtot/).
3. `python3 fit.py <name>` for stage 1, which needs no new simulations.
4. Centre `COMBINED_GRID` on that result, then `./grid.sh <name>` and
   `python3 fit.py <name>` again.

Optional entries in a target file: `INCLINATION` (default 90) and
`SPECTRA_DIR` (default `../<name>/spectra`).

## Running on a server

`grid.sh` needs only `skirt` and `cloudy` on PATH and Python with numpy and
astropy; no HEASoft or SAS. CLOUDY's abundance file is found via
`$CLOUDY_DATA_PATH` or the `data/` directory above the `cloudy` executable.
All paths inside the grids are relative, so the directory can be copied as is:

```bash
rsync -av --exclude 'runs/*/*.tar.gz' --exclude __pycache__ model/ server:model/
# on the server
nohup ./grid.sh 3C445 -j JOBS -t THREADS > 3C445.out 2>&1 &
# -> runs/3C445/combined.fits, copy back
```

To run CLOUDY locally and only SKIRT on the server: `./grid.sh 3C445 --to
skirt-gen` here, copy `runs/3C445/combined_grid/`, and `./grid.sh 3C445 --from
skirt-run` there.

A combined run peaks at ~11–18 GB, so choose JOBS ≈ RAM / 20 GB and
THREADS ≈ cores / JOBS. `fit.py` (pyXspec) runs on the machine with the spectra;
it sets up HEASoft itself (`$HEADAS`, default `../heasoft-6.36/`).

The combined grid runs in SKIRT's `LyaExtinctionOnly` mode so that resonant
line scattering in the shell is active (`config.RESONANT_SCATTERING`).

## Table conventions

Rest-frame energies with redshift as a table parameter (no `zashift`); every
spectrum normalised to an intrinsic rest-frame 2–10 keV photon flux of
1 photon/cm²/s, so `norm` means the same in both tables; the observer
inclination (`config.INCLINATIONS`) is the last table parameter; a grid axis
with a single value is fixed and left out of the table.

`templates/` holds the two SKIRT templates. `scripts/` contains an older
plotting helper written for the first version of the pipeline.
