# SKIRT X-ray pipeline: neutral torus + photoionised polar shell

Neutral torus + photoionised hollow polar bi-cone ("conical shell", n ∝ r⁻²),
simulated together in one SKIRT model and fitted to XMM spectra in stages
(Bert's suggestion). The model (geometry, energy range, torus grid) is in
`config.py`; its docstring describes it. Everything that belongs to one AGN is
in `targets/<name>.py`.

```
model/
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

## Steps

Once, for all targets (already done; redo only if `config.py` changes):

| step | command | what |
|---|---|---|
| 3 | `python3 3_skirt_generator.py torus` | torus-only SKIRT models → `torus_grid/` |
| 4 | `./4_skirt_run.sh torus -j JOBS -t THREADS` | run them (resumable) |
| 5 | `python3 5_build_xspec.py torus` | → `torus.fits` |

Per target (here `3C445`):

| step | command | what |
|---|---|---|
| 6 | `python3 6_xspec.py 3C445` | stage 1 fit (3–10 keV); stops if `runs/3C445/combined.fits` does not exist yet |
| — | edit `targets/3C445.py` | centre the torus axes (and Γ) of `COMBINED_GRID` on the stage-1 result |
| 1 | `python3 1_cloudy_generator.py 3C445` | CLOUDY inputs for every shell configuration → `runs/3C445/shell_cloudy/` |
| 2 | `./2_cloudy_run.sh 3C445 -j JOBS` | run CLOUDY (≈100–400 s per run, resumable) |
| 3 | `python3 3_skirt_generator.py 3C445` | torus + shell SKIRT models → `runs/3C445/combined_grid/`; validates CLOUDY (`shell_cloudy/cloudy_check.csv`) |
| 4 | `./4_skirt_run.sh 3C445 -j JOBS -t THREADS` | run them (`--dry-run` lists what is left, `--force` reruns finished models) |
| 5 | `python3 5_build_xspec.py 3C445` | → `runs/3C445/combined.fits` |
| 6 | `python3 6_xspec.py 3C445` | stages 1–3 → `runs/3C445/fit_stage*.png`, `fit_stage3.xcm` |

## Adding a target

1. Reduce it (`../README.md`, "Adding a target") so that `../<name>/spectra/`
   holds `{pn,mos1,mos2}_grp.fits`.
2. Copy `targets/3C234.py` to `targets/<name>.py`; set `TITLE`, `REDSHIFT` and
   `NH_GALACTIC_1E22` (https://www.swift.ac.uk/analysis/nhtot/).
3. `python3 6_xspec.py <name>` for stage 1, which needs no new simulations.
4. Centre `COMBINED_GRID` on that result and run steps 1–6.

Optional entries in a target file: `INCLINATION` (default 90) and
`SPECTRA_DIR` (default `../<name>/spectra`).

## Running on a server

Steps 1–5 need only `skirt` and `cloudy` on PATH and Python with numpy and
astropy; no HEASoft or SAS. CLOUDY's abundance file is found via
`$CLOUDY_DATA_PATH` or the `data/` directory above the `cloudy` executable.
All paths inside the grids are relative, so the directory can be copied as is:

```bash
rsync -av --exclude 'runs/*/*.tar.gz' --exclude __pycache__ model/ server:model/
# on the server
python3 1_cloudy_generator.py 3C445 && ./2_cloudy_run.sh 3C445 -j $(nproc)
python3 3_skirt_generator.py 3C445
nohup ./4_skirt_run.sh 3C445 -j JOBS -t THREADS > 3C445.out 2>&1 &
python3 5_build_xspec.py 3C445          # -> runs/3C445/combined.fits, copy back
```

A combined run peaks at ~11–18 GB, so choose JOBS ≈ RAM / 20 GB and
THREADS ≈ cores / JOBS. Step 6 (pyXspec) runs on the machine with the spectra.

The combined grid runs in SKIRT's `LyaExtinctionOnly` mode so that resonant
line scattering in the shell is active (`config.RESONANT_SCATTERING`).

## Table conventions

Rest-frame energies with redshift as a table parameter (no `zashift`); every
spectrum normalised to an intrinsic rest-frame 2–10 keV photon flux of
1 photon/cm²/s, so `norm` means the same in both tables; the observer
inclination (`config.INCLINATIONS`) is the last table parameter; a grid axis
with a single value is fixed and left out of the table.

`templates/` holds the two SKIRT templates. `scripts/` contains older helpers
written for the first version of the pipeline.
