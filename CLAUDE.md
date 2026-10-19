# CLAUDE.md

PhD project: XMM-Newton spectra of obscured radio galaxies (3C 234, 3C 445)
fitted with SKIRT X-ray radiative-transfer models (neutral torus + photoionised
polar shell, ionisation from CLOUDY). The application paper is about SKIRT's new
`XRayIonicGasMix`; background, framing and paper status are in
`docs/paper_context.md` (its file paths predate the move to `model/config/`).

Read `README.md` (reduction) and `model/README.md` (modelling pipeline) first;
the docstrings of `model/config/model.py`, `model/fit.py` and `model/steps/*.py`
describe the model and each step in detail.

## Working with the user

- **The user runs the pipeline** (reductions, `grid.sh`, server runs). Do not
  launch reductions or grids in parallel with them. Quick XSPEC fits and single
  test runs are fine, but write their output to the scratchpad / `$TMPDIR`, never
  into `model/runs/` (redirect with `config.model.RUNS_DIR = ...`, see below).
- Grids run on a server (CLOUDY + SKIRT); only `combined.fits`, plots and fit
  results come back. `model/runs/*/combined_grid/` and `torus_grid/` are usually
  NOT present locally.
- Combined SKIRT runs are expensive: up to ~60 GB RAM, so ~2 jobs at a time.
  About 2/3 of the runtime is building the shell's ionic gas mixes (cells), not
  photon packets. Cost per run scales with (#radial shell cells) x (ions/cell);
  `scripts/predict_memory.py` estimates it from a run's `shell.txt`.
- Ask before changing the physics model; explain trade-offs with numbers. The
  user prefers fewer extra spectral components (no phenomenological
  `zpowerlw`/`apec` add-ons for now) and wants to understand choices.

## Layout (model/)

```
config/model.py      the model, shared (geometry, energies, torus grid, helpers)
config/<TARGET>.py   one target: TITLE, REDSHIFT, NH_GALACTIC_1E22, COMBINED_GRID
config/__init__.py   `from .model import *`  -> scripts do `import config as C`
grid.sh              cloudy-gen -> cloudy-run -> skirt-gen -> skirt-run -> table
fit.py               staged XSPEC fit -> runs/TARGET/fit_stage{1,1_full,2,3}.{png,txt}, fit_stage3.xcm
steps/               cloudy_generator.py, cloudy_run.sh, skirt_generator.py, skirt_run.sh, build_xspec.py
scripts/             predict_memory.py, radial_test.py (N_radial 30 vs 60), plot_seds.py (old)
templates/           torus.ski, combined.ski
torus.fits           shared stage-1 table (rest frame, z is a table parameter)
runs/TARGET/         shell_cloudy/, combined_grid/, combined.fits, fit outputs (git-ignored)
```

Target names exclude `__init__` and `model` (`C.target_names()`). `fit.py`
globs `config/*.py` itself because `config` (numpy) must not be imported before
`xspec` (see below).

## Model in one paragraph

Point source (cut-off power law, Γ) + uniform neutral torus wedge (±α from the
equator, R 1e18–1e19 cm, parameterised by equatorial N_H) + hollow bi-cone shell
(wall `SHELL_WALL_DEG` = 20° thick centred at θ_s from the pole, n ∝ r^-2,
R 1e19–1e21 cm, i.e. starting at the torus outer edge). One 1D CLOUDY run per
(logxi, logNHs, θ_s, Γ); ξ = L(1–1000 Ryd)/(n_in R_in²) at the inner face; CLOUDY
zones are rebinned column-conserving onto `SHELL_N_RADIAL` (60) cells, trace
ions < `SHELL_ION_MIN_FRACTION` (1e-6) dropped. Both media in one SKIRT run
(`LyaExtinctionOnly` for resonant line scattering). Tables: every spectrum
normalised to intrinsic rest-frame 2–10 keV photon flux = 1, so `norm` carries
over between tables; single-valued grid axes are fixed and left out of the
table (recorded as `fixed:` header comments).

Constraints enforced in code:
- cone inside torus opening: θ_s + wall/2 ≤ 90 − α (checked in `load_target`;
  this is also where the source ionises the gas directly);
- combined grids record only inclinations strictly inside the narrowest torus
  (`C.combined_inclinations`: 90 − α_min + 5, then every 10° to 90), so XSPEC
  never interpolates between obscured and unobscured views. The torus grid
  keeps `TORUS_INCLINATIONS` = 0..90 step 15 (unchanged, no rerun needed).

## Fit (fit.py)

`constant*TBabs*atable{...}`, C-stat, pn/MOS1/MOS2 with free cross-norms.
Stage 1: torus.fits, 3–10 keV, Incl frozen (90). Stage 2: combined.fits,
0.3–10 keV, torus/Γ/norm frozen at stage 1. Stage 3: all free, Incl limited to
table inclinations inside the narrowest torus. `fit_stage{2,3}.txt` list every
parameter with 90% interval, range, and "AT LIMIT" flags.

## Environment quirks

- HEASoft: `HEADAS=../heasoft-6.36/x86_64-pc-linux-gnu-libc2.43`, then
  `source $HEADAS/headas-init.sh`. In the sandbox HOME is not writable: set
  `HOME=<scratchpad>/home` first. `fit.py` re-execs itself with HEASoft set up.
- **Import `xspec` before numpy/matplotlib/astropy/config**, else segfault.
- Shell is fish; use bash syntax inside the Bash tool as usual. User's own
  setup: `source reduction/setsas.fish [TARGET]`.
- Testing fit.py without touching runs/: a wrapper that does `import xspec`,
  `import config; config.model.RUNS_DIR = Path(<scratch>)`, symlinks
  `combined.fits` into `<scratch>/TARGET/`, then `runpy.run_path("fit.py")`.
  (`Target.dir` reads `config.model.RUNS_DIR`; scripts using `C.RUNS_DIR`
  directly need `config.RUNS_DIR` patched as well.)
- SKIRT resources (for predict_memory.py):
  `~/Documents/PhD/code/SKIRT/resources/SKIRT9_Resources_AtomsMolecules/XRayIonicGas`.
- `cloudy` and `skirt` are on PATH locally. A CLOUDY run at logNHs 24 takes ~10 min.
- Character-device dotfiles (`.bashrc`, `.gitconfig`, ...) in the repo root are
  sandbox mount placeholders, not real files: ignore them.

## State as of 2026-10-19

First full combined grids (old geometry: shell 1e18–1e20 cm inside the torus
radius, 10° wall, logxi −1..3 step 0.5, logNHs 21–23, α ∈ {30,45}) were fitted:
- Hard band fine; soft band underpredicted ~10–20x. Stage 3 drove every shell
  parameter to its limit (logNHs 23, θ_s 40, α 30, Γ top of range, Incl 60 =
  torus edge). Adding a scattered power law showed the data need ~9% of the
  intrinsic 2–10 keV flux scattered.
- Causes found: (1) ~90% of the shell column sat inside the torus radius and was
  hidden by the torus for edge-on views (soft flux 3–10x lower at Incl 90 than
  60); (2) too little covering/column (τ_T ~0.07, 10° wall covers ~11% of 4π);
  (3) soft flux has a sharp ridge along logxi ≈ logNHs − 21 (shell just ionised
  through), so 0.5-dex logxi steps were too coarse and the fit stuck on the grid
  node logxi = 2.
- Partial covering (`zpcfabs`) on stage 1 gave no improvement: dropped.
- Stage-1 torus angle scan (ΔC): 3C 234 flat over 45–60 (best 50); 3C 445 needs
  ~60 (+16 at 50, +26 at 45).
- Luminosity check: L(1–1000 Ryd) ≈ 7.8e44 (3C 234), 1.3e44 erg/s (3C 445);
  ξ·N_H ≈ L/R_in, so the fitted logxi ~2 at logNHs ~23 puts the gas at ~25 pc
  (3C 234) / ~4 pc (3C 445), supporting a shell outside the torus.

Changes made (uncommitted; the user commits):
- shell moved to 1e19–1e21 cm, wall 20°; per-target combined inclinations; cone
  check in `load_target`; spatial grid spans torus + shell (12-digit radial mesh);
- new grids, 405 runs each (was 810): torus logNH and α fixed (3C 234: 23.45,
  α 45; 3C 445: 23.34, α 60), logxi step 0.25 (3C 234: 1.5–3.5, 3C 445:
  1.0–3.0), logNHs 23–24 step 0.25, θ_s 3C 234 {20, 27.5, 35} / 3C 445
  {10, 15, 20}, Γ = stage-1 90% interval (unchanged);
- fit summaries `fit_stage{2,3}.txt`; config moved into `model/config/`;
  READMEs and grid.sh updated.

Next steps / open items:
1. The user runs `scripts/radial_test.py` on the server (3C 445, logNHs 24,
   logxi 2): if n30 differs from n60 no more than the seed-1 rerun does, set
   `SHELL_N_RADIAL = 30` (halves memory and cell-building time).
2. Then `./grid.sh 3C234` / `3C445` on the server with the new grids, and refit.
   Check `fit_stage3.txt` for parameters at limits and steppar logxi/logNHs for
   grid-node artefacts.
3. 3C 445 may still lack soft flux: α = 60 limits the cone to ≤ 30° from the
   pole (covering ≤ 0.12); more N_H saturates (τ_T ~1 at 24). That torus-vs-cone
   tension may itself be the result (3C 445 is a BLRG).
4. Optional ideas the user declined for now: rotating the logxi axis along the
   ridge (Δ = logxi − logNHs + 21) to cut logNHs to 3 values; ion cutoff 1e-4
   (keep 1e-6: at 1e-4 trace ions can reach τ ~10 in strong lines at logNHs 24);
   extra phenomenological components.
5. `docs/paper_context.md` still refers to `model/config.py` and
   `model/targets/`; update when the paper text is revised.
