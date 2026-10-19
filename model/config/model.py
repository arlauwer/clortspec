"""
config/model.py
===============
Single source of truth for the model (scripts do `import config as C`, which
re-exports this file). Geometry, energy range and the torus grid stay in sync
this way. Everything that belongs to one AGN (redshift, Galactic N_H, the
combined grid) lives in config/<name>.py and is read with load_target().

MODEL
-----
One point source (cut-off power law, photon index Gamma) surrounded by

1. TORUS -- neutral gas, uniform-density wedge between TORUS_R_IN and
   TORUS_R_OUT within +-alpha of the equatorial plane (alpha = torus opening
   angle, measured from the EQUATOR). No CLOUDY: all elements neutral.
   Parameterised by its equatorial column N_H (centre -> edge), so every line
   of sight through the wedge sees exactly N_H.

2. SHELL -- photoionised gas in a hollow polar bi-cone, i.e. SKIRT's
   "ConicalShellGeometry" (Stalevski et al. 2012) with radial exponent 2: the
   gas fills the cone WALLS, polar angles theta_s -+ SHELL_WALL_DEG/2 from the
   polar axis (theta_s = shell opening angle, measured from the POLE), from
   SHELL_R_IN to SHELL_R_OUT with n ~ r^-2. The shell starts at the outer
   edge of the torus (narrow-line-region scales), so an observer looking
   through the torus still sees the whole shell; the cone must stay inside the
   torus opening (theta_s + wall/2 <= 90 - alpha), the region the source
   ionises directly.
   Ionisation from CLOUDY. Because every radial ray through the wall is
   identical (same n(r), unobstructed view of the source), one 1D CLOUDY run
   describes the whole shell exactly; with n ~ r^-2 the local ionisation
   parameter is moreover constant along each ray (apart from attenuation).
   Set SHELL_WALL_DEG = None for a filled cone.

Both media live in ONE SKIRT simulation, so photons scattered in the shell can
still be absorbed by the torus on their way to an edge-on observer, and vice
versa. Every simulation records several inclinations, which become a table
parameter (TORUS_INCLINATIONS for the torus grid, combined_inclinations() for a
target's combined grid).

WORKFLOW (Bert's staged fit)
----------------------------
Stage 1: a cheap torus-only grid (TORUS_GRID, ~30 MB / few s per SKIRT run)
         fitted above 3 keV, where the shell is negligible. The tables are in
         rest-frame energies with the redshift as a table parameter, so this
         grid and torus.fits are SHARED by all targets.
Stage 2/3: the combined torus + shell grid (COMBINED_GRID of the target). Its
         torus axes should be centred on the target's stage-1 result -- an axis
         with a single value is held fixed and does not appear in the XSPEC
         table.

FILES
-----
  config/<name>.py             target constants + COMBINED_GRID (in git)
  torus_grid/, torus.fits      shared stage-1 grid and table
  runs/<name>/shell_cloudy/    CLOUDY grid of the target
  runs/<name>/combined_grid/   SKIRT grid of the target
  runs/<name>/combined.fits, fit_stage*.png, fit_stage3.xcm
  ../<name>/spectra/           reduced spectra (reduction/xmm_reduce.sh)
"""

from __future__ import annotations

import math
import os
import runpy
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MODEL_DIR = Path(__file__).resolve().parent.parent
OBS_ROOT = MODEL_DIR.parent          # <target>/spectra/ live here
TEMPLATE_DIR = MODEL_DIR / "templates"
TARGET_DIR = MODEL_DIR / "config"    # config/<name>.py, next to this file
NOT_TARGETS = {"__init__", "model"}
RUNS_DIR = MODEL_DIR / "runs"

# --------------------------------------------------------------------------
# Energy range (REST frame)
# --------------------------------------------------------------------------
# Observed 0.3-10 keV at redshift z is rest 0.3 (1+z) - 10 (1+z) keV, so the
# table must extend beyond 10 keV rest or XSPEC sees zero model flux at the
# top of the band (load_target checks this: fine up to z = 0.5).
E_MIN_KEV = 0.3
E_MAX_KEV = 15.0
FIT_E_MAX_KEV = 10.0     # upper end of the fitted (observed) band
# ~0.8% per bin: finer than the EPIC resolution at Fe K (~2%).
N_ENERGY_BINS = 500

# Every table spectrum is divided by its intrinsic (unabsorbed) photon flux in
# this rest-frame band, so the XSPEC norm = intrinsic 2-10 keV photon flux.
NORM_BAND_KEV = (2.0, 10.0)

# --------------------------------------------------------------------------
# Source
# --------------------------------------------------------------------------

E_CUT_KEV = 200.0        # exponential cutoff of the primary power law
# ExtinctionOnly RT is linear in the luminosity and every table spectrum is
# renormalised, so this value is arbitrary.
SOURCE_LUMINOSITY_ERG_S = 1e43
NUM_PACKETS = "1e7"

ABUNDANCE_SET = "GASS10"   # CLOUDY abundance set, metallicity fixed at solar


def _cloudy_data_dir() -> Path:
	"""CLOUDY's data/ directory: $CLOUDY_DATA_PATH if set (CLOUDY's own
	variable), else the first data/ above the `cloudy` executable
	(c23/source/sys_gcc/cloudy.exe -> c23/data), else the local install."""
	if os.environ.get("CLOUDY_DATA_PATH"):
		return Path(os.environ["CLOUDY_DATA_PATH"].split(":")[0])
	exe = shutil.which("cloudy")
	if exe:
		for parent in Path(exe).resolve().parents:
			if (parent / "data/abundances").is_dir():
				return parent / "data"
	return Path.home() / "Documents/PhD/code/cloudy/data"


CLOUDY_ABUNDANCE_FILE = _cloudy_data_dir() / "abundances/solar_GASS10.abn"

# Observer inclinations (deg from the polar axis), recorded in every SKIRT run
# and used as the last table parameter. XSPEC interpolates linearly between
# them, also across the jump where a line of sight crosses the torus edge
# (90 - alpha), which would blend an obscured and an unobscured spectrum.
# The torus grid (stage 1, Incl frozen at the target's INCLINATION) records:
TORUS_INCLINATIONS = [0, 15, 30, 45, 60, 75, 90]
# A combined grid records only lines of sight through the torus, see
# combined_inclinations(): starting this far inside the narrowest torus of the
# grid, then every COMBINED_INCLINATION_STEP up to 90.
COMBINED_INCLINATION_MARGIN = 5
COMBINED_INCLINATION_STEP = 10


def combined_inclinations(grid: dict) -> list:
	"""Inclinations of a combined grid: strictly inside every torus of the grid
	(> 90 - min(alpha)), so every table point is an obscured view and the table
	never interpolates across the torus edge."""
	first = int(math.floor(90 - min(grid["torus_opening_angle"]))) + COMBINED_INCLINATION_MARGIN
	return sorted(set(range(first, 90, COMBINED_INCLINATION_STEP)) | {90})

# --------------------------------------------------------------------------
# Geometry (fixed)
# --------------------------------------------------------------------------

TORUS_R_IN_CM = 1e18      # ~0.32 pc
TORUS_R_OUT_CM = 1e19     # ~3.2 pc
# The shell starts where the torus ends, so it is never behind the torus for an
# observer looking through the torus (with the shell inside R_out, the torus
# hid most of it at Incl = 90). R_out / R_in = 100 as before.
SHELL_R_IN_CM = TORUS_R_OUT_CM   # ~3.2 pc
SHELL_R_OUT_CM = 1e21            # ~320 pc
SHELL_DENSITY_POWER = -2.0
# Angular thickness of the cone wall; None = filled cone. 20 deg covers about
# twice the solid angle of the earlier 10 deg wall.
SHELL_WALL_DEG = 20.0
# CLOUDY zones are rebinned (column-conserving, per ion) onto this many
# log-spaced radial shells before SKIRT. SKIRT memory scales with the number
# of imported cells (each has its own ionic mix).
SHELL_N_RADIAL = 60
# Ions below this fraction of their element (in a cell) are set to zero and the
# element's other ions scaled up to keep its abundance. Each non-zero ion adds
# its fluorescence and line channels to the cell's SKIRT mix, and SKIRT memory
# grows faster than the number of ions; trace ions add no visible opacity.
# 0 keeps every ion CLOUDY reports. Changing it changes every shell.txt, so
# all combined runs are redone.
SHELL_ION_MIN_FRACTION = 1e-6

# Resonant line scattering by H- and He-like ions in the shell (O VII/VIII,
# Ne IX/X, Mg, Si, ... -- the main soft X-ray signature of photoionised gas).
# SKIRT only honours XRayIonicGasMix(Family).resonantScattering in the
# "LyaExtinctionOnly" simulation mode, and then requires every medium to
# support polarisation, so the torus switches to polarisation-capable
# free-electron scattering (in both templates, so stage 1 and 2 match).
RESONANT_SCATTERING = True
# Core-skipping acceleration of the resonant-line random walk ("None",
# "Constant" or "Variable"). "None" is unbiased but slower in optically
# thick lines.
LYA_ACCELERATION = "None"

# "atomic": XRayAtomicGasMix (neutral reference implementation, bound-electron
#           scattering, fluorescence).
# "ionic":  XRayIonicGasMix with every element in its +0 state -- the same
#           cross sections as the shell, useful when later comparing against
#           a lightly ionised torus.
TORUS_GAS_MIX = "atomic"

# --------------------------------------------------------------------------
# Grids. Axis order = manifest/table row order (last axis fastest).
# --------------------------------------------------------------------------

# Stage 1: torus only, shared by all targets.
TORUS_GRID = {
	# log10 equatorial column density, cm^-2
	"torus_logNH": [21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5],
	# wedge half-opening angle from the EQUATOR (deg)
	"torus_opening_angle": [30.0, 45.0, 60.0],
	"photon_index": [1.0, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7],
}

# Stage 2/3: torus + shell, one grid per target (COMBINED_GRID in
# config/<name>.py) with these axes, in this order.
COMBINED_AXES = ["torus_logNH", "torus_opening_angle", "logxi", "shell_logNH", "shell_opening_angle",
				 "photon_index"]
SHELL_AXES = ["logxi", "shell_logNH", "shell_opening_angle", "photon_index"]  # what CLOUDY depends on

TORUS_GRID_DIR = MODEL_DIR / "torus_grid"
TORUS_TABLE = MODEL_DIR / "torus.fits"

# XSPEC parameter names (OGIP limit: 12 characters).
FITS_NAMES = {
	"torus_logNH": "logNH",
	"torus_opening_angle": "OpenAngle",
	"logxi": "logxi",
	"shell_logNH": "logNHs",
	"shell_opening_angle": "OpenAngS",
	"photon_index": "PhoIndex",
	"inclination": "Incl",
}

# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
	"""One AGN: the contents of config/<name>.py plus where its files live."""
	name: str                 # directory name, e.g. "3C234"
	title: str                # for plots, e.g. "3C 234"
	redshift: float
	nh_galactic_1e22: float   # Galactic N_H for TBabs, 1e22 cm^-2
	combined_grid: dict
	inclination: float        # deg; frozen in stages 1-2, start value in stage 3
	spectra_dir: Path         # {pn,mos1,mos2}_grp.fits

	@property
	def dir(self) -> Path:
		return RUNS_DIR / self.name

	@property
	def cloudy_dir(self) -> Path:
		return self.dir / "shell_cloudy"

	@property
	def combined_dir(self) -> Path:
		return self.dir / "combined_grid"

	@property
	def combined_table(self) -> Path:
		return self.dir / "combined.fits"

	@property
	def inclinations(self) -> list:
		"""Inclinations recorded in the combined grid."""
		return combined_inclinations(self.combined_grid)


def target_names():
	return sorted(p.stem for p in TARGET_DIR.glob("*.py") if p.stem not in NOT_TARGETS)


def load_target(name: str) -> Target:
	"""Read config/<name>.py (see config/3C234.py for the format)."""
	path = TARGET_DIR / f"{name}.py"
	if name in NOT_TARGETS or not path.exists():
		raise SystemExit(f"unknown target '{name}': no {path} (targets: {', '.join(target_names())})")
	conf = runpy.run_path(str(path))
	grid = conf["COMBINED_GRID"]
	if list(grid) != COMBINED_AXES:
		raise SystemExit(f"{path}: COMBINED_GRID must have the axes {COMBINED_AXES}, in this order")
	# every cone must fit in every torus opening: a clipped wall would no longer
	# match the covering factor its CLOUDY run was computed with
	alpha_max = max(grid["torus_opening_angle"])
	theta_max = max(grid["shell_opening_angle"])
	if shell_wall(theta_max)[1] > 90.0 - alpha_max + 1e-9:
		raise SystemExit(f"{path}: shell_opening_angle {theta_max:g} + wall/2 = {shell_wall(theta_max)[1]:g} deg reaches "
						 f"into the torus with alpha = {alpha_max:g} (edge at {90 - alpha_max:g} deg from the pole)")
	if FIT_E_MAX_KEV * (1.0 + conf["REDSHIFT"]) > E_MAX_KEV:
		raise SystemExit(f"{path}: z = {conf['REDSHIFT']} needs tables beyond E_MAX_KEV = {E_MAX_KEV} keV rest")
	return Target(
		name=name,
		title=conf.get("TITLE", name),
		redshift=conf["REDSHIFT"],
		nh_galactic_1e22=conf["NH_GALACTIC_1E22"],
		combined_grid=grid,
		inclination=conf.get("INCLINATION", 90.0),
		spectra_dir=Path(conf.get("SPECTRA_DIR", OBS_ROOT / name / "spectra")),
	)

# --------------------------------------------------------------------------
# Derived quantities and helpers
# --------------------------------------------------------------------------

RYD_KEV = 13.605693122994e-3


def shell_wall(opening_angle_deg: float):
	"""(theta_in, theta_out) of one cone wall, deg from the pole."""
	if SHELL_WALL_DEG is None:
		return 0.0, opening_angle_deg
	half = 0.5 * SHELL_WALL_DEG
	return max(0.0, opening_angle_deg - half), opening_angle_deg + half


def shell_covering_factor(opening_angle_deg: float) -> float:
	"""Fraction of 4 pi covered by the two cone walls."""
	th_in, th_out = shell_wall(opening_angle_deg)
	return math.cos(math.radians(th_in)) - math.cos(math.radians(th_out))


def torus_covering_factor(opening_angle_deg: float) -> float:
	"""Fraction of 4 pi covered by a wedge of half-angle alpha around the equator."""
	return math.sin(math.radians(opening_angle_deg))


def shell_inner_density(log_nh: float) -> float:
	"""n_in such that integral_{R_in}^{R_out} n_in (r/R_in)^p dr = N_H."""
	r_in, r_out, p = SHELL_R_IN_CM, SHELL_R_OUT_CM, SHELL_DENSITY_POWER
	if p == -1.0:
		integral = r_in * math.log(r_out / r_in)
	else:
		integral = r_in * ((r_out / r_in) ** (p + 1.0) - 1.0) / (p + 1.0)
	return 10.0 ** log_nh / integral


def intrinsic_nufnu(e_kev, photon_index: float):
	"""Primary continuum shape as nu F_nu (arbitrary units).

	Cut-off power law N(E) ~ E^-Gamma exp(-E/E_cut) above 1 Ryd; below 1 Ryd it
	is rolled off steeply (nu F_nu ~ E^3) so CLOUDY gets a defined but negligible
	IR/optical field -- the grid is about the X-ray ionisation parameter, not a
	full AGN SED.
	"""
	e = np.asarray(e_kev, dtype=float)
	nufnu = e ** (2.0 - photon_index) * np.exp(-e / E_CUT_KEV)
	e0 = RYD_KEV
	low = e < e0
	nufnu0 = e0 ** (2.0 - photon_index) * math.exp(-e0 / E_CUT_KEV)
	nufnu[low] = nufnu0 * (e[low] / e0) ** 3
	return nufnu


def skirt_sed(photon_index: float, e_lo_kev: float = 0.1, e_hi_kev: float = 100.0, n_points: int = 2000) -> str:
	"""SKIRT FileSED with the intrinsic spectrum (same shape CLOUDY gets)."""
	e = np.logspace(math.log10(e_lo_kev), math.log10(e_hi_kev), n_points)
	nufnu = intrinsic_nufnu(e, photon_index)
	lines = ["# Column 1: Wavelength (keV)", "# Column 2: neutralfluxdensity (erg/s/cm2)"]
	lines += [f"{a:.10e} {b:.10e}" for a, b in zip(e, nufnu)]
	return "\n".join(lines) + "\n"


def sync_run_dir(run_dir: Path, inputs: dict) -> bool:
	"""Write the input files {name: text} of one grid run.

	If the directory already holds exactly these inputs it is left alone, so
	its results stay valid. If an input changed (the grid was edited and this
	run number now means other parameters), everything else in the directory
	-- the results of the old inputs -- is deleted, so the run is redone.
	Returns True if old results were deleted."""
	run_dir.mkdir(parents=True, exist_ok=True)
	same = all((run_dir / name).exists() and (run_dir / name).read_text() == text for name, text in inputs.items())
	if same:
		return False
	stale = [f for f in run_dir.iterdir() if f.is_file() and f.name not in inputs]
	for f in stale:
		f.unlink()
	for name, text in inputs.items():
		(run_dir / name).write_text(text)
	return bool(stale)


def report_sync(root: Path, n_runs: int, n_reset: int) -> None:
	"""Tell the user about results invalidated by a grid change."""
	if n_reset:
		print(f"  {n_reset} run(s) had results for other parameters (the grid changed): results deleted, will be redone")
	extra = sorted(d.name for d in root.glob("run_*") if d.is_dir() and int(d.name[4:]) >= n_runs)
	if extra:
		print(f"  WARNING: {len(extra)} run directories ({extra[0]} .. {extra[-1]}) are no longer part of the grid; "
			  f"delete them, or the run step will still simulate them")


def read_cloudy_abundances(path: Path = CLOUDY_ABUNDANCE_FILE):
	"""Number abundances relative to H for Z = 1..30 from CLOUDY's abundance file,
	so torus and shell share one composition."""
	values = []
	for line in Path(path).read_text().splitlines():
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		if line.startswith("*"):
			break
		values.append(float(line.split()[1]))
	if len(values) != 30:
		raise ValueError(f"expected 30 abundances in {path}, found {len(values)}")
	return values


ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
			"S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]


def grid_size(grid: dict) -> int:
	n = 1
	for v in grid.values():
		n *= len(v)
	return n
