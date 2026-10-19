"""
build_xspec.py
==============
Assemble OGIP XSPEC additive table models (atable)
from the SKIRT SEDs of the torus and combined grids.

  python3 steps/build_xspec.py torus    -> torus.fits (shared by all targets)
  python3 steps/build_xspec.py TARGET   -> runs/TARGET/combined.fits
  (or: ./grid.sh torus|TARGET --only table)

TABLES WRITTEN
--------------
  torus.fits      torus-only grid (stage 1): logNH, OpenAngle, PhoIndex, Incl
  combined.fits   torus + ionised shell, one simulation (stage 2/3):
                  [logNH, OpenAngle,] logxi, logNHs, OpenAngS, PhoIndex, Incl
Both use SKIRT's "total flux" (direct + scattered + fluorescence). The
inclination axis comes for free: every SKIRT run records several inclinations
(config.TORUS_INCLINATIONS, or config.combined_inclinations() of the target). A grid axis with a single value (e.g. a torus parameter
fixed at the stage-1 result) is not a table parameter; its value is written
as a COMMENT in the primary header.

FORMAT
------
OGIP Memo 92-009 ("The File Format for XSPEC Table Models"):
  PRIMARY:    MODLNAME, MODLUNIT, REDSHIFT, ADDMODEL, HDUCLASS/HDUCLAS1.
  PARAMETERS: NAME, METHOD, INITIAL, DELTA, MINIMUM, BOTTOM, TOP, MAXIMUM,
              NUMBVALS, VALUE -- one row per grid axis.
  ENERGIES:   ENERG_LO, ENERG_HI (keV), shared by all spectra.
  SPECTRA:    PARAMVAL, INTPSPEC (photons/cm^2/s per bin); row order with the
              LAST parameter changing fastest = itertools.product order, which
              is how the generators wrote grid_manifest.csv.

CHANGES W.R.T. THE FIRST VERSION
--------------------------------
- REDSHIFT = True: energies are REST frame and XSPEC adds a frozen-able "z"
  parameter, so no separate zashift is needed (as in Bert's tables).
- Every spectrum is divided by its own intrinsic photon flux in the rest-frame
  2-10 keV band (SKIRT "transparent flux" column). Previously the absolute
  flux scaled with xi * n (CLOUDY's luminosity), by orders of magnitude across
  the grid, which made XSPEC interpolate normalisation instead of shape. Now:
    norm = intrinsic rest-frame 2-10 keV photon flux of the source,
  identical in meaning for torus and shell tables, so their norms can be
  linked in XSPEC.

SKIRT gives F_E (photons/s/cm^2/keV) at bin centres; INTPSPEC = F_E * dE with
bin edges reconstructed from the log-spaced centres.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # config.py is one level up
import config as C

# --------------------------------------------------------------------------
# Parsing SKIRT SEDs
# --------------------------------------------------------------------------

_SED_HEADER_RE = re.compile(r"#\s*column\s*(\d+)\s*:\s*([^;]+);", re.IGNORECASE)


def read_model_sed(path: Path):
	"""(E_kev[], {column_name: values[]}) from a SKIRT *_sed.dat file."""
	col_names = {}
	rows = []
	for line in Path(path).read_text().splitlines():
		line = line.strip()
		if not line:
			continue
		if line.startswith("#"):
			m = _SED_HEADER_RE.match(line)
			if m:
				col_names[int(m.group(1))] = m.group(2).strip().lower()
			continue
		rows.append([float(x) for x in line.split()])
	data = np.array(rows)
	columns = {name: data[:, i - 1] for i, name in col_names.items()}
	e_kev = columns.pop("wavelength")
	return e_kev, columns


def reconstruct_log_bin_edges(e_centers):
	"""Bin edges from log-uniform centres (geometric means; outer edges extrapolated)."""
	e = np.asarray(e_centers, dtype=float)
	log_ratios = np.diff(np.log(e))
	mean_ratio = log_ratios.mean()
	if not np.allclose(log_ratios, mean_ratio, rtol=1e-2):
		print("warning: bin centres are not consistently log-spaced -- edge reconstruction may be inaccurate")
	internal = np.sqrt(e[:-1] * e[1:])
	e_lo = np.concatenate(([e[0] * math.exp(-mean_ratio / 2.0)], internal))
	e_hi = np.concatenate((internal, [e[-1] * math.exp(mean_ratio / 2.0)]))
	return e_lo, e_hi


def band_photon_flux(f_e, e_lo, e_hi, band):
	"""Integral of F_E over [band], counting partial bins by overlap fraction."""
	lo = np.clip(e_lo, band[0], band[1])
	hi = np.clip(e_hi, band[0], band[1])
	return float(np.sum(f_e * np.clip(hi - lo, 0.0, None)))

# --------------------------------------------------------------------------
# HDUs
# --------------------------------------------------------------------------

def _ogip(hdu, clas2=None):
	hdu.header["HDUCLASS"] = "OGIP"
	hdu.header["HDUCLAS1"] = "XSPEC TABLE MODEL"
	if clas2:
		hdu.header["HDUCLAS2"] = clas2
	hdu.header["HDUVERS"] = "1.0.0"
	return hdu


def build_parameters_hdu(axes: dict, fits_names: dict):
	max_n = max(len(v) for v in axes.values())
	value_matrix = np.zeros((len(axes), max_n))
	names, initial, delta, vmin, vmax, numbvals = [], [], [], [], [], []
	for i, (name, values) in enumerate(axes.items()):
		fits_name = fits_names[name]
		if len(fits_name) > 12:
			raise ValueError(f"FITS parameter name '{fits_name}' exceeds the OGIP 12-character limit")
		n = len(values)
		names.append(fits_name)
		initial.append(values[n // 2])
		delta.append(max((values[-1] - values[0]) / 50.0, 1e-3))
		vmin.append(values[0])
		vmax.append(values[-1])
		numbvals.append(n)
		value_matrix[i, :n] = values

	f32 = lambda x: np.array(x, dtype=np.float32)
	cols = [
		fits.Column(name="NAME", format="12A", array=np.array(names)),
		# every axis is already linear in the quantity we want to interpolate in
		# (logNH, logxi, angle, Gamma), so METHOD = 0 (linear) throughout
		fits.Column(name="METHOD", format="J", array=np.zeros(len(axes), dtype=np.int32)),
		fits.Column(name="INITIAL", format="E", array=f32(initial)),
		fits.Column(name="DELTA", format="E", array=f32(delta)),
		fits.Column(name="MINIMUM", format="E", array=f32(vmin)),
		fits.Column(name="BOTTOM", format="E", array=f32(vmin)),
		fits.Column(name="TOP", format="E", array=f32(vmax)),
		fits.Column(name="MAXIMUM", format="E", array=f32(vmax)),
		fits.Column(name="NUMBVALS", format="J", array=np.array(numbvals, dtype=np.int32)),
		fits.Column(name="VALUE", format=f"{max_n}E", array=value_matrix.astype(np.float32)),
	]
	hdu = fits.BinTableHDU.from_columns(cols, name="PARAMETERS")
	hdu.header["NINTPARM"] = len(axes)
	hdu.header["NADDPARM"] = 0
	return _ogip(hdu, "PARAMETERS")


def build_energies_hdu(e_lo, e_hi):
	cols = [
		fits.Column(name="ENERG_LO", format="E", unit="keV", array=np.asarray(e_lo, dtype=np.float32)),
		fits.Column(name="ENERG_HI", format="E", unit="keV", array=np.asarray(e_hi, dtype=np.float32)),
	]
	return _ogip(fits.BinTableHDU.from_columns(cols, name="ENERGIES"), "ENERGIES")


def build_spectra_hdu(param_rows, spectra):
	cols = [
		fits.Column(name="PARAMVAL", format=f"{len(param_rows[0])}E", array=np.array(param_rows, dtype=np.float32)),
		fits.Column(name="INTPSPEC", format=f"{len(spectra[0])}E", unit="photons/cm2/s",
					array=np.array(spectra, dtype=np.float32)),
	]
	return _ogip(fits.BinTableHDU.from_columns(cols, name="SPECTRA"), "MODEL SPECTRA")

# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def build_table(*, grid_dir: Path, grid: dict, inclinations, model_name: str, out: Path,
				flux_columns=("total flux",)):
	"""One atable from grid_dir/run_*/model_iNN_sed.dat over the varying axes of
	`grid` plus inclination (last, fastest), normalised by the rest-frame
	intrinsic band flux."""
	with open(grid_dir / "grid_manifest.csv", newline="") as f:
		manifest_rows = list(csv.DictReader(f))
	if len(manifest_rows) != C.grid_size(grid):
		raise ValueError(f"{grid_dir}: manifest has {len(manifest_rows)} rows but the config grid has "
						 f"{C.grid_size(grid)} points -- regenerate the grid (or finish it) first")

	varying = [name for name, values in grid.items() if len(values) > 1]
	fixed = {name: values[0] for name, values in grid.items() if len(values) == 1}

	e_lo = e_hi = None
	param_rows, spectra, missing = [], [], []
	for row in manifest_rows:
		for inc in inclinations:
			sed_path = grid_dir / row["run"] / f"model_i{inc:02d}_sed.dat"
			if not sed_path.exists():
				missing.append(sed_path.parent.name)
				continue
			e, cols = read_model_sed(sed_path)
			order = np.argsort(e)
			lo, hi = reconstruct_log_bin_edges(e[order])
			if e_lo is None:
				e_lo, e_hi = lo, hi
			elif not (np.allclose(lo, e_lo, rtol=1e-4) and np.allclose(hi, e_hi, rtol=1e-4)):
				raise ValueError(f"{sed_path}: energy grid differs from the first run")

			flux = sum(cols[name][order] for name in flux_columns)
			intrinsic = band_photon_flux(cols["transparent flux"][order], lo, hi, C.NORM_BAND_KEV)
			if intrinsic <= 0:
				raise ValueError(f"{sed_path}: zero intrinsic flux in {C.NORM_BAND_KEV} keV")
			param_rows.append([float(row[name]) for name in varying] + [float(inc)])
			spectra.append(flux * (hi - lo) / intrinsic)

	if missing:
		runs = sorted(set(missing))
		raise FileNotFoundError(f"{model_name}: {len(runs)} run(s) without SEDs, e.g. {runs[:5]} -- run the skirt-run step first")

	primary = fits.PrimaryHDU()
	primary.header["MODLNAME"] = model_name
	primary.header["MODLUNIT"] = "photons/cm^2/s"
	primary.header["REDSHIFT"] = True
	primary.header["ADDMODEL"] = True
	primary.header["NORMBAND"] = (f"{C.NORM_BAND_KEV[0]}-{C.NORM_BAND_KEV[1]} keV",
								  "norm = intrinsic rest-frame photon flux in band")
	for name, value in fixed.items():
		primary.header.add_comment(f"fixed: {name} = {value}")
	_ogip(primary)

	axes = {name: sorted(float(v) for v in grid[name]) for name in varying}
	axes["inclination"] = [float(i) for i in inclinations]
	fits.HDUList([
		primary,
		build_parameters_hdu(axes, C.FITS_NAMES),
		build_energies_hdu(e_lo, e_hi),
		build_spectra_hdu(param_rows, spectra),
	]).writeto(out, overwrite=True)

	print(f"Wrote {out}: {len(param_rows)} spectra, {len(e_lo)} bins "
		  f"({e_lo[0]:.3g}-{e_hi[-1]:.3g} keV rest)")
	for name, values in axes.items():
		print(f"  {C.FITS_NAMES[name]:>10}: {len(values)} values, {values[0]:g} .. {values[-1]:g}")
	for name, value in fixed.items():
		print(f"  {C.FITS_NAMES[name]:>10}: fixed at {value:g}")


if __name__ == "__main__":
	which = sys.argv[1] if len(sys.argv) > 1 else ""
	if which == "torus":
		build_table(grid_dir=C.TORUS_GRID_DIR, grid=C.TORUS_GRID, inclinations=C.TORUS_INCLINATIONS,
					model_name="torus", out=C.TORUS_TABLE)
	elif which:
		target = C.load_target(which)
		build_table(grid_dir=target.combined_dir, grid=target.combined_grid, inclinations=target.inclinations,
					model_name="combined",
					out=target.combined_table)
	else:
		sys.exit(f"usage: python3 steps/build_xspec.py torus|TARGET   (targets: {', '.join(C.target_names())})")
