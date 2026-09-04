"""
4_build_xspec_table.py
========================
Step 4 of the pipeline: assemble an XSPEC additive table model (atable) FITS
file from the grid of SKIRT SED outputs (run_*/skirt/model_i00_sed.dat), using
grid_manifest.csv for the parameter grid.

FORMAT: strictly follows OGIP Memo 92-009 ("The File Format for XSPEC Table
Models" -- https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/
ogip_92_009/ogip_92_009.html), confirmed directly against that memo's
PARAMETERS/ENERGIES/SPECTRA extension pages and its worked FITS-header
examples:
  - PRIMARY header: MODLNAME, MODLUNIT='photons/cm^2/s', REDSHIFT (bool),
	ADDMODEL=True, HDUCLASS='OGIP', HDUCLAS1='XSPEC TABLE MODEL'.
  - PARAMETERS extension (BinTableHDU): one row per grid axis --
	NAME, METHOD (0=linear/1=log), INITIAL, DELTA, MINIMUM, BOTTOM, TOP,
	MAXIMUM, NUMBVALS, VALUE (tabulated grid values, must be increasing).
  - ENERGIES extension: ENERG_LO, ENERG_HI (keV) bin edges, shared by every
	spectrum.
  - SPECTRA extension: one row per grid point -- PARAMVAL (this point's
	parameter values) and INTPSPEC (the spectrum, photons/cm2/s PER BIN,
	i.e. already multiplied by bin width -- NOT a flux density). Row order
	must have the LAST parameter changing fastest, e.g. for 2 params each
	with values (1,2,3): (1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1)...
	-- this is exactly the order itertools.product(*axis_values) produces,
	which is how 1_cloudy_generator.py built grid_manifest.csv in the first
	place, so grid_manifest.csv's row order is ALREADY correct as long as
	`axis_names` below is given in the same left-to-right order as that
	script's PARAMETER_SPACE dict.

WHAT'S USED FROM SKIRT'S OUTPUT
----------------------------------
Per the user's SED header:
  # column 1: wavelength; E (keV)
  # column 2: total flux; F_E (1/s/cm2/keV)
  ...
Column 2 ("total flux") is used -- it's the full observable spectrum
(direct + scattered + secondary), matching what an actual XMM observation
would see. F_E is already a PHOTON flux density (1/s/cm2/keV), so getting
to OGIP's required photons/cm2/s-per-bin needs only:
	INTPSPEC[i] = F_E[i] * (E_hi[i] - E_lo[i])
-- no erg/energy-weighting conversion needed at all.

Column 1 only gives bin CENTERS (confirmed log-spaced by the user), so bin
edges are reconstructed via the geometric mean of consecutive centers
(standard for a log-uniform grid); the two outermost edges are extrapolated
using the same log spacing. A uniformity check warns if the grid isn't
actually log-uniform, since that's what the reconstruction assumes.
"""

from __future__ import annotations
import csv
import glob
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits

# --------------------------------------------------------------------------
# Parsing SKIRT's model.sed
# --------------------------------------------------------------------------

_SED_HEADER_RE = re.compile(r"#\s*column\s*(\d+)\s*:\s*([^;]+);", re.IGNORECASE)


def read_model_sed(path: Path):
	"""Return (E_kev[], {column_name: values[]}) from a SKIRT model_i00_sed.dat
	file, using the '# column N: <name>; ...' header comments to identify
	columns by name (robust to column order/count varying between runs)."""
	col_names = {}  # 1-based column index -> name
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

	if 1 not in col_names:
		raise ValueError(f"{path}: could not find '# column 1: ...' header")
	n_cols = max(col_names)
	data = {col_names[i]: [] for i in col_names}
	for row in rows:
		for i in range(1, n_cols + 1):
			if i in col_names:
				data[col_names[i]].append(row[i - 1])

	e_kev = data.pop("wavelength")  # column 1, always present per the header format
	return e_kev, data


def select_flux_column(columns: dict, name: str = "total flux"):
	key = name.strip().lower()
	if key not in columns:
		raise KeyError(f"Column '{name}' not found. Available: {list(columns.keys())}")
	return columns[key]


# --------------------------------------------------------------------------
# Reconstructing bin edges from log-spaced bin centers
# --------------------------------------------------------------------------

def reconstruct_log_bin_edges(e_centers):
	"""Given log-uniformly-spaced bin centers, return (E_lo[], E_hi[]) via
	the geometric mean of consecutive centers for internal edges, and the
	same log spacing extrapolated for the two outermost edges. Warns if the
	centers aren't actually consistently log-spaced (the reconstruction
	assumes they are)."""
	e = np.asarray(e_centers, dtype=float)
	n = len(e)
	if n < 2:
		raise ValueError("Need at least 2 bin centers to reconstruct edges")

	log_ratios = np.diff(np.log(e))
	mean_ratio = log_ratios.mean()
	if not np.allclose(log_ratios, mean_ratio, rtol=1e-2):
		print("warning: bin centers are not consistently log-spaced "
			  f"(log-ratio std/mean = {log_ratios.std()/abs(mean_ratio):.3f}) -- "
			  "edge reconstruction may be inaccurate")

	internal = np.sqrt(e[:-1] * e[1:])  # geometric mean, n-1 internal edges
	e_lo0 = e[0] * math.exp(-mean_ratio / 2.0)
	e_hi_last = e[-1] * math.exp(mean_ratio / 2.0)

	e_lo = np.concatenate(([e_lo0], internal))
	e_hi = np.concatenate((internal, [e_hi_last]))
	return e_lo, e_hi


# --------------------------------------------------------------------------
# Parameter grid
# --------------------------------------------------------------------------

def gather_axes(manifest_rows, axis_names):
	"""Sorted unique values per axis, plus a completeness check (does the
	manifest have exactly one row per point of the full Cartesian grid)."""
	axes = {}
	for name in axis_names:
		values = sorted({float(r[name]) for r in manifest_rows})
		axes[name] = values

	expected = 1
	for v in axes.values():
		expected *= len(v)
	if expected != len(manifest_rows):
		print(f"warning: grid looks incomplete -- expected {expected} rows "
			  f"(product of axis sizes) but manifest has {len(manifest_rows)}. "
			  "The table will still be built from what's present, but this "
			  "is not a valid rectangular grid for XSPEC interpolation "
			  "unless every combination is actually there.")
	return axes


# axis -> ("linear" or "log") interpolation method; adjust to taste.
DEFAULT_METHOD = {
	"xi": "linear",  # already log10(xi) per the manifest (e.g. -1, 0, 1, 2), so linear in this column
	"hden": "log",
	"metallicity": "linear",
	"opening_angle": "linear",
	"photon_index": "linear",
}

# axis -> FITS PARAMETERS.NAME (OGIP limits this to 12 characters -- longer
# names are silently truncated by the FITS writer, which can make two
# different axes collide, so anything long needs an explicit short form).
DEFAULT_FITS_NAME = {
	"opening_angle": "OpenAngle",
	"photon_index": "PhoIndex",
	"metallicity": "Zmetal",
	"xi": "logxi",
	"hden": "hden",
}


def build_parameters_hdu(axes: dict, methods: dict, fits_names: dict):
	names, methods_col, initial, delta = [], [], [], []
	minimum, bottom, top, maximum, numbvals = [], [], [], [], []
	max_n = max(len(v) for v in axes.values())
	value_matrix = np.zeros((len(axes), max_n))

	for i, (name, values) in enumerate(axes.items()):
		n = len(values)
		method = methods.get(name, "linear")
		fits_name = fits_names.get(name, name)
		if len(fits_name) > 12:
			raise ValueError(
				f"Parameter name '{fits_name}' (for axis '{name}') is longer than the "
				f"OGIP 12-character limit -- add a short form to fits_names/DEFAULT_FITS_NAME."
			)
		names.append(fits_name)
		methods_col.append(1 if method == "log" else 0)
		minimum.append(values[0])
		bottom.append(values[0])
		top.append(values[-1])
		maximum.append(values[-1])
		numbvals.append(n)
		initial.append(values[n // 2])
		delta.append(max((values[-1] - values[0]) / 50.0, 1e-6))
		value_matrix[i, :n] = values

	cols = [
		fits.Column(name="NAME", format="12A", array=np.array(names)),
		fits.Column(name="METHOD", format="J", array=np.array(methods_col, dtype=np.int32)),
		fits.Column(name="INITIAL", format="E", array=np.array(initial, dtype=np.float32)),
		fits.Column(name="DELTA", format="E", array=np.array(delta, dtype=np.float32)),
		fits.Column(name="MINIMUM", format="E", array=np.array(minimum, dtype=np.float32)),
		fits.Column(name="BOTTOM", format="E", array=np.array(bottom, dtype=np.float32)),
		fits.Column(name="TOP", format="E", array=np.array(top, dtype=np.float32)),
		fits.Column(name="MAXIMUM", format="E", array=np.array(maximum, dtype=np.float32)),
		fits.Column(name="NUMBVALS", format="J", array=np.array(numbvals, dtype=np.int32)),
		fits.Column(name=f"VALUE", format=f"{max_n}E", array=value_matrix.astype(np.float32)),
	]
	hdu = fits.BinTableHDU.from_columns(cols, name="PARAMETERS")
	hdu.header["NINTPARM"] = len(axes)
	hdu.header["NADDPARM"] = 0
	hdu.header["HDUCLASS"] = "OGIP"
	hdu.header["HDUCLAS1"] = "XSPEC TABLE MODEL"
	hdu.header["HDUCLAS2"] = "PARAMETERS"
	hdu.header["HDUVERS"] = "1.0.0"
	return hdu


def build_energies_hdu(e_lo, e_hi):
	cols = [
		fits.Column(name="ENERG_LO", format="E", unit="keV", array=np.asarray(e_lo, dtype=np.float32)),
		fits.Column(name="ENERG_HI", format="E", unit="keV", array=np.asarray(e_hi, dtype=np.float32)),
	]
	hdu = fits.BinTableHDU.from_columns(cols, name="ENERGIES")
	hdu.header["HDUCLASS"] = "OGIP"
	hdu.header["HDUCLAS1"] = "XSPEC TABLE MODEL"
	hdu.header["HDUCLAS2"] = "ENERGIES"
	hdu.header["HDUVERS"] = "1.0.0"
	return hdu


def build_spectra_hdu(param_rows, spectra, n_energy_bins):
	n_params = len(param_rows[0])
	paramval_col = np.array(param_rows, dtype=np.float32)
	intpspec_col = np.array(spectra, dtype=np.float32)
	cols = [
		fits.Column(name="PARAMVAL", format=f"{n_params}E", array=paramval_col),
		fits.Column(name="INTPSPEC", format=f"{n_energy_bins}E", unit="photons/cm2/s", array=intpspec_col),
	]
	hdu = fits.BinTableHDU.from_columns(cols, name="SPECTRA")
	hdu.header["HDUCLASS"] = "OGIP"
	hdu.header["HDUCLAS1"] = "XSPEC TABLE MODEL"
	hdu.header["HDUCLAS2"] = "MODEL SPECTRA"
	hdu.header["HDUVERS"] = "1.0.0"
	return hdu


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def build_table(cloudy_root: str, *, axis_names,
				 sed_glob: str = "model_i00_sed.dat", flux_column: str = "total flux",
				 model_name: str = "skirtxray", redshift: bool = False,
				 methods: dict = None, fits_names: dict = None):
	root = Path(cloudy_root)
	manifest_path = root / "grid_manifest.csv"
	with open(manifest_path, newline="") as f:
		manifest_rows = list(csv.DictReader(f))

	methods = {**DEFAULT_METHOD, **(methods or {})}
	fits_names = {**DEFAULT_FITS_NAME, **(fits_names or {})}
	axes = gather_axes(manifest_rows, axis_names)

	e_lo = e_hi = None
	param_rows, spectra = [], []

	for row in manifest_rows:
		skirt_dir = root / row["run"] / "skirt"
		matches = sorted(glob.glob(str(skirt_dir / sed_glob)))
		if len(matches) != 1:
			raise FileNotFoundError(
				f"{skirt_dir}: expected exactly 1 file matching '{sed_glob}', found {len(matches)}"
			)
		e_centers, columns = read_model_sed(Path(matches[0]))
		flux = select_flux_column(columns, flux_column)

		# ascending-energy order, guaranteed regardless of the file's native order
		order = np.argsort(e_centers)
		e_centers = np.asarray(e_centers)[order]
		flux = np.asarray(flux)[order]

		lo, hi = reconstruct_log_bin_edges(e_centers)
		if e_lo is None:
			e_lo, e_hi = lo, hi
		elif not (np.allclose(lo, e_lo, rtol=1e-4) and np.allclose(hi, e_hi, rtol=1e-4)):
			raise ValueError(
				f"{row['run']}: energy grid differs from the reference run -- "
				"all runs must share the same wavelength grid"
			)

		intpspec = flux * (hi - lo)  # photons/cm2/s per bin
		param_rows.append([float(row[name]) for name in axis_names])
		spectra.append(intpspec)

	primary = fits.PrimaryHDU()
	primary.header["MODLNAME"] = model_name
	primary.header["MODLUNIT"] = "photons/cm^2/s"
	primary.header["REDSHIFT"] = redshift
	primary.header["ADDMODEL"] = True
	primary.header["HDUCLASS"] = "OGIP"
	primary.header["HDUCLAS1"] = "XSPEC TABLE MODEL"
	primary.header["HDUVERS"] = "1.0.0"

	output_fits = model_name + ".fits"
	hdul = fits.HDUList([
		primary,
		build_parameters_hdu(axes, methods, fits_names),
		build_energies_hdu(e_lo, e_hi),
		build_spectra_hdu(param_rows, spectra, len(e_lo)),
	])
	hdul.writeto(output_fits, overwrite=True)

	print(f"Wrote {output_fits}: {len(manifest_rows)} grid point(s), "
		  f"{len(axis_names)} parameter(s), {len(e_lo)} energy bin(s).")
	for name, values in axes.items():
		print(f"  {fits_names.get(name, name)} ({name}): {len(values)} values, "
			  f"{values[0]:.4g} to {values[-1]:.4g} ({methods.get(name, 'linear')})")


if __name__ == "__main__":
	import sys
	root = sys.argv[1] if len(sys.argv) > 1 else "cloudy_grid"
	# Must match the left-to-right key order of PARAMETER_SPACE in
	# 1_cloudy_generator.py, so SPECTRA row order (already correct in
	# grid_manifest.csv) stays consistent with the "last parameter changes
	# fastest" OGIP requirement.
	axis_names = ["xi", "hden", "metallicity", "opening_angle", "photon_index"]
	build_table(root, axis_names=axis_names, sed_glob="model_i00_sed.dat", model_name="skirt_i00")
	build_table(root, axis_names=axis_names, sed_glob="model_i90_sed.dat", model_name="skirt_i90")