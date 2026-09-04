"""
3_cloudy_skirt_converter.py
=============================
Step 3 of the pipeline: for every run_*/cloudy/ directory produced by
1_cloudy_generator.py and run by 2_cloudy_run.fish, read the Cloudy output
(sim.zones for the radial grid, sim.ovr for temperature, sim.species for
per-ion densities) together with that run's row in grid_manifest.csv, and
write a spatially-resolved SKIRT .ski file (run_*/skirt/model.ski +
run_*/skirt/medium.txt) with one radial shell per Cloudy zone, turned into
a torus wedge using the opening angle.

FILE FORMATS (confirmed against user-supplied examples/headers)
-----------------------------------------------------------------
- sim.zones: tab-separated, header "#NZONE  radius  depth  dr". radius is
  the absolute distance from the origin to the zone CENTER (cm); depth is
  ignored (per the user: "useless"); dr is the zone width (cm). Zone
  boundaries are therefore rmin = radius - dr/2, rmax = radius + dr/2.
  This is the sole source of the radial grid -- sim.ovr's own depth column
  is a zone-center value too and is no longer used for geometry.
- sim.ovr: tab-separated, header line starting with '#', containing a
  "Te" column (K) and a total-hydrogen-density column ("hden" or "Htot")
  used as the reference density -- ion abundances in medium.txt are
  fractions relative to this column (abundance * number_density = actual
  ion number density, per the user's spec). Matched to sim.zones by row
  order (same zone, same order, across all Cloudy save files).
- sim.species: tab-separated, header line starting with '#'; first column
  is depth (ignored, same row-order matching as above), remaining columns
  are species labels in Cloudy's own notation ("H", "H+", "He+2", ...)
  giving that species' absolute number density (cm^-3) in that zone.
- Cloudy species label -> (element, charge) via regex, per the user's
  spec: no '+' = neutral (charge 0); '+' alone = singly ionized (charge
  1); '+N' = charge N.
- medium.txt (SphericalCellMedium import file): tab-separated columns
  rmin(pc) thetamin(deg) phimin(deg) rmax(pc) thetamax(deg) phimax(deg)
  number_density(1/cm3) temperature(K) then one column per ion (in the
  exact order declared in the template's XRayIonicGasMixFamily ions=
  attribute), each an abundance fraction relative to the number-density
  column. phi is always "0 0" (autoRevolve=Azimuth handles the full
  torus). Every Cloudy zone becomes one wedge row at
  theta in [90-openingAngle, 90+openingAngle]; two additional rows (zero
  density) cover the polar caps above/below the wedge so the grid spans
  the full sphere.

ASSUMPTIONS (flagged so you can double check)
------------------------------------------------
- sim.zones' "radius" is the absolute distance from the origin (i.e. it
  already includes the inner_radius_cm set via Cloudy's `radius` command),
  not a depth-into-the-cloud value -- so zone boundaries are used directly
  with no offset added. If your Cloudy version reports radius as a
  depth-from-illuminated-face instead, add inner_radius_cm (from the
  manifest) to rmin/rmax before using them.
- sim.zones' "radius"/"dr" are in cm (Cloudy's internal default unit).
- Reference density column: "hden" if present in sim.ovr, else "Htot".
- Luminosity is derived from the manifest's "intensity" (erg/s/cm^2 at the
  inner face, 0.3-10 keV) and the actual inner radius from sim.zones
  (first zone's rmin), assuming an isotropic point source:
  L = intensity * 4*pi*inner_radius_cm^2 (matches the same assumption used
  when Cloudy's own `intensity` command was set up in script 1).
- Zones are assumed contiguous (zone i's rmax == zone i+1's rmin); the
  shared radial mesh used for the SKIRT spatial grid is built from
  [rmin[0], rmax[0], rmax[1], ..., rmax[-1]].
- Empty polar-cap cells get density=0 and a dummy abundance vector
  (H+0=1, everything else 0) purely so the material mix is well-defined;
  since density is zero there, it has no effect on the radiative transfer.
"""

from __future__ import annotations
import csv
import math
import re
from pathlib import Path

# --------------------------------------------------------------------------
# Cloudy species label parsing
# --------------------------------------------------------------------------

_ION_RE = re.compile(r"^([A-Za-z]+)(\+?)([0-9]*)$")


def parse_cloudy_ion_label(label: str):
	"""'H' -> ('H',0); 'H+' -> ('H',1); 'He+2' -> ('He',2)."""
	m = _ION_RE.match(label.strip())
	if not m:
		return None
	symbol, plus, digits = m.groups()
	if not plus:
		charge = 0
	elif digits:
		charge = int(digits)
	else:
		charge = 1
	return symbol, charge


def canonical_ion_label(symbol: str, charge: int) -> str:
	return f"{symbol}+{charge}"


# --------------------------------------------------------------------------
# Template ion list
# --------------------------------------------------------------------------

def extract_ion_list_from_template(template_path: str):
	"""Pull the full, ordered ion list out of the template's
	XRayIonicGasMixFamily ions="..." attribute (comma-separated)."""
	text = Path(template_path).read_text()
	m = re.search(r'ions="([^"]*)"', text)
	if not m:
		raise ValueError(f"Could not find ions=\"...\" attribute in {template_path}")
	return [tok.strip() for tok in m.group(1).split(",") if tok.strip()]


# --------------------------------------------------------------------------
# Cloudy output parsing
# --------------------------------------------------------------------------

def _read_header_and_rows(path: Path):
	"""Return (header_tokens, list_of_float_row_lists) for a tab-separated
	Cloudy save file whose first non-empty line is the '#'-prefixed
	header."""
	lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
	if not lines:
		raise ValueError(f"{path} is empty")
	header = lines[0].lstrip("#").split("\t")
	header = [h.strip() for h in header]
	rows = []
	for line in lines[1:]:
		rows.append([float(x) for x in line.split("\t")])
	return header, rows


def parse_sim_zones(path: Path):
	"""Return (rmin_cm[], rmax_cm[]) from sim.zones (#NZONE radius depth dr).
	radius = absolute distance to zone center; dr = zone width; depth is
	ignored. rmin/rmax = radius -+ dr/2."""
	header, rows = _read_header_and_rows(path)
	lower = [h.lower() for h in header]
	if "radius" not in lower or "dr" not in lower:
		raise ValueError(f"sim.zones header missing radius/dr columns. Found: {header}")
	i_radius = lower.index("radius")
	i_dr = lower.index("dr")
	rmin, rmax = [], []
	for r in rows:
		radius, dr = r[i_radius], r[i_dr]
		rmin.append(radius - dr / 2.0)
		rmax.append(radius + dr / 2.0)
	return rmin, rmax


def parse_sim_ovr(path: Path):
	"""Return (Te_K[])."""
	header, rows = _read_header_and_rows(path)
	lower = [h.lower() for h in header]

	def find(name):
		if name.lower() in lower:
			return lower.index(name.lower())
		return None

	i_te = find("Te")
	if i_te is None:
		raise ValueError(
			f"sim.ovr header missing required column (Te). "
			f"Found columns: {header}"
		)
	te = [r[i_te] for r in rows]
	return te


def parse_sim_species(path: Path):
	"""Return {canonical_ion_label: [density_cm3, ...]} (row order = zone order)."""
	header, rows = _read_header_and_rows(path)
	# first column is depth (ignored), regardless of its exact header text
	species_cols = {}  # canonical_label -> column index (1-based into row)
	for col_idx, label in enumerate(header[1:], start=1):
		parsed = parse_cloudy_ion_label(label)
		if parsed is None:
			continue  # skip anything that doesn't look like an ion label
		species_cols[canonical_ion_label(*parsed)] = col_idx

	densities = {ion: [r[col] for r in rows] for ion, col in species_cols.items()}
	return densities


# --------------------------------------------------------------------------
# medium.txt construction
# --------------------------------------------------------------------------

def build_medium_txt(*, ion_list, rmin_cm, rmax_cm, te_K, n_ref_cm3, species_density,
					  opening_angle_deg: float) -> str:
	n_zones = len(rmin_cm)
	theta_lo = 90.0 - opening_angle_deg
	theta_hi = 90.0 + opening_angle_deg

	lines = []
	lines.append("# Column 1: box rmin (cm)")
	lines.append("# Column 2: box thetamin (deg)")
	lines.append("# Column 3: box phimin (deg)")
	lines.append("# Column 4: box rmax (cm)")
	lines.append("# Column 5: box thetamax (deg)")
	lines.append("# Column 6: box phimax (deg)")
	lines.append("# Column 7: number density (1/cm3)")
	lines.append("# Column 8: temperature (K)")
	for i, ion in enumerate(ion_list, start=9):
		lines.append(f"# Column {i}: {ion} (1)")

	def fmt_row(rmin_pc, thetamin, phimin, rmax_pc, thetamax, phimax, ndens, temp, abunds):
		cols = [f"{rmin_pc:.10e}", f"{thetamin:.6f}", f"{phimin:.1f}",
				f"{rmax_pc:.10e}", f"{thetamax:.6f}", f"{phimax:.1f}",
				f"{ndens:.6e}", f"{temp:.6e}"] + [f"{a:.6e}" for a in abunds]
		return "\t".join(cols)

	# --- wedge: one row per Cloudy zone ---
	for i in range(n_zones):
		rmin_pc, rmax_pc = rmin_cm[i], rmax_cm[i]

		nref = n_ref_cm3
		abunds = []
		for ion in ion_list:
			abunds.append(species_density.get(ion, [0.0] * n_zones)[i] / nref)

		lines.append(fmt_row(rmin_pc, theta_lo, 0.0, rmax_pc, theta_hi, 0.0,
							  nref, te_K[i], abunds))

	return "\n".join(lines) + "\n"

# --------------------------------------------------------------------------
# Incident continuum (sim.inc) -> SKIRT sed.txt, and luminosity integration
# --------------------------------------------------------------------------
 
def parse_incident_continuum(path: Path):
	"""Read sim.inc (`save incident continuum ... units keV`), columns
	Enr(keV), nFn(erg/s/cm2, = nu*F_nu), Occ Num. Returns (E_kev[], nFn[])
	for every row, in file order (Occ Num is not needed and is dropped)."""
	header, rows = _read_header_and_rows(path)
	lower = [h.lower() for h in header]

	def find(*names):
		for name in names:
			for i, h in enumerate(lower):
				if name.lower() in h:
					return i
		return None

	i_e = find("enr")
	i_nfn = find("nfn")
	if i_e is None or i_nfn is None:
		raise ValueError(f"sim.inc header missing Enr/nFn columns. Found: {header}")
	e_kev = [r[i_e] for r in rows]
	nfn = [r[i_nfn] for r in rows]
	return e_kev, nfn
 
 
def build_sed_txt_from_continuum(e_kev, nfn, sed_txt_path: Path) -> None:
	"""Write a SKIRT FileSED sed.txt straight from Cloudy's incident
	continuum: keep only the physically populated (non-zero) rows -- the
	save file spans Cloudy's full internal energy grid, almost all of
	which is zero outside the band we actually illuminated -- and copy
	nFn through unchanged (it's already nu*F_nu = erg/s/cm2, the same
	convention SKIRT expects; only the shape matters, the absolute scale
	is set separately by IntegratedLuminosityNormalization)."""
	out = ["# Column 1: Wavelength (keV)", "# Column 2: neutralfluxdensity (erg/s/cm2)"]
	for e, val in zip(e_kev, nfn):
		if val > 0:
			out.append(f"{e:.18e} {val:.18e}")
	Path(sed_txt_path).write_text("\n".join(out) + "\n")
 
 
def integrate_band_flux(e_kev, nfn, e_lo_kev: float, e_hi_kev: float) -> float:
	"""Integrate nu*F_nu over [e_lo_kev, e_hi_kev] to get the total flux
	(erg/s/cm2) in that band. nFn = nu*F_nu is a flux PER LOG ENERGY, i.e.
	total flux = integral of F_nu dnu = integral of nFn dlnE -- equivalently,
	converting to F_E = nFn/E first and integrating F_E dE linearly (used
	here) gives the same result and is simpler to trapezoid safely."""
	pairs = sorted(zip(e_kev, nfn))
	band = [(e, val) for e, val in pairs if e_lo_kev <= e <= e_hi_kev]
	if len(band) < 2:
		raise ValueError(
			f"Fewer than 2 sim.inc points fall within [{e_lo_kev}, {e_hi_kev}] keV -- "
			f"cannot integrate the band flux. Check sim.inc's energy coverage."
		)
	total = 0.0
	for (e1, v1), (e2, v2) in zip(band[:-1], band[1:]):
		f1, f2 = v1 / e1, v2 / e2  # F_E = nFn / E
		total += 0.5 * (f1 + f2) * (e2 - e1)
	return total  # erg/s/cm2

# --------------------------------------------------------------------------
# ski file construction
# --------------------------------------------------------------------------

def polar_mesh_points(opening_angle_deg: float) -> str:
	lo = (90.0 - opening_angle_deg) / 180.0
	hi = (90.0 + opening_angle_deg) / 180.0
	return f"0, {lo:.6f}, {hi:.6f}, 1.0"


def build_ski(template_path: str, *, rmin_cm, rmax_cm, luminosity_erg_s: float,
			  opening_angle_deg: float) -> str:
	inner_radius_cm, outer_radius_cm = rmin_cm[0], rmax_cm[-1]

	values = dict(
		integratedLuminosity=f"{luminosity_erg_s:.6e} erg/s",
		minRadius=f"{inner_radius_cm:.10e} cm",
		maxRadius=f"{outer_radius_cm:.10e} cm",
		polarMeshPoints=polar_mesh_points(opening_angle_deg),
	)
	text = Path(template_path).read_text()
	return text.format(**values)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def convert_single(cloudy_root: str = "cloudy_grid", run_name: str = "run_00000",
					template_path: str = "template.ski", ion_list=None) -> dict:
	"""Convert just one run -- useful for testing/inspecting a single case
	before running the whole grid. Reads that run's row from
	grid_manifest.csv for opening_angle/intensity, converts its Cloudy
	output, and writes run_name/skirt/{medium.txt, model.ski}.
 
	Returns a small summary dict (paths written, zone count, radii) so you
	can eyeball the result without having to re-open the files.
	"""
	root = Path(cloudy_root)
	manifest_path = root / "grid_manifest.csv"
	if not manifest_path.exists():
		raise FileNotFoundError(f"{manifest_path} not found -- run 1_cloudy_generator.py first")
 
	with open(manifest_path, newline="") as f:
		manifest_rows = {r["run"]: r for r in csv.DictReader(f)}
	if run_name not in manifest_rows:
		raise KeyError(f"'{run_name}' not found in {manifest_path}")
	row = manifest_rows[run_name]
 
	if ion_list is None:
		ion_list = extract_ion_list_from_template(template_path)
 
	cloudy_dir = root / run_name / "cloudy"
	zones_path = cloudy_dir / "sim.zones"
	ovr_path = cloudy_dir / "sim.ovr"
	species_path = cloudy_dir / "sim.species"
	inc_path = cloudy_dir / "sim.inc"
	for p in (zones_path, ovr_path, species_path, inc_path):
		if not p.exists():
			raise FileNotFoundError(f"{p} not found -- has Cloudy been run for {run_name} yet?")
 
	rmin_cm, rmax_cm = parse_sim_zones(zones_path)
	te_K = parse_sim_ovr(ovr_path)
	species_density = parse_sim_species(species_path)
 
	lengths = {len(rmin_cm), len(te_K), *(len(v) for v in species_density.values())}
	if len(lengths) > 1:
		n = min(lengths)
		print(f"warning {run_name}: zone counts differ across sim.zones/sim.ovr/"
			  f"sim.species ({sorted(lengths)}) -- truncating to {n}")
		rmin_cm, rmax_cm, te_K, n_ref_cm3 = rmin_cm[:n], rmax_cm[:n], te_K[:n], n_ref_cm3[:n]
		species_density = {k: v[:n] for k, v in species_density.items()}
 
	n_ref_cm3 = float(row["hden"])
	opening_angle_deg = float(row["opening_angle"])
 
	skirt_dir = root / run_name / "skirt"
	skirt_dir.mkdir(parents=True, exist_ok=True)
 
	medium_txt = build_medium_txt(
		ion_list=ion_list, rmin_cm=rmin_cm, rmax_cm=rmax_cm, te_K=te_K,
		n_ref_cm3=n_ref_cm3, species_density=species_density,
		opening_angle_deg=opening_angle_deg,
	)
	medium_path = skirt_dir / "medium.txt"
	medium_path.write_text(medium_txt)

	e_kev, nfn = parse_incident_continuum(inc_path)
	sed_txt_path = skirt_dir / "sed.txt"
	build_sed_txt_from_continuum(e_kev, nfn, sed_txt_path)
	
	band_flux = integrate_band_flux(e_kev, nfn, 0.3, 10.0)  # erg/s/cm2
	inner_radius_cm = rmin_cm[0]
	luminosity_erg_s = band_flux * 4.0 * math.pi * inner_radius_cm ** 2

	ski_text = build_ski(
		template_path, rmin_cm=rmin_cm, rmax_cm=rmax_cm,
		luminosity_erg_s=luminosity_erg_s, opening_angle_deg=opening_angle_deg,
	)
	ski_path = skirt_dir / "model.ski"
	ski_path.write_text(ski_text)
 
	summary = dict(
		run=run_name,
		n_zones=len(rmin_cm),
		n_ions=len(ion_list),
		inner_radius_pc=rmin_cm[0],
		outer_radius_pc=rmax_cm[-1],
		opening_angle_deg=opening_angle_deg,
		band_flux_erg_s_cm2=band_flux,
		luminosity_erg_s=luminosity_erg_s,
		medium_txt=str(medium_path),
		model_ski=str(ski_path),
		sed_txt=str(sed_txt_path),
	)
	return summary
 
 
def convert_all(cloudy_root: str = "cloudy_grid", template_path: str = "template.ski"):
	root = Path(cloudy_root)
	manifest_path = root / "grid_manifest.csv"
	if not manifest_path.exists():
		raise FileNotFoundError(f"{manifest_path} not found -- run 1_cloudy_generator.py first")
 
	ion_list = extract_ion_list_from_template(template_path)
 
	with open(manifest_path, newline="") as f:
		manifest_rows = list(csv.DictReader(f))
 
	n_done, n_skipped = 0, 0
	for row in manifest_rows:
		run_name = row["run"]
		cloudy_dir = root / run_name / "cloudy"
		if not all((cloudy_dir / f).exists() for f in ("sim.zones", "sim.ovr", "sim.species")):
			print(f"skip {run_name}: missing sim.zones/sim.ovr/sim.species (not run yet?)")
			n_skipped += 1
			continue
 
		convert_single(cloudy_root, run_name, template_path, ion_list=ion_list)
		n_done += 1
 
	print(f"Converted {n_done} run(s), skipped {n_skipped} (missing Cloudy output).")


if __name__ == "__main__":
	import sys
	# Usage:
	#   python3 3_cloudy_skirt_converter.py [cloudy_grid_dir] [template.ski]
	#   python3 3_cloudy_skirt_converter.py [cloudy_grid_dir] [template.ski] --run run_00000
	root = sys.argv[1] if len(sys.argv) > 1 else "cloudy_grid"
	template = sys.argv[2] if len(sys.argv) > 2 else "template.ski"
	if "--run" in sys.argv:
		run_name = sys.argv[sys.argv.index("--run") + 1]
		summary = convert_single(root, run_name, template)
		for k, v in summary.items():
			print(f"{k}: {v}")
	else:
		convert_all(root, template)
 
