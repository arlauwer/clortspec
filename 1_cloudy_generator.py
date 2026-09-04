"""
generate_cloudy_grid.py
========================
Step 1 of the pipeline: takes a parameter space (Cartesian product of a few
axes) and writes out one directory per grid point, each containing a Cloudy
`cloudy/sim.in` + `cloudy/sed.in`, ready to run with Cloudy to get a 1D
spherical ionization structure. A later script will read Cloudy's output
(sim.ovr / sim.species) plus this same parameter space to build the SKIRT
.ski files with spatially-resolved ion abundances.

PARAMETER SPACE
---------------
Edit PARAMETER_SPACE below. Any key with a list value becomes a grid axis
(Cartesian product across all such keys). Any key with a scalar value is
held fixed for every run. This keeps script #1 and #2 in sync: whatever
axes you add here just need matching handling wherever you build the .ski
per run later.

ASSUMPTIONS (flagged so you can override easily)
--------------------------------------------------
- inner_radius_cm is fixed at 1e18 cm (~0.32 pc), a typical AGN torus
  sublimation-radius scale. Change FIXED["inner_radius_cm"] if you want
  something else, or move it into PARAMETER_SPACE to make it a grid axis.
- intensity is defined at the inner face, integrated over 0.3-10 keV
  (matching the XMM band / the ski file's normalization band), in
  erg/s/cm^2, entered as `intensity <value> range <Elo> to <Ehi> linear`.
- The SED *shape* is generated with a simple power law
  F_lambda(lambda) ~ lambda^(photon_index - 3)  (equivalently F_nu ~
  nu^(1-photon_index)), spanning a broad 1e-6-1e4 Ryd grid so Cloudy has
  enough continuum for its thermal balance calculation. This is a
  placeholder -- if you already have real SED shapes (e.g. from a
  disk+corona model) per grid point, swap `powerlaw_sed_points()` for a
  function that returns your own (Ryd, value) points instead.
- covering factor = sin(opening_angle), where opening_angle is the
  half-angle of the torus wedge measured from the equatorial plane
  (0 deg -> CF=0, 90 deg -> CF=1), consistent with the polar-mesh
  convention used for the .ski files.
- abundances table fixed to GASS10 (matches your example); metallicity is
  a grid axis multiplying it via `metals <value> linear`.
"""

from __future__ import annotations
import csv
import itertools
import math
import numpy as np
from pathlib import Path

# --------------------------------------------------------------------------
# Parameter space -- EDIT THIS
# --------------------------------------------------------------------------

PARAMETER_SPACE = {
	# hydrogen ionization parameter (log)
	"xi": [-1, 0, 1, 2],
	# hydrogen density, cm^-3 (linear)
	"hden": [1e4, 1e6],
	# metallicity relative to GASS10 solar
	"metallicity": [1.0],
	# torus half-opening angle from the equatorial plane, degrees
	"opening_angle": [30],
	# SED power-law photon index
	"photon_index": [1, 2],
}

FIXED = {
	"inner_radius_cm": 1e18,   # ~0.32 pc
	"abundance_set": "GASS10",
}

OUTPUT_ROOT = Path("cloudy_grid")

# --------------------------------------------------------------------------
# Unit conversions
# --------------------------------------------------------------------------

def covering_factor(opening_angle_deg: float) -> float:
	"""Torus wedge covering factor from the half-angle above the equator."""
	return math.sin(math.radians(opening_angle_deg))


# --------------------------------------------------------------------------
# SED generation (placeholder power law -- see module docstring)
# --------------------------------------------------------------------------

def powerlaw_sed_points(photon_index: float, Ryd_min: float = 1.0,
						 Ryd_max: float = 1e3, n_points: int = 200):
	Ryd = np.logspace(math.log10(Ryd_min), math.log10(Ryd_max), n_points)
	Fnu = Ryd ** (1-photon_index) # F_nu ~ F_E ~ E * dN/dE (Energy / time / area / energy)

	return (Ryd, Fnu)


def write_sed_file(path: Path, sed) -> None:
	lines = []
	for i, (Ryd, fnu) in enumerate(zip(*sed)):
		if i == 0:
			lines.append(f"{Ryd:.8e} {fnu:.8e} Fnu")
		else:
			lines.append(f"{Ryd:.8e} {fnu:.8e}")
	path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# sim.in builder
# --------------------------------------------------------------------------

SIM_IN_TEMPLATE = """\
### One Zone ###
iterate until convergence
covering factor {covering_factor:.6f}
radius {inner_radius_log:.6f}
#
######### INPUT #########
### Radiation Field ###
table SED "sed.in"
xi {xi:.6e}
### Medium ###
hden {hden:.6e} linear
abundances {abundance_set}
metals {metallicity:.6f} linear
no molecules
#
######### OUTPUT #########
### Temperature ###
save overview "sim.ovr" last
### Abundances ###
save species densities *temp all "sim.species" last
### Zones ###
save radius "sim.zones" last
### SED ###
save incident continuum "sim.inc" last units keV
"""


def build_sim_in(*, covering_factor: float, inner_radius_cm: float, xi: float,
				  hden: float, metallicity: float, abundance_set: str) -> str:
	return SIM_IN_TEMPLATE.format(
		covering_factor=covering_factor,
		inner_radius_log=math.log10(inner_radius_cm),
		xi=xi,
		ryd_lo=0.3,
		ryd_hi=10.,
		hden=hden,
		metallicity=metallicity,
		abundance_set=abundance_set,
	)


# --------------------------------------------------------------------------
# Grid generation
# --------------------------------------------------------------------------

def generate_grid(output_root: Path = OUTPUT_ROOT):
	output_root.mkdir(parents=True, exist_ok=True)
	axis_names = list(PARAMETER_SPACE.keys())
	axis_values = [PARAMETER_SPACE[name] for name in axis_names]

	manifest_rows = []
	for i, combo in enumerate(itertools.product(*axis_values)):
		params = dict(zip(axis_names, combo))
		run_name = f"run_{i:05d}"
		run_dir = output_root / run_name / "cloudy"
		run_dir.mkdir(parents=True, exist_ok=True)

		cf = covering_factor(params["opening_angle"])

		sim_in_text = build_sim_in(
			covering_factor=cf,
			inner_radius_cm=FIXED["inner_radius_cm"],
			xi=params["xi"],
			hden=params["hden"],
			metallicity=params["metallicity"],
			abundance_set=FIXED["abundance_set"],
		)
		(run_dir / "sim.in").write_text(sim_in_text)

		sed = powerlaw_sed_points(params["photon_index"])
		write_sed_file(run_dir / "sed.in", sed)

		row = {"run": run_name, **params, "covering_factor": cf,
			   "inner_radius_cm": FIXED["inner_radius_cm"],
			   "abundance_set": FIXED["abundance_set"]}
		manifest_rows.append(row)

	manifest_path = output_root / "grid_manifest.csv"
	with open(manifest_path, "w", newline="") as f:
		writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
		writer.writeheader()
		writer.writerows(manifest_rows)

	print(f"Wrote {len(manifest_rows)} runs under {output_root}/")
	print(f"Manifest: {manifest_path}")
	return manifest_rows


if __name__ == "__main__":
	generate_grid()