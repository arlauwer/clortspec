"""
cloudy_generator.py
===================
Write one CLOUDY input for every distinct shell
configuration in the target's COMBINED_GRID, i.e. every combination of the
axes in config.SHELL_AXES (xi, shell N_H, shell opening angle, Gamma). The
torus parameters do not enter: the torus is neutral (no CLOUDY) and does not
shadow the cone walls, so one CLOUDY run is shared by all torus variants.

  python3 steps/cloudy_generator.py TARGET        (or: ./grid.sh TARGET --only cloudy-gen)

  runs/TARGET/shell_cloudy/run_NNNNN/{sim.in, sed.in}
  runs/TARGET/shell_cloudy/grid_manifest.csv

WHAT CLOUDY IS ASKED TO DO
--------------------------
1D spherical, open geometry, from R_in to R_out (fixed for every run):
  radius log(R_in) log(R_out)          -> CLOUDY stops at R_out
  hden log(n_in) power -2              -> n(r) = n_in (r/R_in)^-2
  xi log(xi)                           -> xi = L(1-1000 Ryd)/(n_in R_in^2)
  stop temperature off                 -> do NOT stop at the 4000 K default

n_in is chosen so that the column R_in -> R_out equals the grid N_H exactly:
  N_H = n_in R_in (1 - R_in/R_out)
so geometry and column are controlled by the grid, and xi only changes the
ionisation state. (Previously CLOUDY stopped at the ionisation front, which
made N_H a hidden function of xi and n.)

Because n ~ r^-2, the local ionisation parameter L/(n r^2) is constant with
radius apart from attenuation, so the grid xi describes the whole shell.

The covering factor of the two hollow cone walls, cos(theta_in) -
cos(theta_out), is passed to CLOUDY so its diffuse field matches the geometry.
Every radial ray through the wall is identical, so this one 1D run describes
the whole shell.
"""

from __future__ import annotations

import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # config.py is one level up
import config as C

SIM_IN_TEMPLATE = """\
title shell {run_name}
iterate to convergence
covering factor {covering_factor:.6f}
radius {log_r_in:.6f} {log_r_out:.6f}
#
### Radiation field ###
table SED "sed.in"
xi {logxi:.4f}
#
### Medium ###
hden {log_n_in:.6f} power {density_power:.3f}
abundances {abundance_set}
no molecules
#
### Stopping: outer radius only ###
stop temperature off
set nend 2000
#
### Output ###
save overview "sim.ovr" last
save species densities *temp all "sim.species" last
save radius "sim.zones" last
save incident continuum "sim.inc" last units keV
"""


def cloudy_sed(photon_index: float, ryd_min: float = 1e-4, ryd_max: float = 1e5, n_points: int = 400) -> str:
	"""CLOUDY `table SED` file: (Ryd, F_nu) pairs, F_nu keyword on the first line."""
	ryd = np.logspace(math.log10(ryd_min), math.log10(ryd_max), n_points)
	e_kev = ryd * C.RYD_KEV
	fnu = C.intrinsic_nufnu(e_kev, photon_index) / e_kev
	lines = [f"{ryd[0]:.8e} {fnu[0]:.8e} Fnu"]
	lines += [f"{r:.8e} {f:.8e}" for r, f in zip(ryd[1:], fnu[1:])]
	return "\n".join(lines) + "\n"


def generate_grid(target: C.Target, output_root=None):
	output_root = output_root or target.cloudy_dir
	output_root.mkdir(parents=True, exist_ok=True)
	axis_names = C.SHELL_AXES
	axis_values = [target.combined_grid[name] for name in axis_names]

	manifest_rows, n_reset = [], 0
	for i, combo in enumerate(itertools.product(*axis_values)):
		params = dict(zip(axis_names, combo))
		run_name = f"run_{i:05d}"
		run_dir = output_root / run_name

		cf = C.shell_covering_factor(params["shell_opening_angle"])
		n_in = C.shell_inner_density(params["shell_logNH"])

		sim_in = SIM_IN_TEMPLATE.format(
			run_name=run_name,
			covering_factor=cf,
			log_r_in=math.log10(C.SHELL_R_IN_CM),
			log_r_out=math.log10(C.SHELL_R_OUT_CM),
			logxi=params["logxi"],
			log_n_in=math.log10(n_in),
			density_power=C.SHELL_DENSITY_POWER,
			abundance_set=C.ABUNDANCE_SET,
		)
		n_reset += C.sync_run_dir(run_dir, {"sim.in": sim_in, "sed.in": cloudy_sed(params["photon_index"])})

		manifest_rows.append({
			"run": run_name, **params,
			"covering_factor": cf,
			"n_in_cm3": n_in,
			"r_in_cm": C.SHELL_R_IN_CM,
			"r_out_cm": C.SHELL_R_OUT_CM,
		})

	manifest_path = output_root / "grid_manifest.csv"
	with open(manifest_path, "w", newline="") as f:
		writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
		writer.writeheader()
		writer.writerows(manifest_rows)

	print(f"Wrote {len(manifest_rows)} CLOUDY runs under {output_root}/")
	print(f"Manifest: {manifest_path}")
	C.report_sync(output_root, len(manifest_rows), n_reset)
	return manifest_rows


if __name__ == "__main__":
	if len(sys.argv) != 2:
		sys.exit(f"usage: python3 steps/cloudy_generator.py TARGET   (targets: {', '.join(C.target_names())})")
	generate_grid(C.load_target(sys.argv[1]))
