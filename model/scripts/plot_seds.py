#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SED_COLUMNS = {
	"total": 2,
	"transparent": 3,
	"direct_primary": 4,
	"scattered_primary": 5,
	"direct_secondary": 6,
	"scattered_secondary": 7,
	"transparent_secondary": 8,
}


def read_sed(path):
	data = np.loadtxt(path, comments="#")
	energy = data[:, 0]
	return energy, data


def make_label(row):
	return (
		f"$\\xi={row.xi:g}$,\t"
		f"$n_H={row.hden:.1e}$,\t"
		# f"$Z={row.metallicity:g}$, "
		# f"$\\theta={row.opening_angle:g}^\\circ$, "
		f"$\\Gamma={row.photon_index:g}$,\t"
		# f"$C_f={row.covering_factor:g}$"
	)


def main():

	import argparse

	parser = argparse.ArgumentParser(
		description="Plot SEDs from a cloudy_grid."
	)

	parser.add_argument(
		"--grid",
		default="cloudy_grid",
		help="Grid directory (default: cloudy_grid)",
	)

	parser.add_argument(
		"--runs",
		nargs="+",
		help="Runs to plot, e.g. --runs 0 1 2",
	)

	parser.add_argument(
		"--inclination",
		choices=["i00", "i90", "both"],
		default="both",
		help="Inclination to plot (default: both)",
	)

	parser.add_argument(
		"--component",
		choices=SED_COLUMNS.keys(),
		default="total",
		help="SED component to plot (default: total)",
	)

	parser.add_argument(
		"--separate",
		action="store_true",
		help="Make separate panels for i00 and i90.",
	)

	args = parser.parse_args()

	grid = Path(args.grid)
	manifest_path = grid / "grid_manifest.csv"

	if not manifest_path.exists():
		raise FileNotFoundError(
			f"Could not find {manifest_path}"
		)

	manifest = pd.read_csv(manifest_path)

	# ---------------------------------------------------------
	# If no runs were specified, print available runs
	# ---------------------------------------------------------

	if not args.runs:
		print("\nAvailable runs:\n")

		for _, row in manifest.iterrows():
			print(
				f"{row.run:12s} "
				f"xi={row.xi:g} "
				f"hden={row.hden:g} "
				f"Z={row.metallicity:g} "
				f"OA={row.opening_angle:g} "
				f"Gamma={row.photon_index:g} "
				f"CF={row.covering_factor:g}"
			)

		print(
			"\nExamples:\n"
			"  python plot_seds.py --runs 0 1 2\n"
			"  python plot_seds.py --runs 0 1 2 --component direct_primary\n"
			"  python plot_seds.py --runs 0 1 2 --inclination i00\n"
			"  python plot_seds.py --runs 0 1 2 --separate\n"
		)

		return

	# ---------------------------------------------------------
	# Normalize run names
	#
	# Supports:
	#   --runs 0 1 2
	#   --runs 0-10
	#   --runs 0-5 10 15-20
	#   --runs run_00000 run_00001
	# ---------------------------------------------------------

	run_names = []

	for run in args.runs:

		# Handle ranges such as 0-10
		if "-" in run and not run.startswith("run_"):

			start, end = run.split("-", 1)

			start = int(start)
			end = int(end)

			if end < start:
				raise ValueError(
					f"Invalid run range: {run}"
				)

			for number in range(start, end + 1):
				run_names.append(
					f"run_{number:05d}"
				)

		# Already a full run name
		elif run.startswith("run_"):
			run_names.append(run)

		# Single integer
		else:
			run_names.append(
				f"run_{int(run):05d}"
			)

	# Remove duplicates while preserving order
	run_names = list(dict.fromkeys(run_names))

	selected = manifest[
		manifest["run"].isin(run_names)
	]

	missing = set(run_names) - set(selected["run"])

	if missing:
		raise ValueError(
			f"These runs were not found: {sorted(missing)}"
		)

	# ---------------------------------------------------------
	# Which inclinations?
	# ---------------------------------------------------------

	if args.inclination == "both":
		inclinations = ["i00", "i90"]

	else:
		inclinations = [args.inclination]

	# ---------------------------------------------------------
	# Create figure
	# ---------------------------------------------------------

	if args.separate and len(inclinations) == 2:

		fig, axes = plt.subplots(
			2,
			1,
			figsize=(10, 10),
			sharex=True,
		)

		axes = list(axes)

	else:

		fig, ax = plt.subplots(
			figsize=(11, 7)
		)

		axes = [ax]

	# ---------------------------------------------------------
	# Colors
	#
	# One color = one run
	# Solid = i00
	# Dashed = i90
	# ---------------------------------------------------------

	colors = plt.cm.Set2(
		np.linspace(
			0,
			1,
			max(len(selected), 1),
		)
	)

	linestyles = {
		"i00": "-",
		"i90": "--",
	}

	# ---------------------------------------------------------
	# Plot
	# ---------------------------------------------------------

	for inc_index, inclination in enumerate(inclinations):

		ax = (
			axes[inc_index]
			if args.separate
			else axes[0]
		)

		for color, (_, row) in zip(
			colors,
			selected.iterrows(),
		):

			run_dir = (
				grid
				/ row["run"]
				/ "skirt"
			)

			sed_file = (
				run_dir
				/ f"model_{inclination}_sed.dat"
			)

			if not sed_file.exists():

				print(
					f"WARNING: missing {sed_file}"
				)

				continue

			energy, data = read_sed(
				sed_file
			)

			column = SED_COLUMNS[
				args.component
			]

			flux = data[:, column - 1]

			good = (
				np.isfinite(energy)
				& np.isfinite(flux)
				& (energy > 0)
				& (flux > 0)
			)

			ax.plot(
				energy[good],
				flux[good],
				color=color,
				linestyle=linestyles[inclination],
				lw=1.8,
			)

		ax.set_xscale("log")
		ax.set_yscale("log")

		ax.set_ylabel(
			r"$F_E$ [1/s/cm$^2$/keV]"
		)

		ax.grid(
			True,
			which="both",
			alpha=0.15,
			linestyle=":",
		)

		if args.separate:

			ax.set_title(
				f"Inclination {inclination}"
			)

	axes[-1].set_xlabel(
		"Energy [keV]"
	)

	# ---------------------------------------------------------
	# Build legend
	#
	# First: fixed inclination/style labels
	# Then: one entry for each run/color
	# ---------------------------------------------------------

	inclination_handles = [

		Line2D(
			[0],
			[0],
			color="black",
			lw=2,
			linestyle="-",
			label="i00",
		),

		Line2D(
			[0],
			[0],
			color="black",
			lw=2,
			linestyle="--",
			label="i90",
		),

	]

	dataset_handles = []

	for color, (_, row) in zip(
		colors,
		selected.iterrows(),
	):

		dataset_handles.append(
			Line2D(
				[0],
				[0],
				color=color,
				lw=2,
				linestyle="-",
				label=make_label(row),
			)
		)

	handles = (
		inclination_handles
		+ dataset_handles
	)

	# ---------------------------------------------------------
	# Legend
	# ---------------------------------------------------------

	axes[0].legend(
		handles=handles,
		bbox_to_anchor=(1.02, 1),
		loc="upper left",
		fontsize=8,
	)

	component_name = (
		args.component
		.replace("_", " ")
	)

	fig.suptitle(
		f"SED comparison — {component_name}",
		fontsize=14,
	)

	fig.tight_layout()

	plt.show()


if __name__ == "__main__":
	main()