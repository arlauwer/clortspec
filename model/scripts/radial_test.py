"""
radial_test.py
==============
Does the shell need SHELL_N_RADIAL = 60 radial cells, or is 30 enough?
One shell configuration is simulated three times and the spectra compared:

  run_n60         60 cells (config.SHELL_N_RADIAL now)
  run_n30         30 cells, half the memory and cell-building time
  run_n60_seed1   60 cells with another random seed: the Monte Carlo noise,
                  i.e. how different two equally good runs are

Default configuration: logNHs 24 and logxi 2, where the ionisation front lies
inside the shell, so the radial structure matters most (for 3C445: ~44% of the
hydrogen neutral). Override with --logxi/--logNHs/--angle/--gamma.

On the server, from model/:

  python3 scripts/radial_test.py cloudy 3C445        # CLOUDY input -> runs/3C445/radial_test/cloudy/
  steps/cloudy_run.sh runs/3C445/radial_test/cloudy  # ~10 min
  python3 scripts/radial_test.py skirt 3C445         # the three SKIRT models -> runs/3C445/radial_test/skirt/
  steps/skirt_run.sh runs/3C445/radial_test/skirt -j 2 -t THREADS   # (the seed run is optional, see below)
  python3 scripts/radial_test.py compare 3C445       # differences per inclination (+ compare.png)

`compare` works on the server (numpy only; the plot needs matplotlib) or after
copying runs/3C445/radial_test/ back. 30 cells are good enough if n30 differs
from n60 by no more than n60_seed1 does, in particular in the soft band. The
seed run costs as much as a normal run; without it, `compare` still prints the
n30/n60 differences, but there is nothing to judge them against.
"""

import argparse
import csv
import dataclasses
import re
import shutil
import sys
from pathlib import Path

import numpy as np

MODEL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MODEL_DIR))
sys.path.insert(0, str(MODEL_DIR / "steps"))
import config as C
import cloudy_generator
import skirt_generator
from build_xspec import read_model_sed

VARIANTS = {"run_n60": (60, None), "run_n30": (30, None), "run_n60_seed1": (60, 1)}
BANDS = {"0.3-2 keV": (0.3, 2.0), "2-10 keV": (2.0, 10.0)}
REBIN_DEX = 0.05   # compare in 0.05-dex bins (~12% wide), coarser than the noise of single bins


def test_target(args):
	"""The target with a single-point COMBINED_GRID at the test configuration."""
	T = C.load_target(args.target)
	g = dict(T.combined_grid)
	pick = {"logxi": args.logxi, "shell_logNH": args.logNHs,
			"shell_opening_angle": args.angle if args.angle is not None else max(g["shell_opening_angle"]),
			"photon_index": args.gamma if args.gamma is not None else sorted(g["photon_index"])[len(g["photon_index"]) // 2]}
	g.update({k: [v] for k, v in pick.items()})
	g["torus_logNH"] = [g["torus_logNH"][0]]
	g["torus_opening_angle"] = [min(g["torus_opening_angle"])]
	return dataclasses.replace(T, combined_grid=g)


def test_dir(target_name):
	return C.RUNS_DIR / target_name / "radial_test"


def rel(path):
	"""Path relative to model/ (where the commands are run), if it is inside it."""
	return path.relative_to(MODEL_DIR) if path.is_relative_to(MODEL_DIR) else path


def cmd_cloudy(args):
	T = test_target(args)
	print("test configuration:", {k: v[0] for k, v in T.combined_grid.items()})
	cloudy_generator.generate_grid(T, test_dir(T.name) / "cloudy")
	print(f"next: steps/cloudy_run.sh {rel(test_dir(T.name) / 'cloudy')}")


def cmd_skirt(args):
	T = test_target(args)
	root = test_dir(T.name)
	cloudy_root, skirt_root = root / "cloudy", root / "skirt"
	if not (cloudy_root / "run_00000" / "sim.out").exists():
		sys.exit(f"no CLOUDY output in {cloudy_root}/run_00000: run the cloudy step first")
	skirt_root.mkdir(parents=True, exist_ok=True)
	n_default = C.SHELL_N_RADIAL
	rows = []
	try:
		for name, (n_radial, seed) in VARIANTS.items():
			run = skirt_root / name
			if run.exists():
				print(f"{run} exists, left alone (delete it to regenerate)")
				continue
			C.SHELL_N_RADIAL = n_radial      # skirt_generator reads it from the config package
			tmp = skirt_root / f"_tmp_{name}"
			if skirt_generator.generate_combined(T, tmp, cloudy_root) != 0:
				sys.exit("generating the SKIRT model failed (see above)")
			shutil.move(str(tmp / "run_00000"), str(run))
			shutil.rmtree(tmp)
			if seed is not None:
				ski = run / "model.ski"
				text, n = re.subn(r'<Random seed="\d+"/>', f'<Random seed="{seed}"/>', ski.read_text())
				if n != 1:
					sys.exit(f"{ski}: no Random seed to change")
				ski.write_text(text)
			rows.append({"run": name, "shell_n_radial": n_radial, "seed": seed if seed is not None else 0,
						 **{k: v[0] for k, v in T.combined_grid.items()}})
	finally:
		C.SHELL_N_RADIAL = n_default
	if rows:
		with open(skirt_root / "radial_test.csv", "a", newline="") as f:
			writer = csv.DictWriter(f, fieldnames=list(rows[0]))
			if f.tell() == 0:
				writer.writeheader()
			writer.writerows(rows)
	print(f"next: steps/skirt_run.sh {rel(skirt_root)} -j 2 -t THREADS")


def rebinned(path):
	"""Total flux in REBIN_DEX-wide energy bins (bin centres, flux)."""
	e, cols = read_model_sed(path)
	order = np.argsort(e)
	e, f = e[order], cols["total flux"][order]
	edges = 10 ** np.arange(np.log10(e[0]), np.log10(e[-1]) + REBIN_DEX, REBIN_DEX)
	n = len(edges) - 1
	idx = np.clip(np.digitize(e, edges) - 1, 0, n - 1)
	sums = np.bincount(idx, weights=f, minlength=n)
	counts = np.bincount(idx, minlength=n)
	keep = counts > 0
	return np.sqrt(edges[:-1] * edges[1:])[keep], sums[keep] / counts[keep]


def skirt_log_summary(run):
	"""Peak memory and run time lines from the SKIRT log, if present."""
	lines = []
	for log in sorted(run.glob("*_log.txt")):
		for line in log.read_text(errors="ignore").splitlines():
			if re.search(r"Peak memory|Finished simulation|Finished setup", line):
				lines.append(line.split(" ", 2)[-1].strip())
	return lines


def cmd_compare(args):
	root = test_dir(args.target) / "skirt"
	T = C.load_target(args.target)
	have = {name: root / name for name in VARIANTS if list((root / name).glob("*_sed.dat"))}
	if "run_n60" not in have or "run_n30" not in have:
		sys.exit(f"need finished run_n60 and run_n30 in {root}")
	for name, run in have.items():
		print(f"{name}: " + ("; ".join(skirt_log_summary(run)) or "(no log summary found)"))
	print(f"\nrelative differences of the total flux in {REBIN_DEX}-dex bins: rms / max |ratio - 1| per band"
		  f"\n(n30 vs n60 should be no larger than seed1 vs n60, which is pure Monte Carlo noise)")
	cols = [(band, other) for band in BANDS for other in ("run_n30", "run_n60_seed1")]
	print(f"{'incl':>5}" + "".join(f"{band + (' n30' if other == 'run_n30' else ' seed1'):>22}" for band, other in cols))
	curves = {}
	for inc in T.inclinations:
		e, ref = rebinned(have["run_n60"] / f"model_i{inc:02d}_sed.dat")
		row = f"{inc:>5}"
		for band, other in cols:
			if other not in have:
				row += f"{'-':>22}"
				continue
			lo, hi = BANDS[band]
			sel = (e >= lo) & (e <= hi) & (ref > 0)
			_, f = rebinned(have[other] / f"model_i{inc:02d}_sed.dat")
			d = f[sel] / ref[sel] - 1
			curves[(inc, other)] = (e, f / np.where(ref > 0, ref, np.nan))
			row += f"{f'{100 * np.sqrt(np.mean(d ** 2)):.1f}% / {100 * np.max(np.abs(d)):.1f}%':>22}"
		print(row)
	for band, (lo, hi) in BANDS.items():
		e, ref = rebinned(have["run_n60"] / f"model_i90_sed.dat")
		sel = (e >= lo) & (e <= hi)
		parts = []
		for other in ("run_n30", "run_n60_seed1"):
			if other in have:
				_, f = rebinned(have[other] / "model_i90_sed.dat")
				parts.append(f"{other} {f[sel].sum() / ref[sel].sum():.4f}")
		print(f"band flux ratio to n60 at Incl 90, {band}: " + ", ".join(parts))
	try:
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
	except ImportError:
		print("(no matplotlib: no plot)")
		return
	fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 8), gridspec_kw={"height_ratios": [2, 1]})
	for name, style in (("run_n60", "-"), ("run_n30", "--"), ("run_n60_seed1", ":")):
		if name in have:
			e, f = rebinned(have[name] / "model_i90_sed.dat")
			axes[0].loglog(e, e * f, style, label=f"{name} (Incl 90)")
	for (inc, other), (e, r) in curves.items():
		axes[1].semilogx(e, r, "-" if other == "run_n30" else ":", lw=1,
						 label=f"{other} / n60, Incl {inc}" if inc in (T.inclinations[0], 90) else None)
	axes[0].set_ylabel("E F_E (arbitrary)")
	axes[0].legend()
	axes[1].axhline(1, color="k", lw=0.5)
	axes[1].set_ylim(0.8, 1.2)
	axes[1].set_xlabel("rest-frame energy (keV)")
	axes[1].set_ylabel("ratio to n60")
	axes[1].legend(fontsize=8)
	axes[1].set_xlim(C.E_MIN_KEV, C.E_MAX_KEV)
	out = test_dir(args.target) / "compare.png"
	plt.tight_layout()
	plt.savefig(out, dpi=150)
	print(f"saved {out}")


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("step", choices=["cloudy", "skirt", "compare"])
	parser.add_argument("target", choices=C.target_names())
	parser.add_argument("--logxi", type=float, default=2.0)
	parser.add_argument("--logNHs", type=float, default=24.0)
	parser.add_argument("--angle", type=float, help="shell_opening_angle (default: the largest in the grid)")
	parser.add_argument("--gamma", type=float, help="photon_index (default: the middle value in the grid)")
	args = parser.parse_args()
	{"cloudy": cmd_cloudy, "skirt": cmd_skirt, "compare": cmd_compare}[args.step](args)


if __name__ == "__main__":
	main()
