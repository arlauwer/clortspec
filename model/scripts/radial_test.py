"""
radial_test.py
==============
How many radial shell cells (SHELL_N_RADIAL) are needed, and how should they
be spaced (SHELL_RADIAL_SPACING)? One shell configuration is simulated with
different radial binnings and the spectra compared:

  run_col30       30 equal-column cells
  run_col40       40 equal-column cells
  run_col60       60 equal-column cells: the reference
  run_n60         60 log-spaced cells (the old default)
  run_n60_seed1   as run_n60 with another random seed: the Monte Carlo noise,
                  i.e. how different two equally good runs are
  run_n30         30 log-spaced cells (first test, failed: ~15-20% less soft
                  flux); only compared if it exists, never generated again

The first test (log spacing, 2026-10) showed n30 differing from n60 by 22-28%
rms in 0.3-2 keV against 3-5% seed noise: log spacing puts only ~5 of 60 cells
in the thin ionised layer (front at ~1.6 R_in at logNHs 24), so equal-column
spacing is now the default and this test checks how many cells it needs.

Default configuration: logNHs 24 and logxi 2, where the ionisation front lies
inside the shell, so the radial structure matters most (for 3C445: ~44% of the
hydrogen neutral). Override with --logxi/--logNHs/--angle/--gamma.

On the server, from model/:

  python3 scripts/radial_test.py cloudy 3C234        # CLOUDY input -> runs/3C234/radial_test/cloudy/
  steps/cloudy_run.sh runs/3C234/radial_test/cloudy  # ~10 min
  python3 scripts/radial_test.py skirt 3C234         # the SKIRT models -> runs/3C234/radial_test/skirt/
  steps/skirt_run.sh runs/3C234/radial_test/skirt -j 2 -t THREADS
  python3 scripts/radial_test.py compare 3C234       # differences per inclination (+ compare.png)

`skirt` only writes the variants whose directory does not exist yet, and
skirt_run.sh skips runs that already have their spectra, so after the first
test only the run_col* models are new (the cloudy step can be skipped too).

`compare` works on the server (numpy only; the plot needs matplotlib) or after
copying runs/TARGET/radial_test/ back. N equal-column cells are enough if
colN differs from col60 by no more than n60_seed1 differs from n60, in
particular in the soft band; log60 vs col60 shows how far off the old binning was.
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

# name: (SHELL_N_RADIAL, SHELL_RADIAL_SPACING, random seed or None)
VARIANTS = {
	"run_col30": (30, "column", None),
	"run_col40": (40, "column", None),
	"run_col60": (60, "column", None),
	"run_n60": (60, "log", None),
	"run_n60_seed1": (60, "log", 1),
	"run_n30": (30, "log", None),
}
NOT_GENERATED = {"run_n30"}   # already shown to fail, only compared if present
REFERENCE = "run_col60"
# label: (run, compared with); "noise" is the Monte Carlo yardstick
COMPARISONS = {
	"col30": ("run_col30", REFERENCE),
	"col40": ("run_col40", REFERENCE),
	"log60": ("run_n60", REFERENCE),
	"log30": ("run_n30", REFERENCE),
	"noise": ("run_n60_seed1", "run_n60"),
}
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
	n_default, spacing_default = C.SHELL_N_RADIAL, C.SHELL_RADIAL_SPACING
	rows = []
	try:
		for name, (n_radial, spacing, seed) in VARIANTS.items():
			run = skirt_root / name
			if name in NOT_GENERATED:
				continue
			if run.exists():
				print(f"{run} exists, left alone (delete it to regenerate)")
				continue
			# skirt_generator reads these from the config package
			C.SHELL_N_RADIAL, C.SHELL_RADIAL_SPACING = n_radial, spacing
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
			rows.append({"run": name, "shell_n_radial": n_radial, "spacing": spacing,
						 "seed": seed if seed is not None else 0,
						 **{k: v[0] for k, v in T.combined_grid.items()}})
	finally:
		C.SHELL_N_RADIAL, C.SHELL_RADIAL_SPACING = n_default, spacing_default
	if rows:
		path = skirt_root / "radial_test.csv"
		old = []
		if path.exists():
			with open(path, newline="") as f:
				old = list(csv.DictReader(f))
			for r in old:
				r.setdefault("spacing", "log")   # rows from before the spacing option
		fields = list(rows[0]) + [k for r in old for k in r if k not in rows[0]]
		with open(path, "w", newline="") as f:
			writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(fields)), restval="")
			writer.writeheader()
			writer.writerows(old + rows)
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
	if REFERENCE not in have:
		sys.exit(f"need the finished reference {REFERENCE} in {root}")
	for name, run in have.items():
		print(f"{name}: " + ("; ".join(skirt_log_summary(run)) or "(no log summary found)"))
	comps = {label: pair for label, pair in COMPARISONS.items() if pair[0] in have and pair[1] in have}
	print(f"\nrelative differences of the total flux in {REBIN_DEX}-dex bins: rms / max |ratio - 1| per band"
		  f"\n" + "; ".join(f"{label} = {a}/{b}" for label, (a, b) in comps.items())
		  + "\n(colN vs col60 should be no larger than noise, which is pure Monte Carlo noise)")
	cols = [(band, label) for band in BANDS for label in comps]
	width = 16
	print(f"{'':>5}" + "".join(f"{band:>{width * len(comps)}}" for band in BANDS))
	print(f"{'incl':>5}" + "".join(f"{label:>{width}}" for _, label in cols))
	sed = lambda name, inc: rebinned(have[name] / f"model_i{inc:02d}_sed.dat")
	curves = {}
	for inc in T.inclinations:
		row = f"{inc:>5}"
		for band, label in cols:
			a, b = comps[label]
			e, ref = sed(b, inc)
			_, f = sed(a, inc)
			lo, hi = BANDS[band]
			sel = (e >= lo) & (e <= hi) & (ref > 0)
			d = f[sel] / ref[sel] - 1
			curves[(inc, label)] = (e, f / np.where(ref > 0, ref, np.nan))
			row += f"{f'{100 * np.sqrt(np.mean(d ** 2)):.1f}% / {100 * np.max(np.abs(d)):.1f}%':>{width}}"
		print(row)
	for band, (lo, hi) in BANDS.items():
		parts = []
		for label, (a, b) in comps.items():
			e, ref = sed(b, 90)
			_, f = sed(a, 90)
			sel = (e >= lo) & (e <= hi)
			parts.append(f"{label} {f[sel].sum() / ref[sel].sum():.4f}")
		print(f"band flux ratio at Incl 90, {band}: " + ", ".join(parts))
	try:
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
	except ImportError:
		print("(no matplotlib: no plot)")
		return
	fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 8), gridspec_kw={"height_ratios": [2, 1]})
	for name in have:
		e, f = sed(name, 90)
		axes[0].loglog(e, e * f, "-" if name == REFERENCE else "--", lw=1.5 if name == REFERENCE else 1,
					   label=f"{name} (Incl 90)")
	first = T.inclinations[0]
	for k, label in enumerate(comps):
		colour = f"C{k}"
		for inc, style in ((90, "-"), (first, ":")):
			if (inc, label) in curves:
				e, r = curves[(inc, label)]
				axes[1].semilogx(e, r, style, color=colour, lw=1,
								 label=f"{label}, Incl {inc}" if inc != first or first == 90 else None)
	axes[0].set_ylabel("E F_E (arbitrary)")
	axes[0].legend(fontsize=8)
	axes[1].axhline(1, color="k", lw=0.5)
	axes[1].set_ylim(0.8, 1.2)
	axes[1].set_xlabel("rest-frame energy (keV)")
	axes[1].set_ylabel(f"ratio (solid Incl 90, dotted Incl {first})")
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
