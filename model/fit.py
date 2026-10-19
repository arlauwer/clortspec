"""
fit.py
======
Staged fit of the XMM spectra of one target (pn, MOS1, MOS2).

  python3 fit.py TARGET          plots and fit_stage3.xcm -> runs/TARGET/
  python3 fit.py --help

  stage 1:  constant * TBabs * atable{torus.fits}      3-10 keV observed
  stage 2:  constant * TBabs * atable{combined.fits}   0.3-10 keV,
            torus parameters, photon index and norm frozen at stage 1,
            shell parameters (logxi, logNHs, OpenAngS) free
  stage 3:  as stage 2 with everything free

- constant: cross-normalisation, frozen at 1 for pn, free for MOS1/MOS2.
- TBabs:    Galactic absorption, frozen at the target's NH_GALACTIC_1E22.
- tables:   redshift is a table parameter (frozen at the target's REDSHIFT); norm is
            the intrinsic rest-frame 2-10 keV photon flux in both tables, so
            stage-1 values carry over directly.
- Incl:     frozen at the target's INCLINATION (default 90) in stages 1-2. Every line of sight through
            the uniform torus wedge sees the same N_H, so above 3 keV the
            inclination is unconstrained. Free in stage 3, over the table
            inclinations strictly inside the narrowest torus (> 90 - alpha), so
            the fit never uses or interpolates towards unobscured (type-1) views.

Parameters that are fixed in the combined grid (single-valued axes in the
target's COMBINED_GRID) are not in combined.fits; the script prints them next to
the stage-1 result so you can see whether the grid is centred correctly.
"""

import argparse
from datetime import date
from pathlib import Path
import os
import subprocess
import sys

# config/<name>.py (not config/model.py); config is not imported yet, it must come after xspec
_targets = sorted(p.stem for p in (Path(__file__).resolve().parent / "config").glob("*.py")
				  if p.stem not in ("__init__", "model"))
_parser = argparse.ArgumentParser(
	description="Staged XSPEC fit of a target's XMM spectra (pn, MOS1, MOS2) with the SKIRT tables.",
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog="""\
stages:
  1  torus.fits, 3-10 keV. Always runs; needs no target-specific simulations.
  2  runs/TARGET/combined.fits, 0.3-10 keV, torus and source frozen at stage 1
  3  as stage 2 with everything free -- the result to quote
Stages 2-3 need runs/TARGET/combined.fits (./grid.sh TARGET); without it the
script stops after stage 1 and says so.

output:  runs/TARGET/fit_stage{1,1_full,2,3}.png, fit_stage{1,2,3}.txt, fit_stage3.xcm; parameters on the terminal
         (fit_stage1_full.png: the stage-1 model over 0.3-10 keV, fitted band in white;
          fit_stage1.txt: the stage-1 values to centre COMBINED_GRID on;
          fit_stage{2,3}.txt: every parameter with its 90% interval and range, limits flagged)
needs:   HEASoft (pyXspec; set up automatically from $HEADAS or ../heasoft-6.36) and the
         spectra ../TARGET/spectra/{pn,mos1,mos2}_grp.fits (reduction/xmm_reduce.sh)
targets: """ + ", ".join(_targets) + "  (config/<name>.py)")
_parser.add_argument("target", metavar="TARGET", choices=_targets, help="one of: " + ", ".join(_targets))
_args = _parser.parse_args()

# HEASoft is set up here when it is not already: pyXspec needs LD_LIBRARY_PATH,
# which must be there when Python starts, so the environment is taken from
# headas-init.sh and the script starts again. The location is $HEADAS, default
# heasoft-6.36/ next to this directory (as in reduction/xmm_reduce.sh).
try:
	import xspec
except ImportError:
	if os.environ.get("_FIT_HEADAS_DONE"):
		raise
	_headas = os.environ.get("HEADAS") or str(Path(__file__).resolve().parents[1] / "heasoft-6.36/x86_64-pc-linux-gnu-libc2.43")
	if not Path(_headas, "headas-init.sh").is_file():
		sys.exit(f"error: no HEASoft in {_headas} (set HEADAS)")
	_env = subprocess.run(["bash", "-c", 'source "$HEADAS/headas-init.sh" >/dev/null 2>&1; env -0'],
						  env={**os.environ, "HEADAS": _headas}, capture_output=True, check=True).stdout.decode()
	_env = dict(kv.split("=", 1) for kv in _env.split("\0") if "=" in kv)
	_env["_FIT_HEADAS_DONE"] = "1"
	os.execve(sys.executable, [sys.executable] + sys.argv, _env)

# xspec must be imported before matplotlib.pyplot / astropy: the other order
# segfaults with this HEASoft build.
from xspec import AllData, AllModels, Fit, Model, Plot, Xset
from astropy.io import fits
import matplotlib.pyplot as plt

import config as C

T = C.load_target(_args.target)

TORUS_TABLE = C.TORUS_TABLE
COMBINED_TABLE = T.combined_table
INCLINATION = T.inclination
SPECTRA = ["pn_grp.fits", "mos1_grp.fits", "mos2_grp.fits"]
GROUP_LABELS = ["pn", "MOS1", "MOS2"]
COLORS = ["black", "tab:blue", "tab:red"]

TORUS_PARS = ["logNH", "OpenAngle"]

if not T.spectra_dir.is_dir():
	sys.exit(f"{T.spectra_dir} not found: run reduction/xmm_reduce.sh {T.name} first")
os.chdir(T.spectra_dir)
T.dir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 14})

missing = [s for s in SPECTRA if not Path(s).exists()]
if missing:
	sys.exit(f"missing {missing} in {T.spectra_dir}: run reduction/xmm_reduce.sh {T.name} first")

if not C.TORUS_TABLE.exists():
	sys.exit(f"{C.TORUS_TABLE} not found: build the shared torus table with ./grid.sh torus")

AllData.clear()
AllData(" ".join(f"{i}:{i} {s}" for i, s in enumerate(SPECTRA, start=1)))
n_groups = AllData.nGroups

Fit.statMethod = "cstat"
Fit.nIterations = 10000
Fit.query = "yes"

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def new_model(table: Path):
	"""constant * TBabs * atable, with the common parameters set up; returns the table component."""
	AllModels.clear()
	m = Model(f"constant*TBabs*atable{{{table}}}")
	m(1).values = "1 -1"
	m(2).values = f"{T.nh_galactic_1e22} -1"
	for g in range(2, n_groups + 1):
		p = AllModels(g)(1)
		p.untie()
		p.values = "1 0.01 0.5 0.5 2 2"
		p.frozen = False
	comp = getattr(m, m.componentNames[2])
	comp.z.values = f"{T.redshift} -1"
	comp.Incl.values = INCLINATION
	comp.Incl.frozen = True
	return comp


def notice(lo):
	AllData.notice("all")
	AllData.ignore(f"**-{lo} {C.FIT_E_MAX_KEV}-**")
	AllData.ignore("bad")


def report(title):
	print(f"\n=== {title} ===")
	print(f"C-stat = {Fit.statistic:.1f} for {Fit.dof} d.o.f. (ratio {Fit.statistic / Fit.dof:.3f})")
	m1 = AllModels(1)
	for j in range(1, m1.nParameters + 1):
		p = m1(j)
		state = "frozen" if p.frozen else ("linked" if p.link else "free")
		edge = ""
		if state == "free":
			val, lo, hi = p.values[0], p.values[2], p.values[5]
			if min(abs(val - lo), abs(val - hi)) < 1e-6 * max(1.0, abs(val)):
				edge = "   <-- AT GRID/LIMIT EDGE"
		print(f"  {j:2d} {p.name:>10} = {p.values[0]:.4g} ({state}){edge}")
	for g in range(2, n_groups + 1):
		print(f"     {GROUP_LABELS[g - 1]:>10} constant = {AllModels(g)(1).values[0]:.3f}")


def plot(title, filename, fitted=None):
	"""Data, model and residuals over the noticed channels; `fitted` = (lo, hi) keV shades the rest."""
	Plot.device = "/null"
	Plot.xAxis = "keV"
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 8),
								   gridspec_kw={"height_ratios": [3, 1]})
	fig.suptitle(title)
	for g in range(1, n_groups + 1):
		Plot("ldata")
		x, xerr = Plot.x(g), Plot.xErr(g)
		ax1.errorbar(x, Plot.y(g), xerr=xerr, yerr=Plot.yErr(g), fmt=".", alpha=0.4,
					 color=COLORS[g - 1], label=f"{GROUP_LABELS[g - 1]} data")
		ax1.step(x, Plot.model(g), where="mid", color=COLORS[g - 1])
		Plot("delchi")
		ax2.errorbar(Plot.x(g), Plot.y(g), xerr=Plot.xErr(g), yerr=Plot.yErr(g), fmt=".",
					 alpha=0.4, color=COLORS[g - 1])
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_ylabel("Counts/s/keV")
	if fitted:
		for ax in (ax1, ax2):
			lo, hi = ax.get_xlim()
			ax.axvspan(lo, fitted[0], color="0.9", zorder=0)
			ax.axvspan(fitted[1], hi, color="0.9", zorder=0)
			ax.set_xlim(lo, hi)
	ax1.legend()
	ax2.axhline(0, color="k", lw=0.5)
	ax2.set_xlabel("Energy (keV)")
	ax2.set_ylabel(r"$\Delta\chi$")
	plt.tight_layout()
	out = T.dir / filename
	plt.savefig(out, dpi=150)
	print(f"saved {out}")

PAR_DESCRIPTIONS = {
	"logNH": "torus: log equatorial N_H (cm^-2)",
	"OpenAngle": "torus: half-opening angle from the equator (deg)",
	"logxi": "shell: log ionisation parameter at inner face (erg cm/s)",
	"logNHs": "shell: log radial N_H (cm^-2)",
	"OpenAngS": "shell: cone angle from the pole (deg)",
	"PhoIndex": "source: photon index",
	"Incl": "observer inclination from the pole (deg)",
	"z": "redshift",
	"norm": "intrinsic rest-frame 2-10 keV photon flux (ph/cm2/s)",
}


def fit_summary(comp, title, filename, fixed_axes=None):
	"""Readable overview of a fit: every table parameter with its 90% interval
	and allowed range, parameters at a limit flagged; written to runs/TARGET/filename."""
	free = [getattr(comp, n) for n in comp.parameterNames if not getattr(comp, n).frozen]
	if free:
		Fit.error("2.706 " + " ".join(str(p.index) for p in free))
	lines = [f"{T.title} {title}   ({date.today()})",
			 f"C-stat = {Fit.statistic:.1f} for {Fit.dof} d.o.f. (ratio {Fit.statistic / Fit.dof:.3f})",
			 "cross-normalisation (pn = 1): " + ", ".join(f"{GROUP_LABELS[g - 1]} {AllModels(g)(1).values[0]:.3f}"
														 for g in range(2, n_groups + 1)),
			 "",
			 f"{'parameter':<10}{'value':>10}   {'90% interval':<19}{'allowed range':<17}{'':<8}meaning"]
	for name in comp.parameterNames:
		p = getattr(comp, name)
		val, lo, hi = p.values[0], p.values[2], p.values[5]
		if p.frozen:
			state, interval = "frozen", ""
		else:
			state = "free"
			interval = interval_90(p)
		at_edge = not p.frozen and min(abs(val - lo), abs(val - hi)) < 1e-6 * max(1.0, abs(val))
		rng = "" if name in ("z", "norm") else f"{lo:g} - {hi:g}"
		lines.append(f"{name:<10}{val:>10.4g}   {interval:<19}{rng:<17}{state:<8}{PAR_DESCRIPTIONS.get(name, '')}"
					 + ("   <-- AT LIMIT" if at_edge else ""))
	if fixed_axes:
		lines += ["", "fixed in the combined grid (single-valued axis, not in the table):"]
		lines += [f"{name:<10}{value:>10.4g}   {PAR_DESCRIPTIONS.get(name, '')}" for name, value in fixed_axes.items()]
	out = T.dir / filename
	out.write_text("\n".join(lines) + "\n")
	print("\n" + "\n".join(lines) + f"\nsaved {out}")


def interval_90(par):
	"""90% interval from the last Fit.error; '<lo' / '>hi' where it runs into the limit."""
	lo, hi = par.values[2], par.values[5]
	# error = (lower, upper, flags); flags 4 and 5: ran into the lower / upper limit
	err_lo, err_hi, flags = par.error
	return ((f"<{lo:g}" if flags[3] == "T" else f"{err_lo:.4g}" if err_lo else "?") + " - "
			+ (f">{hi:g}" if flags[4] == "T" else f"{err_hi:.4g}" if err_hi else "?"))


def stage1_summary(comp, filename):
	"""Stage-1 result per COMBINED_GRID axis, next to what config/<name>.py has now."""
	axes = {fits_name: axis for axis, fits_name in C.FITS_NAMES.items()}
	lines = [f"{T.title} stage 1: torus only, 3-{C.FIT_E_MAX_KEV:g} keV, Incl = {INCLINATION:g}   ({date.today()})",
			 f"C-stat = {Fit.statistic:.1f} for {Fit.dof} d.o.f.",
			 "cross-normalisation: " + ", ".join(f"{GROUP_LABELS[g - 1]} {AllModels(g)(1).values[0]:.3f}"
												 for g in range(2, n_groups + 1)),
			 f"norm = {comp.norm.values[0]:.4g} photons/cm2/s (intrinsic rest-frame 2-10 keV)",
			 "",
			 f"{'COMBINED_GRID axis':<21}{'best fit':>9}   {'90% interval':<17}{'torus.fits':<15}config/{T.name}.py now"]
	for name in TORUS_PARS + ["PhoIndex"]:
		par = getattr(comp, name)
		val, lo, hi = par.values[0], par.values[2], par.values[5]
		interval = interval_90(par)
		now = T.combined_grid[axes[name]]
		notes = []
		if min(abs(val - lo), abs(val - hi)) < 1e-6 * max(1.0, abs(val)):
			notes.append("AT EDGE OF torus.fits")
		if not min(now) <= val <= max(now):
			notes.append("fixed at another value" if len(now) == 1 else "OUTSIDE current axis")
		lines.append(f"{axes[name]:<21}{val:>9.4g}   {interval:<17}{f'{lo:g} - {hi:g}':<15}{now}"
					 + ("   <-- " + ", ".join(notes) if notes else ""))
	lines += ["",
			  f"Centre the three axes of COMBINED_GRID in config/{T.name}.py on the best fit (a single",
			  f"value fixes the axis), then run ./grid.sh {T.name}."]
	out = T.dir / filename
	out.write_text("\n".join(lines) + "\n")
	print("\n" + "\n".join(lines) + f"\nsaved {out}")

# --------------------------------------------------------------------------
# Stage 1: torus only, 3-10 keV
# --------------------------------------------------------------------------

torus = new_model(TORUS_TABLE)
notice(3.0)
Fit.perform()
# 90% intervals, for choosing how wide the combined grid has to be
Fit.error("2.706 " + " ".join(str(getattr(torus, name).index) for name in TORUS_PARS + ["PhoIndex"]))
report("Stage 1: torus, 3-10 keV")
stage1_summary(torus, "fit_stage1.txt")
plot(f"{T.title} stage 1: neutral torus (3-10 keV)", "fit_stage1.png")
# the same model over the whole band, to see what the torus alone leaves below 3 keV
notice(C.E_MIN_KEV)
plot(f"{T.title} stage 1: neutral torus, fitted at 3-10 keV (white)", "fit_stage1_full.png",
	 fitted=(3.0, C.FIT_E_MAX_KEV))
notice(3.0)

stage1 = {name: getattr(torus, name).values[0] for name in TORUS_PARS + ["PhoIndex", "norm"]}
cross_norm = [AllModels(g)(1).values[0] for g in range(2, n_groups + 1)]

if not COMBINED_TABLE.exists():
	sys.exit(f"\n{COMBINED_TABLE} not built yet. Centre COMBINED_GRID's torus axes in config/{T.name}.py "
			 f"on the stage-1 result above, then run ./grid.sh {T.name}.")

# --------------------------------------------------------------------------
# Stage 2: torus + shell, 0.3-10 keV, torus / source frozen at stage 1
# --------------------------------------------------------------------------

comb = new_model(COMBINED_TABLE)
for g, c in zip(range(2, n_groups + 1), cross_norm):
	AllModels(g)(1).values = c

# axes fixed in the combined grid are recorded as "fixed: <axis> = <value>" comments
fixed_axes = {}
for card in fits.getheader(COMBINED_TABLE).get("COMMENT", []):
	if str(card).startswith("fixed: "):
		axis, value = str(card)[len("fixed: "):].split(" = ")
		fixed_axes[C.FITS_NAMES[axis]] = float(value)
for name, value in stage1.items():
	if name in comb.parameterNames:
		p = getattr(comb, name)
		lo, hi = p.values[2], p.values[5]
		if not lo <= value <= hi:
			print(f"WARNING: stage-1 {name} = {value:.3g} is outside the combined grid [{lo:g}, {hi:g}]; clipping")
		p.values = min(max(value, lo), hi)
		p.frozen = True
	else:
		print(f"note: {name} fixed at {fixed_axes.get(name, float('nan')):g} in the combined grid "
			  f"(stage 1 gave {value:.3g})")

notice(0.3)
Fit.perform()
report("Stage 2: torus + shell, 0.3-10 keV, torus and source frozen")
fit_summary(comb, "stage 2: torus (frozen) + shell, 0.3-10 keV", "fit_stage2.txt", fixed_axes)
plot(f"{T.title} stage 2: torus (frozen) + ionised shell", "fit_stage2.png")

# --------------------------------------------------------------------------
# Stage 3: everything free (except z and Galactic N_H), Incl within the torus
# --------------------------------------------------------------------------

for name in stage1:
	if name in comb.parameterNames:
		getattr(comb, name).frozen = False

# Incl only over table inclinations strictly inside the NARROWEST torus of the
# table, so XSPEC never interpolates between an obscured and an unobscured view.
# New tables record only such inclinations (config.combined_inclinations);
# older ones also have type-1 views, which this excludes.
alpha_min = comb.OpenAngle.values[2] if "OpenAngle" in comb.parameterNames else fixed_axes["OpenAngle"]
table_incls = [float(v) for v in fits.getdata(COMBINED_TABLE, "PARAMETERS")[-1]["VALUE"]]
incl_min = min(i for i in table_incls if i > 90 - alpha_min)
comb.Incl.values = [max(INCLINATION, incl_min), 1.0, incl_min, incl_min, 90.0, 90.0]
comb.Incl.frozen = False
print(f"note: Incl free in stage 3, limited to [{incl_min:g}, 90] (narrowest torus edge at {90 - alpha_min:g} deg)")
Fit.perform()
report("Stage 3: all free, 0.3-10 keV")
fit_summary(comb, f"stage 3: torus + shell, all free, 0.3-10 keV (Incl >= {incl_min:g})", "fit_stage3.txt",
			fixed_axes)
plot(f"{T.title} stage 3: torus + ionised shell, all free", "fit_stage3.png")

xcm = T.dir / "fit_stage3.xcm"
xcm.unlink(missing_ok=True)  # XSPEC refuses to overwrite
Xset.save(str(xcm), info="m")
print(f"saved model to {xcm}")
plt.show()
