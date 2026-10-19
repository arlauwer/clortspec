"""
skirt_generator.py
==================
Write the SKIRT models.

  python3 steps/skirt_generator.py torus    -> torus_grid/ (stage 1, no CLOUDY,
                                               shared by all targets)
  python3 steps/skirt_generator.py TARGET   -> runs/TARGET/combined_grid/ (torus
                                               + shell, needs the CLOUDY runs)
  (or: ./grid.sh torus|TARGET --only skirt-gen)

Each run gets run_NNNNN/{model.ski, sed.txt[, shell.txt]} and the grid
directory a grid_manifest.csv in itertools.product order of the config grid.
Regenerating is safe: a run whose inputs are unchanged keeps its results, a
run whose inputs changed (because the grid was edited) has its old results
deleted so that the run step redoes it.

TORUS
-----
SKIRT TorusGeometry with p = q = 0: a uniform wedge between TORUS_R_IN and
TORUS_R_OUT within +-alpha of the equator. NumberColumnMaterialNormalization
along X counts the column through the WHOLE model (it crosses the torus
twice), so it is set to 2 N_H; centre -> edge is then exactly N_H.
Abundances: CLOUDY's GASS10, as for the shell.

SHELL (combined grid only)
--------------------------
CLOUDY output (runs/TARGET/shell_cloudy/run_*/sim.{zones,ovr,species,out}) -> shell.txt
for a SphericalCellMedium:
- CLOUDY zones are rebinned onto SHELL_N_RADIAL log-spaced radial shells,
  conserving the column of hydrogen, of every ion, and the n-weighted
  temperature. X-ray absorption along a radial ray depends only on these
  columns, while SKIRT memory scales with the number of imported cells.
- Ions below SHELL_ION_MIN_FRACTION of their element are dropped (the element's
  other ions absorb their share): every non-zero ion enlarges the cell's SKIRT
  mix, and with it the memory (scripts/predict_memory.py estimates it).
- Each radial shell becomes two cells, the northern wall theta_in..theta_out
  and the southern wall 180-theta_out..180-theta_in (deg from the pole). If
  the wall would reach into the torus (theta_out > 90 - alpha), it is clipped
  at the torus edge and the run is flagged in the manifest.
- The spatial grid's radial and polar borders coincide with the shell cells,
  the torus edges, and TORUS_R_OUT, so nothing is resampled.

CLOUDY runs are validated (stop reason "outer radius reached" and the column
within 0.02 dex of the grid value); results go to shell_cloudy/cloudy_check.csv
and combined runs that depend on a failed CLOUDY run are not written.
"""

from __future__ import annotations

import csv
import itertools
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # config.py is one level up
import config as C

TORUS_TEMPLATE = C.TEMPLATE_DIR / "torus.ski"
COMBINED_TEMPLATE = C.TEMPLATE_DIR / "combined.ski"

# --------------------------------------------------------------------------
# Shared .ski pieces
# --------------------------------------------------------------------------

def gas_mix_xml(kind: str = C.TORUS_GAS_MIX) -> str:
	values = ", ".join(f"{a:.4e}" for a in C.read_cloudy_abundances())
	if kind == "atomic":
		bound = "FreeWithPolarization" if C.RESONANT_SCATTERING else "Good"
		return (f'<XRayAtomicGasMix abundancies="{values}" temperature="1e4 K" '
				f'scatterBoundElectrons="{bound}" H2fraction="0"/>')
	if kind == "ionic":
		ions = ", ".join(C.ELEMENTS)  # bare symbol = neutral (+0) ion
		electrons = "FreeWithPolarization" if C.RESONANT_SCATTERING else "FreeBound"
		return (f'<XRayIonicGasMix ions="{ions}" abundances="{values}" temperature="1e4 K" '
				f'electronScattering="{electrons}" resonantScattering="false"/>')
	raise ValueError(f"unknown TORUS_GAS_MIX '{kind}' (use 'atomic' or 'ionic')")


def instruments_xml(inclinations) -> str:
	pad = " " * 20
	return "\n".join(
		f'{pad}<SEDInstrument instrumentName="i{inc:02d}" distance="1 Mpc" inclination="{inc} deg" '
		f'azimuth="0 deg" roll="0 deg" radius="0 pc" recordComponents="true" numScatteringLevels="0" '
		f'recordPolarization="false" recordStatistics="false"/>'
		for inc in inclinations)


def list_mesh(values, lo, hi, digits: int = 8) -> str:
	"""ListMesh points (normalised to [0, 1] over [lo, hi]) for the given borders."""
	pts = sorted({min(1.0, max(0.0, (v - lo) / (hi - lo))) for v in values} | {0.0, 1.0})
	return ", ".join(f"{p:.{digits}f}" for p in pts)


def common_values(inclinations) -> dict:
	lya = ""
	if C.RESONANT_SCATTERING:
		lya = ('                <lyaOptions type="LyaOptions">\n'
			   f'                    <LyaOptions lyaAccelerationScheme="{C.LYA_ACCELERATION}" '
			   'lyaAccelerationStrength="1" includeHubbleFlow="false"/>\n'
			   '                </lyaOptions>\n')
	return dict(
		simulationMode="LyaExtinctionOnly" if C.RESONANT_SCATTERING else "ExtinctionOnly",
		lyaOptions=lya,
		resonantScattering="true" if C.RESONANT_SCATTERING else "false",
		shellElectronScattering="FreeWithPolarization" if C.RESONANT_SCATTERING else "Free",
		numPackets=C.NUM_PACKETS,
		integratedLuminosity=f"{C.SOURCE_LUMINOSITY_ERG_S:.6e} erg/s",
		minWavelength=f"{C.E_MAX_KEV:g} keV",
		maxWavelength=f"{C.E_MIN_KEV:g} keV",
		numWavelengths=C.N_ENERGY_BINS,
		instruments=instruments_xml(inclinations),
	)


def torus_values(log_nh: float, alpha: float, gas_mix: str) -> dict:
	return dict(
		torusOpeningAngle=f"{alpha:g} deg",
		torusMinRadius=f"{C.TORUS_R_IN_CM:.6e} cm",
		torusMaxRadius=f"{C.TORUS_R_OUT_CM:.6e} cm",
		torusGasMix=gas_mix,
		torusAxisColumnDensity=f"{2.0 * 10.0 ** log_nh:.6e} 1/cm2",
	)


def write_manifest(path: Path, rows) -> None:
	with open(path, "w", newline="") as f:
		writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
		writer.writeheader()
		writer.writerows(rows)

# --------------------------------------------------------------------------
# CLOUDY output
# --------------------------------------------------------------------------

_ION_RE = re.compile(r"^([A-Za-z]+)(\+?)([0-9]*)$")


def cloudy_to_skirt_label(label: str):
	"""Cloudy 'H' / 'H+' / 'He+2' -> SKIRT 'H+0' / 'H+1' / 'He+2' (None if not an ion)."""
	m = _ION_RE.match(label.strip())
	if not m:
		return None
	symbol, plus, digits = m.groups()
	charge = 0 if not plus else (int(digits) if digits else 1)
	return f"{symbol}+{charge}"


def template_ions(template_path: Path = COMBINED_TEMPLATE):
	m = re.search(r'ions="([^"]*)"', template_path.read_text())
	return [tok.strip() for tok in m.group(1).split(",") if tok.strip()]


def _read_table(path: Path):
	lines = [l for l in path.read_text().splitlines() if l.strip()]
	header = [h.strip() for h in lines[0].lstrip("#").split("\t")]
	data = np.array([[float(x) for x in l.split("\t")] for l in lines[1:]])
	return header, data


def read_cloudy_run(run_dir: Path):
	"""Zone edges, hydrogen density, temperature and per-ion densities of a CLOUDY run."""
	zh, zones = _read_table(run_dir / "sim.zones")
	zl = [h.lower() for h in zh]
	radius, dr = zones[:, zl.index("radius")], zones[:, zl.index("dr")]

	oh, ovr = _read_table(run_dir / "sim.ovr")
	ol = [h.lower() for h in oh]
	te, hden = ovr[:, ol.index("te")], ovr[:, ol.index("hden")]

	sh, spec = _read_table(run_dir / "sim.species")
	species = {}
	for j, label in enumerate(sh[1:], start=1):  # column 0 is depth
		ion = cloudy_to_skirt_label(label)
		if ion is not None:
			species[ion] = spec[:, j]

	n = min(len(radius), len(te), len(spec))
	if len({len(radius), len(te), len(spec)}) > 1:
		print(f"warning {run_dir.name}: zone counts differ, truncating to {n}")
	reasons = re.findall(r"Calculation stopped because (.*?)\.", (run_dir / "sim.out").read_text(errors="ignore"))
	return dict(
		rmin=(radius - dr / 2)[:n], rmax=(radius + dr / 2)[:n],
		hden=hden[:n], te=te[:n], species={k: v[:n] for k, v in species.items()},
		stop_reason=reasons[-1].strip() if reasons else "unknown",
	)


def rebin(run: dict, edges: np.ndarray):
	"""Column-conserving rebinning of piecewise-constant zone profiles onto `edges`.

	Cumulative columns are piecewise linear in r, so interpolating them at the new
	edges and differencing is exact."""
	zone_edges = np.concatenate(([run["rmin"][0]], run["rmax"]))
	dr = run["rmax"] - run["rmin"]
	width = np.diff(edges)

	def binned_column(density):
		cum = np.concatenate(([0.0], np.cumsum(density * dr)))
		return np.diff(np.interp(edges, zone_edges, cum))

	col_h = binned_column(run["hden"])
	n_h = col_h / width
	temperature = binned_column(run["hden"] * run["te"]) / col_h
	fractions = {ion: binned_column(dens) / col_h for ion, dens in run["species"].items()}
	return n_h, temperature, fractions


def prune_trace_ions(fractions: dict, min_fraction: float = C.SHELL_ION_MIN_FRACTION) -> dict:
	"""Zero the ions below min_fraction of their element in each cell, scaling the
	element's remaining ions so that its total abundance is unchanged."""
	if min_fraction <= 0:
		return fractions
	by_element = {}
	for ion in fractions:
		by_element.setdefault(ion.split("+")[0], []).append(ion)
	pruned = {}
	for ions in by_element.values():
		values = np.array([fractions[ion] for ion in ions])
		total = values.sum(axis=0)
		kept = np.where(values >= min_fraction * total, values, 0.0)
		scale = np.divide(total, kept.sum(axis=0), out=np.zeros_like(total), where=kept.sum(axis=0) > 0)
		pruned.update(zip(ions, kept * scale))
	return pruned


CHECK_FIELDS = ["cloudy_run", "ok", "stop_reason", "logNH_target", "logNH_real", "r_out_real_cm",
				"neutral_H_fraction", "Te_min", "Te_max", "n_zones"]


def check_cloudy(name: str, run: dict, target_log_nh: float) -> dict:
	dr = run["rmax"] - run["rmin"]
	nh = float(np.sum(run["hden"] * dr))
	h0 = run["species"].get("H+0", np.zeros_like(dr))
	return {
		"cloudy_run": name,
		"stop_reason": run["stop_reason"],
		"logNH_target": target_log_nh,
		"logNH_real": math.log10(nh),
		"r_out_real_cm": float(run["rmax"][-1]),
		"neutral_H_fraction": float(np.sum(h0 * dr) / nh),
		"Te_min": float(run["te"].min()),
		"Te_max": float(run["te"].max()),
		"n_zones": len(dr),
		"ok": ("outer radius" in run["stop_reason"]) and abs(math.log10(nh) - target_log_nh) < 0.02,
	}


def shell_txt(edges, n_h, temperature, fractions, ion_list, th_in, th_out) -> str:
	lines = [
		"# Column 1: box rmin (cm)", "# Column 2: box thetamin (deg)", "# Column 3: box phimin (deg)",
		"# Column 4: box rmax (cm)", "# Column 5: box thetamax (deg)", "# Column 6: box phimax (deg)",
		"# Column 7: number density (1/cm3)", "# Column 8: temperature (K)",
	]
	lines += [f"# Column {i}: {ion} (1)" for i, ion in enumerate(ion_list, start=9)]
	zero = np.zeros_like(n_h)
	walls = [(th_in, th_out), (180.0 - th_out, 180.0 - th_in)]
	for i in range(len(n_h)):
		tail = [f"{n_h[i]:.6e}", f"{temperature[i]:.6e}"] + [f"{fractions.get(ion, zero)[i]:.6e}" for ion in ion_list]
		for lo, hi in walls:
			lines.append("\t".join([f"{edges[i]:.10e}", f"{lo:.6f}", "0", f"{edges[i + 1]:.10e}", f"{hi:.6f}", "0"] + tail))
	return "\n".join(lines) + "\n"

# --------------------------------------------------------------------------
# Torus-only grid (stage 1)
# --------------------------------------------------------------------------

def generate_torus(output_root: Path = C.TORUS_GRID_DIR):
	template = TORUS_TEMPLATE.read_text()
	gas_mix = gas_mix_xml()
	axis_names = list(C.TORUS_GRID)
	rows, n_reset = [], 0
	for i, combo in enumerate(itertools.product(*C.TORUS_GRID.values())):
		p = dict(zip(axis_names, combo))
		run_dir = output_root / f"run_{i:05d}"
		alpha = p["torus_opening_angle"]
		values = common_values(C.TORUS_INCLINATIONS)
		# no resonant medium here: plain ExtinctionOnly (the torus gas mix still
		# follows RESONANT_SCATTERING so its physics matches the combined grid)
		values.update(simulationMode="ExtinctionOnly", lyaOptions="")
		values.update(torus_values(p["torus_logNH"], alpha, gas_mix))
		values.update(
			gridMinRadius=f"{C.TORUS_R_IN_CM:.6e} cm",
			gridMaxRadius=f"{C.TORUS_R_OUT_CM:.6e} cm",
			radialMeshPoints=", ".join(f"{x:.4f}" for x in np.linspace(0, 1, 51)),
			polarMeshPoints=list_mesh([90 - alpha, 90 + alpha], 0.0, 180.0),
		)
		n_reset += C.sync_run_dir(run_dir, {"model.ski": template.format(**values),
											 "sed.txt": C.skirt_sed(p["photon_index"])})
		rows.append({"run": run_dir.name, **p})
	write_manifest(output_root / "grid_manifest.csv", rows)
	print(f"Wrote {len(rows)} torus runs under {output_root}/ (gas mix: {C.TORUS_GAS_MIX})")
	C.report_sync(output_root, len(rows), n_reset)

# --------------------------------------------------------------------------
# Combined grid (stage 2/3)
# --------------------------------------------------------------------------

def generate_combined(target: C.Target, output_root: Path = None, cloudy_root: Path = None) -> int:
	output_root = output_root or target.combined_dir
	cloudy_root = cloudy_root or target.cloudy_dir
	grid = target.combined_grid
	if not (cloudy_root / "grid_manifest.csv").exists():
		sys.exit(f"{cloudy_root}/grid_manifest.csv not found -- run the CLOUDY steps for {target.name} first (./grid.sh {target.name})")
	with open(cloudy_root / "grid_manifest.csv", newline="") as f:
		cloudy_rows = list(csv.DictReader(f))
	cloudy_index = {tuple(float(r[a]) for a in C.SHELL_AXES): r for r in cloudy_rows}

	template = COMBINED_TEMPLATE.read_text()
	gas_mix = gas_mix_xml()
	ion_list = template_ions()
	edges = np.logspace(math.log10(C.SHELL_R_IN_CM), math.log10(C.SHELL_R_OUT_CM), C.SHELL_N_RADIAL + 1)
	# the spatial grid spans torus and shell; inside the (uniform) torus a few
	# log-spaced borders are enough, its edges are always cell borders
	r_min = min(C.TORUS_R_IN_CM, C.SHELL_R_IN_CM)
	r_max = max(C.TORUS_R_OUT_CM, C.SHELL_R_OUT_CM)
	torus_edges = np.logspace(math.log10(C.TORUS_R_IN_CM), math.log10(C.TORUS_R_OUT_CM), 11)
	# 12 digits: over three decades, 8 would shift the borders by ~1e-6 of a cell
	radial_mesh = list_mesh(list(edges) + list(torus_edges), r_min, r_max, digits=12)
	inclinations = target.inclinations

	shell_cache, checks = {}, {}
	axis_names = list(grid)
	rows, skipped, n_reset = [], 0, 0
	for i, combo in enumerate(itertools.product(*grid.values())):
		p = dict(zip(axis_names, combo))
		key = tuple(float(p[a]) for a in C.SHELL_AXES)
		if key not in cloudy_index:
			raise KeyError(f"no CLOUDY run for {dict(zip(C.SHELL_AXES, key))} -- rerun ./grid.sh {target.name} --only cloudy-gen")
		cname = cloudy_index[key]["run"]

		if cname not in shell_cache:
			cdir = cloudy_root / cname
			if not (cdir / "sim.out").exists():
				shell_cache[cname] = None
				checks[cname] = {"cloudy_run": cname, "stop_reason": "not run", "ok": False}
			else:
				run = read_cloudy_run(cdir)
				checks[cname] = check_cloudy(cname, run, p["shell_logNH"])
				if checks[cname]["ok"]:
					n_h, temperature, fractions = rebin(run, edges)
					shell_cache[cname] = n_h, temperature, prune_trace_ions(fractions)
				else:
					shell_cache[cname] = None
		if shell_cache[cname] is None:
			skipped += 1
			continue
		n_h, temperature, fractions = shell_cache[cname]

		alpha = p["torus_opening_angle"]
		th_in, th_out = C.shell_wall(p["shell_opening_angle"])
		clipped = th_out > 90.0 - alpha
		th_out = min(th_out, 90.0 - alpha)
		if th_in >= th_out:
			raise ValueError(f"shell wall {p['shell_opening_angle']} deg lies entirely inside a torus "
							 f"with alpha = {alpha} deg -- adjust COMBINED_GRID in config/{target.name}.py")

		run_dir = output_root / f"run_{i:05d}"
		values = common_values(inclinations)
		values.update(torus_values(p["torus_logNH"], alpha, gas_mix))
		values.update(
			gridMinRadius=f"{r_min:.6e} cm",
			gridMaxRadius=f"{r_max:.6e} cm",
			radialMeshPoints=radial_mesh,
			polarMeshPoints=list_mesh([th_in, th_out, 90 - alpha, 90 + alpha, 180 - th_out, 180 - th_in], 0.0, 180.0),
		)
		n_reset += C.sync_run_dir(run_dir, {
			"model.ski": template.format(**values),
			"shell.txt": shell_txt(edges, n_h, temperature, fractions, ion_list, th_in, th_out),
			"sed.txt": C.skirt_sed(p["photon_index"]),
		})
		rows.append({"run": run_dir.name, **p, "cloudy_run": cname,
					 "shell_theta_in": th_in, "shell_theta_out": th_out, "shell_clipped": clipped})

	check_rows = sorted(checks.values(), key=lambda r: r["cloudy_run"])
	with open(cloudy_root / "cloudy_check.csv", "w", newline="") as f:
		writer = csv.DictWriter(f, fieldnames=CHECK_FIELDS, restval="")
		writer.writeheader()
		writer.writerows(check_rows)

	bad = [r for r in check_rows if not r["ok"]]
	for r in bad:
		print(f"PROBLEM {r['cloudy_run']}: {r['stop_reason']}"
			  + (f", logNH {r['logNH_real']:.3f} vs {r['logNH_target']:.2f}" if "logNH_real" in r else ""))
	if rows:
		write_manifest(output_root / "grid_manifest.csv", rows)
	n_clipped = sum(r["shell_clipped"] for r in rows)
	print(f"Wrote {len(rows)} combined runs under {output_root}/ ({n_clipped} with the shell clipped by the torus), "
		  f"inclinations {inclinations}; "
		  f"{len(check_rows) - len(bad)}/{len(check_rows)} CLOUDY runs OK, {skipped} combined runs skipped.")
	C.report_sync(output_root, len(rows) + skipped, n_reset)
	return 1 if (bad or skipped) else 0


if __name__ == "__main__":
	which = sys.argv[1] if len(sys.argv) > 1 else ""
	if which == "torus":
		generate_torus()
	elif which:
		sys.exit(generate_combined(C.load_target(which)))
	else:
		sys.exit(f"usage: python3 steps/skirt_generator.py torus|TARGET   (targets: {', '.join(C.target_names())})")
