"""
plot_reduction.py
=================
Quick-look products for xmm_reduce.sh (its `summary` step):

  reduction/lightcurves.png   >10 keV light curves with flare threshold and GTIs
  reduction/regions.png       sky images with source/background regions
  spectra/summary.txt         exposure, counts, net rates, BACKSCAL per camera

Usage: python3 plot_reduction.py TARGET_DIR
"""

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

CAMERAS = ["pn", "mos1", "mos2"]


def read_conf(path: Path) -> dict:
	conf = {}
	for line in path.read_text().splitlines():
		m = re.match(r'^\s*([A-Z_]+)=(.*)$', line)
		if m:
			conf[m.group(1)] = m.group(2).strip().strip('"')
	return conf


def circle(spec: str):
	x, y, r = (float(v) for v in re.findall(r"[-\d.]+", spec))
	return x, y, r


def plot_lightcurves(red: Path, conf: dict):
	fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
	t0 = None
	for ax, cam in zip(axes, CAMERAS):
		rate = fits.open(red / f"{cam}_rate.fits")["RATE"].data
		gti = fits.open(red / f"{cam}_gti.fits")[1].data
		t0 = rate["TIME"][0] if t0 is None else t0
		thr = float(conf["PN_RATE_MAX" if cam == "pn" else "MOS_RATE_MAX"])
		ax.step((rate["TIME"] - t0) / 1e3, rate["RATE"], where="mid", lw=0.7, color="k")
		ax.axhline(thr, color="tab:red", ls="--", label=f"threshold {thr} cts/s")
		for start, stop in zip(gti["START"], gti["STOP"]):
			ax.axvspan((start - t0) / 1e3, (stop - t0) / 1e3, color="tab:green", alpha=0.15)
		good = float(np.sum(gti["STOP"] - gti["START"]))
		ax.set_ylabel(f"{cam} rate")
		ax.set_yscale("log")
		ax.set_title(f"{cam}: {good / 1e3:.1f} ks in {len(gti)} GTI(s)", fontsize=10, loc="left")
		ax.legend(loc="upper right", fontsize=8)
	axes[-1].set_xlabel("time (ks)")
	fig.suptitle(f"{conf['TARGET']}: >10 keV single-event light curves (green = GTI)")
	fig.tight_layout()
	fig.savefig(red / "lightcurves.png", dpi=120)


def plot_regions(red: Path, conf: dict):
	fig, axes = plt.subplots(1, 3, figsize=(15, 5))
	for ax, cam in zip(axes, CAMERAS):
		hdu = fits.open(red / f"{cam}_image.fits")[0]
		img, h = hdu.data, hdu.header
		# physical X,Y of the image axes (LTM/LTV keywords written by evselect)
		x0 = (1 - h["LTV1"]) / h["LTM1_1"]
		dx = 1 / h["LTM1_1"]
		y0 = (1 - h["LTV2"]) / h["LTM2_2"]
		dy = 1 / h["LTM2_2"]
		extent = [x0 - dx / 2, x0 + dx * (img.shape[1] - 0.5), y0 - dy / 2, y0 + dy * (img.shape[0] - 0.5)]
		ax.imshow(np.log10(img + 1), origin="lower", extent=extent, cmap="magma")
		prefix = "PN" if cam == "pn" else "MOS"
		circles = [circle(conf[f"{prefix}_{kind}"]) for kind in ("SRC", "BKG")]
		for (x, y, r), kind, color in zip(circles, ("src", "bkg"), ("cyan", "lime")):
			ax.add_patch(plt.Circle((x, y), r, fill=False, color=color, lw=1.2, label=kind))
		# zoom on the regions with a margin
		lo_x = min(x - r for x, y, r in circles); hi_x = max(x + r for x, y, r in circles)
		lo_y = min(y - r for x, y, r in circles); hi_y = max(y + r for x, y, r in circles)
		half = 0.5 * max(hi_x - lo_x, hi_y - lo_y) + 3000
		cx, cy = 0.5 * (lo_x + hi_x), 0.5 * (lo_y + hi_y)
		ax.set_xlim(cx - half, cx + half)
		ax.set_ylim(cy - half, cy + half)
		ax.set_title(cam)
		ax.set_xlabel("X (sky px)")
		ax.legend(fontsize=8)
	axes[0].set_ylabel("Y (sky px)")
	fig.suptitle(f"{conf['TARGET']}: cleaned images with extraction regions")
	fig.tight_layout()
	fig.savefig(red / "regions.png", dpi=120)


def summary(spec: Path, conf: dict):
	lines = [f"{conf['TARGET']} (obsid {conf['OBSID']})",
			 f"{'camera':<6} {'exposure(s)':>11} {'src cts':>8} {'bkg cts':>8} {'BACKSCAL s/b':>14} {'net rate':>9} {'bkg frac':>8}"]
	for cam in CAMERAS:
		s = fits.open(spec / f"{cam}_src.fits")["SPECTRUM"]
		b = fits.open(spec / f"{cam}_bkg.fits")["SPECTRUM"]
		ns, nb = s.data["COUNTS"].sum(), b.data["COUNTS"].sum()
		scale = (s.header["BACKSCAL"] / b.header["BACKSCAL"]) * (s.header["EXPOSURE"] / b.header["EXPOSURE"])
		net = (ns - nb * scale) / s.header["EXPOSURE"]
		lines.append(f"{cam:<6} {s.header['EXPOSURE']:>11.0f} {ns:>8d} {nb:>8d} "
					 f"{s.header['BACKSCAL']:>6.3g}/{b.header['BACKSCAL']:<7.3g} {net:>9.4f} {nb * scale / ns:>8.3f}")
	(spec / "summary.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
	target = Path(sys.argv[1]).resolve()
	conf = read_conf(target / "target.conf")
	plot_lightcurves(target / "reduction", conf)
	plot_regions(target / "reduction", conf)
	summary(target / "spectra", conf)
