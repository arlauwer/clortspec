from pathlib import Path
from xspec import *
import os
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent
os.chdir(BASE_DIR)

table_name = 'skirt_i90.fits'

AllData.clear()
AllData.show()
AllData("1:1 pn_grp.fits 2:2 mos_grp.fits")
AllData(1).ignore("**-0.3 10.0-**")
AllData(2).ignore("**-0.3 10.0-**")

Fit.statMethod = "cstat"
Xset.cosmo = "70 0 0.73"

model = Model(f"constant*tbabs*zashift(atable{{{BASE_DIR / 'model' / table_name}}})")

AllModels(1)(1).values = "1 -1"				# PN constant (frozen)
AllModels(1)(2).values = "0.07 -1"        	# Galactic nH
AllModels(1)(3).values = "0.1849 -1"      	# 3C234 redshift
AllModels(1)(6).values = "1.0"            	# metallicity
AllModels(1)(7).values = "30" 			  	# opening angle
AllModels(2)(1).values = "2 0.01 0 0 5 5" 	# MOS constant (not frozen)

Fit.nIterations = 1000
Fit.perform()

print(Fit.statistic)           # fit statistic value
print(Fit.dof)                 # degrees of freedom
print(Fit.statistic / Fit.dof) # should be around ~1

Plot.device = "/null"  # don't open an X window
Plot.xAxis = "keV"
Plot("ldata", "delchi")

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 8),
								gridspec_kw={"height_ratios": [3, 1]})
fig.suptitle(f"3C234 XSPEC using {table_name}")

colors = ['red', 'blue']
groups = ['PN', 'MOS']

for g in range(1, AllData.nGroups + 1):
	Plot("ldata")
	x = Plot.x(g)
	xerr = Plot.xErr(g)
	y = Plot.y(g)
	yerr = Plot.yErr(g)
	model_y = Plot.model(g)

	ax1.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=".", label=f"{groups[g-1]} data", color=colors[g - 1], alpha=0.5)
	ax1.step(x, model_y, where="mid", label=f"{groups[g-1]} model", color=colors[g - 1], alpha=0.5)

	Plot("delchi")
	dx = Plot.x(g)
	dxerr = Plot.xErr(g)
	dy = Plot.y(g)
	dyerr = Plot.yErr(g)
	ax2.errorbar(dx, dy, xerr=dxerr, yerr=dyerr, fmt=".", label=f"{groups[g-1]} residuals", color=colors[g - 1], alpha=0.5)

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylabel("Counts/s/keV")
ax1.legend()

ax2.axhline(0, color="k", lw=0.5)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel(r"$\Delta\chi$")

plt.tight_layout()
plt.savefig("spectrum_fit.png", dpi=150)
plt.show()