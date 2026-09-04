import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Parameters per run: (xi, hden, phot)
params = {
	0: (-2, 10000.0,   0),
	1: (-2, 10000.0,   2),
	2: (-2, 1000000.0, 0),
	3: (-2, 1000000.0, 2),
	4: ( 2, 10000.0,   0),
	5: ( 2, 10000.0,   2),
	6: ( 2, 1000000.0, 0),
	7: ( 2, 1000000.0, 2),
}

radii = np.array([1.0000000119e+18,1.0226601100e+18,1.0000000001e+18,1.0015082180e+18,8.1734240000e+18,7.3640890000e+18,2.9514965000e+18,2.2467895000e+18])
radii -= 1e18

# Encode each variable on its own visual channel
xi_color   = {-2: '#e63946', 2: '#1d3557'}       # red vs dark blue
phot_style = {0: '-', 2: '--'}                    # solid vs dashed
hden_width = {10000.0: 1.5, 1000000.0: 3.0}       # thin vs thick

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

for run_id, (xi, hden, phot) in params.items():
	name = f'run_0000{run_id}'
	with fits.open(f'../cloudy_grid/{name}/skirt/model_opacity_gas_tau.fits') as hdu:
		tau = hdu[0].data[:, 0, 0]
		tau /= radii[run_id]
		E = np.array(hdu[1].data, dtype=float)

	ax.plot(
		E, tau,
		color=xi_color[xi],
		linestyle=phot_style[phot],
		linewidth=hden_width[hden],
		alpha=0.85,
		label=f'xi={xi}, hden={hden:.0e}, phot={phot}'
	)

ax.set_title('Opacity across grid')
ax.set_xlabel('energy (keV)')
ax.set_ylabel(r'$\tau$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=8, ncol=2)
ax.grid(True, which='both', alpha=0.2)

plt.tight_layout()
plt.show()