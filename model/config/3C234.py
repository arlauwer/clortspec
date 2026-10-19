"""
3C 234 (z = 0.1849): target constants and the combined torus + shell grid.
Read by config.load_target("3C234"); the spectra are ../3C234/spectra/.
"""

TITLE = "3C 234"
REDSHIFT = 0.1849
# Galactic N_H,tot (HI + 2 H2, Willingale et al. 2013) towards 3C 234
# (RA 150.45625, Dec 28.78583), weighted value from
# https://www.swift.ac.uk/analysis/nhtot/ -- in units of 1e22 cm^-2 for TBabs.
NH_GALACTIC_1E22 = 0.0187

# Stage 2/3 grid. Axis order = manifest/table row order (last axis fastest),
# fixed by config.COMBINED_AXES; an axis with a single value is held fixed and
# left out of the XSPEC table. Choices (2026-10-19, from the first combined fit):
# - torus: fixed at stage 1. OpenAngle is barely constrained above 3 keV (dC < 1
#   for 45-60, best 50); 45 leaves the cone room up to 45 deg from the pole.
# - logxi / shell_logNH: the soft flux peaks where the shell just becomes
#   ionised through, which moves to higher xi for larger N_H, and changes by up
#   to 5x between 0.5-dex steps -> 0.25 dex. The first fit wanted more scattering
#   (logNHs at the top, 23). With L(1-1000 Ryd) ~ 8e44 erg/s and R_in = 1e19 cm,
#   xi * N_H = L / R_in gives logxi ~ 2.9 / 2.4 / 1.9 at logNHs 23 / 23.5 / 24.
COMBINED_GRID = {
	"torus_logNH": [23.45],
	"torus_opening_angle": [45.0],
	# log10 xi = L(1-1000 Ryd) / (n_in R_in^2), erg cm s^-1, at the inner face
	"logxi": [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5],
	# log10 radial column density of the shell, R_in -> R_out, cm^-2
	"shell_logNH": [23.0, 23.25, 23.5, 23.75, 24.0],
	# cone half-opening angle from the POLE, centre of the wall (deg). Must stay
	# inside the torus opening: theta_s + wall/2 <= 90 - alpha (checked on load).
	"shell_opening_angle": [20.0, 27.5, 35.0],
	# stage-1 best fit and 90% interval
	"photon_index": [1.309, 1.434, 1.524],
}
