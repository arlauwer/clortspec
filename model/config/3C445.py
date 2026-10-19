"""
3C 445 (z = 0.0562): target constants and the combined torus + shell grid.
Read by config.load_target("3C445"); the spectra are ../3C445/spectra/.
"""

TITLE = "3C 445"
REDSHIFT = 0.0562
# Galactic N_H,tot (HI + 2 H2, Willingale et al. 2013) towards 3C 445
# (RA 335.95708, Dec -2.10361), mean = weighted value from
# https://www.swift.ac.uk/analysis/nhtot/ -- in units of 1e22 cm^-2 for TBabs.
NH_GALACTIC_1E22 = 0.0570

# Stage 2/3 grid. Axis order = manifest/table row order (last axis fastest),
# fixed by config.COMBINED_AXES; an axis with a single value is held fixed and
# left out of the XSPEC table. Choices (2026-10-19, from the first combined fit):
# - torus: fixed at stage 1. OpenAngle is well constrained above 3 keV (best 60;
#   dC = +7 at 55, +16 at 50, +26 at 45), which leaves the cone at most 30 deg
#   from the pole. If the shell still falls short, that tension (wide torus for
#   Fe K vs wide cone for the soft band) is itself the result to look at.
# - logxi / shell_logNH: see 3C234.py. With L(1-1000 Ryd) ~ 1.3e44 erg/s and
#   R_in = 1e19 cm, logxi ~ 2.1 / 1.6 / 1.1 at logNHs 23 / 23.5 / 24.
COMBINED_GRID = {
	"torus_logNH": [23.34],
	"torus_opening_angle": [60.0],
	# log10 xi = L(1-1000 Ryd) / (n_in R_in^2), erg cm s^-1, at the inner face
	"logxi": [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0],
	# log10 radial column density of the shell, R_in -> R_out, cm^-2
	"shell_logNH": [23.0, 23.25, 23.5, 23.75, 24.0],
	# cone half-opening angle from the POLE, centre of the wall (deg). Must stay
	# inside the torus opening: theta_s + wall/2 <= 90 - alpha (checked on load).
	"shell_opening_angle": [10.0, 15.0, 20.0],
	# stage-1 best fit and 90% interval
	"photon_index": [1.288, 1.355, 1.407],
}
