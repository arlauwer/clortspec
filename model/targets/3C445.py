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
# fixed by config.COMBINED_AXES. PROVISIONAL: copied from 3C 234; centre the
# torus axes and photon_index on the stage-1 fit (python3 6_xspec.py 3C445)
# before running steps 1-5.
COMBINED_GRID = {
	"torus_logNH": [23.25, 23.5],
	"torus_opening_angle": [45.0],
	# log10 xi = L(1-1000 Ryd) / (n_in R_in^2), erg cm s^-1, at the inner face
	"logxi": [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
	# log10 radial column density of the shell, R_in -> R_out, cm^-2
	"shell_logNH": [21.0, 21.5, 22.0, 22.5, 23.0],
	# cone half-opening angle from the POLE, centre of the wall (deg). Must stay
	# clear of the torus: theta_s + wall/2 <= 90 - alpha (else clipped).
	"shell_opening_angle": [20.0, 30.0, 40.0],
	"photon_index": [1.2, 1.45, 1.7],
}
