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
# fixed by config.COMBINED_AXES. The torus axes are centred on the stage-1 fit
# (logNH = 23.45, OpenAngle = 52 +- 12, PhoIndex = 1.43); an axis with a single
# value is held fixed and left out of the XSPEC table.
COMBINED_GRID = {
	"torus_logNH": [23.45],
	"torus_opening_angle": [30.0, 45.0],
	# log10 xi = L(1-1000 Ryd) / (n_in R_in^2), erg cm s^-1, at the inner face
	"logxi": [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
	# log10 radial column density of the shell, R_in -> R_out, cm^-2
	"shell_logNH": [21.0, 21.5, 22.0, 22.5, 23.0],
	# cone half-opening angle from the POLE, centre of the wall (deg). Must stay
	# clear of the torus: theta_s + wall/2 <= 90 - alpha (else clipped).
	"shell_opening_angle": [20.0, 30.0, 40.0],
	"photon_index": [1.309, 1.434, 1.524],
}
