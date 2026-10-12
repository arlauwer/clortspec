# Interactive counterpart of 6_xspec.py (combined model, everything free).
# Run from a target's spectra directory, e.g. 3C234/spectra/:
#   xspec - ../../model/6_xspec.tcl
# and adapt the three target-specific values below (table path, Galactic N_H, z).
data 1:1 pn_grp.fits 2:2 mos1_grp.fits 3:3 mos2_grp.fits
ignore **:**-0.3 10.0-**
ignore bad
statistic cstat

# constant * TBabs * combined table (redshift is a table parameter).
# Table parameter order follows the varying axes of config.COMBINED_GRID,
# then PhoIndex, Incl, z, norm -- check with `show par` after loading.
model constant*TBabs*atable{../../model/runs/3C234/combined.fits} & /*

newpar 1 1 -1
newpar 2 0.0187 -1
show par
# then e.g.:  newpar <z> 0.1849 -1 ; newpar <Incl> 90 -1
#             untie <MOS1 constant> <MOS2 constant> ; thaw them ; fit 10000
cpd /xw
setplot energy
