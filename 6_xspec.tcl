data 1:1 pn_grp.fits 2:2 mos_grp.fits
ignore 1: **-0.3 10.0-**
ignore 2: **-0.3 10.0-**
statistic cstat
cosmo 70 0 0.73
model constant*tbabs*zashift(atable{model/skirtxray.fits})

newpar 1 1 -1
newpar 2 0.07 -1
newpar 3 0.1849 -1
newpar 6 1.0 -1
newpar 7 30 -1
newpar 10 2 0.01 0 0 5 5

cpd /xw
setplot energy
setplot command r x 0.3 10
plot model ldata delchi