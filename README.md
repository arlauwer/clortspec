# XMM-Newton observations of radio galaxies + SKIRT modelling

```
observations/
├── reduction/            generic EPIC reduction (any target)
│   ├── xmm_reduce.sh     ODF archive -> grouped pn/MOS1/MOS2 spectra
│   └── plot_reduction.py light curves, region check, summary (summary step)
├── calib/ccf/            XMM calibration files (SAS_CCFPATH), shared by all targets
├── docs/                 xmm_spectral_guide.pdf, notes, Bert's heasp table script
├── sas/, heasoft-6.36/   software
├── model/                CLOUDY + SKIRT + XSPEC pipeline for all targets (see model/README.md)
│   ├── grid.sh, fit.py   build the model grids / fit a target
│   ├── targets/3C234.py  redshift, Galactic N_H, model grid of one target
│   └── runs/3C234/       its CLOUDY/SKIRT grids, XSPEC table and fit plots
└── 3C234/, 3C445/        one directory per target
    ├── target.conf       obsid, extraction regions, flare thresholds
    ├── odf/              raw data: the downloaded <obsid>.tar.gz (+ unpacked ODF)
    ├── reduction/        SAS working dir (events, GTIs, images, logs/, *.png)
    └── spectra/          {pn,mos1,mos2}_grp.fits (+ src/bkg/rmf/arf) -- what XSPEC loads
```

## Quick reference

From a downloaded observation to a fit, for a target `X` (every command has
`--help`):

```bash
# 1. reduction (details below; first time: see "Adding a target")
reduction/xmm_reduce.sh X                 # X/odf/*.tar.gz + X/target.conf -> X/spectra/

# 2. model (details in model/README.md); needs model/targets/X.py
source reduction/setsas.fish              # HEASoft for pyXspec (fish)
cd model
python3 fit.py X                          # stage 1: torus only, no new simulations
#   -> centre COMBINED_GRID in targets/X.py on the printed result
./grid.sh X -j JOBS -t THREADS            # CLOUDY + SKIRT + table (hours-days; server)
python3 fit.py X                          # stages 1-3 -> runs/X/fit_stage*.png
```

## Reducing a target

```bash
reduction/xmm_reduce.sh 3C234                  # all steps (~15-30 min)
reduction/xmm_reduce.sh 3C234 --from spectra   # e.g. after changing regions in target.conf
reduction/xmm_reduce.sh 3C234 --only summary
```

The script sets `SAS_CCFPATH`, `SAS_CCF` and `SAS_ODF` itself from the
directories above on every run, so nothing depends on the shell environment
or on where files used to live. (Moving reduction products by hand breaks
SAS because `*SUM.SAS` stores the absolute ODF path; rerunning the `init`
step regenerates it.)

For interactive SAS/XSPEC in fish: `source reduction/setsas.fish [3C234]`
(HEASoft + SAS; with a target also its `SAS_CCF`/`SAS_ODF`).

Check `reduction/lightcurves.png` (flare threshold) and `reduction/regions.png`
(extraction regions) after a run.

## Adding a target

1. Download the ODF from the XSA (obsid → "ODF", a `<obsid>.tar.gz`) into
   `<target>/odf/`.
2. Copy `3C234/target.conf` to `<target>/` and set `TARGET`, `OBSID`, `RA`,
   `DEC`; leave the regions for now.
3. `reduction/xmm_reduce.sh <target> --to images`
4. Open `reduction/{pn,mos1}_image.fits` in ds9, draw the source and background
   circles, and read them off in *physical* coordinates (Region → Save, format
   ds9, coordinate system physical) into `target.conf`.
5. `reduction/xmm_reduce.sh <target> --from spectra`, then check the two PNGs
   and `spectra/summary.txt`.
6. For the modelling, add `model/targets/<target>.py` (see `model/README.md`).

Differences from the guide: MOS1 and MOS2 are reduced and fitted separately
(merging the event lists gives one camera's exposure and response for the
counts of both, i.e. a cross-normalisation of ~2), and flares are filtered per
camera with the standard 0.4 (pn) / 0.35 (MOS) cts/s thresholds.
