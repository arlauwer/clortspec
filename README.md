# XMM-Newton radio galaxies: data reduction

Modelling and fitting: see `model/README.md`.

```
reduction/xmm_reduce.sh    ODF archive -> grouped pn/MOS1/MOS2 spectra
reduction/setsas.fish      HEASoft + SAS in fish, for interactive use
calib/ccf/                 XMM calibration files (shared)
<TARGET>/target.conf       obsid, RA/Dec, extraction regions, flare thresholds
<TARGET>/odf/              downloaded <obsid>.tar.gz
<TARGET>/reduction/        SAS working dir, logs/, check plots
<TARGET>/spectra/          {pn,mos1,mos2}_grp.fits (+ rmf/arf): what XSPEC loads
```

## Reduce

```bash
reduction/xmm_reduce.sh 3C234                   # all steps, ~15-30 min
reduction/xmm_reduce.sh 3C234 --from spectra    # e.g. after changing regions in target.conf
reduction/xmm_reduce.sh 3C234 --only summary    # one step
reduction/xmm_reduce.sh 3C234 --to images       # stop after a step
```

Steps: `odf` (unpack) → `init` (cifbuild, odfingest) → `proc` (epproc, emproc)
→ `flares` (light curves, GTIs) → `images` → `spectra` (extract, rmf/arf,
group) → `summary` (plots, `spectra/summary.txt`).

Check afterwards: `reduction/lightcurves.png` (flare threshold) and
`reduction/regions.png` (extraction regions).

Interactive SAS/XSPEC in fish: `source reduction/setsas.fish [3C234]`.

## New target

1. Download the ODF from the XSA into `<TARGET>/odf/`.
2. Copy `3C234/target.conf`; set `TARGET`, `OBSID`, `RA`, `DEC`.
3. `reduction/xmm_reduce.sh <TARGET> --to images`
4. Draw source/background circles on `reduction/{pn,mos1}_image.fits` in ds9;
   save them in *physical* coordinates into `target.conf`.
5. `reduction/xmm_reduce.sh <TARGET> --from spectra`, check the plots.
6. Add `model/config/<TARGET>.py` (see `model/README.md`).

Notes: MOS1 and MOS2 are reduced separately (merging them gives a
cross-normalisation of ~2). Flare thresholds: 0.4 (pn) / 0.35 (MOS) cts/s.
