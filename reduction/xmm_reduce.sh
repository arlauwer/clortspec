#!/usr/bin/env bash
#
# xmm_reduce.sh
# =============
# XMM-Newton EPIC (pn, MOS1, MOS2) reduction from the ODF archive to grouped
# spectra, following docs/xmm_spectral_guide.pdf but with MOS1 and MOS2 kept
# SEPARATE: merging the event lists (as the guide does) leaves one camera's
# EXPOSURE, RMF and ARF on the counts of two cameras, which showed up as a
# MOS/pn cross-normalisation of ~2.
#
# Usage:
#   reduction/xmm_reduce.sh TARGET_DIR [--from STEP] [--to STEP] [--only STEP]
#
# Steps (in order):
#   odf      unpack TARGET_DIR/odf/*.tar.gz (outer archive + inner .TAR, or the .TAR alone)
#   init     cifbuild + odfingest                        -> reduction/ccf.cif, *SUM.SAS
#   proc     epproc + emproc                             -> reduction/{pn,mos1,mos2}.fits
#   flares   >10 keV light curves, GTIs, clean events    -> reduction/*_{rate,gti,clean}.fits
#   images   sky images for checking the regions         -> reduction/*_image.fits
#   spectra  src/bkg spectra, backscale, rmfgen, arfgen, specgroup -> spectra/
#   summary  plots and a text summary                    -> reduction/*.png, spectra/summary.txt
#
# Layout (everything derived can be regenerated from odf/ and target.conf):
#   TARGET_DIR/target.conf   obsid, regions, flare thresholds
#   TARGET_DIR/odf/          raw ODF (the downloaded .tar.gz is the only original)
#   TARGET_DIR/reduction/    SAS working directory, logs/ per step
#   TARGET_DIR/spectra/      {pn,mos1,mos2}_{src,bkg,grp}.fits + .rmf/.arf -- load *_grp.fits
#   calib/ccf/               calibration files, shared by all targets (SAS_CCFPATH)
#
# All SAS paths are set from these directories on every run, so nothing
# depends on where you call the script from or on stale environment variables.
# Override the software locations with HEADAS or SAS_SETUP, the calibration
# directory with CCF_DIR (an inherited SAS_CCFPATH is deliberately ignored).

set -euo pipefail

OBS_ROOT="$(realpath "$(dirname "$0")/..")"
STEPS=(odf init proc flares images spectra summary)

usage() { sed -n '10,21p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

TARGET_DIR=""
FROM=""
TO=""
ONLY=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --to) TO="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) TARGET_DIR="$1"; shift ;;
    esac
done
[[ -z "$TARGET_DIR" ]] && usage
TARGET_DIR="$(realpath "$TARGET_DIR")"
[[ -f "$TARGET_DIR/target.conf" ]] || { echo "error: $TARGET_DIR/target.conf not found" >&2; exit 1; }
# shellcheck disable=SC1091
source "$TARGET_DIR/target.conf"

ODF_DIR="$TARGET_DIR/odf"
RED_DIR="$TARGET_DIR/reduction"
SPEC_DIR="$TARGET_DIR/spectra"
LOG_DIR="$RED_DIR/logs"
mkdir -p "$RED_DIR" "$SPEC_DIR" "$LOG_DIR"

# --------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------

export HEADAS="${HEADAS:-$OBS_ROOT/heasoft-6.36/x86_64-pc-linux-gnu-libc2.43}"
SAS_SETUP="${SAS_SETUP:-$OBS_ROOT/sas/xmmsas_22.1.0-a8f2c2afa-20250304/setsas.sh}"
export SAS_CCFPATH="${CCF_DIR:-$OBS_ROOT/calib/ccf}"
[[ -d "$SAS_CCFPATH" ]] || { echo "error: calibration directory $SAS_CCFPATH not found" >&2; exit 1; }

if ! command -v evselect >/dev/null 2>&1; then
    set +u  # the setup scripts reference unset variables
    # shellcheck disable=SC1091
    source "$HEADAS/headas-init.sh"
    # shellcheck disable=SC1090
    source "$SAS_SETUP"
    set -u
fi
command -v evselect >/dev/null 2>&1 || { echo "error: SAS not available (check HEADAS / SAS_SETUP)" >&2; exit 1; }

sas_env() {
    # after `init`: calibration index and ODF summary of this target
    export SAS_CCF="$RED_DIR/ccf.cif"
    local sum
    sum=$(ls "$RED_DIR"/*SUM.SAS 2>/dev/null | head -1 || true)
    [[ -n "$sum" ]] || { echo "error: no *SUM.SAS in $RED_DIR -- run the init step" >&2; exit 1; }
    export SAS_ODF="$sum"
}

run() {  # run a command, appending its output to the current step's log
    echo "+ $*" >> "$LOG"
    "$@" >> "$LOG" 2>&1 || { echo "error: '$1' failed, see $LOG" >&2; exit 1; }
}

# --------------------------------------------------------------------------
# Steps
# --------------------------------------------------------------------------

step_odf() {
    cd "$ODF_DIR"
    # the science files (e.g. *_M1S00100IME.FIT); *SUM.ASC alone is not enough,
    # it sits next to the still packed inner .TAR
    if ls ./*_*_[!S]*.FIT >/dev/null 2>&1; then
        echo "  ODF already unpacked"
        return
    fi
    local outer inner
    outer=$(ls ./*.tar.gz ./*.tgz 2>/dev/null | head -1 || true)
    if ! ls ./*.TAR >/dev/null 2>&1; then
        [[ -n "$outer" ]] || { echo "error: no .tar.gz or .TAR in $ODF_DIR" >&2; exit 1; }
        run tar xzf "$outer"
    fi
    for inner in ./*.TAR; do
        [[ -e "$inner" ]] || continue
        run tar xf "$inner"
        # with the outer archive present the inner one is a copy; otherwise it
        # is the only original, so keep it
        [[ -n "$outer" ]] && rm -f "$inner"
    done
    ls ./*_*_[!S]*.FIT >/dev/null 2>&1 || { echo "error: no science files in $ODF_DIR after unpacking" >&2; exit 1; }
}

step_init() {
    cd "$RED_DIR"
    rm -f ./*SUM.SAS ccf.cif
    export SAS_ODF="$ODF_DIR"
    run cifbuild
    export SAS_CCF="$RED_DIR/ccf.cif"
    run odfingest outdir="$RED_DIR"
}

step_proc() {
    cd "$RED_DIR" && sas_env
    run epproc
    run emproc
    local cam pattern files
    for cam in pn mos1 mos2; do
        case "$cam" in
            pn)   pattern="*EPN*ImagingEvts.ds" ;;
            mos1) pattern="*EMOS1*ImagingEvts.ds" ;;
            mos2) pattern="*EMOS2*ImagingEvts.ds" ;;
        esac
        # shellcheck disable=SC2086
        files=$(ls -S $pattern)  # largest first if there are several exposures
        if (( $(wc -w <<< "$files") > 1 )); then
            echo "  warning: several $cam exposures, using the largest: $(head -1 <<< "$files")"
        fi
        mv "$(head -1 <<< "$files")" "$cam.fits"
    done
}

step_flares() {
    cd "$RED_DIR" && sas_env
    local cam flag hard thr
    for cam in pn mos1 mos2; do
        if [[ "$cam" == pn ]]; then
            flag="#XMMEA_EP"; hard="(PI>10000&&PI<12000)"; thr="$PN_RATE_MAX"
        else
            flag="#XMMEA_EM"; hard="(PI>10000)"; thr="$MOS_RATE_MAX"
        fi
        run evselect table="$cam.fits" withrateset=yes rateset="${cam}_rate.fits" \
            maketimecolumn=yes timecolumn=TIME timebinsize=100 makeratecolumn=yes \
            expression="$flag && $hard && (PATTERN==0)"
        run tabgtigen table="${cam}_rate.fits" expression="RATE<=$thr" gtiset="${cam}_gti.fits"
        run evselect table="$cam.fits" withfilteredset=yes filteredset="${cam}_clean.fits" \
            destruct=yes keepfilteroutput=yes \
            expression="$flag && gti(${cam}_gti.fits,TIME) && (PI>150)"
    done
}

step_images() {
    cd "$RED_DIR" && sas_env
    local cam
    for cam in pn mos1 mos2; do
        run evselect table="${cam}_clean.fits" withimageset=yes imageset="${cam}_image.fits" \
            imagebinning=binSize xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80
    done
}

step_spectra() {
    cd "$SPEC_DIR" && sas_env
    local cam clean src bkg pattern chmax
    for cam in pn mos1 mos2; do
        clean="$RED_DIR/${cam}_clean.fits"
        if [[ "$cam" == pn ]]; then
            src="$PN_SRC"; bkg="$PN_BKG"; pattern=4; chmax=20479
        else
            src="$MOS_SRC"; bkg="$MOS_BKG"; pattern=12; chmax=11999
        fi
        local sel="energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=$chmax"
        # shellcheck disable=SC2086
        run evselect table="$clean" withspectrumset=yes spectrumset="${cam}_src.fits" $sel \
            expression="(FLAG==0) && (PATTERN<=$pattern) && ((X,Y) IN $src)"
        # shellcheck disable=SC2086
        run evselect table="$clean" withspectrumset=yes spectrumset="${cam}_bkg.fits" $sel \
            expression="(FLAG==0) && (PATTERN<=$pattern) && ((X,Y) IN $bkg)"
        run backscale spectrumset="${cam}_src.fits" badpixlocation="$clean"
        run backscale spectrumset="${cam}_bkg.fits" badpixlocation="$clean"
        run rmfgen spectrumset="${cam}_src.fits" rmfset="$cam.rmf"
        run arfgen spectrumset="${cam}_src.fits" arfset="$cam.arf" withrmfset=yes rmfset="$cam.rmf" \
            badpixlocation="$clean" detmaptype=psf
        # run in spectra/ so the RESPFILE/ANCRFILE/BACKFILE keywords are bare file
        # names: XSPEC finds them as long as the four files stay together
        run specgroup spectrumset="${cam}_src.fits" groupedset="${cam}_grp.fits" \
            rmfset="$cam.rmf" arfset="$cam.arf" backgndset="${cam}_bkg.fits" mincounts=1 oversample=3
    done
}

step_summary() {
    run python3 "$OBS_ROOT/reduction/plot_reduction.py" "$TARGET_DIR"
    cat "$SPEC_DIR/summary.txt"
}

# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

started=0
[[ -z "$FROM" ]] && started=1
for step in "${STEPS[@]}"; do
    [[ "$step" == "$FROM" ]] && started=1
    [[ "$started" -eq 1 ]] || continue
    [[ -n "$ONLY" && "$step" != "$ONLY" ]] && continue
    LOG="$LOG_DIR/$step.log"
    : > "$LOG"
    echo "[$(date '+%T')] $step"
    "step_$step"
    [[ "$step" == "$TO" ]] && break
done
if [[ -n "$TO" && "$TO" != summary ]]; then
    echo "[$(date '+%T')] stopped after $TO"
else
    echo "[$(date '+%T')] done -- spectra in $SPEC_DIR"
fi
