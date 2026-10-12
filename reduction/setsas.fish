# setsas.fish
# ===========
# HEASoft + SAS for interactive use in fish (xmm_reduce.sh does not need it).
#
#   source reduction/setsas.fish             # HEASoft + SAS, SAS_CCFPATH = calib/ccf
#   source reduction/setsas.fish 3C234       # also SAS_CCF/SAS_ODF of that target
#
# The SAS setup script is bash, so it is sourced in bash and the resulting
# environment is copied into this shell. HEASoft is set up the same way when
# it is not already (via $HEADAS, default: heasoft-6.36/ next to this directory).

set -l obs (realpath (dirname (status filename))/..)
set -q HEADAS; or set -gx HEADAS $obs/heasoft-6.36/x86_64-pc-linux-gnu-libc2.43
set -l setup "set +u; source $HEADAS/headas-init.sh; source $obs/sas/xmmsas_22.1.0-a8f2c2afa-20250304/setsas.sh"

for kv in (bash -c "$setup >/dev/null 2>&1; env -0" | string split0)
    set -l k (string split -m1 = -- $kv)
    switch $k[1]
        case PWD OLDPWD SHLVL _ 'BASH_FUNC_*'
            continue
        case PATH LD_LIBRARY_PATH PYTHONPATH PERL5LIB PERLLIB
            set -gx $k[1] (string split : -- $k[2])
        case '*'
            set -gx $k[1] $k[2]
    end
end

set -gx SAS_CCFPATH $obs/calib/ccf
set -e SAS_CCF SAS_ODF
if set -q argv[1]
    set -l red (realpath $argv[1])/reduction
    set -gx SAS_CCF $red/ccf.cif
    set -gx SAS_ODF (ls $red/*SUM.SAS 2>/dev/null)[1]
end

if command -q evselect
    echo "SAS ready  SAS_CCFPATH=$SAS_CCFPATH"
    set -q SAS_ODF; and echo "           SAS_CCF=$SAS_CCF  SAS_ODF=$SAS_ODF"
else
    echo "setsas.fish: SAS setup failed (check HEADAS and the sas/ directory)" >&2
end
