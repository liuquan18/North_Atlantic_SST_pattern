#!/bin/bash
#SBATCH --job-name=hw_epoc_tmean
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_epoc_tmean.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --account=mh0033
#
# km-scale ICON (EPOC) daily MEAN 2 m temperature, April-October, Europe
# 0.25 deg -- a stand-in for daily Tmax where Tmax is not on disk.
#
# The EPOC control (epoc2_010) wrote an atm_2d_1d_max stream, but those files
# have been removed from /work and survive only on tape
# (slk:///arch/bm1313/b383127/epoc-icon-2024.10/experiments/epoc2_010/work/
# run_*/epoc2_010_atm_2d_1d_max_*.nc). Its daily means are still on disk. To
# see how much the substitution matters, the same Tmean metrics are also made
# for the transient run (epoc2_020), which has both.
#
# Usage: 07_epoc_tmean.sh <exp> [<exp> ...]     (default: epoc2_010 epoc2_020)
# Output: ${HW_WORK}/<exp>_tmean/tmax_<year>.nc -- named "tmax" so the metrics
#         step reads it unchanged, but it holds the daily MEAN temperature.
# Stream layout and end-of-day stamps as in 02_epoc_tmax.sh.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

ATM_GRID=/pool/data/ICON/grids/public/mpim/0054/icon_grid_0054_R02B08_G.nc
declare -A SRC=(
    [epoc2_010]=/work/bm1313/b383127/epoc-icon-2024.10/experiments/epoc2_010/work
    [epoc2_020]=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
)
declare -A LAST=([epoc2_010]=2024 [epoc2_020]=${REF_YEAR})

one_month() {
    ym=$1; work=$2; src=$3; exp=$4; grid=$5; gridfile=$6; wgt=$7
    out="${work}/mon/tmean_${ym}.nc"
    [ -s "$out" ] && return 0
    y=${ym:0:4}; m=${ym:4:2}
    next=$(date -u -d "${y}-${m}-01 +1 month" +%Y%m)
    dir=$(ls -d ${src}/run_${y}${m}01T000000-* 2>/dev/null | head -1)
    f1="${dir}/${exp}_atm_2d_1d_mean_${y}${m}01T000000Z.nc"
    f2="${dir}/${exp}_atm_2d_1d_mean_${next}01T000000Z.nc"
    if [ ! -s "$f1" ] || [ ! -s "$f2" ]; then echo "WARN ${exp} ${ym}: missing input" >&2; return 0; fi
    cdo -s -O -f nc -remap,"$grid","$wgt" -setgrid,"$gridfile" -shifttime,-1day \
        -selname,tas -mergetime "$f1" "$f2" "$out"
}
export -f one_month

for EXP in ${@:-epoc2_010 epoc2_020}; do
    WORK=${HW_WORK}/${EXP}_tmean
    mkdir -p "${WORK}/mon"
    WGT=${HW_WORK}/epoc2_020/wgt_r2b8_to_eu025.nc     # same R02B08 grid for both runs
    if [ ! -s "$WGT" ]; then
        ref=$(ls ${SRC[epoc2_020]}/run_20000701T*/epoc2_020_atm_2d_1d_max_20000701T000000Z.nc)
        cdo -s gencon,"$GRID_EU025" -setgrid,"$ATM_GRID" -seltimestep,1 -selname,tas "$ref" "$WGT"
    fi
    months=()
    for y in $(seq 1990 ${LAST[$EXP]}); do
        [ -s "${WORK}/tmax_${y}.nc" ] && continue
        for m in 04 05 06 07 08 09 10; do months+=("${y}${m}"); done
    done
    echo "=== ${EXP}: daily Tmean for ${#months[@]} months ==="
    [ "${#months[@]}" -gt 0 ] && parallel -j ${HW_JOBS:-32} one_month {} "$WORK" "${SRC[$EXP]}" \
        "$EXP" "$GRID_EU025" "$ATM_GRID" "$WGT" ::: "${months[@]}"
    for y in $(seq 1990 ${LAST[$EXP]}); do
        out="${WORK}/tmax_${y}.nc"
        [ -s "$out" ] && continue
        [ "$(ls ${WORK}/mon/tmean_${y}??.nc 2>/dev/null | wc -l)" -eq 7 ] || { echo "  ${y}: incomplete, skipped"; continue; }
        cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
            -setreftime,1850-01-01,00:00:00,1day -mergetime ${WORK}/mon/tmean_${y}??.nc "$out"
    done
    echo "  ${EXP}: $(ls ${WORK}/tmax_*.nc | wc -l) years"
done
echo "=== done ==="
