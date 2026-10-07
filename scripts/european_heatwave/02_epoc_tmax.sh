#!/bin/bash
#SBATCH --job-name=hw_epoc
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_epoc.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --account=mh0033
#
# km-scale ICON (EPOC, epoc2_020 = hist+ssp585) daily maximum 2 m temperature,
# April-October 1990-2025, Europe 0.25 deg.
#
# Source: `tas` in the epoc2_020_atm_2d_1d_max stream (daily maximum of the
# model-time-step 2 m temperature), on the unstructured R02B08 atmosphere grid
# (grid 0054, 5,242,880 cells), which -- as for the ocean -- is not in the data
# files and has to be attached.
#
# Time stamps are the END of each day: 2002-06-02T00:00 is the 1 June maximum.
# Within run_<YYYYMM>01.../ the month is split over two files: the one stamped
# <YYYYMM>01 holds days 1..n-1, the one stamped with the next month's 1st holds
# day n. Both are merged and shifted back by one day.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

EXP=epoc2_020
SRC=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/${EXP}/work
ATM_GRID=/pool/data/ICON/grids/public/mpim/0054/icon_grid_0054_R02B08_G.nc
WORK=${HW_WORK}/${EXP}
mkdir -p "${WORK}/mon"

FIRST_YEAR=${HW_FIRST_YEAR:-1990}
LAST_YEAR=${HW_LAST_YEAR:-2025}

WGT=${WORK}/wgt_r2b8_to_eu025.nc
if [ ! -s "$WGT" ]; then
    echo "=== conservative weights R02B08 -> Europe 0.25 deg (once) ==="
    ref=$(ls ${SRC}/run_20000701T*/${EXP}_atm_2d_1d_max_20000701T000000Z.nc)
    cdo -s gencon,"$GRID_EU025" -setgrid,"$ATM_GRID" -seltimestep,1 -selname,tas "$ref" "$WGT"
fi

one_month() {
    ym=$1; work=$2; src=$3; exp=$4; grid=$5; gridfile=$6; wgt=$7
    out="${work}/mon/tmax_${ym}.nc"
    [ -s "$out" ] && return 0
    y=${ym:0:4}; m=${ym:4:2}
    next=$(date -u -d "${y}-${m}-01 +1 month" +%Y%m)
    dir=$(ls -d ${src}/run_${y}${m}01T000000-* 2>/dev/null | head -1)
    f1="${dir}/${exp}_atm_2d_1d_max_${y}${m}01T000000Z.nc"
    f2="${dir}/${exp}_atm_2d_1d_max_${next}01T000000Z.nc"
    if [ ! -s "$f1" ] || [ ! -s "$f2" ]; then echo "WARN ${ym}: missing input" >&2; return 0; fi
    cdo -s -O -f nc -remap,"$grid","$wgt" -setgrid,"$gridfile" -shifttime,-1day \
        -selname,tas -mergetime "$f1" "$f2" "$out"
    n=$(cdo -s ntime "$out"); want=$(date -u -d "${y}-${m}-01 +1 month -1 day" +%d)
    [ "$n" -eq "$((10#$want))" ] || echo "WARN ${ym}: ${n} days, expected ${want}" >&2
}
export -f one_month

months=()
for y in $(seq "$FIRST_YEAR" "$LAST_YEAR"); do
    [ -s "${WORK}/tmax_${y}.nc" ] && continue
    for m in 04 05 06 07 08 09 10; do months+=("${y}${m}"); done
done
echo "=== daily Tmax for ${#months[@]} months ==="
# ~5 GB per cdo process while the 5.2M-cell grid is attached
[ "${#months[@]}" -gt 0 ] && parallel -j ${HW_JOBS:-32} one_month {} "$WORK" "$SRC" "$EXP" "$GRID_EU025" "$ATM_GRID" "$WGT" ::: "${months[@]}"

echo "=== merging per year ==="
for y in $(seq "$FIRST_YEAR" "$LAST_YEAR"); do
    out="${WORK}/tmax_${y}.nc"
    [ -s "$out" ] && continue
    [ "$(ls ${WORK}/mon/tmax_${y}??.nc 2>/dev/null | wc -l)" -eq 7 ] || { echo "  ${y}: incomplete, skipped"; continue; }
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -setreftime,1850-01-01,00:00:00,1day -mergetime ${WORK}/mon/tmax_${y}??.nc "$out"
    echo "  ${y}: $(cdo -s ntime "$out") days"
done
echo "=== done ==="
