#!/bin/bash
#SBATCH --job-name=na2026_eerie
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/eerie.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=mh0033
#
# EERIE ICON-ESM-ER (10 km atmosphere / 5 km ocean), CMORized Omon/tos on the
# 0.25 deg regular `gr` grid -- no unstructured-grid handling needed here.
#
#   eerie_hist_ssp245  hist-1950 (1950-2014) + highres-future-ssp245 (2015-2050)
#                      -> the "EERIE historical + scenario" line, covers 2026
#   eerie_control      eerie-control-1950 (1950-2050), constant 1950 forcing
#                      -> unforced reference: how often a 2026-like pattern
#                         turns up with no trend at all
#
# Only r1i1p1f1 is on disk for the forced runs at this location.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"

WORK=${WORK_BASE}/eerie
mkdir -p "$WORK" "${PROJECT_ROOT}/logs"

CMOR=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
CTRL=/pool/data/EERIE/EERIE/MPI-M/ICON-ESM-ER/eerie-control-1950

run_case() {
    label=$1; shift
    files="$*"
    n=$(echo $files | wc -w)
    echo "==================== ${label} (${n} files) ===================="

    cdo -s -O -selmon,6/8 -mergetime $files "${WORK}/${label}_jja_monthly.nc"
    cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${WORK}/${label}_jja_monthly.nc" \
        "${WORK}/${label}_clim.nc"
    cdo -s -O -yearmean -ymonsub "${WORK}/${label}_jja_monthly.nc" "${WORK}/${label}_clim.nc" \
        "${WORK}/${label}_anom_jja.nc"

    # global ocean mean, taken on the global field before regional subsetting
    cdo -s -O -fldmean "${WORK}/${label}_anom_jja.nc" "${OUT_BASE}/${label}_jja_gmsst.nc"

    WGT=${WORK}/wgt_${label}_to_na025.nc
    [ -s "$WGT" ] || cdo -s gencon,"$GRID_NA025" "${WORK}/${label}_anom_jja.nc" "$WGT"
    cdo -s -O -remap,"$GRID_NA025","$WGT" "${WORK}/${label}_anom_jja.nc" \
        "${OUT_BASE}/${label}_jja_anom_na025.nc"

    cdo -s showyear "${OUT_BASE}/${label}_jja_anom_na025.nc" | tr ' ' '\n' | grep -c . \
        | xargs echo "  JJA seasons:"
}

hist=$(ls ${CMOR}/hist-1950/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
ssp=$(ls ${CMOR}/highres-future-ssp245/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
run_case eerie_hist_ssp245 $hist $ssp

ctrl=$(ls ${CTRL}/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
run_case eerie_control $ctrl

echo "=== done ==="
