#!/bin/bash
#SBATCH --job-name=na2026_eerie
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/eerie.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --account=mh0033
#
# EERIE ICON-ESM-ER (10 km atmosphere / 5 km ocean), SST on 0.25 deg regular
# grids -- no unstructured-grid handling needed here.
#
#   eerie_hist_ssp245  three realizations of hist-1950 + highres-future-ssp245
#                      -> the "EERIE historical + scenario" ensemble
#   eerie_control      eerie-control-1950 (1950-2050), constant 1950 forcing
#                      -> unforced reference: how often a 2026-like pattern
#                         turns up with no trend at all
#
# Members of the forced run:
#
#   r1  erc2020  CMORized Omon/tos (`gr`), 1950-2050
#   r2  erc2023  raw ICON output remapped to 0.25 deg by DKRZ, 1975-2014: the
#                ssp245 leg wrote no ocean output past 2015-01
#   r3  erc2024  raw output, 1975-2020; the monthly stream misses 1993-03 ..
#                1994-02 and stops at 2014, so JJA 1993 and 2015-2020 are
#                monthly means of the daily stream
#
# r2 and r3 are only in the "phase2" EERIE catalogue, not in the CMOR tree
# (/work/bm1344/DKRZ/kerchunks_batched/ICON/phase2/hist-1950/v20240618/{2,3}
# reference the files read here). Their output is stamped at the END of the
# averaging period (`..._19800701T000000Z.nc` is the June 1980 mean; a daily
# stamp of 06-02 00:00 is the 1 June mean), hence -shifttime,-1day. Each
# `run_<YYYYMM>...` directory is named after the month it represents.
#
# Every member is an anomaly against its own 1991-2020 JJA climatology --
# except r2, which has SST only to 2014 and so uses 1991-2014.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"

WORK=${WORK_BASE}/eerie
mkdir -p "$WORK" "${PROJECT_ROOT}/logs"

CMOR=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
CTRL=/pool/data/EERIE/EERIE/MPI-M/ICON-ESM-ER/eerie-control-1950
RAW=/work/bm1344/k202193/ICON

# JJA monthly file -> climatology, JJA-mean anomaly, global mean, regional remap.
#   $1 label, $2 JJA monthly input, $3 directory for the two products
anomalies() {
    label=$1; monthly=$2; out=$3
    cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "$monthly" "${WORK}/${label}_clim.nc"
    cdo -s -O -yearmean -ymonsub "$monthly" "${WORK}/${label}_clim.nc" "${WORK}/${label}_anom_jja.nc"

    # global ocean mean, taken on the global field before regional subsetting
    cdo -s -O -fldmean "${WORK}/${label}_anom_jja.nc" "${out}/${label}_jja_gmsst.nc"

    WGT=${WORK}/wgt_${label}_to_na025_${NA025_ID}.nc
    [ -s "$WGT" ] || cdo -s gencon,"$GRID_NA025" "${WORK}/${label}_anom_jja.nc" "$WGT"
    cdo -s -O -remap,"$GRID_NA025","$WGT" "${WORK}/${label}_anom_jja.nc" \
        "${out}/${label}_jja_anom_na025.nc"

    echo "  ${label}: $(cdo -s showyear "${out}/${label}_jja_anom_na025.nc" | wc -w) JJA seasons," \
         "climatology from $(cdo -s showyear -selyear,${CLIM_START}/${CLIM_END} "$monthly" | wc -w) years"
}

run_cmor() {
    label=$1; out=$2; shift 2
    echo "==================== ${label} ($# files) ===================="
    cdo -s -O -selmon,6/8 -mergetime "$@" "${WORK}/${label}_jja_monthly.nc"
    anomalies "$label" "${WORK}/${label}_jja_monthly.nc" "$out"
}

# One JJA month of a raw run as a single time step, stamped inside that month.
#   $1 run, $2 stream base dir (hist or SSP245), $3 YYYY, $4 MM, $5 output
raw_month() {
    run=$1; base=$2; y=$3; m=$4; out=$5
    mon=$(ls -d ${base}/oce_2d_1mth_mean_remap025/run_${y}${m}01T* 2>/dev/null || true)
    if [ -n "$mon" ]; then
        cdo -s -O -shifttime,-1day -selname,to ${mon}/*.nc "$out"
        return
    fi
    day=$(ls -d ${base}/oce_2d_1d_mean_remap025/run_${y}${m}01T* 2>/dev/null || true)
    [ -n "$day" ] || { echo "  ${run} ${y}-${m}: no monthly or daily ocean output" >&2; return 1; }
    cdo -s -O -monmean -selmon,$((10#$m)) -shifttime,-1day -mergetime \
        $(for f in ${day}/*.nc; do echo "-selname,to $f"; done) "$out"
    n=$(cdo -s ntime -selmon,$((10#$m)) -shifttime,-1day -mergetime \
        $(for f in ${day}/*.nc; do echo "-selname,to $f"; done))
    echo "  ${run} ${y}-${m}: from ${n} daily means"
}

run_raw() {
    run=$1; member=$2; first=$3; last=$4
    label=eerie_r${member}
    echo "==================== ${label} (${run}, ${first}-${last}) ===================="
    part=${WORK}/part_${run}; mkdir -p "$part"
    for y in $(seq "$first" "$last"); do
        base=${RAW}/${run}/postprocessing/interpolation
        [ "$y" -ge 2015 ] && base=${base}/SSP245
        for m in 06 07 08; do
            f=${part}/to_${y}${m}.nc
            [ -s "$f" ] || raw_month "$run" "$base" "$y" "$m" "$f"
        done
    done
    cdo -s -O -mergetime ${part}/to_*.nc "${WORK}/${label}_jja_monthly.nc"
    ny=$(cdo -s showyear "${WORK}/${label}_jja_monthly.nc" | wc -w)
    nt=$(cdo -s ntime "${WORK}/${label}_jja_monthly.nc")
    [ "$nt" -eq $((3 * ny)) ] || { echo "  ${label}: ${nt} months for ${ny} years" >&2; exit 1; }
    anomalies "$label" "${WORK}/${label}_jja_monthly.nc" "$WORK"
}

hist=$(ls ${CMOR}/hist-1950/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
ssp=$(ls ${CMOR}/highres-future-ssp245/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
run_cmor eerie_r1 "$WORK" $hist $ssp

run_raw erc2023 2 1975 2014
run_raw erc2024 3 1975 2020

activate_python
python3 "${SCRIPTS}/pre_process/stack_eerie.py"

ctrl=$(ls ${CTRL}/r1i1p1f1/Omon/tos/gr/v*/tos_*.nc)
run_cmor eerie_control "$OUT_BASE" $ctrl

echo "=== done ==="
