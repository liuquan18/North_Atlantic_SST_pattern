#!/bin/bash
#SBATCH --job-name=na2026_sap0006
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/sap0006.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --account=mh0033
#
# ICON Sapphire 100-yr control run sap0006 (10 km atm / 5 km ocean), SST
# pattern products AND daily Tmax for the heatwave metrics, from its HEALPix
# zarr store (see extract_sap0006.py).
#
# The run is a "1950 control run for Sapphire 2.01 tuning" (run script
# README): CO2/CH4/N2O/CFCs at constant 1950 values (irad_*=2), ozone and
# Kinne aerosol files fixed to 1950 for every year, no volcanic aerosol
# (irad_aero=13). Until 1960-07-01 the solar constant was transient
# (isolrad=1), then fixed (isolrad=2); 1950-1960 is also the start of the
# run, so seasons before FIRST_YEAR are not used.
#
# Products (same as every other dataset):
#   data/pattern_2026/sap0006_jja_anom_na025.nc, sap0006_jja_gmsst.nc
#   ${HW_WORK}/sap0006/tmax_<year>.nc   daily Tmax, Apr-Oct, GRID_EU025
#
# Anomalies against the run's own model years 1991-2020, like the other
# control runs. The Tmax is the maximum of 8 three-hourly instantaneous
# values (03..24 UTC) -- slightly below a true daily maximum, but the
# heatwave definition works on anomalies against the run's own percentiles.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

WORK=${WORK_BASE}/sap0006
TMAX=${HW_WORK}/sap0006
mkdir -p "$WORK" "$TMAX"
FIRST_YEAR=${SAP_FIRST_YEAR:-1961}
LAST_YEAR=2049
BASEPY="env PYTHONNOUSERSITE=1 /sw/spack-levante/mambaforge-23.1.0-1-Linux-x86_64-3boc6i/bin/python3"
EXTRACT="${SCRIPTS}/pre_process/extract_sap0006.py"
YEARS=$(seq "$FIRST_YEAR" "$LAST_YEAR")

# zoom-8 cell centres, for selecting the European cells (cdo cannot print them)
[ -s "${WORK}/hp256_centres.nc" ] || cdo -s -f nc -remapnn,hp256_nested \
    -expr,'clon=clon(const);clat=clat(const)' -const,1,r3600x1800 "${WORK}/hp256_centres.nc"

echo "=== extracting from zarr: JJA SST and Apr-Oct Tmax, ${FIRST_YEAR}-${LAST_YEAR} ==="
parallel -j ${SAP_JOBS:-24} $BASEPY "$EXTRACT" {1} "$WORK" {2} ::: sst tmax ::: $YEARS 2>&1 \
    | grep -v -i warn

echo "=== SST anomalies ==="
cdo -s -O -setgrid,hp256_nested -mergetime $(for y in $YEARS; do echo "${WORK}/sst/sst_${y}.nc"; done) \
    "${WORK}/sap0006_jja_monthly.nc"
cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${WORK}/sap0006_jja_monthly.nc" "${WORK}/sap0006_clim.nc"
cdo -s -O -yearmean -ymonsub "${WORK}/sap0006_jja_monthly.nc" "${WORK}/sap0006_clim.nc" \
    "${WORK}/sap0006_anom_jja.nc"
cdo -s -O -fldmean "${WORK}/sap0006_anom_jja.nc" "${OUT_BASE}/sap0006_jja_gmsst.nc"
WGT=${WORK}/wgt_hp256_to_na025_${NA025_ID}.nc
# HEALPix has no stored cell corners, so cdo asks for --force on conservative weights
[ -s "$WGT" ] || cdo -s --force gencon,"$GRID_NA025" "${WORK}/sap0006_anom_jja.nc" "$WGT"
cdo -s -O -remap,"$GRID_NA025","$WGT" "${WORK}/sap0006_anom_jja.nc" "${OUT_BASE}/sap0006_jja_anom_na025.nc"
echo "  JJA seasons: $(cdo -s showyear "${OUT_BASE}/sap0006_jja_anom_na025.nc" | wc -w)"

echo "=== Tmax onto the Europe grid ==="
one_year() {
    y=$1; work=$2; tmax=$3; grid=$4
    out="${tmax}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    cdo -s -O -f nc4 -z zip_1 -setreftime,1850-01-01,00:00:00,1day -remapbil,"$grid" \
        -setgrid,hp256_nested "${work}/tmax_hp/tmax_${y}.nc" "$out"
}
export -f one_year
parallel -j 24 one_year {} "$WORK" "$TMAX" "$GRID_EU025" ::: $YEARS
echo "  $(ls ${TMAX}/tmax_*.nc | wc -l) years of daily Tmax"
echo "=== done ==="
