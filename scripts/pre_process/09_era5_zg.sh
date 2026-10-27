#!/bin/bash
#SBATCH --job-name=na2026_zg
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/era5_zg.%j.out
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --account=mh0033
#
# ERA5 JJA geopotential height anomalies on one pressure level (ZLEVEL hPa,
# default 500; `ZLEVEL=300 bash 09_era5_zg.sh` for Z300), 1940-2026, for the
# contours on the SST maps of fig1_patterns_era5.py. Built like the SST of 01_era5.sh:
# JJA monthly fields, anomalies against the 1991-2020 JJA climatology of the
# same months, JJA mean. The box mean is removed at plot time, as for the SST.
#
# Source: the monthly pressure-level geopotential (code 129, m2 s-2, reduced
# Gaussian N320, 37 levels; GRIB-indexed, so pulling one level is ~1 s/year).
# The archived ERA5 monthly stream ends at 2026-07; the ERA5T monthly stream
# (ET .../pl/an/1M/129) holds all of JJA 2026, so 2026 is taken from ERA5T
# entirely, as the SST is.
#
# Output: ${OUT_BASE}/era5_jja_z<ZLEVEL>_anom_1deg.nc, variable z<ZLEVEL> (gpm),
# on a 1 deg grid over the analysis box plus a 5 deg margin (contours reach the
# map edge). 1 deg is ample for a field as smooth as upper-air geopotential.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load cdo/2.5.0-gcc-11.2.0 2>/dev/null || true

ZLEVEL=${ZLEVEL:-500}
Z=z${ZLEVEL}
WORK=${WORK_BASE}/era5_${Z}
mkdir -p "${WORK}/monthly"
E5_MON=/pool/data/ERA5/E5/pl/an/1M/129
ET_MON=/pool/data/ERA5/ET/pl/an/1M/129
FIRST_YEAR=1940
LAST_ARCHIVED_YEAR=$((REF_YEAR - 1))
G=9.80665

echo "=== 1/4 JJA monthly Z${ZLEVEL}, ${FIRST_YEAR}-${REF_YEAR} ==="
build_year() {
    year=$1; work=$2; src=$3; stem=$4; g=$5; lev=$6
    out="${work}/monthly/z${lev}_jja_${year}.nc"
    [ -s "$out" ] && return 0
    cdo -s -O -f nc -divc,"$g" -setgridtype,regular -sellevel,$((lev * 100)) -selmon,6/8 \
        "${src}/${stem}_1M_${year}_129.grb" "$out"
}
export -f build_year
module load parallel 2>/dev/null || true
parallel -j 8 build_year {} "$WORK" "$E5_MON" E5pl00 "$G" "$ZLEVEL" ::: $(seq $FIRST_YEAR $LAST_ARCHIVED_YEAR)
build_year "$REF_YEAR" "$WORK" "$ET_MON" ETpl00 "$G" "$ZLEVEL"
n=$(cdo -s ntime "${WORK}/monthly/${Z}_jja_${REF_YEAR}.nc")
[ "$n" -eq 3 ] || { echo "ERROR: ERA5T JJA ${REF_YEAR} has ${n} months" >&2; exit 1; }

echo "=== 2/4 merge ==="
cdo -s -O -setattribute,${Z}@units=gpm -setname,${Z} -mergetime \
    ${WORK}/monthly/${Z}_jja_*.nc "${WORK}/${Z}_jja_monthly.nc"
cdo -s showyear "${WORK}/${Z}_jja_monthly.nc" | tr ' ' '\n' | grep -c . | xargs echo "    years:"

echo "=== 3/4 ${CLIM_START}-${CLIM_END} JJA climatology, anomalies, JJA mean ==="
cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${WORK}/${Z}_jja_monthly.nc" "${WORK}/${Z}_clim_jja.nc"
# the file holds only Jun/Jul/Aug, so yearmean == JJA mean
cdo -s -O -yearmean -ymonsub "${WORK}/${Z}_jja_monthly.nc" "${WORK}/${Z}_clim_jja.nc" \
    "${WORK}/${Z}_anom_jja.nc"

echo "=== 4/4 1 deg regional grid ==="
cdo -s -O -sellonlatbox,$((LON_W - 5)),$((LON_E + 5)),$((LAT_S - 5)),$((LAT_N + 5)) \
    -remapcon,r360x180 "${WORK}/${Z}_anom_jja.nc" "${OUT_BASE}/era5_jja_${Z}_anom_1deg.nc"
cdo -s sinfon "${OUT_BASE}/era5_jja_${Z}_anom_1deg.nc" | head -8
