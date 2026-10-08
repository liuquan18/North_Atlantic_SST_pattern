#!/bin/bash
#SBATCH --job-name=na2026_era5
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/era5.%j.out
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# ERA5 observed JJA SST anomalies, 1940-2026, for the North Atlantic + Mediterranean box.
#
# The DKRZ ERA5 pool's *monthly* stream (E5 .../1M/034) currently ends at
# 2026-06, so JJA 2026 cannot be built from it. The near-real-time ERA5T
# stream (ET .../1H/034) runs to 2026-09-16, so July and August 2026 are
# assembled from hourly ERA5T fields. This is not a compromise: for June 2026,
# where both streams exist, the hourly-derived monthly mean agrees with the
# archived monthly mean to 7e-4 K (GRIB packing precision), so 2026 is built
# entirely from ERA5T for internal consistency within the season.

set -euo pipefail
# SLURM copies the batch script into /var/spool, so $0 is not the repo path.
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"

WORK=${WORK_BASE}/era5
mkdir -p "$WORK" "${PROJECT_ROOT}/logs"

ERA5_MON=/pool/data/ERA5/E5/sf/an/1M/034
ERA5T_HOUR=/pool/data/ERA5/ET/sf/an/1H/034

FIRST_YEAR=1940
LAST_ARCHIVED_YEAR=2025

echo "=== 1/6 JJA monthly fields, ${FIRST_YEAR}-${LAST_ARCHIVED_YEAR} (archived ERA5 monthly) ==="
mkdir -p "${WORK}/monthly"
build_year() {
    year=$1; work=$2; src=$3
    out="${work}/monthly/era5_jja_${year}.nc"
    [ -s "$out" ] && return 0
    cdo -s -O -f nc -setgridtype,regular -selmon,6/8 "${src}/E5sf00_1M_${year}_034.grb" "$out"
}
export -f build_year
module load parallel 2>/dev/null || true
parallel -j 16 build_year {} "$WORK" "$ERA5_MON" ::: $(seq $FIRST_YEAR $LAST_ARCHIVED_YEAR)

echo "=== 2/6 JJA ${REF_YEAR} from hourly ERA5T ==="
out2026="${WORK}/monthly/era5_jja_${REF_YEAR}.nc"
if [ ! -s "$out2026" ]; then
    for m in 06 07 08; do
        nfiles=$(ls ${ERA5T_HOUR}/ETsf00_1H_${REF_YEAR}-${m}-*_034.grb 2>/dev/null | wc -l)
        echo "    ${REF_YEAR}-${m}: ${nfiles} hourly files"
        if [ "$nfiles" -lt 28 ]; then
            echo "ERROR: incomplete hourly coverage for ${REF_YEAR}-${m}" >&2; exit 1
        fi
        cdo -s -O -f nc -setgridtype,regular -monmean \
            -mergetime ${ERA5T_HOUR}/ETsf00_1H_${REF_YEAR}-${m}-*_034.grb \
            "${WORK}/monthly/era5t_${REF_YEAR}_${m}.nc"
    done
    cdo -s -O mergetime ${WORK}/monthly/era5t_${REF_YEAR}_0[678].nc "$out2026"
    rm -f ${WORK}/monthly/era5t_${REF_YEAR}_0[678].nc
fi

echo "=== 3/6 merge full JJA record and convert K -> degC ==="
cdo -s -O -subc,273.15 -mergetime ${WORK}/monthly/era5_jja_*.nc "${WORK}/era5_jja_monthly.nc"
cdo -s showyear "${WORK}/era5_jja_monthly.nc" | tr ' ' '\n' | grep -c . | xargs echo "    years:"

echo "=== 4/6 ${CLIM_START}-${CLIM_END} JJA climatology and anomalies (native grid) ==="
cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${WORK}/era5_jja_monthly.nc" "${WORK}/era5_clim_jja.nc"
cdo -s -O -ymonsub "${WORK}/era5_jja_monthly.nc" "${WORK}/era5_clim_jja.nc" "${WORK}/era5_anom_monthly.nc"

# JJA seasonal mean: the file holds only Jun/Jul/Aug, so yearmean == JJA mean
cdo -s -O -yearmean "${WORK}/era5_anom_monthly.nc" "${WORK}/era5_anom_jja.nc"

echo "=== 5/6 global-ocean-mean JJA anomaly (before regional subsetting) ==="
cdo -s -O -fldmean "${WORK}/era5_anom_jja.nc" "${OUT_BASE}/era5_jja_gmsst.nc"

echo "=== 6/6 remap onto the common 0.25 deg regional grid ==="
WGT=${WORK}/wgt_era5_to_na025_${NA025_ID}.nc
if [ ! -s "$WGT" ]; then
    cdo -s gencon,"$GRID_NA025" "${WORK}/era5_anom_jja.nc" "$WGT"
fi
cdo -s -O -remap,"$GRID_NA025","$WGT" "${WORK}/era5_anom_jja.nc" "${OUT_BASE}/era5_jja_anom_na025.nc"

# Also keep the raw (non-anomaly) JJA 2026 field for context maps
cdo -s -O -remap,"$GRID_NA025","$WGT" -yearmean -selyear,${REF_YEAR} "${WORK}/era5_jja_monthly.nc" \
    "${OUT_BASE}/era5_jja_${REF_YEAR}_absolute_na025.nc"

# ... and the individual June / July / August 2026 anomalies, to show whether
# the seasonal pattern is present all summer or is carried by one month
cdo -s -O -remap,"$GRID_NA025","$WGT" -selyear,${REF_YEAR} "${WORK}/era5_anom_monthly.nc" \
    "${OUT_BASE}/era5_${REF_YEAR}_monthly_anom_na025.nc"

echo "=== done ==="
ncdump -h "${OUT_BASE}/era5_jja_anom_na025.nc" | head -12
