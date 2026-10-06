#!/bin/bash
#SBATCH --job-name=hw_era5
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_era5.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# ERA5 daily maximum 2 m temperature, April-October 1991-2026, Europe 0.25 deg.
#
# As in the paper, Tmax is the maximum of the 24 hourly 2 m temperatures of
# each (UTC) day -- not the forecast-stream mx2t (param 201).
#
# The archived hourly stream (E5 .../1H/167) ends 2026-07-31; the rest of 2026
# comes from near-real-time ERA5T (ET .../1H/167, to 2026-09-30). Both are the
# same analysis on the same N320 grid, so a day is taken from E5 when it exists
# and from ET otherwise. 2026 therefore stops at 30 Sep, which is the end of
# the detection window; the reference years 1991-2020 are complete.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

E5=/pool/data/ERA5/E5/sf/an/1H/167
ET=/pool/data/ERA5/ET/sf/an/1H/167
WORK=${HW_WORK}/era5
mkdir -p "${WORK}/day"

FIRST_YEAR=${HW_FIRST_YEAR:-1991}
LAST_YEAR=${HW_LAST_YEAR:-${REF_YEAR}}

# bilinear weights N320 (regularised) -> Europe 0.25 deg, built once
WGT=${WORK}/wgt_n320_to_eu025.nc
if [ ! -s "$WGT" ]; then
    cdo -s genbil,"$GRID_EU025" -setgridtype,regular -seltimestep,1 \
        "${E5}/E5sf00_1H_2000-07-01_167.grb" "$WGT"
fi

one_day() {
    day=$1; work=$2; grid=$3; wgt=$4; e5=$5; et=$6
    out="${work}/day/tmax_${day}.nc"
    [ -s "$out" ] && return 0
    src="${e5}/E5sf00_1H_${day}_167.grb"
    [ -s "$src" ] || src="${et}/ETsf00_1H_${day}_167.grb"
    [ -s "$src" ] || return 0            # beyond the end of ERA5T
    n=$(cdo -s ntime "$src")
    if [ "$n" -ne 24 ]; then echo "WARN ${day}: ${n} hourly steps, skipped" >&2; return 0; fi
    cdo -s -O -f nc -remap,"$grid","$wgt" -daymax -setgridtype,regular "$src" "$out"
}
export -f one_day

days=()
for y in $(seq "$FIRST_YEAR" "$LAST_YEAR"); do
    [ -s "${WORK}/tmax_${y}.nc" ] && continue
    d=$(date -u -d "${y}-04-01" +%F)
    while [ "$d" != "${y}-11-01" ]; do days+=("$d"); d=$(date -u -d "$d +1 day" +%F); done
done
echo "=== daily Tmax for ${#days[@]} days ==="
[ "${#days[@]}" -gt 0 ] && parallel -j ${HW_JOBS:-96} one_day {} "$WORK" "$GRID_EU025" "$WGT" "$E5" "$ET" ::: "${days[@]}"

echo "=== merging per year ==="
merge_year() {
    y=$1; work=$2
    out="${work}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    ls ${work}/day/tmax_${y}-*.nc >/dev/null 2>&1 || return 0
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -mergetime ${work}/day/tmax_${y}-*.nc "$out"
    echo "  ${y}: $(cdo -s ntime "$out") days"
    rm -f ${work}/day/tmax_${y}-*.nc
}
export -f merge_year
parallel -j 16 merge_year {} "$WORK" ::: $(seq "$FIRST_YEAR" "$LAST_YEAR")
echo "=== done ==="
