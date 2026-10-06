#!/bin/bash
#SBATCH --job-name=hw_eerie
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_eerie.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=mh0033
#
# EERIE ICON-ESM-ER (hist-1950 + highres-future-ssp245, r1i1p1f1) daily maximum
# 2 m temperature, April-October 1991-2050, Europe 0.25 deg.
#
# Source: CMORized day/tasmax on the 0.25 deg `gr` grid, one file per year. The
# grid nodes coincide with GRID_EU025, so the bilinear "remap" is an exact
# selection.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

BASE=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
WORK=${HW_WORK}/eerie_hist_ssp245
mkdir -p "$WORK"

FIRST_YEAR=${HW_FIRST_YEAR:-1991}
LAST_YEAR=${HW_LAST_YEAR:-2050}

one_year() {
    y=$1; work=$2; base=$3; grid=$4; months=$5
    out="${work}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    if [ "$y" -le 2014 ]; then exp=hist-1950; else exp=highres-future-ssp245; fi
    src=$(ls ${base}/${exp}/r1i1p1f1/day/tasmax/gr/v*/tasmax_day_ICON-ESM-ER_${exp}_r1i1p1f1_gr_${y}0101-${y}1231.nc)
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -remapbil,"$grid" -selmon,"$months" "$src" "$out"
    echo "  ${y}: $(cdo -s ntime "$out") days"
}
export -f one_year

parallel -j ${HW_JOBS:-30} one_year {} "$WORK" "$BASE" "$GRID_EU025" "$HW_MONTHS" ::: $(seq "$FIRST_YEAR" "$LAST_YEAR")
echo "=== done ==="
