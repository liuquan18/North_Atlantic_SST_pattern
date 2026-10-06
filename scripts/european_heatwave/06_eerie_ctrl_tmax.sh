#!/bin/bash
#SBATCH --job-name=hw_eerie_ctrl
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_eerie_ctrl.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=mh0033
#
# EERIE ICON-ESM-ER control (eerie-control-1950, constant 1950 forcing) daily
# maximum 2 m temperature, April-October 1950-2050, Europe 0.25 deg.
#
# Source: CMORized day/tasmax on the 0.25 deg `gr` grid, one file per year;
# its nodes coincide with GRID_EU025, so the bilinear "remap" is an exact
# selection (as for the forced run, 03_eerie_tmax.sh).
#
# Output: ${HW_WORK}/eerie_control/tmax_<year>.nc

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

BASE=/pool/data/EERIE/EERIE/MPI-M/ICON-ESM-ER/eerie-control-1950/r1i1p1f1/day/tasmax/gr
WORK=${HW_WORK}/eerie_control
mkdir -p "$WORK"

one_year() {
    y=$1; work=$2; base=$3; grid=$4; months=$5
    out="${work}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    src=$(ls ${base}/v*/tasmax_day_ICON-ESM-ER_eerie-control-1950_r1i1p1f1_gr_${y}0101-${y}1231.nc)
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -remapbil,"$grid" -selmon,"$months" "$src" "$out"
}
export -f one_year

parallel -j ${HW_JOBS:-30} one_year {} "$WORK" "$BASE" "$GRID_EU025" "$HW_MONTHS" ::: $(seq 1950 2050)
echo "  $(ls ${WORK}/tmax_*.nc | wc -l) years of daily Tmax"
echo "=== done ==="
