#!/bin/bash
#SBATCH --job-name=hw_mpige
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_mpige.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# MPI-ESM1-2-LR grand ensemble (MPI-GE), 50 members, historical + ssp245:
# daily maximum 2 m temperature (CMIP6 day/tasmax), April-October 1850-2100.
#
# Unlike 01-03, MPI-GE stays on its NATIVE T63 Gaussian grid (1.875 deg),
# cropped to Europe. Regridding 50 members x 251 years of daily fields to
# 0.25 deg would be ~440 GB of scratch and only interpolate a ~200 km model;
# heatwaves are detected in the model's own grid cells instead. Land is the
# model's own: sftlf > 50 %.
#
# Output: ${HW_WORK}/mpige/tmax_r<e>.nc (degC), one file per member.
# The archive version directory differs between members (v20190710,
# v20190815, v20210901, ...), so it is globbed, never hard-coded.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

CMIP=/pool/data/CMIP6/data
HIST=${CMIP}/CMIP/MPI-M/MPI-ESM1-2-LR/historical
SSP=${CMIP}/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245
WORK=${HW_WORK}/mpige
mkdir -p "${WORK}/part"

# a little wider than the 0.25 deg domain so the 1.875 deg cells cover its edges
BOX=-17,47,28,74
MEMBERS=$(seq ${HW_FIRST_MEMBER:-1} ${HW_LAST_MEMBER:-50})

LANDMASK_MPIGE=${HW_OUT}/landmask_mpige_native.nc
if [ ! -s "$LANDMASK_MPIGE" ]; then
    cdo -s -O -f nc -setctomiss,0 -gtc,50 -sellonlatbox,$BOX \
        $(ls ${HIST}/r1i1p1f1/fx/sftlf/gn/v*/sftlf_fx_*.nc | head -1) "$LANDMASK_MPIGE"
fi

one_file() {
    src=$1; work=$2; box=$3; months=$4
    mem=$(basename "$src" | grep -o 'r[0-9]*i1p1f1' | sed 's/i1p1f1//')
    span=$(basename "$src" .nc | awk -F_ '{print $NF}')
    out="${work}/part/tmax_${mem}_${span}.nc"
    [ -s "$out" ] && return 0
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -selmon,"$months" -sellonlatbox,"$box" "$src" "$out"
}
export -f one_file

files=()
for e in $MEMBERS; do
    [ -s "${WORK}/tmax_r${e}.nc" ] && continue
    for base in "$HIST" "$SSP"; do
        files+=($(ls ${base}/r${e}i1p1f1/day/tasmax/gn/v*/tasmax_day_*.nc))
    done
done
echo "=== cropping ${#files[@]} source files ==="
[ "${#files[@]}" -gt 0 ] && parallel -j ${HW_JOBS:-64} one_file {} "${WORK}" "$BOX" "$HW_MONTHS" ::: "${files[@]}"

echo "=== merging per member ==="
merge_member() {
    e=$1; work=$2
    out="${work}/tmax_r${e}.nc"
    [ -s "$out" ] && return 0
    n=$(ls ${work}/part/tmax_r${e}_*.nc | wc -l)
    [ "$n" -eq 14 ] || { echo "  r${e}: ${n} of 14 parts, skipped" >&2; return 0; }
    cdo -s -O -f nc4 -z zip_1 -mergetime ${work}/part/tmax_r${e}_*.nc "$out"
    echo "  r${e}: $(cdo -s ntime "$out") days"
    rm -f ${work}/part/tmax_r${e}_*.nc
}
export -f merge_member
parallel -j 25 merge_member {} "$WORK" ::: $MEMBERS
echo "=== done ==="
