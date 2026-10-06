#!/bin/bash
#SBATCH --job-name=hw_eerie
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_eerie.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# EERIE ICON-ESM-ER (hist-1950 + highres-future-ssp245) daily maximum 2 m
# temperature, April-October, Europe 0.25 deg, for the three realizations:
#
#   r1  erc2020  CMORized day/tasmax (`gr`), one file per year, 1991-2050
#   r2  erc2023  raw atm_2d_1d_max_remap025 (variable `tas`), 1991-2020
#   r3  erc2024  same, 1991-2020
#
# r2/r3 are the raw output remapped to 0.25 deg by DKRZ (see
# pre_process/04_eerie.sh): one `run_<YYYYMM>...` directory per month, daily
# values stamped at the END of the day (06-02 00:00 is the 1 June maximum),
# hence -shifttime,-1day before selecting the months.
#
# All three grids have nodes on multiples of 0.25 deg, so the bilinear
# "remap" onto GRID_EU025 is an exact selection.
#
# Output: ${HW_WORK}/eerie_hist_ssp245/r<m>/tmax_<year>.nc

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/european_heatwave/hw_config.sh"
module load parallel 2>/dev/null || true

CMOR=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
RAW=/work/bm1344/k202193/ICON
WORK=${HW_WORK}/eerie_hist_ssp245
mkdir -p "$WORK"/r{1,2,3}

FIRST_YEAR=${HW_FIRST_YEAR:-1991}

cmor_year() {
    y=$1; work=$2; base=$3; grid=$4; months=$5
    out="${work}/r1/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    if [ "$y" -le 2014 ]; then exp=hist-1950; else exp=highres-future-ssp245; fi
    src=$(ls ${base}/${exp}/r1i1p1f1/day/tasmax/gr/v*/tasmax_day_ICON-ESM-ER_${exp}_r1i1p1f1_gr_${y}0101-${y}1231.nc)
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -remapbil,"$grid" -selmon,"$months" "$src" "$out"
    echo "  r1 ${y}: $(cdo -s ntime "$out") days"
}

raw_year() {
    run=$1; m=$2; y=$3; work=$4; raw=$5; grid=$6; months=$7
    out="${work}/r${m}/tmax_${y}.nc"
    [ -s "$out" ] && return 0
    base=${raw}/${run}/postprocessing/interpolation
    [ "$y" -ge 2015 ] && base=${base}/SSP245
    m0=${months%/*}; m1=${months#*/}
    files=""
    for mm in $(seq "$m0" "$m1"); do
        files+=" $(ls ${base}/atm_2d_1d_max_remap025/run_${y}$(printf %02d $mm)01T*/*.nc)"
    done
    cdo -s -O -f nc4 -z zip_1 -setname,tmax -setunit,degC -subc,273.15 \
        -remapbil,"$grid" -selmon,"$months" -shifttime,-1day -mergetime \
        $(for f in $files; do echo "-selname,tas $f"; done) "$out"
    echo "  r${m} ${y}: $(cdo -s ntime "$out") days"
}
export -f cmor_year raw_year

parallel -j ${HW_JOBS:-30} cmor_year {} "$WORK" "$CMOR" "$GRID_EU025" "$HW_MONTHS" \
    ::: $(seq "$FIRST_YEAR" 2050)
parallel -j ${HW_JOBS:-30} raw_year {1} {2} {3} "$WORK" "$RAW" "$GRID_EU025" "$HW_MONTHS" \
    ::: erc2023 ::: 2 ::: $(seq "$FIRST_YEAR" 2020)
parallel -j ${HW_JOBS:-30} raw_year {1} {2} {3} "$WORK" "$RAW" "$GRID_EU025" "$HW_MONTHS" \
    ::: erc2024 ::: 3 ::: $(seq "$FIRST_YEAR" 2020)
echo "=== done ==="
