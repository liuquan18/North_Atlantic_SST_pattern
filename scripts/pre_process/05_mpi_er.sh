#!/bin/bash
#SBATCH --job-name=na2026_mpier
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/mpier.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --account=mh0033
#
# MPI-ESM1.2-ER: T127 atmosphere (~100 km) coupled to the eddy-resolving MPIOM
# TP6M ocean (3602 x 2394 curvilinear, ~0.1 deg / ~10 km).
#
# Three realizations, each historical (1950-2014) + ssp585 (2015-2099):
#   ER   -> member 1 (branched from the control at 1980)
#   ER3  -> member 2 (branched at 1996)
#   ER5  -> member 3 (branched at 2040)
# Membership and branch years are from /work/uo0122/u241089/MPIESM/readme.
#
# Notes on this output, compared with the other sources here:
#  - The SST variable is lowercase `tos` (the readme says TOS), in
#    `*_mpiom_data_2d_mm_*.nc`, with a proper _FillValue -- no land-zero trap
#    like the EPOC ICON output, and the mask is bit-identical between 1950 and
#    2090, so conservative remapping stays linear in the data.
#  - Monthly means are stamped at the END of the month they represent
#    (1950-06-30 23:55 is the June mean), so unlike EPOC no time-axis shift is
#    needed -- `-selmon,6/8` already selects JJA correctly.
#  - Files are one year each and ~2.2 GB (deflate-compressed, ~14 2-D fields),
#    so the cost is decompression: ~10 s to pull JJA `tos` out of a year. The
#    extracted field is written once and then used for both the remap and the
#    global mean rather than decompressing the source twice.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

WORK=${WORK_BASE}/mpier
mkdir -p "$WORK/reg" "$WORK/gm" "${PROJECT_ROOT}/logs"

SRC=/work/uo0122/u241089/MPIESM
declare -A MEMBER=( [ER]=1 [ER3]=2 [ER5]=3 )

WGT=${WORK}/wgt_mpier_to_na025.nc
if [ ! -s "$WGT" ]; then
    echo "=== generating remap weights once (grid shared by all members) ==="
    sample=$(ls ${SRC}/ER-hist/outdata/mpiom/ER-hist_mpiom_data_2d_mm_19500101_19501231.nc)
    cdo -s gencon,"$GRID_NA025" -seltimestep,1 -selname,tos "$sample" "$WGT"
fi

process_file() {
    infile=$1; run=$2; work=$3; grid=$4; wgt=$5
    year=$(basename "$infile" | sed -E 's/.*_2d_mm_([0-9]{4})[0-9]{4}_.*/\1/')
    reg="${work}/reg/${run}_${year}.nc"
    gm="${work}/gm/${run}_${year}.nc"
    [ -s "$reg" ] && [ -s "$gm" ] && return 0

    tmp=$(mktemp -d "${work}/tmp_XXXXXX")
    # decompress once; both products are built from this
    cdo -s -O -selmon,6/8 -selname,tos "$infile" "${tmp}/jja.nc"
    cdo -s -O -remap,"$grid","$wgt" "${tmp}/jja.nc" "$reg"
    cdo -s -O -fldmean "${tmp}/jja.nc" "$gm"
    rm -rf "$tmp"
}
export -f process_file

for run in ER ER3 ER5; do
    files=$(ls ${SRC}/${run}-hist/outdata/mpiom/${run}-hist_mpiom_data_2d_mm_*.nc \
               ${SRC}/${run}-ssp585/outdata/mpiom/${run}-ssp585_mpiom_data_2d_mm_*.nc 2>/dev/null)
    n=$(echo "$files" | wc -w)
    echo "==================== ${run} (member ${MEMBER[$run]}): ${n} year-files ===================="
    parallel -j 12 process_file {} "$run" "$WORK" "$GRID_NA025" "$WGT" ::: $files

    echo "  merging and building anomalies"
    for kind in reg gm; do
        out=${WORK}/${run}_jja_monthly_${kind}.nc
        cdo -s -O -setreftime,1850-01-01,00:00:00,1day -mergetime ${WORK}/${kind}/${run}_*.nc "$out"
        cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "$out" "${WORK}/${run}_clim_${kind}.nc"
        cdo -s -O -yearmean -ymonsub "$out" "${WORK}/${run}_clim_${kind}.nc" \
            "${WORK}/${run}_anom_${kind}.nc"
    done
done

echo "=== stacking the three realizations ==="
activate_python
python3 "${SCRIPTS}/pre_process/stack_mpi_er.py"
echo "=== done ==="
