#!/bin/bash
#SBATCH --job-name=na2026_slp
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/slp.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --account=mh0033
#
# JJA mean sea-level pressure (hPa) over the North Atlantic-European sector
# (SLP_BOX, 20-80N, 90W-40E) on a 2.5 deg grid, for the NAO / EA EOFs of
# scripts/circulation/01_nao_ea_eof.py. Absolute JJA means, no anomalies: the
# ERA5 detrending and the MPI-GE ensemble-mean removal are done there.
#
#   ERA5    monthly MSL (code 151, reduced Gaussian N320), 1940-2026. The
#           archived monthly stream holds 2026 only to July, so 2026 comes from
#           the ERA5T monthly stream entirely (as in 09_era5_zg.sh).
#   MPI-GE  CMIP6 Amon/psl, 50 members, historical + ssp245, 1935-2030: the
#           10-year windows (y-5 to y+4) of 1940 and 2026 need the years around.
#
# Output: data/circulation/era5_jja_slp_2p5deg.nc (year, lat, lon) and
# mpige_jja_slp_2p5deg.nc (member, year, lat, lon), variable slp, from stack_slp.py.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

export SLP_BOX=-90,40,20,80
export MPIGE_Y0=1935 MPIGE_Y1=2030
export WORK=${WORK_BASE}/slp
export OUT_SLP=${PROJECT_ROOT}/data/circulation
mkdir -p "${WORK}/era5" "${WORK}/mpige" "$OUT_SLP" "${PROJECT_ROOT}/logs"

echo "=== 1/3 ERA5 JJA MSL, 1940-${REF_YEAR} ==="
era5_year() {
    year=$1
    out=${WORK}/era5/slp_jja_${year}.nc
    [ -s "$out" ] && return 0
    if [ "$year" -eq "$REF_YEAR" ]; then
        src=/pool/data/ERA5/ET/sf/an/1M/151/ETsf00_1M_${year}_151.grb
    else
        src=/pool/data/ERA5/E5/sf/an/1M/151/E5sf00_1M_${year}_151.grb
    fi
    [ "$(cdo -s ntime -selmon,6/8 "$src")" -eq 3 ] || { echo "ERROR: ${src} lacks JJA" >&2; return 1; }
    cdo -s -O -f nc -divc,100 -timmean -sellonlatbox,${SLP_BOX} -remapcon,r144x73 \
        -setgridtype,regular -selmon,6/8 "$src" "$out"
}
export -f era5_year
parallel -j 16 era5_year ::: $(seq 1940 "$REF_YEAR")
cdo -s -O -setattribute,slp@units=hPa -setname,slp -mergetime ${WORK}/era5/slp_jja_*.nc \
    "${WORK}/era5_jja_slp.nc"

echo "=== 2/3 MPI-GE JJA psl, 50 members, ${MPIGE_Y0}-${MPIGE_Y1} ==="
mpige_member() {
    mem=r$1i1p1f1
    out=${WORK}/mpige/slp_jja_${mem}.nc
    [ -s "$out" ] && return 0
    files=""
    for exp in CMIP/MPI-M/MPI-ESM1-2-LR/historical ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245; do
        for f in /pool/data/CMIP6/data/${exp}/${mem}/Amon/psl/gn/v*/psl_*.nc; do
            # keep the files that overlap the years (name ends _YYYYMM-YYYYMM.nc)
            span=${f##*_}; a=${span:0:4}; b=${span:7:4}
            [ "$b" -ge "$MPIGE_Y0" ] && [ "$a" -le "$MPIGE_Y1" ] && files="$files $f"
        done
    done
    # the file holds only Jun/Jul/Aug, so yearmean == JJA mean
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -yearmean \
        -sellonlatbox,${SLP_BOX} -remapcon,r144x73 -selyear,${MPIGE_Y0}/${MPIGE_Y1} \
        -selmon,6/8 -mergetime $files "$out"
}
export -f mpige_member
parallel -j 25 mpige_member ::: $(seq 1 50)

echo "=== 3/3 stack ==="
activate_python
python3 "${SCRIPTS}/pre_process/stack_slp.py"
echo "=== done ==="
