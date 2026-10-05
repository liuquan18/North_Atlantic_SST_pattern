#!/bin/bash
#SBATCH --job-name=na2026_mpige
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/mpige.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --account=mh0033
#
# MPI-ESM1-2-LR grand ensemble (50 members), Omon/tos, historical + ssp245,
# JJA 1850-2100, reduced to the two common products (see ../config.sh).
#
# One node with GNU parallel, not multi-node MPI: the work is I/O-bound.
# Bilinear/conservative remap weights are generated once from the shared
# tripolar ocean grid and reused by all 50 members (identical grid), which is
# the main lever on runtime.

set -euo pipefail
# SLURM copies the batch script into /var/spool, so $0 is not the repo path.
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

WORK=${WORK_BASE}/mpige
mkdir -p "$WORK" "${PROJECT_ROOT}/logs"

CMIP=/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical
SCEN=/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245

WGT=${WORK}/wgt_mpige_to_na025.nc
if [ ! -s "$WGT" ]; then
    echo "=== generating remap weights once (shared tripolar grid) ==="
    sample=$(ls ${CMIP}/r1i1p1f1/Omon/tos/gn/v*/tos_*_185001-186912.nc | head -1)
    cdo -s gencon,"$GRID_NA025" -seltimestep,1 "$sample" "$WGT"
fi

process_member() {
    ens=$1
    source "$2/config.sh"
    WORK=${WORK_BASE}/mpige
    WGT=${WORK}/wgt_mpige_to_na025.nc
    CMIP=/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical
    SCEN=/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245

    mem=r${ens}i1p1f1
    reg_out=${WORK}/anom_na025_${mem}.nc
    gm_out=${WORK}/gmsst_${mem}.nc
    [ -s "$reg_out" ] && [ -s "$gm_out" ] && { echo "  ${mem} already done"; return 0; }

    hist=$(ls ${CMIP}/${mem}/Omon/tos/gn/v*/tos_*.nc 2>/dev/null)
    ssp=$(ls ${SCEN}/${mem}/Omon/tos/gn/v*/tos_*.nc 2>/dev/null)
    if [ -z "$hist" ] || [ -z "$ssp" ]; then echo "  ${mem}: MISSING INPUT" >&2; return 1; fi

    tmp=${WORK}/tmp_${mem}
    mkdir -p "$tmp"

    # JJA-only monthly series on the native ocean grid
    cdo -s -O -selmon,6/8 -mergetime $hist $ssp "${tmp}/jja_monthly.nc"
    # 1991-2020 JJA monthly climatology, native grid
    cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${tmp}/jja_monthly.nc" "${tmp}/clim.nc"
    # anomaly, then JJA seasonal mean (file holds only Jun/Jul/Aug, so yearmean == JJA mean)
    cdo -s -O -yearmean -ymonsub "${tmp}/jja_monthly.nc" "${tmp}/clim.nc" "${tmp}/anom_jja.nc"

    # global-ocean-mean JJA anomaly, taken on the native global field
    cdo -s -O -fldmean "${tmp}/anom_jja.nc" "$gm_out"
    # regional field on the common grid
    cdo -s -O -remap,"$GRID_NA025","$WGT" "${tmp}/anom_jja.nc" "$reg_out"

    rm -rf "$tmp"
    echo "  ${mem} done"
}
export -f process_member

echo "=== 50 members, 20-way parallel ==="
parallel -j 20 process_member {} "$SCRIPTS" ::: $(seq 1 50)

echo "=== stacking members ==="
activate_python
python3 "${SCRIPTS}/pre_process/stack_mpige.py"

echo "=== done ==="
