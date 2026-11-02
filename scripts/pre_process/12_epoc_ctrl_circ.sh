#!/bin/bash
#SBATCH --job-name=na2026_ectrl
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/epoc_ctrl_circ.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --account=mh0033
#
# JJA Z300 and SLP of the EPOC control epoc2_010 (constant 1990 forcing,
# 1990-2024), for fig 28. Same grids and products as 10_model_zg.sh and 11_slp.sh.
#
# The control has no pressure-level output, so Z300 comes from the monthly
# model-level stream atm_ml_1mth_mean: gpsm and pfull are remapped on the levels
# around 300 hPa (MLEVS; 18.9-42 kPa in the box in July 2000) and interpolated
# in log pressure (zg_ml2pl.py). A monthly mean interpolated this way is not the
# monthly mean of pressure-level values, so the same is done for the transient
# epoc2_020 over 1991-2020 as a check against its 6-hourly pl stream
# (ml_check_epoc2_020_anom.nc, not stacked).
#
# Anomalies against 1991-2020 of the run (Z300; SLP at plot time, as 11_slp.sh).
# ~5 GB per cdo process while the grid is attached: JOBS sets the parallelism.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load cdo/2.5.0-gcc-11.2.0 parallel 2>/dev/null || true
activate_python

export Z=z300 ZLEVEL=300 MLEVS=52/62 SCRIPTS
export WORK=${WORK_BASE}/model_${Z} SLPWORK=${WORK_BASE}/slp
export GRID_EPOC=/pool/data/ICON/grids/public/mpim/0054/icon_grid_0054_R02B08_G.nc
export GRID_Z=${WORK}/grid_${Z}_1deg.txt WGT_EPOC=${WORK}/wgt_epoc_${Z}_1deg.nc
export WGT_SLP=${SLPWORK}/wgt_epoc_r144x73.nc SLP_BOX=-90,40,20,80
CTRL=/work/bm1313/b383127/epoc-icon-2024.10/experiments/epoc2_010/work
HIST=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
mkdir -p "${WORK}/ml" "${SLPWORK}/epocctrl"
JOBS=${JOBS:-24}
for f in "$GRID_Z" "$WGT_EPOC" "$WGT_SLP"; do
    [ -s "$f" ] || { echo "ERROR: run 10_model_zg.sh and 11_slp.sh first ($f)" >&2; exit 1; }
done

ml_month() {
    exp=$1; src=$2; ym=$3
    out=${WORK}/ml/${exp}_${ym}.nc
    [ -s "$out" ] && return 0
    cdo -s -O -remap,"$GRID_Z","$WGT_EPOC" -setgrid,"$GRID_EPOC" -sellevidx,"$MLEVS" \
        -selname,gpsm,pfull ${src}/run_${ym}01T*/${exp}_atm_ml_1mth_mean_*.nc "${out}.tmp"
    python3 "${SCRIPTS}/pre_process/zg_ml2pl.py" "${out}.tmp" "$out" "$ZLEVEL"
    rm -f "${out}.tmp"
}
export -f ml_month

anomaly() {
    cdo -s -O -mergetime ${WORK}/ml/$1_*.nc "${WORK}/ml/$1.mon"
    cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "${WORK}/ml/$1.mon" "${WORK}/ml/$1.clim"
    cdo -s -O -setattribute,${Z}@units=gpm -setname,${Z} -yearmean \
        -ymonsub "${WORK}/ml/$1.mon" "${WORK}/ml/$1.clim" "$2"
    rm -f "${WORK}/ml/$1.mon" "${WORK}/ml/$1.clim"
}

months() { for y in $(seq "$1" "$2"); do echo ${y}06 ${y}07 ${y}08; done; }

echo "=== 1/3 Z300 from model levels: control 1990-2024, transient 1991-2020 ==="
parallel -j "$JOBS" ml_month epoc2_010 "$CTRL" {} ::: $(months 1990 2024)
parallel -j "$JOBS" ml_month epoc2_020 "$HIST" {} ::: $(months 1991 2020)
anomaly epoc2_010 "${WORK}/epocctrl_r1_anom.nc"
anomaly epoc2_020 "${WORK}/ml_check_epoc2_020_anom.nc"

echo "=== 2/3 control SLP ==="
slp_year() {
    y=$1; src=$2
    out=${SLPWORK}/epocctrl/slp_jja_${y}.nc
    [ -s "$out" ] && return 0
    cdo -s -O -setattribute,slp@units=hPa -setname,slp -divc,100 -sellonlatbox,${SLP_BOX} \
        -remap,r144x73,"$WGT_SLP" -setgrid,"$GRID_EPOC" -timmean -selname,psl \
        -mergetime ${src}/run_${y}0{6,7,8}01T*/epoc2_010_atm_2d_1mth_mean_*.nc "$out"
}
export -f slp_year
parallel -j "$JOBS" slp_year {} "$CTRL" ::: $(seq 1990 2024)
cdo -s -O -mergetime ${SLPWORK}/epocctrl/slp_jja_[12]*.nc "${SLPWORK}/epocctrl/slp_jja_r1.nc"

echo "=== 3/3 stack ==="
python3 "${SCRIPTS}/pre_process/stack_zg.py" epocctrl
python3 "${SCRIPTS}/pre_process/stack_slp.py" epocctrl
echo "=== done ==="
