#!/bin/bash
#SBATCH --job-name=na2026_mzg
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/model_zg.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# JJA geopotential height anomalies at ZLEVEL hPa (default 300) of the runs that
# have upper-air output, YEARS (1991-2030), for the contours of
# fig25_composite_maps.py. Built like 09_era5_zg.sh: JJA monthly fields on the
# same 1 deg grid as data/sst_anomaly/era5_jja_z<ZLEVEL>_anom_1deg.nc (analysis
# box + 5 deg margin), anomalies against the 1991-2020 JJA climatology of the
# same months of the same member, JJA mean. The box mean is removed at plot time.
#
#   MPI-GE      CMIP6 Amon/zg, 50 members, historical + ssp245
#   ICON-EPOC   epoc2_020 atm_pl_6h_inst `gpsm` (geopotential above surface,
#               m2 s-2) on the native R02B08 grid 0054, not attached in the
#               files; the 6-hourly fields are averaged to months. The surface
#               geopotential is constant in time, so it drops out of the anomaly.
#               Years to 2025 (last complete JJA).
#   EERIE r1    CMOR Amon/zg (0.25 deg), hist-1950 + highres-future-ssp245
#   MPI-ER r2   ECHAM Amon3d code 156 (geopotential height, m) of ER3,
#               hist + ssp585 (the ssp585 files are also named *-hist_*).
#               Months stamped at their end, so -selmon,6/8 is JJA.
#
# No upper-air output exists on /work for EERIE r2/r3 after 1993 (6-hourly
# pressure levels end 1982 / 1993-09) or for MPI-ER r3 after 1989 (Amon3d of
# ER5), nor any for MPI-ER r1 beyond the 2014 Amon3d; those runs get SST only.
#
# Output: ${OUT_BASE}/<dataset>_jja_z<ZLEVEL>_anom_1deg.nc, variable z<ZLEVEL>
# (gpm) over (member, year, lat, lon), stacked by stack_zg.py.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

export ZLEVEL=${ZLEVEL:-300}
export Z=z${ZLEVEL}
export Y0=1991 Y1=2030
export G=9.80665
export WORK=${WORK_BASE}/model_${Z}
mkdir -p "${WORK}/epoc_month" "${PROJECT_ROOT}/logs"

ERA5_Z=${OUT_BASE}/era5_jja_${Z}_anom_1deg.nc
[ -s "$ERA5_Z" ] || { echo "ERROR: run ZLEVEL=${ZLEVEL} 09_era5_zg.sh first ($ERA5_Z)" >&2; exit 1; }
export GRID_Z=${WORK}/grid_${Z}_1deg.txt
cdo -s griddes "$ERA5_Z" > "$GRID_Z"

# monthly JJA series on GRID_Z ($1, any years) -> yearly JJA anomaly ($2)
anomaly() {
    cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "$1" "$1.clim"
    cdo -s -O -setattribute,${Z}@units=gpm -setname,${Z} -yearmean -ymonsub "$1" "$1.clim" "$2"
    rm -f "$1.clim"
}
export -f anomaly

echo "=== 1/5 MPI-GE, 50 members ==="
mpige_member() {
    mem=r$1i1p1f1
    out=${WORK}/mpige_${mem}_anom.nc
    [ -s "$out" ] && return 0
    files=""
    for exp in CMIP/MPI-M/MPI-ESM1-2-LR/historical ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245; do
        for f in /pool/data/CMIP6/data/${exp}/${mem}/Amon/zg/gn/v*/zg_*.nc; do
            # keep the files that overlap Y0..Y1 (name ends _YYYYMM-YYYYMM.nc)
            span=${f##*_}; a=${span:0:4}; b=${span:7:4}
            [ "$b" -ge "$Y0" ] && [ "$a" -le "$Y1" ] && files="$files $f"
        done
    done
    cdo -s -O -remapcon,"$GRID_Z" -sellevel,$((ZLEVEL * 100)) -selyear,${Y0}/${Y1} -selmon,6/8 \
        -mergetime $files "${WORK}/mpige_${mem}_mon.nc"
    anomaly "${WORK}/mpige_${mem}_mon.nc" "$out"
    rm -f "${WORK}/mpige_${mem}_mon.nc"
}
export -f mpige_member
parallel -j 25 mpige_member ::: $(seq 1 50)

echo "=== 2/5 EERIE r1 ==="
CMOR=/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER
if [ ! -s "${WORK}/eerie_r1_anom.nc" ]; then
    files=$(for y in $(seq $Y0 $Y1); do
        ls ${CMOR}/{hist-1950,highres-future-ssp245}/r1i1p1f1/Amon/zg/gr/v20240618/zg_*_${y}01-${y}12.nc 2>/dev/null || true
    done)
    eerie_year() { cdo -s -O -remapcon,"$GRID_Z" -sellevel,$((ZLEVEL * 100)) -selmon,6/8 "$1" "${WORK}/eerie_$(basename $1)"; }
    export -f eerie_year
    parallel -j 20 eerie_year ::: $files
    cdo -s -O -mergetime ${WORK}/eerie_zg_*.nc "${WORK}/eerie_r1_mon.nc"
    anomaly "${WORK}/eerie_r1_mon.nc" "${WORK}/eerie_r1_anom.nc"
    rm -f ${WORK}/eerie_zg_*.nc "${WORK}/eerie_r1_mon.nc"
fi

echo "=== 3/5 MPI-ER r2 (ER3) ==="
if [ ! -s "${WORK}/mpier_r2_anom.nc" ]; then
    ER=/work/uo0122/u241089/MPIESM
    files=$(for y in $(seq $Y0 $Y1); do for m in 06 07 08; do
        ls ${ER}/ER3-{hist,ssp585}/outdata/echam6/Amon3d_ER3-hist_${y}${m}.grb 2>/dev/null || true
    done; done)
    cdo -s -O -f nc -remapcon,"$GRID_Z" -sellevel,$((ZLEVEL * 100)) -selcode,156 \
        -mergetime $files "${WORK}/mpier_r2_mon.nc"
    anomaly "${WORK}/mpier_r2_mon.nc" "${WORK}/mpier_r2_anom.nc"
    rm -f "${WORK}/mpier_r2_mon.nc"
fi

echo "=== 4/5 ICON-EPOC, monthly means of the 6-hourly fields ==="
EPOC=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
export EPOC
export GRID_EPOC=/pool/data/ICON/grids/public/mpim/0054/icon_grid_0054_R02B08_G.nc
# last complete JJA of the still-running transient
EPOC_Y1=2025
export WGT_EPOC=${WORK}/wgt_epoc_${Z}_1deg.nc
if [ ! -s "${WORK}/epoc_r1_anom.nc" ]; then
    sample=$(ls ${EPOC}/run_20000601T*/epoc2_020_atm_pl_6h_inst_20000602T000000Z.nc)
    [ -s "$WGT_EPOC" ] || cdo -s gencon,"$GRID_Z" -setgrid,"$GRID_EPOC" \
        -seltimestep,1 -sellevel,$((ZLEVEL * 100)) -selname,gpsm "$sample" "$WGT_EPOC"
    epoc_month() {
        ym=$1
        out=${WORK}/epoc_month/${ym}.nc
        [ -s "$out" ] && return 0
        # a run_<YYYYMM> directory holds that month's daily 6-hourly files;
        # -selmon drops any step stamped outside the month
        cdo -s -O -divc,$G -remap,"$GRID_Z","$WGT_EPOC" -setgrid,"$GRID_EPOC" -monmean \
            -selmon,${ym:4:2} -sellevel,$((ZLEVEL * 100)) -selname,gpsm \
            -mergetime ${EPOC}/run_${ym}01T*/epoc2_020_atm_pl_6h_inst_*.nc "$out"
    }
    export -f epoc_month
    months=$(for y in $(seq $Y0 $EPOC_Y1); do echo ${y}06 ${y}07 ${y}08; done)
    parallel -j 24 epoc_month ::: $months
    cdo -s -O -mergetime ${WORK}/epoc_month/*.nc "${WORK}/epoc_r1_mon.nc"
    anomaly "${WORK}/epoc_r1_mon.nc" "${WORK}/epoc_r1_anom.nc"
fi

echo "=== 5/5 stack ==="
activate_python
python3 "${SCRIPTS}/pre_process/stack_zg.py"
echo "=== done ==="
