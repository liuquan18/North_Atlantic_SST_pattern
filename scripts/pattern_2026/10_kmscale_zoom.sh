#!/bin/bash
#SBATCH --job-name=na2026_zoom
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts/pattern_2026/logs/zoom.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# Supporting data for the km-scale figure: the 1991-2020 JJA mean SST of each
# dataset on a 0.05 deg zoom grid over two regions where ocean resolution is
# known to matter -- the Gulf Stream / North Atlantic Current separation, and
# the Mediterranean.
#
# Every dataset is put on the *same* fine grid, so what the reader sees is each
# model's own resolution rather than a plotting artefact: conservative
# remapping from a ~1 deg ocean grid onto 0.05 deg is piecewise constant and
# looks blocky, which is an honest rendering of what that model resolves.
#
# The climatologies are rebuilt here rather than reused, because the per-member
# MPI-GE intermediates are deleted by 02_mpige.sh once its anomalies are out.

set -euo pipefail
SCRIPT_DIR=${NA2026_SCRIPT_DIR:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts/pattern_2026}
source "${SCRIPT_DIR}/00_config.sh"

WORK=${WORK_BASE}/zoom
mkdir -p "$WORK" "${SCRIPT_DIR}/logs"
ICON_OCE_GRID=/pool/data/ICON/grids/public/mpim/0045/icon_grid_0045_R02B09_O.nc

# --- the two zoom grids ---------------------------------------------------
cat > "${WORK}/grid_gulfstream.txt" <<GEOF
gridtype = lonlat
xsize    = 800
ysize    = 300
xfirst   = -74.975
xinc     = 0.05
yfirst   = 33.025
yinc     = 0.05
GEOF
cat > "${WORK}/grid_medsea.txt" <<GEOF
gridtype = lonlat
xsize    = 640
ysize    = 240
xfirst   = -1.975
xinc     = 0.05
yfirst   = 33.025
yinc     = 0.05
GEOF

remap_to_zooms() {
    name=$1; infile=$2; setgrid_arg=${3:-}
    for z in gulfstream medsea; do
        out="${OUT_BASE}/zoom_${z}_${name}.nc"
        [ -s "$out" ] && { echo "    ${name}/${z} already done"; continue; }
        wgt="${WORK}/wgt_${name}_${z}.nc"
        if [ -n "$setgrid_arg" ]; then
            [ -s "$wgt" ] || cdo -s gencon,"${WORK}/grid_${z}.txt" -setgrid,"$setgrid_arg" "$infile" "$wgt"
            cdo -s -O -remap,"${WORK}/grid_${z}.txt","$wgt" -setgrid,"$setgrid_arg" "$infile" "$out"
        else
            [ -s "$wgt" ] || cdo -s gencon,"${WORK}/grid_${z}.txt" "$infile" "$wgt"
            cdo -s -O -remap,"${WORK}/grid_${z}.txt","$wgt" "$infile" "$out"
        fi
        echo "    ${name}/${z} -> $(basename "$out")"
    done
}

echo "=== ERA5 (0.25 deg) ==="
# JJA mean of the monthly climatology already built by 01_era5.sh
cdo -s -O -timmean "${WORK_BASE}/era5/era5_clim_jja.nc" "${WORK}/clim_era5.nc"
remap_to_zooms ERA5 "${WORK}/clim_era5.nc"

echo "=== MPI-ESM1-2-LR (~1 deg ocean, member r1) ==="
CMIP=/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/tos/gn
SCEN=/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r1i1p1f1/Omon/tos/gn
cdo -s -O -timmean -selmon,6/8 -selyear,${CLIM_START}/${CLIM_END} \
    -mergetime $(ls ${CMIP}/v*/tos_*.nc) $(ls ${SCEN}/v*/tos_*_201501-203412.nc) \
    "${WORK}/clim_mpige.nc"
remap_to_zooms MPI-GE "${WORK}/clim_mpige.nc"

echo "=== EERIE ICON-ESM-ER (5 km ocean, CMORized to 0.25 deg) ==="
cdo -s -O -timmean "${WORK_BASE}/eerie/eerie_hist_ssp245_clim.nc" "${WORK}/clim_eerie.nc"
remap_to_zooms EERIE "${WORK}/clim_eerie.nc"

echo "=== km-scale ICON / EPOC (5 km ocean, native unstructured) ==="
# built from the native-grid monthly files so the 5 km structure survives
EPOC=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
MASK=${WORK_BASE}/epoc/landmask_epoc2_020.nc
if [ ! -s "$MASK" ]; then echo "  ERROR: run 03_epoc_icon.sh first (need the land mask)" >&2; exit 1; fi
if [ ! -s "${WORK}/clim_epoc.nc" ]; then
    files=""
    for y in $(seq ${CLIM_START} ${CLIM_END}); do
        for m in 07 08 09; do   # stamps for June / July / August of year y
            f=$(ls ${EPOC}/run_*/epoc2_020_oce_2d_1mth_mean_${y}${m}01T000000Z.nc 2>/dev/null | head -1)
            [ -n "$f" ] && files="$files $f"
        done
    done
    # Two passes on purpose: chaining ~90 nested operator sub-chains into one
    # -ensmean call runs into cdo's operator limit, so extract the surface field
    # once per file and average the small results.
    echo "  extracting $(echo $files | wc -w) native monthly files"
    mkdir -p "${WORK}/epoc_parts"
    extract_one() {
        f=$1; mask=$2; outdir=$3
        b=$(basename "$f" .nc)
        [ -s "${outdir}/${b}.nc" ] && return 0
        cdo -s -O -ifthen "$mask" -sellevidx,1 -selname,to "$f" "${outdir}/${b}.nc"
    }
    export -f extract_one
    module load parallel 2>/dev/null || true
    parallel -j 8 extract_one {} "$MASK" "${WORK}/epoc_parts" ::: $files
    echo "  averaging"
    cdo -s -O -ensmean ${WORK}/epoc_parts/*.nc "${WORK}/clim_epoc.nc"
    rm -rf "${WORK}/epoc_parts"
fi
remap_to_zooms ICON-EPOC "${WORK}/clim_epoc.nc" "$ICON_OCE_GRID"

echo "=== done ==="
