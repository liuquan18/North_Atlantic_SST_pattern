#!/bin/bash
#SBATCH --job-name=na2026_epoc
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/epoc.%j.out
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --account=mh0033
#
# km-scale ICON (EPOC): 10 km atmosphere / 5 km ocean coupled runs.
#   epoc2_010  control,              1990-2024
#   epoc2_020  historical + ssp585,  1990-2026-03 (still running)  <- the "ICON-historical" run
#
# Three things about this output are not like the CMIP archives:
#
#  1. The ocean field is `to` (sea water potential temperature) in the
#     oce_2d_1mth_mean stream, on the *unstructured* ICON R02B09 ocean grid
#     (14,914,033 cells, grid 0045, uuid b7ad8f68-...). The grid is not in the
#     data files; it has to be attached from the grid repository.
#  2. Land cells on that ocean grid are written as exactly 0.0 with no
#     _FillValue declared. Left alone they bleed into every coastal target
#     cell during remapping -- which would wreck the Mediterranean, half of
#     this study's region. The zero set is bit-identical between 1990 and 2020,
#     so a static mask built once from a reference file is both correct and
#     exact.
#  3. Monthly-mean files are stamped with the END of their averaging period:
#     ..._19900801T000000Z.nc holds the July 1990 mean. The time axis is reset
#     to the middle of the month actually represented.
#
# Because the land mask is static, conservative remapping is linear in the
# data, so the climatology and anomalies are computed on the small remapped
# regional fields rather than on 6 GB of native-grid intermediates.

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

WORK=${WORK_BASE}/epoc
mkdir -p "$WORK" "${PROJECT_ROOT}/logs"

ICON_OCE_GRID=/pool/data/ICON/grids/public/mpim/0045/icon_grid_0045_R02B09_O.nc

declare -A EXP_DIR=(
  [epoc2_010]=/work/bm1313/b383127/epoc-icon-2024.10/experiments/epoc2_010/work
  [epoc2_020]=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
)

# ---- one monthly file -> (regional 0.25 deg field, global ocean mean) -------
process_file() {
    infile=$1; exp=$2; work=$3; grid=$4; gridfile=$5; wgt=$6; mask=$7

    stamp=$(basename "$infile" | sed -E 's/.*_([0-9]{8})T[0-9]+Z\.nc/\1/')
    # monthly means are stamped at the END of the period: subtract a day to get
    # the month actually represented
    rep=$(date -u -d "${stamp} -1 day" +%Y-%m)
    year=${rep%-*}; month=${rep#*-}

    reg="${work}/reg/${exp}_${year}${month}.nc"
    gm="${work}/gm/${exp}_${year}${month}.nc"
    [ -s "$reg" ] && [ -s "$gm" ] && return 0

    tmp=$(mktemp -d "${work}/tmp_XXXXXX")
    # surface layer, land zeros removed via the static mask, grid attached
    cdo -s -O -setgrid,"$gridfile" -ifthen "$mask" -sellevidx,1 -selname,to "$infile" "${tmp}/native.nc"
    cdo -s -O -settaxis,"${year}-${month}-15",12:00:00,1day -remap,"$grid","$wgt" "${tmp}/native.nc" "$reg"
    cdo -s -O -settaxis,"${year}-${month}-15",12:00:00,1day -fldmean "${tmp}/native.nc" "$gm"
    rm -rf "$tmp"
}
export -f process_file

for exp in epoc2_010 epoc2_020; do
    src=${EXP_DIR[$exp]}
    echo "==================== ${exp} ===================="
    mkdir -p "${WORK}/reg" "${WORK}/gm"

    # JJA monthly means live in files stamped July / August / September
    mapfile -t files < <(ls ${src}/run_*/${exp}_oce_2d_1mth_mean_*.nc 2>/dev/null \
        | awk -F_ '{d=$NF; sub(/T.*/,"",d); m=substr(d,5,2); if (m=="07"||m=="08"||m=="09") print}' \
        | sort -u)
    echo "  JJA monthly files found: ${#files[@]}"
    [ "${#files[@]}" -eq 0 ] && { echo "  no files, skipping" >&2; continue; }

    # static land mask + remap weights, built once per experiment
    MASK=${WORK}/landmask_${exp}.nc
    WGT=${WORK}/wgt_${exp}_to_na025.nc
    ref=${files[0]}
    if [ ! -s "$MASK" ]; then
        echo "  building static land mask (exact zeros -> missing)"
        cdo -s -O -gtc,-100 -setctomiss,0 -sellevidx,1 -selname,to "$ref" "$MASK"
    fi
    if [ ! -s "$WGT" ]; then
        echo "  generating conservative remap weights (once, reused for all months)"
        cdo -s gencon,"$GRID_NA025" -setgrid,"$ICON_OCE_GRID" -ifthen "$MASK" \
            -sellevidx,1 -selname,to "$ref" "$WGT"
    fi
    # measured peak is ~4.3 GB per cdo process while the 14.9M-cell grid is
    # attached, so 8 workers fit comfortably in the 48 GB requested above
    parallel -j ${NA2026_JOBS:-8} process_file {} "$exp" "$WORK" "$GRID_NA025" "$ICON_OCE_GRID" "$WGT" "$MASK" ::: "${files[@]}"

    echo "  merging and building anomalies"
    # setreftime puts the merged axis on a single day-based reference, which
    # also repairs per-file axes written by an earlier "1mon" run
    cdo -s -O -setreftime,1850-01-01,00:00:00,1day \
        -mergetime ${WORK}/reg/${exp}_*.nc "${WORK}/${exp}_jja_monthly_na025.nc"
    cdo -s -O -setreftime,1850-01-01,00:00:00,1day \
        -mergetime ${WORK}/gm/${exp}_*.nc  "${WORK}/${exp}_jja_monthly_gm.nc"

    for kind in na025 gm; do
        in=${WORK}/${exp}_jja_monthly_${kind}.nc
        cdo -s -O -ymonmean -selyear,${CLIM_START}/${CLIM_END} "$in" "${WORK}/${exp}_clim_${kind}.nc"
        # drop incomplete seasons: keep only years with all three JJA months
        cdo -s -O -yearmean -ymonsub "$in" "${WORK}/${exp}_clim_${kind}.nc" "${WORK}/${exp}_anom_${kind}_allyears.nc"
    done

    activate_python
    EXP=$exp WORK=$WORK OUT=$OUT_BASE python3 - <<'PYEOF'
import os, xarray as xr, numpy as np
exp, work, out = os.environ["EXP"], os.environ["WORK"], os.environ["OUT"]
# yearmean silently averages whatever months are present; keep only complete JJA
mon = xr.open_dataset(f"{work}/{exp}_jja_monthly_na025.nc")
# cdo's yearmean averages whatever months happen to be present, so a year that
# is missing e.g. August would silently become a Jun-Jul mean. Keep only the
# years where all three JJA months are there.
counts = mon["time"].groupby("time.year").count()
complete = counts.where(counts == 3, drop=True)["year"].values
print(f"  {exp}: complete JJA seasons: {len(complete)} ({complete.min()}-{complete.max()})")
for kind, fname in [("na025", f"{exp}_jja_anom_na025.nc"), ("gm", f"{exp}_jja_gmsst.nc")]:
    ds = xr.open_dataset(f"{work}/{exp}_anom_{kind}_allyears.nc")
    ds = ds.sel(time=ds["time.year"].isin(complete))
    ds.to_netcdf(f"{out}/{fname}")
    print(f"    -> {out}/{fname}", dict(ds.sizes))
PYEOF
done
echo "=== done ==="
