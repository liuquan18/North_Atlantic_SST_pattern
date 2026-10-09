#!/bin/bash
#SBATCH --job-name=na2026_nino34
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/nino34.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --account=mh0033
#
# JJA ENSO index of every dataset in the heatwave scatter figures (fig9):
#
#   n34   Nino-3.4 box mean SST anomaly, 5S-5N, 170W-120W
#   trop  tropical-mean SST anomaly, 20S-20N, all longitudes
#   rn34  relative Nino-3.4 = n34 - trop (van Oldenborgh et al. 2021; the
#         index NOAA CPC's relative ONI is built on). Subtracting the tropical
#         mean removes the forced warming, so in the scenario runs the index
#         does not share a trend with the heatwave days.
#
# Anomalies are JJA means relative to the 1991-2020 JJA monthly climatology of
# the same dataset (as everywhere in this study); the box means are area
# weighted on each dataset's own grid -- except MPI-GE, whose tripolar grid
# cdo cannot cut a lon/lat box out of (sellonlatbox keeps an index rectangle),
# so it goes through a conservative remap to 1 deg first.
#
# Sources:
#   ERA5, EERIE r1-r3, EERIE control, sap0006: the native-grid JJA anomaly
#     files 01_era5.sh, 04_eerie.sh and 07_sap0006.sh leave in $WORK_BASE.
#     Scratch is purged -- rerun those stages if they are gone.
#   MPI-GE, EPOC hist + control: re-read from the permanent archives.
#
# Output: data/sst_anomaly/nino34_jja.nc, laid out like sst_pattern_corr/corr_*.nc
# (one variable per dataset key and index, `year` plus `member_<key>`).

set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
module load parallel 2>/dev/null || true

WORK=${WORK_BASE}/nino34
mkdir -p "$WORK"/{box,mpige,epoc}

export N34_BOX="-sellonlatbox,-170,-120,-5,5"
export TROP_BOX="-sellonlatbox,-180,180,-20,20"

# $1 input (file, or "-operator file"), $2 output prefix -> <prefix>_{n34,trop}.nc
boxes() {
    cdo -s -O -fldmean $N34_BOX $1 "${2}_n34.nc"
    cdo -s -O -fldmean $TROP_BOX $1 "${2}_trop.nc"
}
export -f boxes

# $1 JJA monthly box-mean series (file, or "-operator file"),
# $2 output: JJA-mean anomaly against the 1991-2020 JJA monthly climatology
jja_anomaly() {
    cdo -s -O -yearmean -ymonsub $1 -ymonmean -selyear,${CLIM_START}/${CLIM_END} $1 "$2"
}
export -f jja_anomaly

echo "=== 1/3 datasets with native JJA anomalies on scratch ==="
declare -A ANOM=(
    [era5]=${WORK_BASE}/era5/era5_anom_jja.nc
    [eerie_r1]=${WORK_BASE}/eerie/eerie_r1_anom_jja.nc
    [eerie_r2]=${WORK_BASE}/eerie/eerie_r2_anom_jja.nc
    [eerie_r3]=${WORK_BASE}/eerie/eerie_r3_anom_jja.nc
    [eerie_control]=${WORK_BASE}/eerie/eerie_control_anom_jja.nc
    [sap0006]=${WORK_BASE}/sap0006/sap0006_anom_jja.nc
)
for name in "${!ANOM[@]}"; do
    f=${ANOM[$name]}
    [ -s "$f" ] || { echo "  missing $f -- rerun its pre_process stage" >&2; exit 1; }
    # the HEALPix file carries no usable grid; equal-area cells, so the
    # constant weights fldmean falls back to are the right ones
    [ "$name" = sap0006 ] && src="-setgrid,hp256_nested $f" || src=$f
    for kind in n34 trop; do rm -f "${WORK}/box/${name}_${kind}.nc"; done
    boxes "$src" "${WORK}/box/${name}"
    echo "  ${name}: $(cdo -s ntime "${WORK}/box/${name}_n34.nc") seasons"
done

echo "=== 2/3 MPI-GE, 50 members (1 deg conservative remap) ==="
CMIP=/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical
SCEN=/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245
WGT1=${WORK}/wgt_mpige_r360x180.nc
[ -s "$WGT1" ] || cdo -s gencon,r360x180 -seltimestep,1 \
    $(ls ${CMIP}/r1i1p1f1/Omon/tos/gn/v*/tos_*_185001-186912.nc) "$WGT1"

mpige_member() {
    ens=$1; work=$2; wgt=$3; cmip=$4; scen=$5
    out=${work}/box/mpige_r${ens}
    [ -s "${out}_n34.nc" ] && [ -s "${out}_trop.nc" ] && return 0
    tmp=${work}/mpige/r${ens}
    mkdir -p "$tmp"
    cdo -s -O -remap,r360x180,"$wgt" -selmon,6/8 -mergetime \
        ${cmip}/r${ens}i1p1f1/Omon/tos/gn/v*/tos_*.nc \
        ${scen}/r${ens}i1p1f1/Omon/tos/gn/v*/tos_*.nc "${tmp}/jja.nc"
    boxes "${tmp}/jja.nc" "${tmp}/mon"
    for kind in n34 trop; do jja_anomaly "${tmp}/mon_${kind}.nc" "${out}_${kind}.nc"; done
    rm -rf "$tmp"
}
export -f mpige_member
parallel -j 25 mpige_member {} "$WORK" "$WGT1" "$CMIP" "$SCEN" ::: $(seq 1 50)
echo "  done: $(ls ${WORK}/box/mpige_r*_n34.nc | wc -l) members"

echo "=== 3/3 EPOC transient + control, native 5 km ocean ==="
ICON_OCE_GRID=/pool/data/ICON/grids/public/mpim/0045/icon_grid_0045_R02B09_O.nc
declare -A EXP_DIR=(
  [epoc2_010]=/work/bm1313/b383127/epoc-icon-2024.10/experiments/epoc2_010/work
  [epoc2_020]=/work/bm1313/b383127/epoc-icon-2024.10_aerosols/experiments/epoc2_020/work
)

# one monthly file -> two box means; stamps are the END of the averaging month
epoc_file() {
    infile=$1; exp=$2; work=$3; gridfile=$4; mask=$5
    stamp=$(basename "$infile" | sed -E 's/.*_([0-9]{8})T[0-9]+Z\.nc/\1/')
    rep=$(date -u -d "${stamp} -1 day" +%Y-%m)
    year=${rep%-*}; month=${rep#*-}
    out=${work}/epoc/${exp}_${year}${month}
    [ -s "${out}_n34.nc" ] && [ -s "${out}_trop.nc" ] && return 0
    native=$(mktemp "${work}/epoc/native_XXXXXX.nc")
    cdo -s -O -settaxis,"${year}-${month}-15",12:00:00,1day -setgrid,"$gridfile" \
        -ifthen "$mask" -sellevidx,1 -selname,to "$infile" "$native"
    boxes "$native" "$out"
    rm -f "$native"
}
export -f epoc_file

for exp in epoc2_010 epoc2_020; do
    mapfile -t files < <(ls ${EXP_DIR[$exp]}/run_*/${exp}_oce_2d_1mth_mean_*.nc 2>/dev/null \
        | awk -F_ '{d=$NF; sub(/T.*/,"",d); m=substr(d,5,2); if (m=="07"||m=="08"||m=="09") print}' \
        | sort -u)
    MASK=${WORK_BASE}/epoc/landmask_${exp}.nc
    [ -s "$MASK" ] || cdo -s -O -gtc,-100 -setctomiss,0 -sellevidx,1 -selname,to "${files[0]}" "$MASK"
    echo "  ${exp}: ${#files[@]} JJA monthly files"
    # ~4.3 GB per process with the 14.9M-cell grid attached
    parallel -j ${NA2026_JOBS:-32} epoc_file {} "$exp" "$WORK" "$ICON_OCE_GRID" "$MASK" ::: "${files[@]}"
    for kind in n34 trop; do
        cdo -s -O -setreftime,1850-01-01,00:00:00,1day -mergetime \
            ${WORK}/epoc/${exp}_??????_${kind}.nc "${WORK}/epoc/${exp}_mon_${kind}.nc"
        # keep complete seasons only (the transient run is still in progress)
        yrs=$(cdo -s showdate "${WORK}/epoc/${exp}_mon_${kind}.nc" | tr -s ' ' '\n' \
              | cut -c1-4 | grep . | uniq -c | awk '$1==3{printf "%s,", $2}' | sed 's/,$//')
        jja_anomaly "-selyear,${yrs} ${WORK}/epoc/${exp}_mon_${kind}.nc" "${WORK}/box/${exp}_${kind}.nc"
    done
    echo "  ${exp}: $(cdo -s ntime "${WORK}/box/${exp}_n34.nc") seasons"
done

echo "=== assemble ${OUT_BASE}/nino34_jja.nc ==="
activate_python
WORK=$WORK OUT=$OUT_BASE python3 - <<'PYEOF'
import os
import numpy as np
import xarray as xr

work, out = os.environ["WORK"] + "/box", os.environ["OUT"]


def series(name, kind):
    da = xr.open_dataset(f"{work}/{name}_{kind}.nc", decode_times=True)
    da = da[next(v for v in da.data_vars if "bnds" not in v and "bounds" not in v)].squeeze(drop=True)
    da = da.drop_vars([c for c in da.coords if c != "time"])
    return da.assign_coords(year=("time", da.time.dt.year.values)).swap_dims(
        time="year").drop_vars("time").astype("f8")


def members(names, ids, dim):
    return lambda kind: xr.concat([series(n, kind) for n in names], dim, join="outer") \
        .assign_coords({dim: ids})


SOURCES = {
    "ERA5": lambda k: series("era5", k),
    "MPI-GE": members([f"mpige_r{m}" for m in range(1, 51)], np.arange(1, 51), "member_MPI-GE"),
    "ICON-EPOC-hist": lambda k: series("epoc2_020", k),
    "ICON-EPOC-ctrl": lambda k: series("epoc2_010", k),
    "EERIE": members([f"eerie_r{m}" for m in (1, 2, 3)], [1, 2, 3], "member_EERIE"),
    "EERIE-ctrl": lambda k: series("eerie_control", k),
    "ICON-sap-ctrl": lambda k: series("sap0006", k),
}

parts = {}
for key, get in SOURCES.items():
    n34, trop = get("n34"), get("trop")
    parts[f"{key}"] = (n34 - trop).assign_attrs(
        long_name="relative Nino-3.4: Nino-3.4 minus 20S-20N mean SST anomaly", units="K")
    parts[f"{key}_n34"] = n34.assign_attrs(long_name="Nino-3.4 SST anomaly", units="K")
    parts[f"{key}_trop"] = trop.assign_attrs(long_name="20S-20N mean SST anomaly", units="K")
    print(f"  {key:15s} {int(n34.year.min())}-{int(n34.year.max())}  sd(rn34) = "
          f"{float((n34 - trop).std()):.2f} K  sd(n34) = {float(n34.std()):.2f} K")
# outer join: assigning into one Dataset would cut every record to the first one's years
ds = xr.merge([v.rename(k) for k, v in parts.items()], join="outer", combine_attrs="drop")
ds.attrs.update(
    season="JJA", climatology="1991-2020 JJA of the same dataset",
    nino34_box="5S-5N, 170W-120W", tropics_box="20S-20N",
    note="variable <key> is the relative index; <key>_n34 and <key>_trop its parts")
ds.to_netcdf(f"{out}/nino34_jja.nc")
print(f"  -> {out}/nino34_jja.nc")
PYEOF
