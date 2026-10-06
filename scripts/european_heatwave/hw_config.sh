#!/bin/bash
# Shared settings for the European summer heatwave calculation.
#
# Every dataset is reduced to one common product, so the heatwave detection in
# 04_heatwave_metrics.py is dataset-agnostic:
#
#   ${HW_WORK}/<ds>/tmax_<year>.nc   daily maximum 2 m temperature (degC),
#                                    April-October of <year>, on GRID_EU025
#
# April-October rather than May-September: heatwaves are detected in May-Sep
# (the paper's NH warm-season window), but the 15-day percentile window around
# 1 May and 30 Sep reaches into 24 Apr and 7 Oct.
#
# The daily fields are scratch intermediates (regenerable from the archives);
# only the per-year heatwave metrics go to data/heatwave_2026/.

source /work/mh0033/m300883/North_Atlantic_SST_pattern/scripts/config.sh

export HW_MONTHS=4/10
export HW_WORK=${WORK_BASE}/heatwave
export HW_OUT=${PROJECT_ROOT}/data/heatwave_2026
mkdir -p "$HW_WORK" "$HW_OUT" "${PROJECT_ROOT}/logs"

# --- common Europe grid: 30-72N, 15W-45E, 0.25 deg -------------------------
# Nodes sit on multiples of 0.25 deg, so the EERIE `gr` grid maps onto it
# exactly. Must match REGION in src/heatwave.py.
export GRID_EU025=${HW_WORK}/grid_eu025.txt
if [ ! -f "$GRID_EU025" ]; then
cat > "$GRID_EU025" <<GRIDEOF
gridtype  = lonlat
xsize     = 241
ysize     = 169
xfirst    = -15
xinc      = 0.25
yfirst    = 30
yinc      = 0.25
GRIDEOF
fi

# --- land mask (ERA5 land-sea mask > 0.5), applied to every dataset ----------
# The paper detects heatwaves over land only. One common mask keeps the
# panels comparable; it is cheap, so it is built on first use.
export LANDMASK_EU025=${HW_OUT}/landmask_eu025.nc
if [ ! -s "$LANDMASK_EU025" ]; then
    cdo -s -O -f nc -setctomiss,0 -gtc,0.5 -remapbil,"$GRID_EU025" -setgridtype,regular \
        /pool/data/ERA5/E5/sf/an/IV/172/E5sf00_IV_INVARIANT_172.grb "$LANDMASK_EU025"
fi
