#!/bin/bash
# Shared configuration for the 2026 JJA North Atlantic + Mediterranean SST pattern analysis.
#
# Every dataset in this study (ERA5, MPI-GE, km-scale ICON/EPOC, EERIE ICON-ESM-ER)
# is reduced to the same two products so the downstream Python is dataset-agnostic:
#
#   <ds>_jja_anom_na025.nc   JJA-mean SST anomaly, one time step per year,
#                            on the common 0.25 deg regional grid (GRID_NA025)
#   <ds>_jja_gmsst.nc        global-ocean-mean JJA SST anomaly, one value per year
#
# Anomalies are always relative to the 1991-2020 JJA monthly climatology of the
# *same* dataset, computed on that dataset's native grid; remapping happens once,
# on the small anomaly field, never on raw monthly data.

set -euo pipefail

module load cdo/2.5.0-gcc-11.2.0 2>/dev/null || true
module load nco 2>/dev/null || true

# --- analysis region: North Atlantic + Mediterranean ---------------------
export LON_W=-80
export LON_E=40
export LAT_S=30
export LAT_N=60
export CLIM_START=1991
export CLIM_END=2020
export REF_YEAR=2026

export PROJECT_ROOT=/work/mh0033/m300883/North_Atlantic_SST_pattern
export WORK_BASE=/scratch/m/m300883/nalt2026
export OUT_BASE=${PROJECT_ROOT}/data/pattern_2026

mkdir -p "$WORK_BASE" "$OUT_BASE"

# --- common target grid --------------------------------------------------
# 0.25 deg matches the native resolution of ERA5 and of the EERIE `gr` output,
# and is fine enough that the km-scale runs are not smoothed away before the
# maps are drawn. The pattern correlation itself is computed on a 1 deg
# coarsening of this grid (done in Python) so that MPI-GE, whose ocean grid is
# ~1 deg, is not scored on structure it cannot resolve.
export GRID_NA025=${WORK_BASE}/grid_na025.txt
if [ ! -f "$GRID_NA025" ]; then
cat > "$GRID_NA025" <<GRIDEOF
gridtype  = lonlat
xsize     = 480
ysize     = 120
xfirst    = -79.875
xinc      = 0.25
yfirst    = 30.125
yinc      = 0.25
GRIDEOF
fi

# Conda env for the Python steps (see CLAUDE.md for the PYTHONNOUSERSITE gotcha)
activate_python() {
    # conda's magics activation hook dereferences MAGPLUS_HOME unguarded, which
    # aborts under `set -u`; relax it just for the activation.
    set +u
    source /sw/spack-levante/mambaforge-23.1.0-1-Linux-x86_64-3boc6i/etc/profile.d/conda.sh
    conda activate north_atlantic_sst
    set -u
    export PYTHONNOUSERSITE=1
}
