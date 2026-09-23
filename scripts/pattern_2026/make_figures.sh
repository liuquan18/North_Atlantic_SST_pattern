#!/bin/bash
# Regenerate the analysis and every figure of the JJA 2026 study.
# Assumes the cdo stages 01-04 (and 10, for figure 5) have already run.
set -euo pipefail

SCRIPT_DIR=${NA2026_SCRIPT_DIR:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts/pattern_2026}
source "${SCRIPT_DIR}/00_config.sh"
activate_python
cd "$PROJECT_ROOT"

python3 "${SCRIPT_DIR}/05_pattern_corr.py"
python3 "${SCRIPT_DIR}/06_plot1_patterns.py"
python3 "${SCRIPT_DIR}/07_plot2_corr_timeseries.py"
python3 "${SCRIPT_DIR}/08_plot3_distribution.py"
python3 "${SCRIPT_DIR}/09_plot4_era5_context.py"
python3 "${SCRIPT_DIR}/11_plot5_kmscale.py"
python3 "${SCRIPT_DIR}/12_region_means.py"
python3 "${SCRIPT_DIR}/13_plot6_region_means.py"

echo
echo "figures in ${PROJECT_ROOT}/figures/pattern_2026:"
ls -1 "${PROJECT_ROOT}/figures/pattern_2026"/*.png
