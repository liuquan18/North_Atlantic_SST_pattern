#!/bin/bash
# Regenerate the analysis and every figure of the JJA 2026 study.
# Assumes the cdo stages in pre_process/ (01-05, plus 06 for figure 5) have run.
set -euo pipefail

SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
activate_python
cd "$PROJECT_ROOT"

python3 "${SCRIPTS}/analysis/01_pattern_corr.py"
python3 "${SCRIPTS}/analysis/02_region_means.py"

python3 "${SCRIPTS}/plotting/fig1_patterns_era5.py"
python3 "${SCRIPTS}/plotting/fig1_patterns_simulations.py"
python3 "${SCRIPTS}/plotting/fig1_heatwave_analogues.py"
python3 "${SCRIPTS}/plotting/fig2_corr_timeseries.py"
python3 "${SCRIPTS}/plotting/fig2_corr_timeseries_short.py"
python3 "${SCRIPTS}/plotting/fig3_distribution.py"
python3 "${SCRIPTS}/plotting/fig3_hist.py"
python3 "${SCRIPTS}/plotting/fig4_era5_context.py"
python3 "${SCRIPTS}/plotting/fig5_kmscale.py"
python3 "${SCRIPTS}/plotting/fig6_region_means.py"
python3 "${SCRIPTS}/plotting/fig7_event_composite.py"
python3 "${SCRIPTS}/plotting/fig7_event_composite_select.py"

echo
echo "figures in ${PROJECT_ROOT}/figures/pattern_2026:"
ls -1 "${PROJECT_ROOT}/figures/pattern_2026"/*.png
