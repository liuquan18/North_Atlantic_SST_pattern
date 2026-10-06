#!/bin/bash
# Submit the daily-Tmax stages (01-04) and the heatwave metrics (05), which
# runs only after all four have succeeded. Run from a login node.
set -euo pipefail
HW=/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts/european_heatwave

ids=()
for s in 01_era5_tmax.sh 02_epoc_tmax.sh 03_eerie_tmax.sh 04_mpige_tmax.sh; do
    ids+=("$(sbatch --parsable "${HW}/${s}")")
done
dep=$(IFS=:; echo "${ids[*]}")

sbatch --dependency=afterok:${dep} --job-name=hw_metrics --account=mh0033 \
    --partition=compute --nodes=1 --time=01:00:00 \
    --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_metrics.%j.out \
    --wrap "source ${HW}/hw_config.sh && activate_python && python3 ${HW}/05_heatwave_metrics.py"
