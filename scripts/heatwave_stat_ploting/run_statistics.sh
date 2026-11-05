#!/bin/bash
#SBATCH --job-name=hw_stats
#SBATCH --output=/work/mh0033/m300883/North_Atlantic_SST_pattern/logs/hw_stats.%j.out
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=mh0033
#
# Heatwave statistics for ERA5 2003 / 2026 (see 02_statistics.py) and their
# figures 13-16 and 18-24. Needs data/heatwave/ (european_heatwave/) and
# data/sst_anomaly/*_jja_gmsst.nc (pre_process/). The bootstrap refits use the
# node's cores; on a login node `bash run_statistics.sh` also works, slower.
set -euo pipefail
SCRIPTS=${NA2026_SCRIPTS:-/work/mh0033/m300883/North_Atlantic_SST_pattern/scripts}
source "${SCRIPTS}/config.sh"
activate_python
cd "$PROJECT_ROOT"
HS=${SCRIPTS}/heatwave_statistics

python3 "${HS}/01_season_table.py"
python3 "${HS}/02_statistics.py"
python3 "${HS}/03_mpige_window_rp.py" hwd
python3 "${HS}/04_txx_table.py"
python3 "${HS}/03_mpige_window_rp.py" txx
python3 "${HS}/05_mpige_window_rank.py" hwd
python3 "${HS}/05_mpige_window_rank.py" txx
for f in "${HS}"/fig1[3-9]_*.py "${HS}"/fig2[0-4]_*.py; do
    python3 "$f"
done
