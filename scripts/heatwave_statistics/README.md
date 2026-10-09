# heatwave_statistics

How extreme and how rare the European heatwave summers of 2003 and 2026 were,
in ERA5 (1940–2026) and in the simulations. Methods and results:
[doc/heatwave_2026_stats.md](../../doc/heatwave_2026_stats.md). Extreme-value
code: `src/extremes.py` (tests `test/test_extremes.py`).

Needs `data/heatwave/` (from `european_heatwave/`, ERA5 from 1940) and
`data/sst_anomaly/*_jja_gmsst.nc` (from `pre_process/`).

| script | output |
|---|---|
| `01_season_table.py` | `data/heatwave_stats/season_table.csv`: land-mean metrics + warming covariate per dataset/member/season |
| `02_statistics.py` | A: ERA5 nonstationary GEV variants (`era5_*.csv`); B: anomalies from the forced mean (`anomaly_exceedance.csv`, `season_anomalies.csv`); C: return periods in the simulations (`model_*.csv`, `reference_1991_2020.csv`) |
| `fig13_era5_heatwave_extremes.py` | ERA5 record with fitted return levels; robustness matrix of 2003/2026 |
| `fig14_anomaly_distributions.py` | simulated vs observed anomalies, absolute and standardised |
| `fig15_return_levels.py` | return-level curves: ERA5 in three climates; simulations |
| `fig16_return_period_summary.py` | every return-period estimate of 2003 and 2026 |
| `run_statistics.sh` | all of the above; `sbatch run_statistics.sh` (~4 min on one node) |

`HWS_N_BOOT=20 python3 02_statistics.py` gives a quick test run.
