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
| `03_mpige_window_rp.py` | stationary GEV across the 50 MPI-GE members in a 10-year window around each ERA5 year (500 seasons); return period of each ERA5 summer in that climate (`mpige_window_gev.csv`, `era5_rp_in_mpige.csv`) |
| `04_txx_table.py` | `txx_table.csv`: TXx, the hottest day of JJA in land-mean daily Tmax over the box, and its anomaly from the dataset's own 1991–2020 mean (ERA5 1940–2026, MPI-GE 50 × 1850–2100; reads daily Tmax on scratch). `03_mpige_window_rp.py txx` then writes `*_txx.csv` |
| `05_mpige_window_rank.py` | as 03 without a fitted distribution: return period 501 / (k + 1) from the ERA5 value's rank among the 500 window summers, percentiles as return levels (`mpige_window_rank*.csv`, `era5_rank_in_mpige*.csv`); argument `hwd` / `txx` |
| `fig13_era5_heatwave_extremes.py` | ERA5 record with fitted return levels; robustness matrix of 2003/2026 |
| `fig14_anomaly_distributions.py` | simulated vs observed anomalies, absolute and standardised |
| `fig15_return_levels.py` | return-level curves: ERA5 in three climates; simulations |
| `fig16_return_period_summary.py` | every return-period estimate of 2003 and 2026 |
| `fig18_era5_return_period_mpige.py` | (a) MPI-GE window return levels with ERA5 on top; (b) return period of each ERA5 summer vs year |
| `fig19_era5_txx_return_period_mpige.py` | the same three panels for TXx |
| `fig20_…_rank.py`, `fig21_…_rank.py` | figs 18 / 19 from the ranks (05) instead of the GEV |
| `run_statistics.sh` | all of the above; `sbatch run_statistics.sh` (~35 min on one node; the two 03 runs take most of it) |

`HWS_N_BOOT=20 python3 02_statistics.py` gives a quick test run.
