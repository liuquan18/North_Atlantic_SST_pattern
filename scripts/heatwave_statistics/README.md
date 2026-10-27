# heatwave_statistics

How extreme and how rare the European heatwave summers of 2003 and 2026 were,
in ERA5 (1940–2026) and in the simulations. Methods and results:
[doc/heatwave_2026_stats.md](../../doc/heatwave_2026_stats.md). Extreme-value
code: `src/extremes.py` (tests `test/test_extremes.py`).

Needs `data/heatwave/` (from `european_heatwave/`, ERA5 from 1940) and
`data/sst_anomaly/*_jja_gmsst.nc` (from `pre_process/`); figs 22–23 also
`data/sst_pattern_corr/` (from `sst_pattern_correlation/`).

| script | output |
|---|---|
| `01_season_table.py` | `data/heatwave_stats/season_table.csv`: land-mean metrics + warming covariate per dataset/member/season |
| `02_statistics.py` | A: ERA5 nonstationary GEV variants (`era5_*.csv`); B: anomalies from the forced mean (`anomaly_exceedance.csv`, `season_anomalies.csv`); C: return periods in the simulations (`model_*.csv`, `reference_1991_2020.csv`) |
| `03_mpige_window_rp.py` | stationary GEV across the 50 MPI-GE members in a 10-year window around each ERA5 year (500 seasons); return period of each ERA5 summer in that climate (`mpige_window_gev.csv`, `era5_rp_in_mpige.csv`) |
| `04_txx_table.py` | `txx_table.csv`: TXx, the hottest day of JJA in land-mean daily Tmax over the box, and its anomaly from the dataset's own 1991–2020 mean (ERA5 1940–2026, MPI-GE 50 × 1850–2100; reads daily Tmax on scratch). `03_mpige_window_rp.py txx` then writes `*_txx.csv` |
| `05_mpige_window_rank.py` | as 03 without a fitted distribution or return periods: percentile rank (0–100%) of each ERA5 value among the 500 window summers, and the window's 50th / 90th / 99th percentiles (`mpige_window_rank*.csv`, `era5_rank_in_mpige*.csv`); argument `hwd` / `txx` |
| `fig13_era5_heatwave_extremes.py` | ERA5 record with fitted return levels; robustness matrix of 2003/2026 |
| `fig14_anomaly_distributions.py` | simulated vs observed anomalies, absolute and standardised |
| `fig15_return_levels.py` | return-level curves: ERA5 in three climates; simulations |
| `fig16_return_period_summary.py` | every return-period estimate of 2003 and 2026 |
| `fig18_era5_return_period_mpige.py` | (a) MPI-GE window return levels with ERA5 on top; (b) return period of each ERA5 summer vs year |
| `fig19_era5_txx_return_period_mpige.py` | the same three panels for TXx |
| `fig20_…_rank.py`, `fig21_…_rank.py` | figs 18 / 19 from the ranks (05) instead of the GEV, in two panels: (a) ERA5 on the window percentiles, (b) percentile rank of each ERA5 summer (50% line = MPI-GE median) |
| `fig22_era5_hwd_corr_mpige.py` | (a) as fig 20 (a); (b) the same for the SST pattern: ERA5 pattern correlation with JJA 2026 on the MPI-GE window percentiles and mean (reads `data/sst_pattern_corr/`); right column: box plots of both quantities for ERA5 (1991–2026) and ICON-EPOC, ICON-EERIE and MPI-ER (historical + scenario 1991–2030 where available) |
| `fig22_boxdiff_mpige.py` | (a) the box-difference index MED − NAL (`box_difference.nc`, from `sst_pattern_correlation/03`; fig 22's former panel c) on the MPI-GE window percentiles and mean; (b) the Mediterranean and (c) the North Atlantic box mean alone, on one y axis, where the warming that cancels in (a) shows; right column as fig 22 |
| `fig23_epoc_hwd_corr_mpige.py` | as fig 22 (a, b), 1991–2030, with one line per higher-resolution run in place of ERA5: ICON-EPOC, ICON-EERIE r1–r3, MPI-ER r2–r3 (heatwave days once `08_mpi_er_tmax.sh` has filled the season table); MPI-GE windows computed from `season_table.csv`; no marked summers |
| `fig23_boxdiff_mpige.py` | as `fig22_boxdiff_mpige.py` with the same runs as fig 23: MED − NAL, Mediterranean and North Atlantic box means, 1991–2030 |
| `fig24_med_nal_correlation.py` | co-variability of the MED and NAL JJA SST box means, 1991–2030 (ERA5 also 1940–2026), linearly detrended: (a) correlation, (b) standard deviations, (c) var(MED − NAL) = σ²_MED + σ²_NAL − 2 cov; bars for ERA5 and the fig 23 runs, MPI-GE members as box plots / member means |
| `run_statistics.sh` | all of the above; `sbatch run_statistics.sh` (~35 min on one node; the two 03 runs take most of it) |

`HWS_N_BOOT=20 python3 02_statistics.py` gives a quick test run.
