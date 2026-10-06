# scripts

Pipeline for the JJA 2026 North Atlantic + Mediterranean SST pattern study.
Methodology, data traps and results: [doc/pattern_2026.md](../doc/pattern_2026.md).

```
scripts/
├── config.sh          shared settings (region, climatology, paths, common grid, conda env)
├── run_analysis.sh    runs analysis/ then plotting/ in order
├── european_heatwave/ daily Tmax -> JJA heatwave metrics -> data/heatwave_2026/ (see its README)
├── pre_process/       cdo stages, submitted with sbatch; raw model/obs data -> data/pattern_2026/
├── analysis/          pattern correlation and region means -> data/pattern_2026/results/
└── plotting/          one script per figure -> figures/pattern_2026/
```

## pre_process/

Each stage reduces one dataset to `<ds>_jja_anom_na025.nc` and `<ds>_jja_gmsst.nc`
in `data/pattern_2026/`. Stages 01–05 are independent and idempotent.

| script | dataset |
|---|---|
| `01_era5.sh` | ERA5 + ERA5T, 1940–2026 |
| `02_mpige.sh` | MPI-ESM1-2-LR grand ensemble, 50 members (calls `stack_mpige.py`) |
| `03_epoc_icon.sh` | km-scale ICON (EPOC), control + transient |
| `04_eerie.sh` | EERIE ICON-ESM-ER, control + hist/ssp245 |
| `05_mpi_er.sh` | MPI-ESM1.2-ER, 3 realizations (calls `stack_mpi_er.py`) |
| `06_kmscale_zoom.sh` | native-resolution zoom fields for figure 5 (needs 01 and 03) |

## analysis/

| script | output |
|---|---|
| `01_pattern_corr.py` | correlations, best analogues, bootstrap → `results/` |
| `02_region_means.py` | area-mean JJA anomalies per region → `results/region_means.nc` |

## plotting/

`fig1_patterns_era5.py`, `fig1_patterns_simulations.py` … `fig7_event_composite.py`, one per figure. Figures 1–4 and 7 need
`analysis/01`; figure 1 also needs `european_heatwave/` (its top row), figure 5 needs `pre_process/06`, figure 6 needs `analysis/02`.

SLURM logs go to `logs/` at the project root.
