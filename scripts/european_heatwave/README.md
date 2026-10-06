# european_heatwave

European summer (JJA) heatwaves for the seasons shown in figure 1, using the
definition of Xu et al. (2026, *Nature Climate Change*,
doi:10.1038/s41558-026-02762-2, `.reference/s41558-026-02762-2.pdf`):

> at least three consecutive days on which the daily maximum temperature
> anomaly, relative to the multi-year mean for each calendar day, exceeds the
> calendar-day 90th-percentile threshold computed with a 15-day moving window.

| | paper | here |
|---|---|---|
| reference period (climatology and threshold) | 1979–2023 | **1991–2020**, each dataset its own |
| Tmax | daily max of hourly ERA5 2 m T | same for ERA5; model daily-max output |
| detection window | May–Sep (NH) | May–Sep; JJA metrics from the JJA days of those events |
| domain | global land | European land, 30–72°N, 15°W–45°E, 0.25° (ERA5 land-sea mask > 0.5); MPI-GE on its native T63 grid (land fraction > 50 %) |

The algorithm is in `src/heatwave.py` (tests: `test/test_heatwave.py`).

## Stages

| script | what |
|---|---|
| `hw_config.sh` | shared settings, common grid, land mask (sources `../config.sh`) |
| `01_era5_tmax.sh` | ERA5 hourly 2 m T (`E5/sf/an/1H/167`, ERA5T after 2026-07) → daily max, 1991–2026 |
| `02_epoc_tmax.sh` | EPOC `epoc2_020` `atm_2d_1d_max` `tas` (R02B08) → 0.25°, 1990–2025 |
| `03_eerie_tmax.sh` | EERIE ICON-ESM-ER: r1 `day/tasmax` (`gr`) 1991–2050; r2, r3 raw `atm_2d_1d_max_remap025` 1991–2020 |
| `04_mpige_tmax.sh` | MPI-GE (MPI-ESM1-2-LR) CMIP6 `day/tasmax`, 50 members, historical + ssp245, 1850–2100, native grid |
| `05_heatwave_metrics.py` | heatwave detection → `data/heatwave_2026/<ds>_heatwave_<eu025\|native>.nc`; takes dataset keys as arguments (default: all) |
| `run_heatwave.sh` | submits 01–04, then 05 once they finish |

01–04 are SLURM jobs, independent and idempotent. They write daily Tmax for
April–October to scratch (`/scratch/m/m300883/nalt2026/heatwave/<ds>/tmax_<year>.nc`;
MPI-GE: `mpige/tmax_r<member>.nc`).
April and October are only there so the 15-day window around 1 May and 30 Sep
is complete. If scratch has been purged, rerun them.

MPI-GE is kept on its native ~1.9° grid, cropped to Europe. Regridding
50 members × 251 years of daily fields to 0.25° would take ~440 GB of scratch
and would only interpolate a ~200 km model. Its metrics file therefore has
dims (member, year, lat, lon) on the T63 grid, with the model's own land mask
(`landmask_mpige_native.nc`). Each member is referenced to its own 1991–2020.

Output variables, per year: `hwd` (JJA heatwave days), `hwn` (events touching
JJA), `hwcum` (cumulative JJA heatwave anomaly, °C·day), `hwmax` (longest event
touching JJA), `hwpeak` (anomaly of the hottest JJA heatwave day), `hwexcess`
(excess heat above the threshold, °C·day, the Perkins-Kirkpatrick & Lewis 2020
cumulative-heat definition), `onset` (first heatwave day, May–Sep). Each file also stores
`threshold` and `climatology` per calendar day.

## Not covered

- **MPI-ESM1.2-ER**: there is no daily atmosphere output on disk.
  `/work/uo0122/u241089/MPIESM/*/outdata/echam6/` holds monthly means only
  (`Amon2d_*`, `*_ATM_mm_*`, `*_BOT_mm_*`), and `/work/mh1421/data/mpiesm/`
  holds ocean output only.
- **ERA5 2026-09-27**: missing from the ERA5T pool. It is treated as "not hot";
  this does not affect JJA.
