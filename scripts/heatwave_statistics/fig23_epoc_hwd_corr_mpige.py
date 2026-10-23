"""
Figure 23 -- as figure 22, with the km-scale ICON-EPOC run (hist + SSP5-8.5,
10 km atmosphere / 5 km ocean, 1990-2025) in place of ERA5 as the black line.

(a) ICON-EPOC heatwave days on the 50th / 90th / 99th percentiles and mean of
    the 500 MPI-GE summers in the 10 years around each year.
(b) the pattern correlation of each ICON-EPOC JJA SST anomaly with ERA5 JJA
    2026, on the same windows of the MPI-GE pattern correlations.

The axis covers the EPOC years only. Marked: the EPOC summers ranked at or
above MARK_PCT (05_mpige_window_rank.py) of their window in heatwave days --
EPOC runs free, so its weather years have nothing to do with the observed 2003
or 2026.
"""
import importlib
import sys
from pathlib import Path

import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
w03 = importlib.import_module("03_mpige_window_rp")
w05 = importlib.import_module("05_mpige_window_rank")
import src.heatwave as hw
import src.pattern_2026 as p2

KEY = "ICON-EPOC-hist"
NAME = "ICON-EPOC"
OUT = p2.FIG_DIR / "fig23_epoc_hwd_corr_mpige.png"
CAPTION_TITLE = "Heatwave days and the SST pattern of each km-scale ICON summer against the MPI-GE summers of its climate"
CAPTION = (
    "As Fig. 22, with the km-scale ICON-EPOC run (historical + SSP5-8.5 from 2015, 10 km "
    "atmosphere / 5 km ocean, one realisation, 1990–2025) in black in place of ERA5. (a) "
    f"ICON-EPOC June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL}, against its "
    "own 1991–2020 threshold; blue: heatwave days of MPI-ESM1.2-LR: the 50th, 90th and 99th "
    "percentiles (solid) and the mean (dashed) of its 50 members in the 10 years around each "
    "year (y−5 to y+4, 500 summers). (b) The same for the SST pattern: centred pattern "
    "correlation of each JJA SST anomaly with that of ERA5 JJA 2026 over 20–60°N, 80°W–40°E "
    "(area-weighted box mean removed, anomalies against 1991–2020 of the same dataset). "
    f"Marked: the ICON-EPOC summers ranked at or above {w05.MARK_PCT}% of their window in "
    "heatwave days; ICON-EPOC runs free, so its years do not follow the observed weather.")


def main():
    meth = fig18.METHODS["rank"]
    t = pd.read_csv(fig18.DIR / "season_table.csv")
    mpige = t[t.dataset == "MPI-GE"]
    epoc = t[t.dataset == KEY].sort_values("year")
    years, hwd = epoc.year.values, epoc.hwd.values

    win_hwd = pd.read_csv(fig18.DIR / f"{meth['windows']}.csv")
    win_hwd = win_hwd[win_hwd.year.isin(years)].reset_index(drop=True)
    pct = {int(y): w05.percentile_rank(w03.window(mpige, y).hwd.values, x)["pct"]
           for y, x in zip(years, hwd)}
    marked = sorted(y for y, v in pct.items() if v >= w05.MARK_PCT)

    with xr.open_dataset(p2.corr_file(p2.PLOT_VARIANTS[0])) as ds:
        r = ds[KEY].load().dropna("year")
    fig22.figure(win_hwd, years, hwd, r, marked, NAME, CAPTION_TITLE, CAPTION, OUT)
    for y in marked:
        print(f"  {y}: {float(epoc.hwd[epoc.year == y].iloc[0]):5.1f} heatwave days, "
              f"rank {pct[y]:.1f}% of its MPI-GE window")


if __name__ == "__main__":
    main()
