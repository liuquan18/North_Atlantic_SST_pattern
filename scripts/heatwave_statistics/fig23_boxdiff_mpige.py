"""
Figure 23 box difference -- as fig22_boxdiff_mpige.py, with every
higher-resolution historical + scenario run (fig 23's RUNS, one line per
realization) in place of ERA5, over 1991-2030.

(a) the box-difference index MED − NAL of each run, (b) its Mediterranean and
(c) its North Atlantic box mean alone (data/sst_pattern_corr/box_difference.nc),
each on the 50th / 90th / 99th percentiles and mean of the 500 MPI-GE summers
in the 10 years around each year. No spatial mean is removed; (b) and (c)
share one y axis. No summers are marked.
"""
import importlib
import sys
from pathlib import Path

import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
fig22b = importlib.import_module("fig22_boxdiff_mpige")
fig23 = importlib.import_module("fig23_epoc_hwd_corr_mpige")
import src.pattern_2026 as p2

OUT = p2.FIG_DIR / "fig23_boxdiff_mpige.png"
CAPTION_TITLE = "The Mediterranean–North Atlantic SST contrast of the higher-resolution runs against the MPI-GE summers of their climate"
CAPTION = (
    f"As Fig. 22 box difference, {fig23.YEARS[0]}–{fig23.YEARS[1]}, with the higher-resolution "
    "historical + scenario runs in place of ERA5, one line per realisation: ICON-EPOC "
    "(SSP5-8.5 from 2015, 10 km atmosphere / 5 km ocean, 1 run, to 2025), ICON-EERIE "
    "(SSP2-4.5 from 2015, 10 km / 5 km, r1 to 2050, r2 to 2014, r3 to 2020) and MPI-ESM1.2-ER "
    "(SSP5-8.5 from 2015, ~100 km atmosphere / ~10 km ocean, r2 and r3); line style = member. "
    f"(a) Box-difference index {p2.DIPOLE_LABEL}: JJA SST anomaly averaged over the "
    f"Mediterranean ({p2.DIPOLE_BOXES['med']['detail']}) minus that over the North Atlantic "
    f"({p2.DIPOLE_BOXES['nal']['detail']}), ocean cells on the common 1° grid, anomalies "
    "against 1991–2020 of the same dataset; blue: the same index in MPI-ESM1.2-LR: the 50th, "
    "90th and 99th percentiles (solid) and the mean (dashed) of its 50 members in the 10 years "
    "around each year (y−5 to y+4, 500 summers). (b) The same for the Mediterranean box mean "
    "alone and (c) for the North Atlantic box mean alone (same y axis): no spatial mean is "
    "removed, so the warming that cancels in (a) shows here.")


def main():
    years = list(range(fig23.YEARS[0], fig23.YEARS[1] + 1))
    rows = []
    with xr.open_dataset(p2.dipole_file()) as ds:
        for suffix, ylabel in fig22b.ROWS:
            runs = {(key, m): fig23.run_series(ds[key + suffix].load(), m)
                    for key, (_, members) in fig23.RUNS.items() for m in members}
            rows.append((fig22.window_levels(ds["MPI-GE" + suffix].load(), years), runs, ylabel))
    fig23.figure(rows, CAPTION_TITLE, CAPTION, OUT, share=(1, 2))


if __name__ == "__main__":
    main()
