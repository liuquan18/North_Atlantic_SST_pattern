"""
Figure 1 (ERA5) -- European summer heatwaves and the North Atlantic +
Mediterranean SST pattern in observed summers: JJA 2026 and the two earlier
benchmark heatwave summers, JJA 2015 and JJA 2003.

  top row     JJA heatwave days over European land (Xu et al. 2026 definition,
              1991-2020 reference; scripts/european_heatwave/)
  bottom row  the SST pattern of the same season (box mean removed), with its
              pattern correlation against JJA 2026

Same layout and colour scales as fig1_patterns_simulations.py, whose
plot_figure() draws it. With an argument `atl` or `med`, r is scored over
that half of the box alone (file suffix _atl / _med).
"""
import sys
from pathlib import Path

import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
from fig1_patterns_simulations import VARIANT, column, plot_figure, scored_over, sst_box

YEARS = [p2.REF_YEAR, 2015, 2003]


def main(half=None):
    with xr.open_dataset(p2.corr_file(VARIANT, half)) as ds:
        corr = ds["ERA5"].load()
    label = p2.DATASETS["ERA5"].label
    cols = []
    for year in YEARS:
        ref = year == p2.REF_YEAR
        r = None if ref else float(corr.sel(year=year))
        title = f"{label}\nJJA {year}" + ("  ·  reference" if ref else "")
        cols.append(column("ERA5", year, None, r, title))
    plot_figure(cols, "European summer heatwaves and the North Atlantic + Mediterranean SST pattern "
                      f"in ERA5: JJA {', '.join(map(str, YEARS[:-1]))} and {YEARS[-1]}{scored_over(half)}",
                p2.FIG_DIR / f"fig1_patterns_2026_era5{p2.half_suffix(half)}.png",
                sst_region_label=p2.half_label(half), sst_box=sst_box(half))


if __name__ == "__main__":
    main(p2.half_from_argv(sys.argv))
