"""
Figure 1 (ERA5) -- European summer heatwaves and the North Atlantic +
Mediterranean SST pattern in observed summers: JJA 2026 and the two earlier
summers that, like it, rank at or above 90% of the MPI-GE summers of their own
climate in heatwave days (figure 20): JJA 2003 and JJA 1947.

  top row     JJA heatwave days over European land (Xu et al. 2026 definition,
              1991-2020 reference; scripts/european_heatwave/)
  bottom row  the SST pattern of the same season (box mean removed), with its
              pattern correlation against JJA 2026, and the geopotential height
              anomaly of the season as contours (500 hPa by default), treated
              like the SST: 1991-2020 anomaly, area-weighted box mean removed
              (data/sst_anomaly/era5_jja_z<level>_anom_1deg.nc, pre_process/09_era5_zg.sh)

Same layout and colour scales as fig1_patterns_simulations.py, whose
plot_figure() draws it. With an argument `atl` or `med`, r is scored over
that half of the box alone (file suffix _atl / _med); with `z<level>`, e.g.
`z300`, the contours are of that level (file suffix _z300; 500 is the
default and gets no suffix).
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
from fig1_patterns_simulations import VARIANT, column, plot_figure, scored_over, sst_box

YEARS = [p2.REF_YEAR, 2003, 1947]
ZG_DEFAULT = 500


def zg_field(year, level=ZG_DEFAULT):
    """
    JJA geopotential height anomaly (gpm) at `level` hPa of `year` minus its
    area-weighted mean over the analysis box, as the SST pattern; kept on the
    5 deg margin around the box.
    """
    with xr.open_dataset(p2.DATA_DIR / f"era5_jja_z{level}_anom_1deg.nc") as ds:
        z = ds[f"z{level}"].squeeze("plev", drop=True).load()
    z = z.assign_coords(year=z.time.dt.year).swap_dims(time="year").drop_vars("time")
    z = z.sel(year=year)
    box = z.sel(lon=slice(*p2.REGION["lon"]), lat=slice(*p2.REGION["lat"]))
    w = np.cos(np.deg2rad(box.lat)) * xr.ones_like(box.lon)
    return z - float((box * w).sum() / w.sum())


def main(half=None, level=ZG_DEFAULT):
    with xr.open_dataset(p2.corr_file(VARIANT, half)) as ds:
        corr = ds["ERA5"].load()
    label = p2.DATASETS["ERA5"].label
    cols = []
    for year in YEARS:
        ref = year == p2.REF_YEAR
        r = None if ref else float(corr.sel(year=year))
        title = f"{label}\nJJA {year}" + ("  ·  reference" if ref else "")
        cols.append(dict(column("ERA5", year, None, r, title), zg=zg_field(year, level)))
    plot_figure(cols, "European summer heatwaves and the North Atlantic + Mediterranean SST pattern "
                      f"in ERA5: JJA {', '.join(map(str, YEARS[:-1]))} and {YEARS[-1]}{scored_over(half)}",
                p2.FIG_DIR / f"fig1_patterns_2026_era5{p2.half_suffix(half)}"
                             f"{'' if level == ZG_DEFAULT else f'_z{level}'}.png",
                sst_region_label=p2.half_label(half), sst_box=sst_box(half), zg_level=level)


if __name__ == "__main__":
    levels = [int(a[1:]) for a in sys.argv[1:] if a.startswith("z")]
    main(p2.half_from_argv([sys.argv[0]] + [a for a in sys.argv[1:] if not a.startswith("z")]),
         *levels[:1])
