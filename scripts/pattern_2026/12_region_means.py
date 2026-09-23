"""
Area-mean JJA SST anomaly time series for each dataset, over four regions.

Regions (see src/pattern_2026.REGIONS):
  global   each dataset's own global ocean mean (already produced by 01-04)
  box      the analysis box, 30-60N 80W-40E
  natl     the box with the Mediterranean box removed
  med      30-50N, 10W-40E

box = natl + med exactly, so the three regional curves decompose rather than
duplicate one another.

The three regional means are taken on the common 1 deg ocean mask -- the same
set of cells the pattern correlations use -- so differences between datasets
are differences in SST, not in where each model thinks the coastline is. The
global mean cannot be put on a common mask (the datasets have different global
grids) and is each dataset's own, which for a global ocean mean is robust.

Output: data/pattern_2026/results/region_means.nc  (+ region_means_2026.csv)
"""
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2

OUT = p2.RESULT_DIR / "region_means.nc"


def main():
    mask = xr.open_dataset(p2.RESULT_DIR / "ocean_mask_1deg.nc")["ocean_mask"]
    masks = p2.region_masks(mask)

    keys = [k for k in p2.DATASETS if p2.DATASETS[k].anom_file.exists()]
    print(f"datasets: {keys}")

    series = {}
    for k in keys:
        print(f"  {k} ...", flush=True)
        gm = p2.load_gmsst(k)
        series[f"{k}|global"] = gm

        # one streaming pass over the 0.25 deg field, then the three regions
        coarse = p2.coarsen_to_1deg(p2.load_anomaly(k)).compute()
        for region, m in masks.items():
            series[f"{k}|{region}"] = p2.area_mean(coarse, m)

    ds = xr.Dataset(series)
    ds.attrs["description"] = ("area-weighted JJA SST anomaly (degC, 1991-2020 base); "
                               "variable names are '<dataset>|<region>'")
    ds.to_netcdf(OUT)
    print(f"\nwrote {OUT}")

    rows = []
    for name, da in series.items():
        k, region = name.split("|")
        v = da.sel(year=p2.REF_YEAR) if p2.REF_YEAR in da.year.values else None
        rows.append(dict(
            dataset=k, region=region,
            year_first=int(da.year.min()), year_last=int(da.year.max()),
            value_2026=None if v is None else round(float(v.mean()), 3),
            p05_1991_2020=round(float(da.sel(year=slice(*p2.CLIM_PERIOD)).quantile(0.05)), 3),
            p95_1991_2020=round(float(da.sel(year=slice(*p2.CLIM_PERIOD)).quantile(0.95)), 3),
            trend_per_decade=round(float(
                np.polyfit(da.year.values,
                           da.mean("member").values if "member" in da.dims else da.values,
                           1)[0] * 10), 4),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(p2.RESULT_DIR / "region_means_2026.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
