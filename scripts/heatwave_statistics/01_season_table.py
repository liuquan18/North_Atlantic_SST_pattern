"""
Season table for the heatwave statistics: one row per dataset, member and JJA.

Inputs : data/heatwave/<ds>_heatwave_<grid>.nc      (european_heatwave/05)
         data/sst_anomaly/<ds>_jja_gmsst.nc          (pre_process/)
Output : data/heatwave_stats/season_table.csv

Columns
  dataset, kind, member, year
  hwd, hwcum, hwexcess, hwpeak   land means over src.heatwave.MEAN_BOX
                                 (southern Europe, 35-60N, 10W-30E), as in fig 8-11
  gmsst      JJA global-mean SST anomaly of that season (1991-2020 reference
             of the same dataset); the member's own where it exists, else the
             mean of the dataset's members
  T          the warming covariate: gmsst averaged over the dataset's members
             (MPI-GE: the 50-member mean, i.e. the forced signal) and smoothed
             with a trailing 4-year mean (src.extremes.smooth_covariate)

MPI-ESM1.2-ER is absent (no daily Tmax). EPOC control heatwaves come from
daily mean temperature (its Tmax is on tape only).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import src.extremes as ex
import src.heatwave as hw
import src.pattern_2026 as p2
from fig1_patterns_simulations import land_mean_series

OUT_DIR = ROOT / "data" / "heatwave_stats"
METRICS = ["hwd", "hwcum", "hwexcess", "hwpeak"]

#: dataset key -> kind; order is the order of every table and figure
KINDS = {
    "ERA5": "observed",
    "MPI-GE": "large ensemble",
    "ICON-EPOC-hist": "km-scale transient",
    "EERIE": "km-scale transient",
    "EERIE-ctrl": "km-scale control",
    "ICON-sap-ctrl": "km-scale control",
    "ICON-EPOC-ctrl": "km-scale control",
}


def heatwave_table(key):
    with xr.open_dataset(hw.metrics_file(key)) as ds:
        df = xr.Dataset({m: land_mean_series(ds[m].load()) for m in METRICS}).to_dataframe()
    df = df.reset_index()
    if "member" not in df:
        df["member"] = 1
    return df[["member", "year"] + METRICS].dropna(subset=["hwd"])


def gmsst_table(key):
    with xr.open_dataset(p2.DATASETS[key].gm_file) as g:
        name = [v for v in g.data_vars if "bnds" not in v][0]
        da = g[name].squeeze(drop=True).load()
    da = da.groupby("time.year").mean()
    if "member" not in da.dims:
        da = da.expand_dims(member=[1])
    df = da.rename("gmsst").to_dataframe().reset_index()[["member", "year", "gmsst"]]
    df["member"] = df.member.astype(int)
    return df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for key, kind in KINDS.items():
        h = heatwave_table(key)
        g = gmsst_table(key)
        forced = g.groupby("year").gmsst.mean()
        T = pd.Series(ex.smooth_covariate(forced.values, forced.index.values), forced.index)
        h["member"] = h.member.astype(int)
        h = h.merge(g, on=["member", "year"], how="left")
        # members without their own SST in some years (EERIE r2 after 2014)
        h["gmsst"] = h.gmsst.fillna(h.year.map(forced))
        h["T"] = h.year.map(T)
        h.insert(0, "kind", kind)
        h.insert(0, "dataset", key)
        rows.append(h)
        yrs = f"{h.year.min()}-{h.year.max()}"
        print(f"  {key:16s} {kind:20s} {h.member.nunique():2d} member(s)  {yrs}  "
              f"{len(h):5d} seasons  T missing: {int(h['T'].isna().sum())}")
    t = pd.concat(rows, ignore_index=True)
    out = OUT_DIR / "season_table.csv"
    t.to_csv(out, index=False, float_format="%.5g")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
