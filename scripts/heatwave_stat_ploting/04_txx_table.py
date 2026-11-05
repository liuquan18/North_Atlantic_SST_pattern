"""
TXx of southern Europe: the hottest day of each JJA in the land-mean daily
maximum temperature over src.heatwave.MEAN_BOX -- the block maximum (block =
one summer) that the GEV is meant for.

  series   daily Tmax averaged over the land cells inside MEAN_BOX
           (cos-latitude weights; the same land masks as the heatwave
           metrics, data/heatwave/landmask_*.nc)
  txx      its maximum over 1 June - 31 August, degC
  txx_anom txx minus the dataset's own 1991-2020 mean of txx (MPI-GE: over
           all 50 members), which takes out the model's mean Tmax bias in the
           same way the heatwave threshold does
  doy      day of year of that maximum (365-day calendar)

Input : daily Tmax written by scripts/european_heatwave/01_era5_tmax.sh and
        04_mpige_tmax.sh to HW_WORK (/scratch/m/m300883/nalt2026/heatwave/):
        era5/tmax_<year>.nc (0.25 deg), mpige/tmax_r<member>.nc (native T63)
Output: data/heatwave_stats/txx_table.csv  (dataset, member, year, txx, txx_anom, doy)
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
import src.heatwave as hw

DIR = ROOT / "data" / "heatwave_stats"
HW_WORK = Path(os.environ.get("HW_WORK", "/scratch/m/m300883/nalt2026/heatwave"))
N_MEMBERS = 50


def box_weights(key, lat, lon):
    """cos(lat) weights on the land cells inside MEAN_BOX, zero elsewhere."""
    land = xr.open_dataset(hw.landmask_file(key)).squeeze(drop=True)
    name = next(v for v in land.data_vars if set(land[v].dims) == {"lat", "lon"})
    land = land[name].notnull()
    if not (np.allclose(land.lat, lat) and np.allclose(land.lon, lon)):
        raise ValueError(f"{key}: land mask and Tmax are not on the same grid")
    (lon0, lon1), (lat0, lat1) = hw.MEAN_BOX["lon"], hw.MEAN_BOX["lat"]
    inside = (lat[:, None] >= lat0) & (lat[:, None] <= lat1) & (lon >= lon0) & (lon <= lon1)
    w = np.cos(np.deg2rad(lat))[:, None] * (land.values & inside)
    return w / w.sum()


def jja_max(series):
    """(year, txx, doy) of each JJA in a daily pandas Series indexed by time."""
    s = series[series.index.month.isin([6, 7, 8])]
    rows = []
    for y, g in s.groupby(s.index.year):
        if len(g) < 92:                    # incomplete summer
            continue
        t = g.idxmax()
        rows.append((int(y), float(g.max()), int(t.dayofyear)))
    return rows


def box_series(da, w):
    """Weighted land-box mean of (time, lat, lon) Tmax as a daily Series."""
    v = (da.values * w).sum(axis=(1, 2))
    return pd.Series(v, index=pd.DatetimeIndex(da.time.values).normalize())


def era5():
    files = sorted((HW_WORK / "era5").glob("tmax_[12][0-9][0-9][0-9].nc"))
    w = None
    rows = []
    for f in files:
        with xr.open_dataset(f) as ds:
            da = ds.tmax.load()
        if w is None:
            w = box_weights("ERA5", da.lat.values, da.lon.values)
        rows += jja_max(box_series(da, w))
    return pd.DataFrame(rows, columns=["year", "txx", "doy"]).assign(dataset="ERA5", member=1)


def mpige():
    out = []
    for m in range(1, N_MEMBERS + 1):
        with xr.open_dataset(HW_WORK / "mpige" / f"tmax_r{m}.nc") as ds:
            da = ds.tmax.load()
        w = box_weights("MPI-GE", da.lat.values, da.lon.values)
        out.append(pd.DataFrame(jja_max(box_series(da, w)), columns=["year", "txx", "doy"])
                   .assign(dataset="MPI-GE", member=m))
    return pd.concat(out, ignore_index=True)


def main():
    DIR.mkdir(parents=True, exist_ok=True)
    t = pd.concat([era5(), mpige()], ignore_index=True)
    for key, g in t.groupby("dataset"):
        ref = g[g.year.between(*hw.REF_PERIOD)].txx.mean()
        t.loc[g.index, "txx_anom"] = g.txx - ref
        print(f"  {key:7s} {g.member.nunique():2d} member(s) {g.year.min()}-{g.year.max()}  "
              f"{len(g):5d} summers  1991-2020 mean TXx {ref:.2f} degC")
    t = t[["dataset", "member", "year", "txx", "txx_anom", "doy"]]
    out = DIR / "txx_table.csv"
    t.to_csv(out, index=False, float_format="%.4f")
    e = t[t.dataset == "ERA5"].set_index("year")
    print("  ERA5 hottest summers:", ", ".join(
        f"{y} ({r.txx:.1f} degC, {r.txx_anom:+.1f})" for y, r in
        e.sort_values("txx", ascending=False).head(6).iterrows()))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
