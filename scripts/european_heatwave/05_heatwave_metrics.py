"""
European summer heatwave metrics per year, for every dataset whose daily Tmax
has been extracted by 01-04.

Definition (Xu et al. 2026, Nat. Clim. Change; see src/heatwave.py): at least
three consecutive days on which the daily Tmax anomaly -- against the
calendar-day mean -- exceeds the calendar-day 90th percentile of anomalies
pooled over a centred 15-day window; detected in May-September. Reference
period 1991-2020, each dataset (and each ensemble member) against its own
climatology. Land only.

Output, data/heatwave_2026/<ds>_heatwave_<grid>.nc:
  hwd, hwn, hwcum, hwmax,
  hwpeak, hwexcess, onset         ([member,] year, lat, lon)  -- see heatwave.season_metrics
  threshold                       ([member,] doy, lat, lon)   90th-percentile anomaly threshold
  climatology                     ([member,] doy, lat, lon)   calendar-day mean Tmax

<grid> is eu025 (common 0.25 deg Europe grid, ERA5 land mask) or native
(MPI-GE: T63 cropped to Europe, the model's own land mask).

Usage: python3 05_heatwave_metrics.py [dataset key ...]   (default: all)
"""
import os
import sys
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.heatwave as hw

HW_WORK = Path(os.environ.get("HW_WORK", "/scratch/m/m300883/nalt2026/heatwave"))
N_MEMBERS = 50

ATTRS = {
    "hwd": ("number of JJA days belonging to a heatwave", "1"),
    "hwn": ("heatwave events with at least one day in JJA", "1"),
    "hwcum": ("cumulative Tmax anomaly over JJA heatwave days", "degC day"),
    "hwmax": ("full length in days of the longest heatwave touching JJA", "1"),
    "hwpeak": ("Tmax anomaly of the hottest JJA heatwave day", "degC"),
    "hwexcess": ("excess heat: sum over JJA heatwave days of the anomaly above the threshold "
                 "(Perkins-Kirkpatrick & Lewis 2020)", "degC day"),
    "onset": ("first heatwave day in the May-Sep window (365-day calendar)", "day of year"),
    "threshold": ("90th percentile of Tmax anomaly, 15-day window", "degC"),
    "climatology": ("calendar-day mean Tmax", "degC"),
}


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def open_tmax(files):
    da = xr.open_mfdataset(files, combine="by_coords")["tmax"].squeeze(drop=True)
    # daily stamps differ between sources (00:00, 11:30, 12:00): keep the date
    idx = da.indexes["time"]
    t = pd.DatetimeIndex(idx.to_datetimeindex() if hasattr(idx, "to_datetimeindex") else idx)
    return da.assign_coords(time=t.normalize()).sortby("time")


def fill_gaps(da, label):
    """
    Reindex to a gap-free daily axis within each year, missing days NaN. A
    missing day can then not be "hot" and splits a spell running across it
    (ERA5T, e.g., has no file for 2026-09-27).
    """
    t = da.indexes["time"]
    if t.duplicated().any():
        raise ValueError(f"{label}: duplicated days in the daily Tmax series")
    full = pd.DatetimeIndex(np.concatenate([
        pd.date_range(g.min(), g.max(), freq="D") for _, g in t.to_series().groupby(t.year)]))
    if len(full) != len(t):
        missing = full.difference(t)
        print(f"   {label}: {len(missing)} missing days set to NaN: "
              + ", ".join(f"{d:%Y-%m-%d}" for d in missing[:10]))
        da = da.reindex(time=full)
    return da


def load_land(key):
    land = xr.open_dataset(hw.landmask_file(key)).squeeze(drop=True)
    # the CMIP-derived mask also carries lat_bnds/lon_bnds: take the 2-D field
    name = next(v for v in land.data_vars if set(land[v].dims) == {"lat", "lon"})
    return land[name].notnull()


# --------------------------------------------------------------------------
# one realisation
# --------------------------------------------------------------------------
def compute(da, land, label):
    """Heatwave metrics for one daily Tmax series (time, lat, lon) -> Dataset."""
    da = fill_gaps(da, label)
    t = da.indexes["time"]
    if land.shape != da.shape[1:] or not (np.allclose(land.lat, da.lat) and np.allclose(land.lon, da.lon)):
        raise ValueError(f"{label}: land mask and Tmax are not on the same grid")
    land = land.assign_coords(lat=da.lat, lon=da.lon)

    iy, ix = np.nonzero(land.values)
    x = da.values[:, iy, ix].astype(np.float64)        # (day, land point)
    doy = hw.noleap_doy(t)
    year = t.year.values

    ref = (year >= hw.REF_PERIOD[0]) & (year <= hw.REF_PERIOD[1])
    ref_years = np.unique(year[ref])
    if len(ref_years) != hw.REF_PERIOD[1] - hw.REF_PERIOD[0] + 1:
        raise ValueError(f"{label}: reference period incomplete, have {ref_years}")

    clim = hw.daily_climatology(x, doy, ref)
    anom = x - clim[doy]
    detect = np.arange(hw.DETECT_DOY[0], hw.DETECT_DOY[1] + 1)
    thr = hw.percentile_threshold(anom, doy, ref, detect)

    years, out = [], {}
    for y in np.unique(year):
        sel = year == y
        d = doy[sel]
        if d.min() > hw.DETECT_DOY[0] or d.max() < hw.DETECT_DOY[1]:
            print(f"   {label} {y}: May-Sep not complete ({t[sel][0]:%m-%d}..{t[sel][-1]:%m-%d}), skipped")
            continue
        m = hw.season_metrics(anom[sel], thr[d], d)
        years.append(int(y))
        for k, v in m.items():
            out.setdefault(k, []).append(v)

    def unpack(arr, lead_name, lead):
        full = np.full((len(lead), land.sizes["lat"], land.sizes["lon"]), np.nan, np.float32)
        full[:, iy, ix] = arr
        return xr.DataArray(full, coords={lead_name: lead, "lat": land.lat, "lon": land.lon},
                            dims=(lead_name, "lat", "lon"))

    ds = xr.Dataset({k: unpack(np.stack(v), "year", years) for k, v in out.items()})
    ds["threshold"] = unpack(thr[detect], "doy", detect)
    ds["climatology"] = unpack(clim[detect], "doy", detect)
    print(f"{label}: {len(t)} days, {t[0]:%Y-%m-%d} .. {t[-1]:%Y-%m-%d}, "
          f"{x.shape[1]} land points, {len(years)} seasons", flush=True)
    return ds


def _member(args):
    e, land = args
    f = HW_WORK / hw.DATASETS["MPI-GE"] / f"tmax_r{e}.nc"
    if not f.exists():
        return None
    return compute(open_tmax([f]).load(), land, f"MPI-GE r{e}").expand_dims(member=[e])


# --------------------------------------------------------------------------
# per dataset
# --------------------------------------------------------------------------
def process(key):
    land = load_land(key)
    if key == "MPI-GE":
        # spawn, not fork: by the time MPI-GE runs, earlier datasets have opened
        # netCDF/HDF5 files in this process, and forked workers inherit HDF5's
        # held lock and hang forever
        ctx = mp.get_context("spawn")
        nproc = int(os.environ.get("HW_PROCS", min(N_MEMBERS, os.cpu_count())))
        with ctx.Pool(nproc) as pool:
            parts = [p for p in pool.map(_member, [(e, land) for e in range(1, N_MEMBERS + 1)])
                     if p is not None]
        if not parts:
            print(f"!! {key}: no daily Tmax, skipped")
            return
        ds = xr.concat(parts, "member")
        land_desc = "MPI-ESM1-2-LR sftlf > 50 %, native T63 grid"
    elif key in hw.MEMBERS:
        parts = []
        for e in hw.MEMBERS[key]:
            files = sorted((HW_WORK / hw.DATASETS[key] / f"r{e}").glob("tmax_[12][0-9][0-9][0-9].nc"))
            if not files:
                print(f"!! {key} r{e}: no daily Tmax, skipped")
                continue
            parts.append(compute(open_tmax(files), land, f"{key} r{e}").expand_dims(member=[e]))
        if not parts:
            return
        # members cover different years: outer-join, NaN where a member was not run
        ds = xr.concat(parts, "member", join="outer")
        land_desc = "ERA5 land-sea mask > 0.5"
    else:
        files = sorted((HW_WORK / hw.DATASETS[key]).glob("tmax_[12][0-9][0-9][0-9].nc"))
        if not files:
            print(f"!! {key}: no daily Tmax in {HW_WORK / hw.DATASETS[key]}, skipped")
            return
        ds = compute(open_tmax(files), land, key)
        land_desc = "ERA5 land-sea mask > 0.5"

    for k, (ln, u) in ATTRS.items():
        ds[k].attrs.update(long_name=ln, units=u)
    ds.attrs.update(
        title=f"European summer heatwaves, {key}",
        definition="Xu et al. 2026, Nat. Clim. Change, doi:10.1038/s41558-026-02762-2: >=3 consecutive "
                   "days with Tmax anomaly (vs calendar-day mean) above the calendar-day 90th "
                   "percentile of anomalies (15-day centred window); detection window May-Sep",
        reference_period=f"{hw.REF_PERIOD[0]}-{hw.REF_PERIOD[1]}",
        land_mask=land_desc,
    )

    hw.DATA_DIR.mkdir(parents=True, exist_ok=True)
    f = hw.metrics_file(key)
    ds.to_netcdf(f, encoding={v: {"zlib": True, "complevel": 4} for v in ds.data_vars})

    wts = np.cos(np.deg2rad(ds.lat))
    mean_hwd = ds.hwd.weighted(wts).mean(("lat", "lon"))
    if "member" in mean_hwd.dims:
        mean_hwd = mean_hwd.mean("member")
    clim_mean = float(mean_hwd.sel(year=slice(*hw.REF_PERIOD)).mean())
    top = mean_hwd.to_series().sort_values(ascending=False).head(5)
    print(f"   -> {f}  {dict(ds.hwd.sizes)}")
    print(f"   land-mean JJA heatwave days: {hw.REF_PERIOD[0]}-{hw.REF_PERIOD[1]} mean {clim_mean:.1f}; top 5:",
          ", ".join(f"{k}: {v:.1f}" for k, v in top.items()))


def main():
    keys = sys.argv[1:] or list(hw.DATASETS)
    for key in keys:
        process(key)


if __name__ == "__main__":
    main()
