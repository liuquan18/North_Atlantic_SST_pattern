"""
Stack the three MPI-ESM1.2-ER realizations into one array with a `member` dim.

The three runs do not share an identical time axis: MPIOM stamps its monthly
means at the end of the month, and ER lands on 23:55:00 while ER3 and ER5 land
on 23:56:00. Those one-minute offsets survive `yearmean`, so concatenating on
the raw axis silently outer-joins into a 450-step union that is two-thirds
missing instead of a 150-step, 3-member array. Each member's axis is therefore
rewritten to a canonical mid-JJA date from its own year before concatenating,
and the years are checked to be identical across members.
"""
import os

import numpy as np
import pandas as pd
import xarray as xr

WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + "/mpier"
OUT = os.environ.get("OUT_BASE",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_2026")

# readme: ER = member 1, ER3 = member 2, ER5 = member 3
RUNS = [("ER", 1), ("ER3", 2), ("ER5", 3)]


def load_canonical(path):
    ds = xr.open_dataset(path)
    ds = ds.drop_vars([v for v in ds.data_vars if v.endswith("_bnds")], errors="ignore")
    years = ds["time"].dt.year.values
    ds = ds.assign_coords(time=pd.to_datetime([f"{y}-07-16" for y in years]))
    for c in ds.coords:
        ds[c].attrs.pop("bounds", None)
        ds[c].encoding.pop("bounds", None)
    return ds, years


def stack(kind, fname):
    parts, all_years = [], []
    for run, _ in RUNS:
        path = f"{WORK}/{run}_anom_{kind}.nc"
        if not os.path.exists(path):
            raise SystemExit(f"missing: {path}")
        ds, years = load_canonical(path)
        parts.append(ds)
        all_years.append(years)

    for run_years, (run, _) in zip(all_years[1:], RUNS[1:]):
        if not np.array_equal(run_years, all_years[0]):
            raise SystemExit(f"{run} covers different years than {RUNS[0][0]}")

    ds = xr.concat(parts, dim="member", join="exact")
    ds = ds.assign_coords(member=("member", [m for _, m in RUNS]))
    ds["time"].encoding = {"units": "days since 1850-01-01",
                           "calendar": "standard", "dtype": "float64"}
    ds.to_netcdf(f"{OUT}/{fname}")
    print(f"  {kind}: {len(parts)} members, {len(all_years[0])} seasons "
          f"({all_years[0][0]}-{all_years[0][-1]}) -> {OUT}/{fname}  {dict(ds.sizes)}")


if __name__ == "__main__":
    stack("reg", "mpier_jja_anom_na025.nc")
    stack("gm", "mpier_jja_gmsst.nc")
