"""
Stack the three EERIE ICON-ESM-ER hist+ssp245 realizations into one array
with a `member` dim.

Unlike MPI-ESM1.2-ER, the members do not cover the same years: r1 runs
1950-2050, r2 1975-2014 and r3 1975-2020 (see 04_eerie.sh). They are put on
one canonical mid-JJA axis and outer-joined, so a member is NaN in the years
it was not run -- the analysis already skips NaN seasons. The SST variable is
also named differently (`tos` from CMOR, `to` with a depth axis from the raw
output), so each member is reduced to a bare 3-D field first.
"""
import os

import pandas as pd
import xarray as xr

WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + "/eerie"
OUT = os.environ.get("OUT_BASE",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_2026")

MEMBERS = [1, 2, 3]


def load_canonical(path):
    ds = xr.open_dataset(path)
    name = next(v for v in ds.data_vars if not v.endswith("_bnds") and ds[v].ndim >= 1)
    da = ds[name].squeeze(drop=True)
    da = da.drop_vars([c for c in da.coords if c not in da.dims])
    years = da["time"].dt.year.values
    da = da.assign_coords(time=pd.to_datetime([f"{y}-07-16" for y in years]))
    da.attrs = {k: v for k, v in da.attrs.items() if k in ("units", "standard_name")}
    return da.rename("tos"), years


def stack(kind, fname):
    parts = []
    for m in MEMBERS:
        path = f"{WORK}/eerie_r{m}_jja_{kind}.nc"
        if not os.path.exists(path):
            raise SystemExit(f"missing: {path}")
        da, years = load_canonical(path)
        parts.append(da)
        print(f"  {kind} r{m}: {years[0]}-{years[-1]} ({len(years)} seasons)")

    da = xr.concat(parts, dim="member", join="outer")
    da = da.assign_coords(member=("member", MEMBERS))
    ds = da.to_dataset()
    ds["time"].encoding = {"units": "days since 1850-01-01",
                           "calendar": "standard", "dtype": "float64"}
    ds.attrs["note"] = ("r1 erc2020 (CMOR, 1950-2050), r2 erc2023 (1975-2014), "
                        "r3 erc2024 (1975-2020); NaN where a member was not run")
    ds.to_netcdf(f"{OUT}/{fname}")
    print(f"  -> {OUT}/{fname}  {dict(ds.sizes)}")


if __name__ == "__main__":
    stack("anom_na025", "eerie_hist_ssp245_jja_anom_na025.nc")
    stack("gmsst", "eerie_hist_ssp245_jja_gmsst.nc")
