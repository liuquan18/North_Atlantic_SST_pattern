"""
Stack the 50 per-member MPI-GE files into one array with a `member` dimension.

Kept out of the sbatch heredoc so it can be run and re-run on its own once the
per-member cdo stage is done.
"""
import glob
import os
import re
import sys

import dask
import xarray as xr

WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + "/mpige"
OUT = os.environ.get("OUT_BASE",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_2026")


def stack(pattern, fname, limit=None):
    files = sorted(glob.glob(f"{WORK}/{pattern}"),
                   key=lambda f: int(re.search(r"_r(\d+)i1p1f1", f).group(1)))
    if limit:
        files = files[:limit]
    if not files:
        raise SystemExit(f"no files matching {WORK}/{pattern}")
    members = [int(re.search(r"_r(\d+)i1p1f1", f).group(1)) for f in files]

    # NB: xarray >= 2025.x rejects data_vars="minimal" together with
    # coords="minimal" when concatenating over a *new* dimension.
    ds = xr.open_mfdataset(files, combine="nested", concat_dim="member",
                           coords="minimal", compat="override")
    ds = ds.assign_coords(member=("member", members))

    # cdo carries time_bnds along as a data variable of datetime type. It is of
    # no use downstream and it is what breaks the write: compressing it sets a
    # dtype-less encoding on a chunked datetime array, which xarray refuses.
    ds = ds.drop_vars([v for v in ds.data_vars if v.endswith("_bnds")], errors="ignore")
    for c in ds.coords:
        ds[c].attrs.pop("bounds", None)
        ds[c].encoding.pop("bounds", None)

    # No zlib here: the payload is float32 geophysical noise, so level-4 zlib
    # buys ~15% at ~30 MB/min single-threaded -- an hour to save 400 MB. The
    # uncompressed write finishes in under a minute.
    enc = {}
    if "time" in ds.coords:
        ds["time"].encoding = {
            "units": ds["time"].encoding.get("units", "days since 1850-01-01"),
            "calendar": ds["time"].encoding.get("calendar", "standard"),
            "dtype": "float64",
        }
    out = f"{OUT}/{fname}"
    # Single-threaded: on a 256-core compute node dask's default thread pool
    # reads the 50 netCDF files concurrently and the write deadlocks in the
    # (not thread-safe) netCDF/HDF5 library -- it hung for 35 min after the
    # header on 2026-10-07. Plain sequential I/O takes about a minute.
    with dask.config.set(scheduler="synchronous"):
        ds.to_netcdf(out, encoding=enc)
    print(f"  {len(files)} members -> {out}  {dict(ds.sizes)}")


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    suffix = f"_test{limit}" if limit else ""
    # per-member regional files are named after the latitude band (config.sh NA025_ID)
    stack(f"anom_na025_{os.environ['NA025_ID']}_r*i1p1f1.nc", f"mpige_jja_anom_na025{suffix}.nc", limit)
    stack("gmsst_r*i1p1f1.nc", f"mpige_jja_gmsst{suffix}.nc", limit)
