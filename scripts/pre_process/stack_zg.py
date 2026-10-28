"""
Stack the per-member JJA geopotential height anomalies of 10_model_zg.sh into
one file per dataset, variable z<ZLEVEL> (gpm) over (member, year, lat, lon).
"""
import glob
import os
import re

import xarray as xr

Z = "z" + os.environ.get("ZLEVEL", "300")
WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + f"/model_{Z}"
OUT = os.environ.get("OUT_BASE",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/sst_anomaly")
#: work-file prefix -> output stem (as the SST anomaly files of the dataset)
DATASETS = {"mpige": "mpige", "epoc": "epoc2_020", "eerie": "eerie_hist_ssp245",
            "mpier": "mpier"}


def main():
    for prefix, stem in DATASETS.items():
        files = glob.glob(f"{WORK}/{prefix}_r*_anom.nc")
        if not files:
            print(f"  {prefix}: no files, skipped")
            continue
        members = [int(re.search(r"_r(\d+)(?:i1p1f1)?_anom", f).group(1)) for f in files]
        fields = []
        for f in files:
            with xr.open_dataset(f) as ds:
                da = ds[Z].squeeze(drop=True).load()
            # the GRIB-born MPI-ER file keeps cdo's absolute axis ("day as
            # %Y%m%d.%f"), which xarray leaves undecoded
            years = (da.time.dt.year if da.time.dtype.kind in "MO"
                     else (da.time // 10000).astype(int))
            fields.append(da.assign_coords(year=years).swap_dims(time="year")
                          .drop_vars("time"))
        da = xr.concat(fields, "member", join="outer").assign_coords(member=members)
        da = da.sortby("member").astype("float32")
        da.attrs = dict(units="gpm", long_name=f"JJA {Z[1:]} hPa geopotential height anomaly "
                                               "against 1991-2020 of the same member")
        out = f"{OUT}/{stem}_jja_{Z}_anom_1deg.nc"
        da.to_dataset(name=Z).to_netcdf(out)
        print(f"  {prefix} -> {out}  {dict(da.sizes)}")


if __name__ == "__main__":
    main()
