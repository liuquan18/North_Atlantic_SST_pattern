"""
The JJA SLP of 11_slp.sh on a plain integer year axis: ERA5 as (year, lat,
lon), every model as (member, year, lat, lon), its members stacked.
Arguments: the models to stack (default: ERA5 and all of MODELS).
"""
import glob
import os
import re
import sys

import xarray as xr

WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + "/slp"
OUT = os.environ.get("OUT_SLP",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/circulation")
ATTRS = dict(units="hPa", long_name="JJA mean sea-level pressure")
#: work sub-directory -> output stem (as the SST anomaly files of the dataset)
MODELS = {"mpige": "mpige", "epoc": "epoc2_020", "eerie": "eerie_hist_ssp245", "mpier": "mpier",
          "epocctrl": "epoc2_010"}


def by_year(f):
    with xr.open_dataset(f) as ds:
        da = ds["slp"].squeeze(drop=True).load()
    return (da.assign_coords(year=da.time.dt.year).swap_dims(time="year")
            .drop_vars("time").astype("float32"))


def main(keys=None):
    if keys is None:
        era5 = by_year(f"{WORK}/era5_jja_slp.nc").assign_attrs(ATTRS)
        era5.to_dataset(name="slp").to_netcdf(f"{OUT}/era5_jja_slp_2p5deg.nc")
        print(f"  ERA5 {dict(era5.sizes)}")

    for sub, stem in MODELS.items():
        if keys is not None and sub not in keys:
            continue
        files = glob.glob(f"{WORK}/{sub}/slp_jja_r*.nc")
        if not files:
            print(f"  {sub}: no files, skipped")
            continue
        files.sort(key=lambda f: int(re.search(r"slp_jja_r(\d+)", f).group(1)))
        members = [int(re.search(r"slp_jja_r(\d+)", f).group(1)) for f in files]
        da = xr.concat([by_year(f) for f in files], "member", join="outer")
        da = da.assign_coords(member=members).assign_attrs(ATTRS)
        da.to_dataset(name="slp").to_netcdf(f"{OUT}/{stem}_jja_slp_2p5deg.nc")
        print(f"  {sub} {dict(da.sizes)}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
