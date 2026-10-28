"""
The JJA SLP of 11_slp.sh on a plain integer year axis: ERA5 as (year, lat,
lon), the 50 MPI-GE members stacked as (member, year, lat, lon).
"""
import glob
import os
import re

import xarray as xr

WORK = os.environ.get("WORK_BASE", "/scratch/m/m300883/nalt2026") + "/slp"
OUT = os.environ.get("OUT_SLP",
                     "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/circulation")
ATTRS = dict(units="hPa", long_name="JJA mean sea-level pressure")


def by_year(f):
    with xr.open_dataset(f) as ds:
        da = ds["slp"].squeeze(drop=True).load()
    return (da.assign_coords(year=da.time.dt.year).swap_dims(time="year")
            .drop_vars("time").astype("float32"))


def main():
    era5 = by_year(f"{WORK}/era5_jja_slp.nc").assign_attrs(ATTRS)
    era5.to_dataset(name="slp").to_netcdf(f"{OUT}/era5_jja_slp_2p5deg.nc")
    print(f"  ERA5 {dict(era5.sizes)}")

    files = sorted(glob.glob(f"{WORK}/mpige/slp_jja_r*i1p1f1.nc"),
                   key=lambda f: int(re.search(r"_r(\d+)i1p1f1", f).group(1)))
    members = [int(re.search(r"_r(\d+)i1p1f1", f).group(1)) for f in files]
    ge = xr.concat([by_year(f) for f in files], "member").assign_coords(member=members)
    ge = ge.assign_attrs(ATTRS)
    ge.to_dataset(name="slp").to_netcdf(f"{OUT}/mpige_jja_slp_2p5deg.nc")
    print(f"  MPI-GE {dict(ge.sizes)}")


if __name__ == "__main__":
    main()
