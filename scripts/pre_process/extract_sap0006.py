"""
Pull SST and daily Tmax out of the 100-yr ICON Sapphire control run sap0006
(10 km atm / 5 km ocean, constant 1950 forcing, 1950-2049), which is only
available as a HEALPix zarr store written by hiopy:

  /work/mh1570/k203123/prj/dolpung_coupled/icon-mpim/experiments/sap0006/outdata/sap0006.zarr

Zoom 8 (nside 256, nested, ~0.23 deg) is used, the closest to the 0.25 deg
grids of the other datasets.

  sst   P1D_mean_z8/to, top ocean level, averaged to JJA monthly means
        -> <work>/sst/sst_<year>.nc        (time=3, cell=786432)
  tmax  PT3H_point_z8/tas, maximum of the 8 three-hourly instantaneous values
        of each day (stamps 03:00 .. 24:00 UTC; the run has no daily-maximum
        output), April-October, cells inside EU_BOX only (NaN elsewhere)
        -> <work>/tmax_hp/tmax_<year>.nc   (time=214, cell=786432)

Daily means are stamped at the END of the day (time[0] = 1950-01-02 00:00 is
the 1 January mean), so array index i is day 1950-01-01 + i.

Needs zarr, which the project env lacks: run with the shared mambaforge base
python (see 07_sap0006.sh). Usage: extract_sap0006.py {sst|tmax} <work> <year>...
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import zarr

STORE = ("/work/mh1570/k203123/prj/dolpung_coupled/icon-mpim/experiments/"
         "sap0006/outdata/sap0006.zarr")
T0 = pd.Timestamp("1950-01-01")
NCELL = 786432
#: a little wider than the 0.25 deg Europe grid (15W-45E, 30-72N)
EU_BOX = (-18.0, 48.0, 27.0, 75.0)


def day_index(date):
    return (pd.Timestamp(date) - T0).days


def sst_year(g, work, year):
    out = work / "sst" / f"sst_{year}.nc"
    if out.exists():
        return
    to = g["P1D_mean_z8/to"]
    fields, times = [], []
    for m in (6, 7, 8):
        d0 = day_index(f"{year}-{m:02d}-01")
        d1 = day_index(pd.Timestamp(f"{year}-{m:02d}-01") + pd.offsets.MonthBegin(1))
        fields.append(np.nanmean(to[d0:d1, 0, :], axis=0))
        times.append(pd.Timestamp(f"{year}-{m:02d}-15"))
    ds = xr.Dataset({"to": (("time", "cell"), np.stack(fields).astype("f4"))},
                    coords={"time": times})
    ds["to"].attrs.update(units="degC", long_name="sea water potential temperature, top level")
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, encoding={"to": {"zlib": True, "complevel": 1, "_FillValue": np.nan},
                                "time": {"units": "days since 1850-01-01", "dtype": "f8"}})


def eu_cells(work):
    """
    Indices of the zoom-8 cells inside EU_BOX. Cell centres come from
    hp256_centres.nc, which 07_sap0006.sh makes by nearest-neighbour remapping
    0.1 deg lon/lat fields onto hp256_nested (cdo cannot print HEALPix
    coordinates directly); good to ~0.05 deg, ample for a box selection.
    """
    with xr.open_dataset(work / "hp256_centres.nc") as c:
        lon = (c.clon.values + 180) % 360 - 180
        lat = c.clat.values
    lon0, lon1, lat0, lat1 = EU_BOX
    return np.flatnonzero((lon >= lon0) & (lon <= lon1) & (lat >= lat0) & (lat <= lat1))


def tmax_year(g, work, year, cells):
    out = work / "tmax_hp" / f"tmax_{year}.nc"
    if out.exists():
        return
    tas = g["PT3H_point_z8/tas"]
    d0, d1 = day_index(f"{year}-04-01"), day_index(f"{year}-11-01")
    c0, c1 = int(cells.min()), int(cells.max()) + 1
    block = tas[8 * d0:8 * d1, c0:c1]                     # (3-hourly, cell range)
    tmax = block.reshape(d1 - d0, 8, c1 - c0).max(axis=1)[:, cells - c0]
    full = np.full((d1 - d0, NCELL), np.nan, dtype="f4")
    full[:, cells] = tmax - 273.15
    times = pd.date_range(f"{year}-04-01", periods=d1 - d0, freq="D")
    ds = xr.Dataset({"tmax": (("time", "cell"), full)}, coords={"time": times})
    ds["tmax"].attrs.update(units="degC",
                            long_name="daily maximum of 3-hourly instantaneous 2 m temperature")
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, encoding={"tmax": {"zlib": True, "complevel": 1, "_FillValue": np.nan},
                                "time": {"units": "days since 1850-01-01", "dtype": "f8"}})


def main():
    what, work, years = sys.argv[1], Path(sys.argv[2]), [int(y) for y in sys.argv[3:]]
    g = zarr.open_group(STORE, mode="r")
    cells = eu_cells(work) if what == "tmax" else None
    for y in years:
        sst_year(g, work, y) if what == "sst" else tmax_year(g, work, y, cells)
        print(f"  {what} {y}", flush=True)


if __name__ == "__main__":
    main()
