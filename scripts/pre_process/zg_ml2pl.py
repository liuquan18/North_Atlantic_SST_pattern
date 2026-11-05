"""
Geopotential height on one pressure level from monthly model-level fields:
log-pressure linear interpolation between the two levels around it, column by
column. Used by 12_epoc_ctrl_circ.sh, for runs without pressure-level output.

Usage: zg_ml2pl.py <in.nc with gpsm, pfull> <out.nc> <level in hPa>
"""
import sys

import numpy as np
import pandas as pd
import xarray as xr

G = 9.80665


def to_level(gpsm, pfull, p0):
    """gpsm / g (gpm) at pressure p0 (Pa); both (height, lat, lon), pressure rising with index."""
    p, z = pfull.values, gpsm.values / G
    k = (p < p0).sum(axis=0) - 1
    if k.min() < 0 or k.max() + 1 >= p.shape[0]:
        raise ValueError(f"{p0} Pa is outside the levels read: {p.min():.0f}..{p.max():.0f} Pa")
    lat, lon = np.ogrid[:p.shape[1], :p.shape[2]]
    w = (np.log(p0) - np.log(p[k, lat, lon])) / (np.log(p[k + 1, lat, lon]) - np.log(p[k, lat, lon]))
    return z[k, lat, lon] + w * (z[k + 1, lat, lon] - z[k, lat, lon])


def main(src, out, level):
    with xr.open_dataset(src) as ds:
        ds = ds.load()
    fields = [to_level(ds.gpsm.isel(time=i), ds.pfull.isel(time=i), level * 100.0)
              for i in range(ds.sizes["time"])]
    # monthly means are stamped at the END of the month
    time = pd.to_datetime(ds.time.values) - pd.Timedelta(days=1)
    da = xr.DataArray(np.stack(fields).astype("f4"), name=f"z{level}", attrs=dict(units="gpm"),
                      coords=dict(time=time, lat=ds.lat, lon=ds.lon), dims=("time", "lat", "lon"))
    da.to_dataset().to_netcdf(out, encoding={"time": {"units": "days since 1850-01-01",
                                                      "dtype": "f8"}})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
