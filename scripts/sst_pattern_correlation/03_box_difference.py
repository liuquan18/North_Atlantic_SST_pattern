"""
Box-difference index of the JJA 2026 SST pattern for each dataset:
area-mean JJA SST anomaly over the Mediterranean box (30-45N, 0-20E) minus
that over the mid-latitude North Atlantic box (40-50N, 60-20W), see
src/pattern_2026.DIPOLE_BOXES.

A simpler measure of 2026-likeness than the pattern correlation: it asks only
whether a summer has the basin-scale contrast of 2026, not whether every
eddy-scale feature sits in the same place, so it does not penalise the coarse
models for their resolution.

Both box means are taken on the common 1 deg ocean mask of 01_pattern_corr.py,
as in 02_region_means.py. Anomalies are against each dataset's own 1991-2020.

Output: data/sst_pattern_corr/box_difference.nc
        variables '<dataset>' (the index, degC) and '<dataset>|nal', '<dataset>|med'
"""
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2


def main():
    mask = xr.open_dataset(p2.RESULT_DIR / "ocean_mask_1deg.nc")["ocean_mask"]
    # read only the 0.25 deg block that covers both boxes; its edges are whole
    # degrees, so the 4x4 coarsening lands on the same 1 deg cells as the mask
    lat = (min(b["lat"][0] for b in p2.DIPOLE_BOXES.values()),
           max(b["lat"][1] for b in p2.DIPOLE_BOXES.values()))
    lon = (min(b["lon"][0] for b in p2.DIPOLE_BOXES.values()),
           max(b["lon"][1] for b in p2.DIPOLE_BOXES.values()))
    mask = mask.sel(lat=slice(*lat), lon=slice(*lon))

    keys = [k for k in p2.DATASETS if p2.DATASETS[k].anom_file.exists()]
    series = {}
    for k in keys:
        print(f"  {k} ...", flush=True)
        coarse = p2.coarsen_to_1deg(
            p2.load_anomaly(k).sel(lat=slice(*lat), lon=slice(*lon))).compute()
        np.testing.assert_allclose(coarse.lat, mask.lat)
        np.testing.assert_allclose(coarse.lon, mask.lon)
        coarse = coarse.assign_coords(lat=mask.lat, lon=mask.lon)
        parts = p2.dipole_index(coarse, mask)
        series[k] = parts["index"]
        series[f"{k}|nal"] = parts["nal"]
        series[f"{k}|med"] = parts["med"]

    # per-ensemble member dimension; see src.pattern_2026.member_dim
    ds = xr.Dataset({
        name: (da.rename({"member": f"member_{name.split('|')[0]}"})
               if "member" in da.dims else da)
        for name, da in series.items()
    })
    ds.attrs.update(
        description=(f"JJA SST anomaly box difference {p2.DIPOLE_LABEL} (degC, 1991-2020 base); "
                     "'<dataset>' is the index, '<dataset>|nal' / '|med' the box means"),
        **{f"box_{k}": f"{b['label']}, {b['detail']}" for k, b in p2.DIPOLE_BOXES.items()})
    ds.to_netcdf(p2.dipole_file())
    print(f"\nwrote {p2.dipole_file()}")

    era5 = ds["ERA5"].to_series()
    clim = era5.loc[slice(*p2.CLIM_PERIOD)]
    print(f"ERA5 {p2.REF_YEAR}: {era5[p2.REF_YEAR]:+.2f} degC "
          f"(nal {float(ds['ERA5|nal'].sel(year=p2.REF_YEAR)):+.2f}, "
          f"med {float(ds['ERA5|med'].sel(year=p2.REF_YEAR)):+.2f}); "
          f"{p2.CLIM_PERIOD[0]}-{p2.CLIM_PERIOD[1]} sd {clim.std():.2f}")
    print("ERA5 top 5:", ", ".join(f"{y} {v:+.2f}" for y, v in era5.nlargest(5).items()))


if __name__ == "__main__":
    main()
