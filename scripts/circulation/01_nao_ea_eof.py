"""
The summer NAO and East Atlantic (EA) patterns and their PC indices, from an
EOF analysis of JJA mean sea-level pressure over the North Atlantic-European
sector (20-80N, 90W-40E, 2.5 deg; data/circulation/, pre_process/11_slp.sh),
separately for ERA5 and MPI-GE, over FIT_YEARS (1940-2026).

Anomalies
  ERA5    the JJA SLP with a linear trend over 1940-2026 removed at every grid
          point (which also removes the mean).
  MPI-GE  each member minus the 50-member mean of the same year: the forced
          signal (and the climatology) goes, the internal variability stays.
          The 4350 member-summers of 1940-2026 are pooled into one EOF analysis.

EOFs: SVD of the anomalies weighted by sqrt(cos lat), so each grid point
counts by its area. The PC of a mode is the projection of the weighted
anomalies onto its EOF, standardised by its standard deviation over the
FIT_YEARS samples; MPI-GE summers outside them (1935-1939, 2027-2030, for the
10-year windows of the index figure) are projected onto the same EOFs. The
pattern is the regression of the unweighted SLP anomaly on the standardised
PC: hPa per one standard deviation of the index.

Which mode is which (MODES) and its sign (SIGN_BOX) were set by inspecting the
leading four patterns of both datasets: a mode is flipped so that the
regression averaged over its SIGN_BOX has the sign given there, so a positive
index always goes with the positive pattern.

Output: data/circulation/nao_ea.nc with the PCs (<dataset>|<mode>), patterns
(pattern_<dataset>|<mode>) and explained-variance fractions (attributes) of
NAO and EA, and of the leading N_EOF modes (eof<k>) for reference.
"""
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2

DIR = p2.PROJECT_ROOT / "data" / "circulation"
OUT = DIR / "nao_ea.nc"
FIT_YEARS = (1940, p2.REF_YEAR)
N_EOF = 4
DOMAIN_LABEL = "20–80°N, 90°W–40°E"
#: dataset -> {mode name: EOF number (0-based)}
MODES = {"ERA5": {"NAO": 0, "EA": 1}, "MPI-GE": {"NAO": 0, "EA": 1}}
#: mode -> (box, sign of the regression averaged over it in the positive phase)
SIGN_BOX = {
    # positive NAO: low pressure over Greenland / Iceland
    "NAO": (dict(lat=(60.0, 80.0), lon=(-60.0, -10.0)), -1),
    # positive EA: low pressure at its main centre (the CPC sign convention of the
    # EA); in JJA the centre sits over the British Isles / North Sea, ~55N 0-5E,
    # east of its winter position west of Ireland
    "EA": (dict(lat=(45.0, 65.0), lon=(-20.0, 20.0)), -1),
}


def load():
    """(ERA5 anomalies (year, lat, lon), MPI-GE anomalies (member, year, lat, lon))."""
    with xr.open_dataset(DIR / "era5_jja_slp_2p5deg.nc") as ds:
        era5 = ds.slp.sel(year=slice(*FIT_YEARS)).load()
    t = era5.year.values.astype(float)
    coef = np.polyfit(t, era5.values.reshape(len(t), -1), 1)
    trend = (np.outer(t, coef[0]) + coef[1]).reshape(era5.shape)
    era5 = era5 - trend
    with xr.open_dataset(DIR / "mpige_jja_slp_2p5deg.nc") as ds:
        ge = ds.slp.load()
    ge = ge - ge.mean("member")
    return era5, ge


def sample_dims(da):
    return [d for d in da.dims if d not in ("lat", "lon")]


def eof(anom, n=N_EOF):
    """
    The leading `n` modes of `anom` fitted over FIT_YEARS: standardised PCs over
    all its samples, regression patterns (hPa per std), explained variance.
    """
    dims = sample_dims(anom)
    w = np.sqrt(np.cos(np.deg2rad(anom.lat)))
    xw = (anom * w).stack(sample=dims).transpose("sample", "lat", "lon")
    fit = xw.sel(sample=xw.year.isin(range(FIT_YEARS[0], FIT_YEARS[1] + 1)))
    x = fit.values.reshape(fit.sizes["sample"], -1)
    _, s, vt = np.linalg.svd(x, full_matrices=False)
    var = s ** 2 / np.sum(s ** 2)
    pcs = xw.values.reshape(xw.sizes["sample"], -1) @ vt[:n].T
    in_fit = xw.year.isin(range(FIT_YEARS[0], FIT_YEARS[1] + 1)).values
    pcs = pcs / pcs[in_fit].std(axis=0)
    pc = xr.DataArray(pcs, dims=("sample", "mode"),
                      coords=dict(sample=xw.sample, mode=np.arange(n))).unstack("sample")
    a = anom.sel(year=slice(*FIT_YEARS))
    p = pc.sel(year=slice(*FIT_YEARS))
    pattern = (a * p).mean(dims) / (p * p).mean(dims)
    # North et al. (1982) sampling error of the eigenvalues, N = number of fitted samples
    err = var[:n] * np.sqrt(2 / x.shape[0])
    return pc.transpose("mode", *dims), pattern.transpose("mode", "lat", "lon"), var[:n], err


def box_mean(da, box):
    sub = da.sel(lat=slice(*box["lat"]), lon=slice(*box["lon"]))
    w = np.cos(np.deg2rad(sub.lat)) * xr.ones_like(sub.lon)
    return (sub * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))


def main():
    era5, ge = load()
    out, attrs = {}, {}
    for key, anom in (("ERA5", era5), ("MPI-GE", ge)):
        pc, pattern, var, err = eof(anom)
        for k in range(N_EOF):
            out[f"pattern_{key}|eof{k + 1}"] = pattern.sel(mode=k, drop=True)
            out[f"{key}|eof{k + 1}"] = pc.sel(mode=k, drop=True)
            attrs[f"var_{key}|eof{k + 1}"] = float(var[k])
        print(f"{key}: explained variance "
              + ", ".join(f"EOF{k + 1} {v:.1%} ± {e:.1%}" for k, (v, e) in enumerate(zip(var, err))))
        for name, k in MODES[key].items():
            box, sign = SIGN_BOX[name]
            s = sign * np.sign(float(box_mean(pattern.sel(mode=k), box)))
            out[f"pattern_{key}|{name}"] = s * pattern.sel(mode=k, drop=True)
            out[f"{key}|{name}"] = s * pc.sel(mode=k, drop=True)
            attrs[f"var_{key}|{name}"] = float(var[k])
            check = float(box_mean(out[f"pattern_{key}|{name}"], box))
            print(f"  {name} = EOF{k + 1} ({var[k]:.1%}), sign {'kept' if s > 0 else 'flipped'}; "
                  f"regression in its sign box {check:+.2f} hPa per std")
    ds = xr.Dataset({k: v.rename({"member": "member_MPI-GE"}) if "member" in v.dims else v
                     for k, v in out.items()})
    ds.attrs = dict(description="JJA SLP EOFs over " + DOMAIN_LABEL + f", {FIT_YEARS[0]}-"
                    f"{FIT_YEARS[1]}; ERA5 linearly detrended, MPI-GE ensemble mean removed; "
                    "PCs standardised, patterns = regression in hPa per std", **attrs)
    DIR.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(OUT)
    print(f"wrote {OUT}")
    # the projection onto the other dataset's pattern: do ERA5 and MPI-GE agree?
    for name in ("NAO", "EA"):
        r = p2.pattern_corr(ds[f"pattern_ERA5|{name}"], ds[f"pattern_MPI-GE|{name}"],
                            centered=False)
        print(f"  {name}: ERA5 vs MPI-GE pattern (uncentred) r = {float(r):+.2f}")


if __name__ == "__main__":
    main()
