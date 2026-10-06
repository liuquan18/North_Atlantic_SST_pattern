"""
Figure 1 (heatwave analogues) -- the same layout as fig1_patterns_simulations.py,
but each model's season is chosen by its heat instead of by its SST pattern.

For every simulation with daily Tmax, the season (and member) with the most
JJA heatwave days, averaged over the European land south of SOUTH_LAT, is
selected. Only seasons that also have an SST pattern are eligible (EERIE r2
has Tmax to 2020 but SST only to 2014). The columns then show

  top row     heatwave days of that season over the whole of Europe, with the
              southern land mean that selected it (dashed line: SOUTH_LAT)
  bottom row  the SST pattern of the same season (box mean removed), with its
              pattern correlation against ERA5 JJA 2026

So it asks the converse of fig1_patterns_simulations.py: when a model produces
its hottest southern-European summer, does the ocean look like 2026?
MPI-ESM1.2-ER has no daily Tmax, so it has no season to select and is left out.
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, column, land_mean, land_mean_series, plot_figure

KEYS = ["MPI-GE", "ICON-EPOC-hist", "EERIE"]
SOUTH_LAT = 60.0
SOUTH_EXTENT = vz.EUROPE_EXTENT[:3] + [SOUTH_LAT]


def south_mean(hwd):
    """Area-weighted land-mean heatwave days south of SOUTH_LAT, per [member,] year."""
    return land_mean_series(hwd, SOUTH_EXTENT)


def hottest_season(key, corr):
    """(year, member or None, southern heatwave days, r) of the hottest scored season."""
    with xr.open_dataset(hw.metrics_file(key)) as ds:
        s = south_mean(ds.hwd.load())
    r = corr[key]
    mdim = p2.member_dim(r)
    if mdim:
        r = r.rename({mdim: "member"})
    s, r = xr.align(s, r, join="inner")
    s = s.where(r.notnull())            # only seasons that have an SST pattern
    if "member" in s.dims:
        z = s.stack(z=("member", "year"))
        best = z.isel(z=int(z.fillna(-1).argmax("z")))
        year, mem = int(best["year"]), int(best["member"])
        return year, mem, float(best), float(r.sel(member=mem, year=year))
    best = s.isel(year=int(s.fillna(-1).argmax("year")))
    year = int(best["year"])
    return year, None, float(best), float(r.sel(year=year))


def main():
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    with xr.open_dataset(hw.metrics_file("ERA5")) as ds:
        obs = float(south_mean(ds.hwd.sel(year=p2.REF_YEAR).load()))

    cols = []
    for key in KEYS:
        year, mem, days, r = hottest_season(key, corr)
        d = p2.DATASETS[key]
        when = f"JJA {year}" + (f"  ·  member r{mem}" if mem else "")
        print(f"{key:16s} hottest south-European season: {when}  {days:.1f} days  r = {r:+.2f}")
        c = column(key, year, mem, r, f"{d.label}\n{d.resolution}\n{when}")
        c["hw_tag"] = f"<{SOUTH_LAT:.0f}°N: {days:.1f} days\nall: {land_mean(c['hwd']):.1f} days"
        cols.append(c)

    plot_figure(
        cols,
        "European summer heatwaves and the North Atlantic + Mediterranean SST pattern: "
        f"each model's hottest southern-European summer",
        p2.FIG_DIR / "fig1_heatwave_analogues_2026.png",
        hw_note=f"\nseason selected by land-mean heatwave days south of {SOUTH_LAT:.0f}°N "
                f"(dashed line; ERA5 JJA {p2.REF_YEAR}: {obs:.1f} days), over all members and "
                "seasons with an SST pattern  ·  top-row values: that mean, and over the whole map",
        lat_line=SOUTH_LAT)


if __name__ == "__main__":
    main()
