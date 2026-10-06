"""
Figure 1 (simulations) -- European summer heatwaves and the North Atlantic +
Mediterranean SST pattern in the season of each model that best matches
ERA5 JJA 2026.

Columns: the best-matching season of MPI-GE, MPI-ESM1.2-ER, km-scale ICON
(EPOC) and EERIE. The ERA5 seasons are in fig1_patterns_era5.py, which draws
with plot_figure() from here. MPI-GE's heatwave panel is on the model's native ~1.9 deg grid, so its cells
are drawn at the resolution the model actually has.

  top row     JJA heatwave days over European land (Xu et al. 2026 definition,
              1991-2020 reference; scripts/european_heatwave/)
  bottom row  the SST pattern of the same season (box mean removed), with its
              pattern correlation against ERA5 JJA 2026

So each column asks: in the season whose SST pattern looks most like 2026,
what did European summer heat look like? MPI-ESM1.2-ER has no daily
atmosphere output, so its heatwave panel is left empty.

The SST panels share one symmetric diverging scale, capped below the
Mediterranean extremes: letting +3 degC Mediterranean cells set the range
washes the Atlantic part of the pattern out to near-white.
"""
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

VARIANT = p2.PLOT_VARIANTS[0]
VMAX_CAP = 2.0
#: (dataset, which season): "ref" = JJA 2026 itself, "best" = best analogue
COLUMNS = [("MPI-GE", "best"), ("MPI-ER", "best"), ("ICON-EPOC-hist", "best"), ("EERIE", "best")]
HWD_LEVELS = np.arange(0, 46, 5)


def season(key, which, best):
    """(year, member or None, r or None) of the season shown in one column."""
    if which == "ref":
        return p2.REF_YEAR, None, None
    row = best[(best.variant == VARIANT) & (best.dataset == key)].iloc[0]
    mem = int(row.best_member) if pd.notna(row.best_member) else None
    return int(row.best_year), mem, float(row.best_r)


def sst_field(key, year, member):
    sel = {"year": year}
    if member is not None:
        sel["member"] = member
    anom, gm = p2.load_anomaly(key), p2.load_gmsst(key)
    return p2.make_pattern(anom.sel(sel), gm.sel(sel), VARIANT).load()


def heatwave_season(key, year, member):
    """All heatwave metrics of one season (Dataset), or None without daily Tmax."""
    if key not in hw.DATASETS or not hw.metrics_file(key).exists():
        return None
    with xr.open_dataset(hw.metrics_file(key)) as ds:
        if year not in ds.year:
            return None
        ds = ds.drop_vars(["threshold", "climatology"]).sel(year=year)
        if "member" in ds.dims:
            ds = ds.sel(member=member)
        return ds.load()


def land_mean_series(field, extent=None):
    """
    Area-weighted mean over the land cells inside the drawn map extent, per
    season ([member,] year). Seasons with no data (a member not run that
    year) come out NaN, not 0.
    """
    lon0, lon1, lat0, lat1 = extent or vz.EUROPE_EXTENT
    inside = ((field.lat >= lat0) & (field.lat <= lat1) & (field.lon >= lon0) & (field.lon <= lon1))
    w = np.cos(np.deg2rad(field.lat)) * (field.notnull() & inside)
    wsum = w.sum(("lat", "lon"))
    return (field.fillna(0) * w).sum(("lat", "lon")) / wsum.where(wsum > 0)


def land_mean(field, extent=None):
    """Area-weighted mean over the land cells inside the drawn map extent."""
    return float(land_mean_series(field, extent))


def column(key, year, member, r, title):
    """One figure column: its title, SST pattern, heatwave days and r."""
    hws = heatwave_season(key, year, member)
    if hws is None:
        print(f"!! {key} JJA {year}: no heatwave field")
    return dict(title=title, sst=sst_field(key, year, member), r=r,
                hwd=None if hws is None else hws.hwd)


def plot_figure(cols, suptitle, out, *, hw_note="top-row value: land-mean heatwave days",
                lat_line=None):
    """
    Heatwave days (top) over SST pattern (bottom), one column per season.

    A column may carry its own "hw_tag" for the heatwave panel (default: the
    land-mean heatwave days) and "r_tag" for the SST panel (default: r); `lat_line` draws a dashed parallel on the
    heatwave maps, e.g. to mark a sub-region the tag refers to.
    """
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)

    vmax = min(VMAX_CAP, max(float(np.nanpercentile(np.abs(c["sst"].values), 98)) for c in cols))
    vmax = np.round(vmax * 2) / 2
    sst_levels = np.linspace(-vmax, vmax, 21)
    sst_cmap = vz.anomaly_cmap()
    hw_cmap = vz.heatwave_cmap()
    hw_norm = BoundaryNorm(HWD_LEVELS, hw_cmap.N, extend="max")

    n = len(cols)
    # Rows are sized by height, not width: the columns are as wide as the SST
    # maps, and the top row is 1.3x their height, so the near-square Europe
    # map fills its row height and sits centred, narrower than the SST map.
    panel_w = 3.1
    h_bot = panel_w / vz.panel_aspect()
    h_top = 1.3 * h_bot
    fig_w = panel_w * n + 2.0
    # head room for the column titles: ~0.18 in per title line beyond the first
    title_h = 0.18 * max(c["title"].count("\n") for c in cols)
    # ...and for caption lines beyond the default two (one is in hw_note's place)
    caption_h = 0.17 * hw_note.count("\n")
    top_h = 1.39 + title_h + caption_h
    fig_h = h_top + h_bot + top_h + 0.65

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(2, n + 1, figure=fig,
                  width_ratios=[1] * n + [0.035], height_ratios=[h_top, h_bot],
                  hspace=0.32, wspace=0.07,
                  left=0.075, right=0.94, top=1 - top_h / fig_h, bottom=0.42 / fig_h)

    im_hw = None
    for j, c in enumerate(cols):
        # --- top: heatwave days -------------------------------------------
        ax = fig.add_subplot(gs[0, j], projection=vz.europe_projection())
        im = vz.draw_land_map(ax, c["hwd"], hw_cmap, hw_norm, lon_step=20,
                              labels_bottom=True, labels_left=(j == 0))
        im_hw = im or im_hw
        # the season is shared by both rows, so it heads the column; the
        # in-map tags then carry only each row's own number
        ax.set_title(c["title"], pad=5, fontsize=8.6, color=vz.INK, linespacing=1.3)
        if c["hwd"] is None:
            ax.text(0.5, 0.5, "no daily\nTmax output", transform=ax.transAxes,
                    ha="center", va="center", fontsize=7.5, color=vz.INK_SOFT, zorder=6,
                    bbox=dict(facecolor="#ffffff", alpha=0.9, edgecolor="none",
                              boxstyle="round,pad=0.4"))
        else:
            vz.panel_tag(ax, c.get("hw_tag") or f"{land_mean(c['hwd']):.1f} days",
                         loc="upper right")
        if lat_line is not None:
            ax.plot([-180, 180], [lat_line, lat_line], transform=ccrs.PlateCarree(),
                    color=vz.INK, lw=0.9, ls=(0, (4, 3)), zorder=5)

        # --- bottom: SST pattern ------------------------------------------
        ax = fig.add_subplot(gs[1, j], projection=vz.map_projection())
        im_sst = vz.draw_map(ax, c["sst"], sst_levels, sst_cmap,
                             labels_bottom=True, labels_left=(j == 0))
        if c["r"] is not None:
            vz.panel_tag(ax, c.get("r_tag") or f"r = {c['r']:+.2f}", loc="upper right")

    for i, label in enumerate(["JJA heatwave days", "JJA SST pattern"]):
        fig.text(0.012, gs[i, 0].get_position(fig).y0 + gs[i, 0].get_position(fig).height / 2,
                 label, rotation=90, va="center", ha="center", fontsize=9.5, color=vz.INK)

    cb = fig.colorbar(im_hw, cax=fig.add_subplot(gs[0, -1]), extend="max", ticks=HWD_LEVELS)
    cb.set_label("heatwave days in JJA", fontsize=8.2, color=vz.INK_SOFT)
    cb = fig.colorbar(im_sst, cax=fig.add_subplot(gs[1, -1]), extend="both",
                      ticks=np.linspace(-vmax, vmax, 5))
    cb.set_label("SST anomaly (°C)", fontsize=8.2, color=vz.INK_SOFT)
    for ax in fig.axes[-2:]:
        ax.tick_params(labelsize=7.5)

    fig.suptitle(suptitle, fontsize=12, color=vz.INK, y=1 - 0.30 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             "heatwave: ≥3 consecutive days with Tmax anomaly above the calendar-day 90th percentile "
             "(15-day window), detected May–Sep, land only (Xu et al. 2026)  ·  "
             f"{hw_note}\n"
             f"SST: {p2.VARIANTS[VARIANT]['long']}, 30–60°N, 80°W–40°E, r = pattern correlation with ERA5 JJA 2026  ·  "
             "both against 1991–2020 of the same dataset",
             ha="center", va="top", fontsize=8.3, color=vz.INK_SOFT, linespacing=1.4)

    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}  (SST scale ±{vmax} °C)")


def main():
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")
    cols = []
    for key, which in COLUMNS:
        year, mem, r = season(key, which, best)
        d = p2.DATASETS[key]
        when = f"JJA {year}" + (f"  ·  member r{mem}" if mem else "")
        cols.append(column(key, year, mem, r, f"{d.label}\n{d.resolution}\n{when}"))
    plot_figure(cols, "European summer heatwaves and the North Atlantic + Mediterranean SST pattern: "
                      f"the season of each model closest to ERA5 JJA {p2.REF_YEAR}",
                p2.FIG_DIR / "fig1_patterns_2026_simulations.png")


if __name__ == "__main__":
    main()
