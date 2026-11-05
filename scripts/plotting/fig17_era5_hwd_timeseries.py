"""
Figure 17 -- communication figure: ERA5 summer heatwave days over southern
Europe and the resemblance of each summer's SST to JJA 2026, 1940-2026, two
panels on a shared time axis.

  a  heatwave days. The series as recorded is drawn very light; the series
     with the long-term trend removed is the solid, highlighted one. JJA 2003
     and 2026 are labelled with their year only. An inset map shows the JJA
     2026 heatwave days per grid cell, with the averaging box outlined.
  b  pattern correlation of each ERA5 summer's SST against JJA 2026 (the
     "spatial" variant of figure 2: box mean removed, common 1 deg ocean grid).
     The line runs through 2026, which against itself is 1 by construction.
     An inset map shows the JJA 2026 SST pattern being matched.

Two versions of the same figure (panel b is the same in both):
  _line  both heatwave series as lines
  _bar   trend-removed series as bars, recorded series as a light line

  series  JJA heatwave days, area mean over southern-European land
          (src.heatwave.MEAN_BOX), Xu et al. (2026) definition against the
          1991-2020 calendar-day 90th percentile of daily Tmax
  trend   quadratic least-squares fit in time over 1940-2026. The record is
          flat (~2 days) until about 1990 and then rises steeply, so a straight
          line would sit too high early on and too low recently, and would rank
          the 1940s as unusual as 2003. A quadratic follows that bend with one
          extra parameter (residual sd 2.1 days vs 2.9 for linear); LOWESS or a
          regression on 4-yr smoothed global-mean SST move each season by
          ~1.5 days at most. Which of 2026 and 2003 ranks first above the trend
          depends on that choice (LOWESS with a 30% span flips them).

Publication layout: 180 mm wide, caption set below the plot, white
background, no fill over the sea in the heatwave map, 300 dpi PNG plus vector
PDF.
"""
import sys
import textwrap
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.feature_artist import FeatureArtist
from cartopy.mpl.gridliner import Gridliner
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, VMAX_CAP, land_mean_series, sst_field

MARKED = [2003, p2.REF_YEAR]
TREND_DEGREE = 2

HOT = vz.heatwave_cmap()(0.80)    # the highlighted, trend-removed series
LIGHT = "#c9c7c0"                 # the recorded series, context only
CORR = vz.COLORS["ERA5"]          # the pattern-correlation series, as in fig 2
#: inset maps: same colours and steps as the maps of fig 1
MAP_LEVELS = np.arange(0, 46, 5)
SST_LEVELS = np.linspace(-VMAX_CAP, VMAX_CAP, 21)
#: heatwave map extent: the averaging box plus a margin, so its dashed outline
#: is not clipped at the map edge
MAP_EXTENT = [-13, 34, 33.5, 62]
#: inset positions in their time-series axes' fraction, both top left, above
#: the low 1940-1965 values. The SST map is ~2.2:1, so it is drawn wider and
#: flatter than the near-square Europe map.
INSET = (0.005, 0.50, 0.30, 0.50)
SST_INSET = (0.005, 0.64, 0.26, 0.36)
#: correlation axis: room above r = 1 for the 2026 label, and the 1940s values
#: (at most +0.39) kept below the inset
CORR_YLIM = (-0.75, 1.2)

#: publication layout: two-column width (180 mm), 7-8 pt type; panel heights,
#: gaps and the caption block in inches
FS = 7.5
FIG_W = 7.09
PANEL_H = (2.35, 2.1)
GAP_H = 0.15
TOP_H = 0.1
XTICK_H = 0.3
CAPTION_H = 1.45
FIG_H = TOP_H + sum(PANEL_H) + GAP_H + XTICK_H + CAPTION_H
CAPTION_TITLE = f"Summer heatwave days in southern Europe and the JJA {p2.REF_YEAR} SST pattern, 1940–2026."
CAPTION = (
    f"(a) June–August heatwave days in ERA5, averaged over land in {hw.MEAN_BOX_LABEL} "
    "(dashed box in the inset). Grey: observed. Red: departure from a quadratic long-term "
    "trend fitted to 1940–2026. A heatwave is a period of at least three consecutive days "
    "with daily maximum temperature above the calendar-day 90th percentile of 1991–2020 "
    f"(Xu et al. 2026). Inset: heatwave days in JJA {p2.REF_YEAR} per grid cell. "
    "(b) Area-weighted pattern correlation between each summer's JJA SST anomaly in ERA5 "
    f"(against 1991–2020, mean over the region removed) and that of JJA {p2.REF_YEAR}, "
    f"over {p2.REGION_LABEL} on a 1° grid; {p2.REF_YEAR} itself is 1 by construction. "
    f"Inset: the JJA {p2.REF_YEAR} SST pattern.")


def load():
    with xr.open_dataset(hw.metrics_file("ERA5")) as ds:
        s = land_mean_series(ds.hwd.load()).to_series()
    year = s.index.values.astype(float)
    trend = np.polyval(np.polyfit(year, s.values, TREND_DEGREE), year)
    return s.index.values, s.values, s.values - trend


def load_map(year=p2.REF_YEAR):
    """Heatwave days of one ERA5 season per land cell (NaN over sea)."""
    with xr.open_dataset(hw.metrics_file("ERA5")) as ds:
        return ds.hwd.sel(year=year).load()


def load_corr():
    """ERA5 pattern correlation with JJA 2026, per season (2026 itself is 1)."""
    with xr.open_dataset(p2.corr_file(VARIANT)) as ds:
        r = ds["ERA5"].dropna("year").to_series()
    return r.index.values, r.values


def inset_axes(fig, ax, rect, projection):
    """Map axes at `rect` (fraction of `ax`), with no graticule or frame."""
    x0, y0, w, h = ax.get_position().bounds
    ix, iy, iw, ih = rect
    return fig.add_axes([x0 + ix * w, y0 + iy * h, iw * w, ih * h], projection=projection)


def finish_inset(fig, m, im, label, ticks):
    """Strip graticule and frame from inset map `m`; colour bar to its right."""
    for art in list(m.artists):      # no graticule at this size
        if isinstance(art, Gridliner):
            art.remove()
    m.spines["geo"].set_visible(False)
    # the map keeps its aspect inside the box it was given: hang the colour
    # bar on the extent it actually drew
    fig.canvas.draw()
    mx0, my0, mw, mh = m.get_position(original=False).bounds
    cax = fig.add_axes([mx0 + mw + 0.008, my0 + 0.1 * mh, 0.007, 0.8 * mh])
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    cb.set_label(label, fontsize=FS - 1, labelpad=7, color=vz.INK)
    cb.ax.tick_params(labelsize=FS - 1.5, length=1.5, width=0.5, pad=1.5)
    cb.outline.set_linewidth(0)


def inset_map(fig, ax, field):
    """Heatwave days of `field` at INSET in `ax`."""
    m = inset_axes(fig, ax, INSET, ccrs.PlateCarree())
    cmap = vz.heatwave_cmap()
    norm = BoundaryNorm(MAP_LEVELS, cmap.N, extend="max")
    im = vz.draw_land_map(m, field, cmap, norm, MAP_EXTENT)
    # no fill over the sea: the coastline alone outlines the land
    for art in list(m.collections) + list(m.artists):
        if isinstance(art, FeatureArtist) and art._kwargs.get("facecolor") == vz.OCEAN:
            art.remove()
    vz.draw_box(m, hw.MEAN_BOX, lw=0.7)
    finish_inset(fig, m, im, f"heatwave days\nJJA {p2.REF_YEAR}", MAP_LEVELS[::2])


def inset_sst(fig, ax, field):
    """SST pattern `field` at SST_INSET in `ax`, over the whole analysis box."""
    m = inset_axes(fig, ax, SST_INSET, vz.map_projection())
    im = vz.draw_map(m, field, SST_LEVELS, vz.anomaly_cmap())
    # no colour over land: the coastline alone outlines it, as in panel a's
    # map. Land is painted the page white rather than left unfilled, so it still
    # covers the ERA5 lake cells (Great Lakes, Hudson Bay) under it.
    m.set_facecolor("white")
    for art in list(m.collections) + list(m.artists):
        if isinstance(art, FeatureArtist) and art._kwargs.get("facecolor") == vz.LAND:
            art.remove()
    m.add_feature(cfeature.LAND, facecolor="white", edgecolor="none", zorder=2)
    finish_inset(fig, m, im, f"SST pattern (°C)\nJJA {p2.REF_YEAR}",
                 np.linspace(-VMAX_CAP, VMAX_CAP, 5))


def caption(fig):
    """Bold title sentence, then the caption body, wrapped to the plot width."""
    x = 0.012
    y = (CAPTION_H - 0.1) / FIG_H
    fig.text(x, y, CAPTION_TITLE, ha="left", va="top", fontsize=FS,
             fontweight="bold", color=vz.INK)
    fig.text(x, y - 0.2 / FIG_H, textwrap.fill(CAPTION, 122), ha="left", va="top",
             fontsize=FS, color=vz.INK, linespacing=1.45)


def style_axes(ax):
    ax.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_color(vz.INK)
    ax.tick_params(labelsize=FS, direction="out", length=2.5, width=0.6, pad=2,
                   colors=vz.INK)


def plot_hwd(ax, style, years, raw, anom):
    """Panel a: heatwave days, recorded and trend removed."""
    ax.plot(years, raw, color=LIGHT, lw=1.0, zorder=2, label="observed")
    if style == "bar":
        ax.bar(years, anom, width=0.7, color=HOT, lw=0, zorder=3, label="trend removed")
    else:
        ax.plot(years, anom, color=HOT, lw=1.3, zorder=3, label="trend removed")
    ax.axhline(0, color=vz.INK_SOFT, lw=0.5, zorder=1)

    for y in MARKED:
        k = int(np.flatnonzero(years == y)[0])
        if style == "line":
            ax.plot(y, anom[k], "o", ms=3.5, color=HOT, mec="white", mew=0.6, zorder=4)
        ax.annotate(str(y), (y, anom[k]), xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=FS, fontweight="bold", color=vz.INK,
                    zorder=5)
        ax.annotate(str(y), (y, raw[k]), xytext=(0, 2.5), textcoords="offset points",
                    ha="center", va="bottom", fontsize=FS - 1, color=vz.INK_MUTED, zorder=5)

    ax.set_ylim(min(anom.min(), 0) * 1.2, raw.max() * 1.12)
    ax.set_ylabel("Heatwave days (JJA)", fontsize=FS, color=vz.INK)
    style_axes(ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.53, 1.0), frameon=False, fontsize=FS,
              handlelength=1.6, borderaxespad=0.2)


def plot_corr(ax, years, r):
    """Panel b: pattern correlation with JJA 2026, through 2026 itself (r = 1)."""
    ax.plot(years, r, color=CORR, lw=1.1, zorder=3)
    ax.axhline(0, color=vz.INK_SOFT, lw=0.5, zorder=1)

    for y in MARKED:
        k = int(np.flatnonzero(years == y)[0])
        ax.plot(y, r[k], "o", ms=3.5, color=CORR, mec="white", mew=0.6, zorder=4)
        label = str(y) if y == p2.REF_YEAR else f"{y}  r = {r[k]:+.2f}"
        ax.annotate(label, (y, r[k]), xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=FS, fontweight="bold", color=vz.INK,
                    zorder=5)

    ax.set_ylim(*CORR_YLIM)
    ax.set_yticks([-0.5, 0, 0.5, 1])
    ax.set_ylabel(f"Pattern correlation\nwith JJA {p2.REF_YEAR} SST", fontsize=FS, color=vz.INK)
    style_axes(ax)


def plot(style, years, raw, anom, field, ryears, r, sst):
    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, sharex=True, figsize=(FIG_W, FIG_H),
        gridspec_kw=dict(height_ratios=PANEL_H, hspace=GAP_H / np.mean(PANEL_H)))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.985, top=1 - TOP_H / FIG_H,
                        bottom=(CAPTION_H + XTICK_H) / FIG_H)
    for ax in (ax_a, ax_b):
        ax.set_facecolor("white")

    plot_hwd(ax_a, style, years, raw, anom)
    plot_corr(ax_b, ryears, r)
    ax_b.set_xlim(years[0] - 1, years[-1] + 2)
    ax_a.tick_params(axis="x", length=0)
    fig.align_ylabels([ax_a, ax_b])

    for ax, letter in ((ax_a, "a"), (ax_b, "b")):
        ax.text(-0.075, 1.0, letter, transform=ax.transAxes, ha="left", va="top",
                fontsize=FS + 1.5, fontweight="bold", color=vz.INK)

    inset_map(fig, ax_a, field)
    inset_sst(fig, ax_b, sst)
    caption(fig)

    out = p2.FIG_DIR / f"fig17_era5_heatwave_days_timeseries_{style}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    vz.use_style()
    years, raw, anom = load()
    field = load_map()
    ryears, r = load_corr()
    sst = sst_field("ERA5", p2.REF_YEAR, None)
    for style in ("line", "bar"):
        plot(style, years, raw, anom, field, ryears, r, sst)


if __name__ == "__main__":
    main()
