"""
Figure 17 -- communication figure: ERA5 summer heatwave days over southern
Europe, 1940-2026, one panel. The series as recorded is drawn very light; the
series with the long-term trend removed is the solid, highlighted one. JJA
2003 and 2026 are labelled with their year only. An inset map shows the
JJA 2026 heatwave days per grid cell, with the averaging box outlined.

Two versions of the same figure:
  _line  both series as lines
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
background, no fill over the sea in the map, 300 dpi PNG plus vector PDF.
"""
import sys
from pathlib import Path

import cartopy.crs as ccrs
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
from fig1_patterns_simulations import land_mean_series

MARKED = [2003, p2.REF_YEAR]
TREND_DEGREE = 2

HOT = vz.heatwave_cmap()(0.80)    # the highlighted, trend-removed series
LIGHT = "#c9c7c0"                 # the recorded series, context only
#: inset map: same colours and 5-day steps as the heatwave-day maps of fig 1
MAP_LEVELS = np.arange(0, 46, 5)
#: map extent: the averaging box plus a margin, so its dashed outline is not
#: clipped at the map edge
MAP_EXTENT = [-13, 34, 33.5, 62]
#: inset position in the time-series axes' fraction: top left, above the low
#: 1940-1965 values
INSET = (0.005, 0.50, 0.30, 0.50)

#: publication layout: two-column width (180 mm), 7-8 pt type
FIG_SIZE = (7.09, 3.75)
FS = 7.5
#: fraction of the figure height kept below the plot for the caption
CAPTION_SPACE = 0.25
CAPTION_TITLE = "Summer heatwave days in southern Europe, 1940–2026."
CAPTION = (
    f"June–August heatwave days in ERA5, averaged over land in {hw.MEAN_BOX_LABEL} "
    "(dashed box in the inset). Grey: observed. Red: departure from a quadratic long-term "
    "trend fitted to 1940–2026. A heatwave is a period of at least three consecutive days "
    "with daily maximum temperature above the calendar-day 90th percentile of 1991–2020 "
    f"(Xu et al. 2026). Inset: heatwave days in JJA {p2.REF_YEAR} per grid cell.")


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


def inset_map(fig, ax, field):
    """Frameless map of `field` at INSET in `ax`, colour bar to its right."""
    x0, y0, w, h = ax.get_position().bounds
    ix, iy, iw, ih = INSET
    m = fig.add_axes([x0 + ix * w, y0 + iy * h, iw * w, ih * h], projection=ccrs.PlateCarree())
    cmap = vz.heatwave_cmap()
    norm = BoundaryNorm(MAP_LEVELS, cmap.N, extend="max")
    im = vz.draw_land_map(m, field, cmap, norm, MAP_EXTENT)
    for art in list(m.artists):      # no graticule at this size
        if isinstance(art, Gridliner):
            art.remove()
    # no fill over the sea: the coastline alone outlines the land
    for art in list(m.collections) + list(m.artists):
        if isinstance(art, FeatureArtist) and art._kwargs.get("facecolor") == vz.OCEAN:
            art.remove()
    m.spines["geo"].set_visible(False)
    vz.draw_box(m, hw.MEAN_BOX, lw=0.7)
    # the map keeps its aspect inside the box it was given: hang the colour
    # bar on the extent it actually drew
    fig.canvas.draw()
    mx0, my0, mw, mh = m.get_position(original=False).bounds
    cax = fig.add_axes([mx0 + mw + 0.008, my0 + 0.1 * mh, 0.007, 0.8 * mh])
    cb = fig.colorbar(im, cax=cax, extend="max", ticks=MAP_LEVELS[::2])
    cb.set_label(f"heatwave days\nJJA {p2.REF_YEAR}", fontsize=FS - 1, labelpad=7, color=vz.INK)
    cb.ax.tick_params(labelsize=FS - 1.5, length=1.5, width=0.5, pad=1.5)
    cb.outline.set_linewidth(0)


def caption(fig):
    """Bold title sentence, then the caption body, wrapped to the plot width."""
    import textwrap
    x = 0.012
    fig.text(x, CAPTION_SPACE - 0.03, CAPTION_TITLE, ha="left", va="top", fontsize=FS,
             fontweight="bold", color=vz.INK)
    fig.text(x, CAPTION_SPACE - 0.03 - 0.052, textwrap.fill(CAPTION, 122), ha="left", va="top",
             fontsize=FS, color=vz.INK, linespacing=1.45)


def plot(style, years, raw, anom, field):
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.975, bottom=CAPTION_SPACE + 0.08)

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

    ax.set_xlim(years[0] - 1, years[-1] + 2)
    ax.set_ylim(min(anom.min(), 0) * 1.2, raw.max() * 1.12)
    ax.set_ylabel("Heatwave days (JJA)", fontsize=FS, color=vz.INK)
    ax.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_color(vz.INK)
    ax.tick_params(labelsize=FS, direction="out", length=2.5, width=0.6, pad=2,
                   colors=vz.INK)
    ax.legend(loc="upper center", bbox_to_anchor=(0.53, 1.0), frameon=False, fontsize=FS,
              handlelength=1.6, borderaxespad=0.2)
    inset_map(fig, ax, field)
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
    for style in ("line", "bar"):
        plot(style, years, raw, anom, field)


if __name__ == "__main__":
    main()
