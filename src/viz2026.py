"""
Shared figure styling for the 2026 pattern study.

Colour assignments follow one rule: **hue identifies the model family, line
style identifies the experiment** (solid = forced historical/scenario,
dashed = unforced control). That keeps the number of hues on screen down to
three plus the observations, so every pair stays distinguishable for
colour-vision-deficient readers rather than relying on an eight-way cycle.

The three model hues are the first three slots of the validated categorical
palette, which pass the all-pairs CVD and normal-vision separation checks in
both light and dark modes. Observations are drawn in near-black ink, not a
hue: ERA5 is the reference every other line is measured against, so it should
read as the axis of the figure rather than as "series 4".
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- categorical: model family -------------------------------------------
INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#8a8985"
SURFACE = "#fcfcfb"

COLORS = {
    "ERA5": INK,
    "MPI-GE": "#2a78d6",           # slot 1, blue
    "ICON-EPOC-hist": "#eb6834",   # slot 2, orange
    "ICON-EPOC-ctrl": "#eb6834",
    "EERIE": "#1baf7a",            # slot 3, aqua
    "EERIE-ctrl": "#1baf7a",
}
LINESTYLES = {
    "ERA5": "-",
    "MPI-GE": "-",
    "ICON-EPOC-hist": "-",
    "ICON-EPOC-ctrl": (0, (4, 2)),
    "EERIE": "-",
    "EERIE-ctrl": (0, (4, 2)),
}

# --- diverging: SST anomaly ----------------------------------------------
# cmocean "balance" is the oceanographic standard: two hues, perceptually
# uniform, neutral (not a hue) at the zero midpoint.
def anomaly_cmap():
    try:
        import cmocean
        return cmocean.cm.balance
    except ImportError:
        return plt.get_cmap("RdBu_r")


LAND = "#e8e6e1"
COAST = "#6f6e6a"


def use_style():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "axes.edgecolor": INK_MUTED,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": INK_SOFT,
        "text.color": INK,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "grid.color": "#dedcd6",
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 1.6,
        "figure.dpi": 130,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })


def season_label(row):
    """'JJA 2016 (r12)' style label from a row of best_analogues.csv."""
    import pandas as pd
    year = int(row["best_year"])
    mem = row["best_member"]
    if pd.notna(mem):
        return f"JJA {year}  ·  member r{int(mem)}"
    return f"JJA {year}"


# --- maps -----------------------------------------------------------------
# The analysis box (30-60N, 80W-40E) is 120 deg wide by 30 deg tall. Drawn in
# PlateCarree that is a 4:1 letterbox, and a grid of such panels collapses to
# unreadable strips. Mercator is the standard choice for a mid-latitude SST map
# and brings the panel aspect down to ~2.73:1, which is what the figure
# geometry below is sized against.
MAP_EXTENT = [-80, 40, 30, 60]


def map_projection():
    import cartopy.crs as ccrs
    return ccrs.Mercator(central_longitude=-20, min_latitude=25, max_latitude=65)


def panel_aspect(extent=None):
    """Width / height of one map panel in the Mercator projection above."""
    import numpy as np
    lon0, lon1, lat0, lat1 = extent or MAP_EXTENT

    def merc_y(lat):
        return np.log(np.tan(np.pi / 4 + np.deg2rad(lat) / 2))

    return (lon1 - lon0) * np.pi / 180 / (merc_y(lat1) - merc_y(lat0))


def draw_map(ax, field, levels, cmap, extent=None, *, labels_bottom=False,
             labels_left=False, lon_step=20, lat_step=10):
    """Filled contours of one SST field, with land, coastline and graticule."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    pc = ccrs.PlateCarree()
    im = ax.contourf(field.lon, field.lat, field.values, levels=levels,
                     cmap=cmap, extend="both", transform=pc)
    ax.add_feature(cfeature.LAND, facecolor=LAND, zorder=2)
    ax.coastlines(resolution="50m", linewidth=0.4, color=COAST, zorder=3)
    ax.set_extent(extent or MAP_EXTENT, crs=pc)

    lon0, lon1, lat0, lat1 = extent or MAP_EXTENT
    gl = ax.gridlines(draw_labels=labels_bottom or labels_left, linewidth=0.4,
                      color="#ffffff", alpha=0.55, zorder=4,
                      xlocs=range(int(lon0), int(lon1) + 1, lon_step),
                      ylocs=range(int(lat0), int(lat1) + 1, lat_step))
    gl.top_labels = gl.right_labels = False
    gl.bottom_labels = labels_bottom
    gl.left_labels = labels_left
    gl.xlabel_style = gl.ylabel_style = {"size": 7, "color": INK_SOFT}
    return im


def panel_tag(ax, text, loc="lower left"):
    """Small annotation box inside a map panel."""
    x, ha = (0.022, "left") if "left" in loc else (0.978, "right")
    y, va = (0.05, "bottom") if "lower" in loc else (0.95, "top")
    ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha, fontsize=7.6,
            color=INK, zorder=6, linespacing=1.35,
            bbox=dict(facecolor="#ffffff", alpha=0.84, edgecolor="none",
                      boxstyle="round,pad=0.28"))
