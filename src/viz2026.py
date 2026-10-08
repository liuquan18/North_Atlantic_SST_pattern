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

#: Neutral grey for context marks that carry no identity of their own -- e.g.
#: the non-highlighted members of a small ensemble, which are there to show
#: spread, not to be told apart from one another.
CONTEXT = "#a8a7a1"

#: Below this many members, a 5-95% band is a quantile estimate of something
#: that isn't there; draw the members instead.
BAND_MIN_MEMBERS = 10

COLORS = {
    "ERA5": INK,
    "MPI-GE": "#2a78d6",           # slot 1, blue
    "ICON-EPOC-hist": "#eb6834",   # slot 2, orange
    "ICON-EPOC-ctrl": "#eb6834",
    "EERIE": "#1baf7a",            # slot 3, aqua
    "EERIE-ctrl": "#1baf7a",
    # Violet, not the palette's nominal 4th slot (yellow). With four model
    # hues on one axis the documented order puts yellow beside orange, which
    # fails the all-pairs normal-vision floor (dE 13.7 < 15); the palette notes
    # that this trade should be undone when a 4th series appears, and undoing
    # it is a pure re-order, not a re-step. Violet clears every all-pairs gate
    # against the three hues already in use (worst normal-vision dE 16.3,
    # worst CVD dE 9.2) and repaints nothing.
    "MPI-ER": "#4a3aa7",           # violet
}
LINESTYLES = {
    "ERA5": "-",
    "MPI-GE": "-",
    "MPI-ER": "-",
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
# The analysis box (20-60N, 80W-40E) is 120 deg wide by 40 deg tall. Drawn in
# PlateCarree that is a 3:1 letterbox, and a grid of such panels collapses to
# thin strips. Mercator is the standard choice for a mid-latitude SST map
# and brings the panel aspect down to ~2.2:1, which is what the figure
# geometry below is sized against (via panel_aspect).
MAP_EXTENT = [-80, 40, 20, 60]


def map_projection():
    import cartopy.crs as ccrs
    return ccrs.Mercator(central_longitude=-20, min_latitude=15, max_latitude=65)


def panel_aspect(extent=None):
    """Width / height of one map panel in a Mercator projection."""
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


#: Europe view for land fields (heatwaves). The data grid is wider (30-72N,
#: 15W-45E); this crops most of North Africa so the panel reads as Europe,
#: and still holds the whole heatwave averaging box (src.heatwave.MEAN_BOX).
EUROPE_EXTENT = [-12, 42, 35, 70]
OCEAN = "#eef0f2"


def europe_projection():
    import cartopy.crs as ccrs
    return ccrs.Mercator(central_longitude=15, min_latitude=25, max_latitude=75)


def draw_land_map(ax, field, cmap, norm, extent=None, *, labels_bottom=False,
                  labels_left=False, lon_step=10, lat_step=10, land_under=False):
    """
    A land-only field (NaN over sea) drawn cell by cell, with the sea filled
    flat. Unlike draw_map, nothing is painted over land, since land is where
    the data is. Pass field=None for an empty panel (sea and coastline only).

    land_under=True paints land grey *beneath* the field, so land cells that
    are NaN (e.g. "no heatwave" for intensity or onset) read as grey rather
    than blank.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    pc = ccrs.PlateCarree()
    extent = extent or EUROPE_EXTENT
    ax.add_feature(cfeature.OCEAN, facecolor=OCEAN, zorder=0)
    if land_under:
        ax.add_feature(cfeature.LAND, facecolor=LAND, zorder=0.5)
    im = None
    if field is not None:
        im = ax.pcolormesh(field.lon, field.lat, field.values, cmap=cmap, norm=norm,
                           shading="nearest", transform=pc, zorder=1, rasterized=True)
    ax.coastlines(resolution="50m", linewidth=0.4, color=COAST, zorder=3)
    ax.add_feature(cfeature.BORDERS, linewidth=0.25, edgecolor=COAST, alpha=0.6, zorder=3)
    ax.set_extent(extent, crs=pc)

    lon0, lon1, lat0, lat1 = extent
    gl = ax.gridlines(draw_labels=labels_bottom or labels_left, linewidth=0.4,
                      color="#ffffff", alpha=0.55, zorder=4,
                      xlocs=range(-20, 61, lon_step), ylocs=range(30, 81, lat_step))
    gl.top_labels = gl.right_labels = False
    gl.bottom_labels = labels_bottom
    gl.left_labels = labels_left
    gl.xlabel_style = gl.ylabel_style = {"size": 7, "color": INK_SOFT}
    return im


def draw_box(ax, box, **kw):
    """Dashed outline of a {"lon": (w, e), "lat": (s, n)} box, e.g. the averaging region."""
    import cartopy.crs as ccrs
    import numpy as np
    (w, e), (s, n) = box["lon"], box["lat"]
    # densified edges so parallels follow the projection's curvature
    lon = np.r_[np.linspace(w, e, 50), np.full(50, e), np.linspace(e, w, 50), np.full(50, w)]
    lat = np.r_[np.full(50, s), np.linspace(s, n, 50), np.full(50, n), np.linspace(n, s, 50)]
    style = dict(color=INK, lw=0.9, ls=(0, (4, 3)), zorder=5)
    style.update(kw)
    ax.plot(lon, lat, transform=ccrs.PlateCarree(), **style)


def heatwave_cmap():
    """Sequential, light at zero: cmocean 'amp' (perceptually uniform), else YlOrRd."""
    try:
        import cmocean
        return cmocean.cm.amp
    except ImportError:
        return plt.get_cmap("YlOrRd")


def panel_tag(ax, text, loc="lower left"):
    """Small annotation box inside a map panel."""
    x, ha = (0.022, "left") if "left" in loc else (0.978, "right")
    y, va = (0.05, "bottom") if "lower" in loc else (0.95, "top")
    ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha, fontsize=7.6,
            color=INK, zorder=6, linespacing=1.35,
            bbox=dict(facecolor="#ffffff", alpha=0.84, edgecolor="none",
                      boxstyle="round,pad=0.28"))


def valid_members(da, member_dim="member"):
    """
    The member labels that actually carry data.

    Not simply `da[member_dim].values`: a results file written before the
    per-dataset member naming pads a 3-member ensemble out to the 50 of the
    largest one with NaN, so the dimension length overstates the ensemble.
    """
    other = [d for d in da.dims if d != member_dim]
    has = da.notnull().any(other) if other else da.notnull()
    return [m for m, ok in zip(da[member_dim].values, has.values) if bool(ok)]


def member_band(da, member_dim="member"):
    """
    (lower, upper, label) for an ensemble's spread.

    A 5-95% quantile band is only meaningful with enough members to have tails.
    MPI-GE has 50 and gets one; MPI-ESM1.2-ER has 3, where the honest summary is
    the range the three realizations actually span.

    The member count is the number of members carrying data, not the length of
    the dimension: a file written before the per-dataset member naming can pad
    a 3-member ensemble out to 50 with NaN, and reporting "50 members" then
    would be wrong.
    """
    n = len(valid_members(da, member_dim))
    if n >= BAND_MIN_MEMBERS:
        return da.quantile(0.05, member_dim), da.quantile(0.95, member_dim), \
            f"5–95% of {n} members"
    return da.min(member_dim), da.max(member_dim), f"range of {n} members"
