"""
Figure 26 -- the summer NAO and EA patterns of ERA5 and MPI-GE: regression of
the JJA SLP anomaly on each standardised PC (hPa per standard deviation), from
01_nao_ea_eof.py. Rows NAO, EA; columns ERA5, MPI-GE. The sign box of each
mode, the box whose regression fixes the positive phase, is outlined; the
method is written out in the framed box under the maps (METHOD, shared with
fig27_nao_ea_index_mpige.py).
"""
import importlib
import sys
import textwrap
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
eof = importlib.import_module("01_nao_ea_eof")
import src.pattern_2026 as p2
import src.viz2026 as vz

OUT = p2.FIG_DIR / "fig26_nao_ea_patterns.png"
EXTENT = [-90, 40, 20, 80]
LEVELS = np.arange(-2.75, 2.76, 0.25)
FS = 7.5
METHOD_TITLE = "How the NAO and EA patterns and indices are computed"
METHOD = (
    "1. Data: JJA mean sea-level pressure, ERA5 (1940–2026; 2026 from ERA5T) and the 50 members "
    "of MPI-ESM1.2-LR (historical + SSP2-4.5), conservatively remapped to 2.5°, over "
    f"{eof.DOMAIN_LABEL}.  "
    "2. Anomalies: ERA5 — a linear trend over 1940–2026 removed at every grid point; MPI-GE — "
    "the 50-member mean of each year removed from every member (forced signal out, internal "
    "variability kept).  "
    "3. EOFs: SVD of the anomalies weighted by √cos(latitude), 1940–2026, separately for ERA5 "
    "(87 summers) and MPI-GE (50 × 87 member-summers pooled). EOF1 is the summer NAO and EOF2 "
    "the summer EA in both, checked against the patterns.  "
    "4. Index: the PC (projection of the weighted anomalies onto the EOF), divided by its standard "
    "deviation over 1940–2026, so it has mean 0 and standard deviation 1; MPI-GE summers of "
    "1935–1939 and 2027–2030 (for the 10-year windows) are projected onto the same EOF.  "
    "5. Pattern: regression of the unweighted SLP anomaly on the index, hPa per standard deviation.  "
    "6. Sign: positive NAO = low pressure over Greenland/Iceland (60–80°N, 60–10°W); positive EA = "
    "low pressure at its main centre over the British Isles/North Sea (45–65°N, 20°W–20°E; the "
    "CPC convention of low pressure at the EA centre). Each mode is flipped where needed so that "
    "its regression averaged over that box (dashed on the maps of Fig. 26) is negative: a positive index always goes "
    "with the positive pattern drawn here. Checked: composites of the summers with index > +1 "
    "correlate at r ≥ 0.95 with the positive pattern.")


def method_box(fig, y_top, height_in, fig_h, width=118):
    """The METHOD text in a framed box whose top is at figure fraction `y_top`."""
    text = textwrap.fill(METHOD, width)
    fig.text(0.02, y_top, METHOD_TITLE + "\n" + text, ha="left", va="top", fontsize=FS - 0.5,
             color=vz.INK, linespacing=1.45,
             bbox=dict(facecolor="#f4f4f2", edgecolor=vz.INK_MUTED, lw=0.6,
                       boxstyle="round,pad=0.6"))
    return text.count("\n") + 2


def method_lines(width=118):
    return textwrap.fill(METHOD, width).count("\n") + 2


def main():
    vz.use_style()
    ds = xr.open_dataset(eof.OUT)
    proj = ccrs.LambertConformal(central_longitude=-25, central_latitude=50,
                                 standard_parallels=(35, 65))
    n_lines = method_lines()
    box_in = 0.13 * n_lines + 0.35
    map_in = 2.25
    cb_in = 0.8
    fig_h = 0.35 + 2 * map_in + cb_in + box_in
    fig = plt.figure(figsize=(7.09, fig_h))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, left=0.05, right=0.97, top=1 - 0.3 / fig_h,
                          bottom=(box_in + cb_in) / fig_h, hspace=0.18, wspace=0.06)
    cmap = vz.anomaly_cmap()
    letters = iter("abcd")
    for i, mode in enumerate(("NAO", "EA")):
        box, _ = eof.SIGN_BOX[mode]
        for j, key in enumerate(("ERA5", "MPI-GE")):
            ax = fig.add_subplot(gs[i, j], projection=proj)
            p = ds[f"pattern_{key}|{mode}"]
            im = ax.contourf(p.lon, p.lat, p.values, levels=LEVELS, cmap=cmap, extend="both",
                             transform=ccrs.PlateCarree())
            cs = ax.contour(p.lon, p.lat, p.values, levels=[v for v in np.arange(-3, 3.1, 0.5) if v],
                            colors=vz.INK, linewidths=0.5,
                            linestyles=["--" if v < 0 else "-" for v in np.arange(-3, 3.1, 0.5) if v],
                            transform=ccrs.PlateCarree())
            ax.clabel(cs, fmt="%.1f", fontsize=5.5, inline_spacing=1)
            ax.coastlines(resolution="50m", linewidth=0.4, color=vz.COAST)
            ax.add_feature(cfeature.LAND, facecolor="none", edgecolor="none")
            ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
            gl = ax.gridlines(linewidth=0.3, color=vz.INK_MUTED, alpha=0.5,
                              xlocs=range(-90, 41, 30), ylocs=range(20, 81, 20))
            vz.draw_box(ax, box, lw=0.9)
            n = "1940–2026" if key == "ERA5" else "50 × 1940–2026"
            sub = "linearly detrended" if key == "ERA5" else "ensemble mean removed"
            var = ds.attrs[f"var_{key}|{mode}"]
            k = eof.MODES[key][mode] + 1
            ax.set_title(f"({next(letters)}) {key} {mode}: EOF{k}, {var:.0%} of variance\n"
                         f"{n}, {sub}", loc="left", fontsize=FS, color=vz.INK, pad=3,
                         linespacing=1.3)
    cax = fig.add_axes([0.3, (box_in + 0.5) / fig_h, 0.4, 0.1 / fig_h])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=np.arange(-2, 2.1, 1))
    cb.set_label("JJA SLP regression on the standardised index (hPa per std)", fontsize=FS,
                 color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=FS - 0.5)
    method_box(fig, (box_in - 0.12) / fig_h, box_in, fig_h)
    fig.savefig(OUT, dpi=300, facecolor="white")
    fig.savefig(OUT.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
