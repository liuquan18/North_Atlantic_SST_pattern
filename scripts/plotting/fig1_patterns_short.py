"""
Figure 1 (short) -- ERA5 and km-scale ICON (EPOC) only, as a 2x2 grid.

Top row: ERA5's closest earlier summer (JJA 2003) and JJA 2026 itself.
Bottom row: EPOC's best-matching season, and an empty slot where EPOC's own
JJA 2026 will go -- the transient run has not reached summer 2026 yet.

One symmetric diverging scale across all panels, capped below the
Mediterranean extremes as in figure 1.
"""
import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

VARIANT = p2.PLOT_VARIANTS[0]
VMAX_CAP = 2.0
EPOC = "ICON-EPOC-hist"


def season_field(key, year):
    anom, gm = p2.load_anomaly(key), p2.load_gmsst(key)
    sel = {"year": year}
    return p2.make_pattern(anom.sel(sel), gm.sel(sel), VARIANT).load()


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")
    row = lambda k: best[(best.variant == VARIANT) & (best.dataset == k)].iloc[0]

    era5_best, epoc_best = row("ERA5"), row(EPOC)
    # (row, col) -> (dataset, year, tag)
    panels = {
        (0, 0): ("ERA5", int(era5_best.best_year),
                 f"JJA {int(era5_best.best_year)}\nr = {era5_best.best_r:+.2f}"),
        (0, 1): ("ERA5", p2.REF_YEAR, f"JJA {p2.REF_YEAR}  ·  observed"),
        (1, 0): (EPOC, int(epoc_best.best_year),
                 f"JJA {int(epoc_best.best_year)}\nr = {epoc_best.best_r:+.2f}"),
    }
    fields = {pos: season_field(k, y) for pos, (k, y, _) in panels.items()}

    vmax = min(VMAX_CAP,
               max(float(np.nanpercentile(np.abs(f.values), 98)) for f in fields.values()))
    vmax = np.round(vmax * 2) / 2
    levels = np.linspace(-vmax, vmax, 21)
    cmap = vz.anomaly_cmap()
    proj = vz.map_projection()

    panel_w = 3.9
    panel_h = panel_w / vz.panel_aspect()
    fig_w = panel_w * 2 + 2.6
    fig_h = panel_h * 2 + 1.9

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 0.035],
                  hspace=0.22, wspace=0.06,
                  left=0.16, right=0.93, top=1 - 1.3 / fig_h, bottom=0.42 / fig_h)

    col_titles = ["best-matching earlier season", f"JJA {p2.REF_YEAR}"]
    for i, key in enumerate(["ERA5", EPOC]):
        for j in range(2):
            ax = fig.add_subplot(gs[i, j], projection=proj)
            if i == 0:
                ax.set_title(col_titles[j], pad=5, fontsize=9, color=vz.INK)
            if j == 0:
                d = p2.DATASETS[key]
                ax.text(-0.2, 0.5, f"{d.label.replace(' (', chr(10) + '(')}\n{d.resolution}",
                        transform=ax.transAxes, rotation=90, va="center", ha="center",
                        fontsize=8.5, color=vz.INK, linespacing=1.35)

            if (i, j) not in panels:
                # EPOC has not reached JJA 2026 yet: keep the slot, draw nothing
                ax.set_extent(vz.MAP_EXTENT, crs=ccrs.PlateCarree())
                ax.spines["geo"].set_edgecolor(vz.INK_MUTED)
                ax.spines["geo"].set_linestyle((0, (3, 3)))
                ax.text(0.5, 0.5, f"JJA {p2.REF_YEAR}\nnot yet simulated",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=8.5, color=vz.INK_MUTED)
                continue

            im = vz.draw_map(ax, fields[(i, j)], levels, cmap,
                             labels_bottom=(i == 1), labels_left=(j == 0))
            vz.panel_tag(ax, panels[(i, j)][2])

    cax = fig.add_subplot(gs[:, -1])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.set_label("JJA SST anomaly (°C)", fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    fig.suptitle("The JJA 2026 SST pattern in ERA5 and km-scale ICON (EPOC)",
                 fontsize=12, color=vz.INK, y=1 - 0.30 / fig_h)
    fig.text(0.5, 1 - 0.66 / fig_h,
             f"{p2.VARIANTS[VARIANT]['long']}  ·  30–60°N, 80°W–40°E  ·  "
             "anomalies against 1991–2020  ·  r = pattern correlation with ERA5 JJA 2026",
             ha="center", fontsize=8.0, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig1_patterns_2026_short.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}  (colour scale ±{vmax} °C)")


if __name__ == "__main__":
    main()
