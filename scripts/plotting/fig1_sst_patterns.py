"""
Figure 1 (publication) -- the JJA 2026 North Atlantic + Mediterranean SST
pattern in ERA5 and, from three models, the season that matches it best.

  a  ERA5, JJA 2026 (the reference)
  b  MPI-ESM1.2-ER, its best analogue
  c  km-scale ICON (EPOC, hist+ssp585), its best analogue
  d  MPI-ESM1-2-LR grand ensemble, its best analogue

2x2, panel letters only and no titles; which season each panel shows goes in
the caption, which the script prints. Panels b-d carry their pattern
correlation r with panel a. Fields are the same as in the bottom row of
fig1_patterns_simulations.py (box mean removed, 1991-2020 climatology of
each dataset) and share its colour scale.

With an argument `atl` or `med`, the analogues are those scored over that
half of the box alone, the half is outlined, and the file gets a _atl / _med
suffix.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, VMAX_CAP, season, sst_box, sst_field

#: (panel letter, dataset, which season); "ref" = JJA 2026, "best" = best analogue
PANELS = [("a", "ERA5", "ref"), ("b", "MPI-ER", "best"),
          ("c", "ICON-EPOC-hist", "best"), ("d", "MPI-GE", "best")]
PAPER = "#ffffff"


def main(half=None):
    best = pd.read_csv(p2.best_file(half))
    panels = []
    for letter, key, which in PANELS:
        year, mem, r = season(key, which, best)
        panels.append(dict(letter=letter, key=key, year=year, member=mem, r=r,
                           sst=sst_field(key, year, mem)))

    vz.use_style()
    plt.rcParams.update({"figure.facecolor": PAPER, "savefig.facecolor": PAPER,
                         "axes.facecolor": PAPER})
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)

    vmax = min(VMAX_CAP, max(float(np.nanpercentile(np.abs(p["sst"].values), 98)) for p in panels))
    vmax = np.round(vmax * 2) / 2
    levels = np.linspace(-vmax, vmax, 21)

    # double-column width (180 mm); panel height follows the Mercator aspect
    fig_w = 7.1
    left, right, wspace = 0.62, 0.12, 0.14
    panel_w = (fig_w - left - right - wspace) / 2
    panel_h = panel_w / vz.panel_aspect()
    top, hspace, cb_gap, cb_h, bottom = 0.18, 0.30, 0.36, 0.11, 0.42
    fig_h = top + 2 * panel_h + hspace + cb_gap + cb_h + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))

    box = sst_box(half)
    for i, p in enumerate(panels):
        row, col = divmod(i, 2)
        x0 = left + col * (panel_w + wspace)
        y0 = fig_h - top - (row + 1) * panel_h - row * hspace
        ax = fig.add_axes([x0 / fig_w, y0 / fig_h, panel_w / fig_w, panel_h / fig_h],
                          projection=vz.map_projection())
        im = vz.draw_map(ax, p["sst"], levels, vz.anomaly_cmap(),
                         labels_bottom=(row == 1), labels_left=(col == 0))
        if box is not None:
            vz.draw_box(ax, box, lw=1.0)
        ax.text(-0.012, 1.025, p["letter"], transform=ax.transAxes, ha="right", va="bottom",
                fontsize=10, fontweight="bold", color=vz.INK)
        if p["r"] is not None:
            vz.panel_tag(ax, f"r = {p['r']:.2f}", loc="upper right")

    cb_w = 0.5 * fig_w
    cax = fig.add_axes([(fig_w - cb_w) / 2 / fig_w, bottom / fig_h, cb_w / fig_w, cb_h / fig_h])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both",
                      ticks=np.linspace(-vmax, vmax, 9))
    cb.set_label("SST anomaly (°C)", fontsize=8, color=vz.INK)
    cb.ax.tick_params(labelsize=7, length=2.5)
    cb.outline.set_linewidth(0.5)

    out = p2.FIG_DIR / f"fig1_sst_patterns_2026{p2.half_suffix(half)}.png"
    fig.savefig(out, dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}  (SST scale ±{vmax} °C)")

    print("\ncaption panels:")
    for p in panels:
        d = p2.DATASETS[p["key"]]
        mem = f", member r{p['member']}" if p["member"] else ""
        r = "" if p["r"] is None else f", r = {p['r']:.2f}"
        print(f"  ({p['letter']}) {d.label}, JJA {p['year']}{mem}{r}")


if __name__ == "__main__":
    main(p2.half_from_argv(sys.argv))
