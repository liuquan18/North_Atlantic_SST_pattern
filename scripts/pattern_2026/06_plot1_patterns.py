"""
Figure 1 -- the observed JJA 2026 North Atlantic + Mediterranean SST pattern
and the closest season each simulation can produce.

Layout: columns are datasets (ERA5, MPI-GE, km-scale ICON, EERIE), rows are
the two ways of defining "the pattern" (box mean removed / global-ocean mean
removed). Observations show JJA 2026 itself; each model shows its best-matching
season, so the question the figure answers is "what is the closest thing this
simulation has to 2026, and how close is it?"

All panels share one symmetric diverging scale, so amplitudes are comparable
across panels and not just shapes. The scale is deliberately capped below the
Mediterranean extremes: letting +3 degC Mediterranean cells set the range
washes the Atlantic part of the pattern out to near-white.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

COLS = p2.MAIN_KEYS                      # ERA5, MPI-GE, ICON-EPOC-hist, EERIE
ROWS = ["spatial", "global"]
VMAX_CAP = 2.0


def panel_field(key, variant, best):
    """
    The 0.25 deg pattern to draw for one dataset under one variant.

    The season is selected *before* the pattern is formed. Both variants are
    per-season operations (subtract this field's box mean / subtract this
    year's global mean), so the result is identical either way -- but on the
    2.9 GB 50-member MPI-GE array it is the difference between touching one
    120x480 slice and touching all 12,550 of them.
    """
    anom, gm = p2.load_anomaly(key), p2.load_gmsst(key)

    if key == "ERA5":
        sel, label, r = {"year": p2.REF_YEAR}, f"JJA {p2.REF_YEAR}  ·  observed", None
    else:
        row = best[(best.variant == variant) & (best.dataset == key)].iloc[0]
        sel = {"year": int(row.best_year)}
        if pd.notna(row.best_member):
            sel["member"] = int(row.best_member)
        label, r = vz.season_label(row), float(row.best_r)

    field = p2.make_pattern(anom.sel(sel), gm.sel(sel), variant).load()
    return field, label, r


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    keys = [k for k in COLS if p2.DATASETS[k].anom_file.exists()]
    if len(keys) < len(COLS):
        print(f"!! missing {[k for k in COLS if k not in keys]}; plotting {keys}")

    fields = {(k, v): panel_field(k, v, best) for k in keys for v in ROWS}

    vmax = min(VMAX_CAP,
               max(float(np.nanpercentile(np.abs(f[0].values), 98)) for f in fields.values()))
    vmax = np.round(vmax * 2) / 2
    levels = np.linspace(-vmax, vmax, 21)
    cmap = vz.anomaly_cmap()
    proj = vz.map_projection()

    # size the grid so the fixed map aspect leaves no dead space
    panel_w = 3.55
    panel_h = panel_w / vz.panel_aspect()
    fig_w = panel_w * len(keys) + 1.6
    fig_h = panel_h * len(ROWS) + 2.05

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(len(ROWS), len(keys) + 1, figure=fig,
                  width_ratios=[1] * len(keys) + [0.035],
                  height_ratios=[1] * len(ROWS),
                  hspace=0.16, wspace=0.06,
                  left=0.082, right=0.955,
                  top=1 - 1.35 / fig_h, bottom=0.42 / fig_h)

    for i, variant in enumerate(ROWS):
        for j, key in enumerate(keys):
            ax = fig.add_subplot(gs[i, j], projection=proj)
            field, season, r = fields[(key, variant)]
            im = vz.draw_map(ax, field, levels, cmap,
                             labels_bottom=(i == len(ROWS) - 1), labels_left=(j == 0))

            if i == 0:
                d = p2.DATASETS[key]
                ax.set_title(f"{d.label}\n{d.resolution}", pad=5, fontsize=9, color=vz.INK)

            vz.panel_tag(ax, season if r is None else f"{season}\nr = {r:+.2f}")

            if j == 0:
                ax.text(-0.105, 0.5, p2.VARIANTS[variant]["title"],
                        transform=ax.transAxes, rotation=90, va="center", ha="center",
                        fontsize=9, color=vz.INK)

    cax = fig.add_subplot(gs[:, -1])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.set_label("JJA SST anomaly (°C, 1991–2020 base)", fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    fig.suptitle("The JJA 2026 North Atlantic + Mediterranean SST pattern, "
                 "and the closest season each model produces",
                 fontsize=12, color=vz.INK, y=1 - 0.30 / fig_h)
    fig.text(0.5, 1 - 0.66 / fig_h,
             "observations show JJA 2026; each model shows its best-matching season out of every "
             "year and member simulated  ·  30–60°N, 80°W–40°E  ·  Mercator",
             ha="center", fontsize=8.3, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig1_patterns_2026.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}  (colour scale ±{vmax} °C)")


if __name__ == "__main__":
    main()
