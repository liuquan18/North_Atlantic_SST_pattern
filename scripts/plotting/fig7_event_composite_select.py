"""
Figure 7 (select) -- the event composites of figure 7 for the observations and
the forced runs only.

Same method, event selection, smoothing and null band as fig7_event_composite.py
(everything is imported from there); this version drops the two control runs
and the autocorrelation panel and lays the rest out as observations and the
MPI-ESM pair on top, the two km-scale ICON runs below. The empty sixth slot
holds the legend.
"""
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
import src.viz2026 as vz
import fig7_event_composite as f7

#: None marks the slot used for the legend
LAYOUT = ["ERA5", "MPI-GE", "MPI-ER",
          "ICON-EPOC-hist", "EERIE", None]
NCOL = 3


def main():
    warnings.filterwarnings("ignore", "Mean of empty slice", RuntimeWarning)
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{f7.VARIANT}.nc")
    keys = [k for k in LAYOUT if k is not None and k in corr.data_vars]

    nrow = -(-len(LAYOUT) // NCOL)
    fig_h = 3.05 * nrow + 1.25
    fig, axes = plt.subplots(nrow, NCOL, figsize=(12.0, fig_h), squeeze=False)
    fig.subplots_adjust(hspace=0.42, wspace=0.17, top=1 - 1.1 / fig_h,
                        bottom=0.55 / fig_h, left=0.07, right=0.99)
    flat = axes.ravel()

    for key in keys:
        ax = flat[LAYOUT.index(key)]
        c = f7.composite(key, corr[key])
        f7.draw_composite(ax, key, c)
        if LAYOUT.index(key) % NCOL == 0:
            ax.set_ylabel(f"r, forced part removed,\n{f7.SMOOTH}-yr running mean")
        print(f"{key:15s} n={c['n']:3d} ({c['rule']})")

    # the legend gets the empty slot rather than covering a panel's data;
    # its entries are generic, so the first panel's handles serve for all
    handles, labels = flat[LAYOUT.index(keys[0])].get_legend_handles_labels()
    for i, name in enumerate(LAYOUT):
        if name is None or name not in keys:
            flat[i].axis("off")
    slot = flat[LAYOUT.index(None)]
    slot.legend(handles, labels, loc="center left", fontsize=8.5,
                handlelength=2.4, labelspacing=0.9)

    fig.suptitle("Do seasons that resemble JJA 2026 come in decadal spells?",
                 fontsize=12, color=vz.INK, y=1 - 0.18 / fig_h)
    fig.text(0.5, 1 - 0.5 / fig_h,
             f"pattern correlation with ERA5 JJA 2026 (box mean removed), {f7.event_text()}  ·  "
             f"{f7.SMOOTH}-yr running mean before compositing\n"
             "forced part removed: 50-member mean (MPI-GE), linear/quadratic trend (other runs)  ·  "
             f"dots: outside the random-event band, |lag| > {f7.SMOOTH // 2}",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.5)

    out = p2.FIG_DIR / "fig7_event_composite_select.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
