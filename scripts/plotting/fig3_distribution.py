"""
Figure 3 -- not just "can a model make this pattern once", but how often.

Figure 1 shows each simulation's single best season, which by construction is
its most flattering one. This figure shows the whole distribution behind that
number, so a model that lands near 2026 routinely is visibly different from one
that got there once in 250 years.

The benchmark line is taken from the observations themselves: the highest
correlation any *other* observed summer since 1940 reaches against JJA 2026.
Reading "fraction of simulated seasons at least as 2026-like as the closest
observed analogue" then needs no arbitrary threshold, and the benchmark moves
with the data rather than being pinned to a number in a docstring.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

ROW_ORDER = ["MPI-GE", "MPI-ER", "ICON-EPOC-hist", "ICON-EPOC-ctrl",
             "EERIE", "EERIE-ctrl"]


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    corr = {v: xr.open_dataset(p2.RESULT_DIR / f"corr_{v}.nc") for v in p2.VARIANTS}
    keys = [k for k in ROW_ORDER if k in corr["spatial"].data_vars]

    variants = p2.PLOT_VARIANTS
    # margins in absolute inches, so the chrome does not grow with the number
    # of datasets the way a fractional top/bottom would
    fig_h = 0.92 * len(keys) + 1.65
    fig, axes = plt.subplots(1, len(variants),
                             figsize=(5.6 * len(variants) + 3.4, fig_h),
                             sharey=True, squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(wspace=0.06, top=1 - 1.05 / fig_h, bottom=0.55 / fig_h,
                        left=0.205, right=0.985)

    for j, variant in enumerate(variants):
        ax = axes[j]
        obs_ref = float(best[(best.variant == variant) & (best.dataset == "ERA5")].best_r.iloc[0])
        obs_year = int(best[(best.variant == variant) & (best.dataset == "ERA5")].best_year.iloc[0])

        ax.axvline(obs_ref, color=vz.INK, lw=1.0, ls=":", zorder=4)
        ax.axvline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
        ax.grid(axis="x", zorder=0)

        for i, key in enumerate(keys):
            y = len(keys) - 1 - i
            vals = corr[variant][key].values.ravel()
            vals = vals[np.isfinite(vals)]
            color = vz.COLORS[key]

            # `vert` is deprecated from matplotlib 3.10 in favour of `orientation`
            parts = ax.violinplot([vals], positions=[y], orientation="horizontal",
                                  widths=0.72, showextrema=False, showmedians=False)
            for b in parts["bodies"]:
                b.set_facecolor(color)
                b.set_alpha(0.30)
                b.set_edgecolor(color)
                b.set_linewidth(1.0)

            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.plot([q1, q3], [y, y], color=color, lw=3.0, solid_capstyle="butt", zorder=5)
            ax.plot([med], [y], marker="|", ms=9, color=vz.SURFACE, mew=1.8, zorder=6)
            ax.plot([vals.max()], [y], marker="o", ms=5, mfc=color,
                    mec=vz.SURFACE, mew=1.1, zorder=6)

            # As a count, not a rounded percentage: MPI-GE clears the benchmark
            # in a handful of its 12,550 seasons, which "0.0%" would hide.
            n_reach = int((vals >= obs_ref).sum())
            pct = 100.0 * n_reach / vals.size
            reach = (f"{n_reach} of {vals.size:,} seasons reach the observed analogue"
                     + (f"  ({pct:.2g}%)" if n_reach else ""))
            ax.text(0.012, y + 0.30, reach, transform=ax.get_yaxis_transform(),
                    ha="left", va="bottom", fontsize=7.4, color=vz.INK_SOFT)

        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels([p2.DATASETS[k].label.replace(" (", "\n(")
                            for k in reversed(keys)], fontsize=8.3)
        if len(variants) > 1:
            ax.set_title(p2.VARIANTS[variant]["title"], fontsize=9.5, color=vz.INK, pad=7)
        ax.set_xlabel("pattern correlation with observed JJA 2026")
        ax.set_xlim(-0.9, 0.9)
        ax.set_ylim(-0.6, len(keys) - 0.25)
        # left of the benchmark line, in the empty strip above the top violin
        ax.text(obs_ref, len(keys) - 0.30, f"closest observed\nanalogue: JJA {obs_year} ",
                fontsize=7.4, color=vz.INK, va="top", ha="right", linespacing=1.3)

    fig.suptitle("How often does each simulation produce a 2026-like summer SST pattern?",
                 fontsize=12, color=vz.INK, y=1 - 0.26 / fig_h)
    fig.text(0.5, 1 - 0.60 / fig_h,
             "distribution over every simulated JJA season and ensemble member  ·  "
             f"{p2.VARIANTS[variants[0]]['long']}  ·  "
             "bar = interquartile range, tick = median, dot = best season",
             ha="center", fontsize=8.2, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig3_correlation_distribution.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
