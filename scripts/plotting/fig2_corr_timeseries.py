"""
Figure 2 -- pattern correlation against observed JJA 2026, season by season,
through every simulation (and through the observed record itself).

Small multiples rather than six lines on one axis: the records cover very
different periods (MPI-GE 1850-2100, EERIE 1950-2050, km-scale ICON 1990-2025,
ERA5 1940-2026), so overlaying them would make the x-axis unreadable and force
a six-way colour cycle. One row per simulation on a shared time axis keeps each
record legible and still lets the eye compare heights down the column.

The two ensembles are drawn differently, by member count. MPI-GE's 50 members
support a 5-95% band with the member-maximum on top -- the ensemble's best
attempt in each year, which is what the best-analogue panel of figure 1 picks
from. MPI-ESM1.2-ER has 3, where a quantile band would be estimating spread
that is not there, so its members are drawn individually: the one holding the
best analogue in the dataset's colour, the other two as neutral context.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

# the same four datasets as figure 1's columns; the two control runs belong in
# figure 3, where comparing whole distributions is the point
ROW_ORDER = p2.MAIN_KEYS


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    corr = {v: xr.open_dataset(p2.RESULT_DIR / f"corr_{v}.nc") for v in p2.VARIANTS}
    keys = [k for k in ROW_ORDER if k in corr["spatial"].data_vars]

    xmin = min(int(corr["spatial"][k].year.min()) for k in keys)
    xmax = max(int(corr["spatial"][k].year.max()) for k in keys)

    # wide, so that the 35-year km-scale record is not a sliver on an axis that
    # has to stretch from 1850 to 2100 for MPI-GE
    variants = p2.PLOT_VARIANTS
    fig_h = 1.5 * len(keys) + 1.6
    fig, axes = plt.subplots(len(keys), len(variants),
                             figsize=(6.9 * len(variants) + 1.6, fig_h),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(hspace=0.24, wspace=0.085,
                        top=1 - 1.15 / fig_h, bottom=0.62 / fig_h,
                        left=0.185 if len(variants) == 1 else 0.105, right=0.988)

    for j, variant in enumerate(variants):
        for i, key in enumerate(keys):
            ax = axes[i][j]
            r = corr[variant][key]
            color = vz.COLORS[key]
            ls = vz.LINESTYLES[key]

            ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
            ax.grid(axis="y", zorder=0)

            if key == "ERA5":
                # ERA5 vs itself in 2026 is 1 by construction and would run off
                # the top of a shared axis; the dotted line marks the year
                r = r.where(r["year"] != p2.REF_YEAR)

            row = best[(best.variant == variant) & (best.dataset == key)].iloc[0]
            mdim = p2.member_dim(r)
            members = vz.valid_members(r, mdim) if mdim else []

            if mdim and len(members) >= vz.BAND_MIN_MEMBERS:
                lo, hi, band_label = vz.member_band(r, mdim)
                ax.fill_between(r.year, lo, hi, color=color, alpha=0.22, lw=0, zorder=2,
                                label=band_label)
                ax.plot(r.year, r.max(mdim), color=color, lw=1.1, ls=ls, zorder=3,
                        label="member maximum")
            elif mdim:
                # Too few members for a band to mean anything: draw them. The
                # one holding the best analogue is in the dataset's colour so
                # the marked season sits on a line the reader can follow; the
                # rest are neutral context.
                bm = int(row.best_member) if pd.notna(row.best_member) else None
                others = [m for m in members if m != bm]
                for n, m in enumerate(others):
                    ax.plot(r.year, r.sel({mdim: m}), color=vz.CONTEXT, lw=0.8,
                            zorder=2, label="other members" if n == 0 else None)
                if bm is not None:
                    ax.plot(r.year, r.sel({mdim: bm}), color=color, lw=1.2, ls=ls,
                            zorder=3, label=f"member r{bm} (best analogue)")
            else:
                ax.plot(r.year, r, color=color, lw=1.2, ls=ls, zorder=3)

            # mark the best-matching season shown in figure 1
            ax.plot([row.best_year], [row.best_r], marker="o", ms=5.5,
                    mfc=color, mec=vz.SURFACE, mew=1.2, zorder=5, clip_on=False)
            # label on whichever side of the marker has room
            right_side = row.best_year < xmin + 0.72 * (xmax - xmin)
            ax.annotate(f"{int(row.best_year)}  r={row.best_r:+.2f}",
                        (row.best_year, row.best_r), textcoords="offset points",
                        xytext=(8 if right_side else -8, 1), fontsize=7.4,
                        color=vz.INK, zorder=6, va="center",
                        ha="left" if right_side else "right")

            if key == "ERA5":
                ax.axvline(p2.REF_YEAR, color=vz.INK, lw=0.8, ls=":", zorder=4)
                ax.text(p2.REF_YEAR, -0.78, " 2026", fontsize=7, color=vz.INK_SOFT,
                        ha="left", va="bottom", zorder=6)

            # set_ylabel keeps the caption clear of the tick labels; the
            # correlation axis needs no separate "r" label, the subtitle says it
            if j == 0:
                d = p2.DATASETS[key]
                ax.set_ylabel(d.label.replace(" (", "\n("), rotation=0,
                              ha="right", va="center", labelpad=12,
                              fontsize=8.5, color=vz.INK, linespacing=1.4)
            # with a single variant its name belongs in the subtitle, not
            # above every column
            if i == 0 and len(variants) > 1:
                ax.set_title(p2.VARIANTS[variant]["title"], fontsize=9.5,
                             color=vz.INK, pad=7)
            if i == len(keys) - 1:
                ax.set_xlabel("year")

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-0.85, 0.85)
            ax.set_yticks([-0.5, 0, 0.5])

    # every ensemble row explains its own encoding: a band and a maximum for
    # MPI-GE, individual members for the 3-member MPI-ESM1.2-ER
    for i in range(len(keys)):
        handles, labels = axes[i][0].get_legend_handles_labels()
        if handles:
            axes[i][-1].legend(handles, labels, loc="lower right", fontsize=7.6,
                               ncol=2, borderaxespad=0.3)

    fig.suptitle("How close does each simulated summer come to the observed JJA 2026 SST pattern?",
                 fontsize=12, color=vz.INK, y=1 - 0.28 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             f"area-weighted pattern correlation against ERA5 JJA 2026 over {p2.REGION_LABEL}  ·  "
             f"{p2.VARIANTS[variants[0]]['long']}  ·  "
             "1° common ocean grid  ·  dotted line marks 2026 in the observations",
             ha="center", fontsize=8.2, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig2_pattern_correlation_timeseries.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
