"""
Figure 2 (short) -- pattern correlation against observed JJA 2026, restricted
to ERA5 and the km-scale ICON (EPOC) transient run over the years both cover.

The full figure 2 stretches its time axis from 1850 to 2100 for MPI-GE, which
leaves the 35-year EPOC record a sliver. Here the axis is the intersection of
the two records (1990 onward), so the observed and km-scale summers can be
compared year by year. The axis still runs to 2026 so the dotted reference
line marks where the observed pattern comes from.
"""

import sys

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

ROW_ORDER = ["ERA5", "ICON-EPOC-hist"]


def common_years(corr, keys):
    """First and last year in which every key has a correlation."""
    valid = [corr[k].dropna("year").year for k in keys]
    return max(int(y.min()) for y in valid), min(int(y.max()) for y in valid)


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    corr = {v: xr.open_dataset(p2.RESULT_DIR / f"corr_{v}.nc") for v in p2.VARIANTS}
    keys = [k for k in ROW_ORDER if k in corr["spatial"].data_vars]

    y0, y1 = common_years(corr["spatial"], keys)
    xmin, xmax = y0 - 0.5, max(y1, p2.REF_YEAR) + 0.5

    variants = p2.PLOT_VARIANTS
    fig_h = 3 * len(keys) + 1.6
    fig, axes = plt.subplots(len(keys), len(variants),
                             figsize=(5.2 * len(variants) + 1.6, fig_h),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(hspace=0.24, wspace=0.085,
                        top=1 - 1.3 / fig_h, bottom=0.62 / fig_h,
                        left=0.215 if len(variants) == 1 else 0.12, right=0.985)

    for j, variant in enumerate(variants):
        for i, key in enumerate(keys):
            ax = axes[i][j]
            r = corr[variant][key].sel(year=slice(y0, y1))
            color = vz.COLORS[key]

            ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
            ax.grid(axis="y", zorder=0)

            # ERA5 vs itself in 2026 is 1 by construction; 2026 is outside the
            # common window anyway, the dotted line marks it
            ax.plot(r.year, r, color=color, lw=1.2, ls=vz.LINESTYLES[key],
                    marker="o", ms=2.6, zorder=3, label="JJA season")
            # # centred 3-year running mean; the end years have no full window
            # ax.plot(r.year, r.rolling(year=5, center=True).mean(), color=color,
            #         alpha=0.5, lw=1.6, ls="--", zorder=2, label="3-year running mean")

            # best analogue within the common window
            rw = r.where(r["year"] != p2.REF_YEAR)
            by = int(rw.idxmax("year"))
            br = float(rw.sel(year=by))
            row = best[(best.variant == variant) & (best.dataset == key)].iloc[0]
            if int(row.best_year) != by:
                print(f"  {key}: best season in {y0}-{y1} is {by}, "
                      f"over the full record {int(row.best_year)}")
            ax.plot([by], [br], marker="o", ms=5.5, mfc=color, mec=vz.SURFACE,
                    mew=1.2, zorder=5, clip_on=False)
            right_side = by < xmin + 0.72 * (xmax - xmin)
            ax.annotate(f"{by}  r={br:+.2f}", (by, br), textcoords="offset points",
                        xytext=(8 if right_side else -8, 1), fontsize=7.4,
                        color=vz.INK, zorder=6, va="center",
                        ha="left" if right_side else "right")

            ax.axvline(p2.REF_YEAR, color=vz.INK, lw=0.8, ls=":", zorder=4)
            if i == 0:
                ax.text(p2.REF_YEAR, -0.78, "2026 ", fontsize=7, color=vz.INK_SOFT,
                        ha="right", va="bottom", zorder=6)

            if j == 0:
                d = p2.DATASETS[key]
                ax.set_ylabel(d.label.replace(" (", "\n("), rotation=0,
                              ha="right", va="center", labelpad=12,
                              fontsize=8.5, color=vz.INK, linespacing=1.4)
            if i == 0 and len(variants) > 1:
                ax.set_title(p2.VARIANTS[variant]["title"], fontsize=9.5,
                             color=vz.INK, pad=7)
            if i == len(keys) - 1:
                ax.set_xlabel("year")

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-0.85, 0.85)
            ax.set_yticks([-0.5, 0, 0.5])

    handles, labels = axes[0][-1].get_legend_handles_labels()
    axes[0][-1].legend(handles, labels, loc="lower left", fontsize=7.6, ncol=2,
                       borderaxespad=0.3)

    fig.suptitle("How close do the observed and km-scale summers come to JJA 2026?",
                 fontsize=12, color=vz.INK, y=1 - 0.28 / fig_h)
    fig.text(0.5, 1 - 0.5 / fig_h,
             f"pattern correlation against ERA5 JJA 2026, 30–60°N, 80°W–40°E  ·  {y0}–{y1}, years both records cover\n"
             f"{p2.VARIANTS[variants[0]]['long']}  ·  dotted line marks 2026",
             ha="center", va="top", linespacing=1.5, fontsize=8.2, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig2_pattern_correlation_timeseries_short.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
