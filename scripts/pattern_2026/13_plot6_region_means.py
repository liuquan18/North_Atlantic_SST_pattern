"""
Figure 6 -- area-mean JJA SST anomaly through each simulation, for four regions.

One panel per region, stacked on a shared time axis and a shared temperature
axis. Sharing the y-axis is the point of the figure: it is what makes the
Mediterranean's amplitude visibly larger than the global ocean's rather than
something the reader has to infer by comparing tick labels between panels.

box = natl + med exactly, so the middle two panels decompose the first rather
than repeating it.

MPI-GE is a 50-member ensemble: the band is the 5-95% spread across members and
the line is the ensemble mean, so its forced response can be told apart from the
internal variability that the single-realisation records also contain.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

ORDER = ["global", "box", "natl", "med"]


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    ds = xr.open_dataset(p2.RESULT_DIR / "region_means.nc")

    keys = [k for k in p2.MAIN_KEYS if f"{k}|global" in ds.data_vars]
    xmin = min(int(ds[f"{k}|global"].year.min()) for k in keys)
    xmax = max(int(ds[f"{k}|global"].year.max()) for k in keys)

    fig_h = 1.55 * len(ORDER) + 1.75
    fig, axes = plt.subplots(len(ORDER), 1, figsize=(12.4, fig_h),
                             sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.22, top=1 - 1.62 / fig_h, bottom=0.55 / fig_h,
                        left=0.163, right=0.985)

    for ax, region in zip(axes, ORDER):
        ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
        ax.grid(axis="y", zorder=0)
        ax.axvline(p2.REF_YEAR, color=vz.INK_MUTED, lw=0.7, ls=":", zorder=1)

        for k in keys:
            da = ds[f"{k}|{region}"]
            color, ls = vz.COLORS[k], vz.LINESTYLES[k]
            label = p2.DATASETS[k].label

            if "member" in da.dims:
                ax.fill_between(da.year, da.quantile(0.05, "member"),
                                da.quantile(0.95, "member"), color=color,
                                alpha=0.20, lw=0, zorder=2,
                                label=f"{label} (5–95% of members)")
                ax.plot(da.year, da.mean("member"), color=color, lw=1.5, ls=ls,
                        zorder=3, label=f"{label} (ensemble mean)")
            else:
                ax.plot(da.year, da, color=color, lw=1.4, ls=ls, zorder=4,
                        label=label)

        # the observed 2026 value, which is what the study is about
        obs = ds[f"ERA5|{region}"]
        if p2.REF_YEAR in obs.year.values:
            v = float(obs.sel(year=p2.REF_YEAR))
            ax.plot([p2.REF_YEAR], [v], marker="o", ms=6, mfc=vz.INK,
                    mec=vz.SURFACE, mew=1.3, zorder=6, clip_on=False)
            ax.annotate(f"2026: {v:+.2f} °C", (p2.REF_YEAR, v),
                        textcoords="offset points", xytext=(11, 0), ha="left",
                        va="center", fontsize=8.2, color=vz.INK, zorder=7,
                        bbox=dict(facecolor="#ffffff", alpha=0.8,
                                  edgecolor="none", boxstyle="round,pad=0.22"))

        # set_ylabel places the caption clear of the tick labels for us
        info = p2.REGIONS[region]
        ax.set_ylabel(f"{info['label']}\n{info['detail']}", rotation=0,
                      ha="right", va="center", labelpad=12, fontsize=8.6,
                      color=vz.INK, linespacing=1.5)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(-1.45, 2.05)

    axes[-1].set_xlabel("year")
    # a figure-level legend: inside the top panel it sat on top of the data
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=8.3,
               bbox_to_anchor=(0.5, 1 - 0.80 / fig_h), columnspacing=2.4)

    fig.suptitle("Area-mean JJA sea surface temperature anomaly (°C), by region",
                 fontsize=12, color=vz.INK, y=1 - 0.28 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             "anomalies against each dataset's own 1991–2020 JJA climatology  ·  "
             "regional means on the common 1° ocean mask  ·  shared vertical scale  ·  "
             "dotted line marks 2026",
             ha="center", fontsize=8.3, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig6_region_mean_timeseries.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
