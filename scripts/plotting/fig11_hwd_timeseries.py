"""
Figure 11 -- JJA heatwave days over southern Europe, season by season, in the
observations and in every simulation with daily Tmax; the heatwave
counterpart of figure 2's pattern-correlation time series.

  y   land-mean JJA heatwave days in src.heatwave.MEAN_BOX (35-60N, 10W-30E;
      Xu et al. 2026 definition, 1991-2020 reference of the same dataset)

Two figures, one row per dataset (MPI-ESM1.2-ER has no daily Tmax):

  fig11_heatwave_days_timeseries.png        each record's full length, on the
                                            shared axis of figure 2
  fig11_heatwave_days_timeseries_short.png  1990-2026, the years the
                                            observations cover, so ERA5 and
                                            the models compare year by year

Ensembles are drawn as in figure 2, by member count: MPI-GE's 50 members as a
5-95% band, with the ensemble mean on top (the forced part of the heatwave
count -- figure 2 draws the member maximum instead, because that is where its
best analogue comes from); EERIE's 3 members individually, r1 (the only one
running past 2020) in the dataset's colour. ERA5's 2003, 2015 and 2026 are
marked, and the short figure adds ERA5's least-squares linear trend.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import land_mean_series
from fig8_hwd_vs_corr import ERA5_MARKED

ROW_ORDER = ["ERA5", "MPI-GE", "ICON-EPOC-hist", "EERIE"]
SHORT = (1990, p2.REF_YEAR)
#: the EERIE member drawn in colour: r1 is the only one past 2020
HIGHLIGHT_MEMBER = {"EERIE": 1}


def load_series():
    out = {}
    for key in ROW_ORDER:
        if not hw.metrics_file(key).exists():
            print(f"!! {key}: no heatwave metrics, row skipped")
            continue
        with xr.open_dataset(hw.metrics_file(key)) as ds:
            out[key] = land_mean_series(ds.hwd.load())
    return out


#: ERA5 label offsets (points, ha) on the compressed full-length axis, where
#: 2003 and 2015 sit too close for side-by-side labels
FULL_OFFSET = {2003: (-3, 9, "right"), 2015: (-8, 28, "right"), p2.REF_YEAR: (7, 2, "left")}


def draw_row(ax, key, s, xlim, *, compact=False):
    color, ls = vz.COLORS[key], vz.LINESTYLES[key]
    ax.grid(axis="y", zorder=0)
    members = vz.valid_members(s) if "member" in s.dims else []

    if len(members) >= vz.BAND_MIN_MEMBERS:
        lo, hi, band_label = vz.member_band(s)
        ax.fill_between(s.year, lo, hi, color=color, alpha=0.22, lw=0, zorder=2,
                        label=band_label)
        ax.plot(s.year, s.mean("member"), color=color, lw=1.2, ls=ls, zorder=3,
                label="ensemble mean")
    elif members:
        hl = HIGHLIGHT_MEMBER.get(key)
        others = [m for m in members if m != hl]
        for n, m in enumerate(others):
            ax.plot(s.year, s.sel(member=m), color=vz.CONTEXT, lw=0.8, zorder=2,
                    label="other members" if n == 0 else None)
        if hl in members:
            ax.plot(s.year, s.sel(member=hl), color=color, lw=1.2, ls=ls, zorder=3,
                    label=f"member r{hl}")
    else:
        ax.plot(s.year, s, color=color, lw=1.2, ls=ls, zorder=3, marker="o", ms=2.2)

    if key == "ERA5":
        for y in ERA5_MARKED:
            v = float(s.sel(year=y))
            ax.plot([y], [v], marker="o", ms=5.5, mfc=color, mec=vz.SURFACE, mew=1.2,
                    zorder=5, clip_on=False)
            if compact:
                dx, dy, ha = FULL_OFFSET.get(y, (7, 2, "left"))
                text = f"{y}" if y != p2.REF_YEAR else f"{y}  {v:.1f} d"
            else:
                right = y < xlim[0] + 0.85 * (xlim[1] - xlim[0])
                dx, dy, ha = (7, 2, "left") if right else (-7, 2, "right")
                text = f"{y}  {v:.1f} d"
            leader = dict(arrowstyle="-", color=vz.INK_SOFT, lw=0.6, shrinkA=1, shrinkB=3) \
                if abs(dy) > 15 else None
            ax.annotate(text, (y, v), textcoords="offset points", xytext=(dx, dy),
                        fontsize=7.4, color=vz.INK, zorder=6, va="center", ha=ha,
                        arrowprops=leader)


def era5_trend(ax, s):
    """Least-squares linear trend of the ERA5 series, drawn and returned (days/decade)."""
    s = s.dropna("year")
    b = np.polyfit(s.year, s, 1)
    xs = np.array([int(s.year.min()), int(s.year.max())])
    ax.plot(xs, np.polyval(b, xs), color=vz.INK, lw=0.9, ls=(0, (4, 3)), zorder=4,
            label=f"linear trend, {b[0] * 10:+.1f} days/decade")
    return b[0] * 10


def plot(series, xlim, out, title, note, *, trend=False):
    vz.use_style()
    keys = list(series)
    fig_h = 1.5 * len(keys) + 1.6
    fig, axes = plt.subplots(len(keys), 1, figsize=(8.5, fig_h), sharex=True, sharey=True,
                             squeeze=False)
    axes = axes[:, 0]
    fig.subplots_adjust(hspace=0.24, top=1 - 1.15 / fig_h, bottom=0.62 / fig_h,
                        left=0.2, right=0.985)

    ymax = max(float(s.max()) for s in
               (v.sel(year=slice(*xlim)) for v in series.values()))
    for ax, key in zip(axes, keys):
        s = series[key].sel(year=slice(*xlim))
        draw_row(ax, key, s, xlim, compact=not trend)
        if key == "ERA5" and trend:
            era5_trend(ax, s)
        d = p2.DATASETS[key]
        ax.set_ylabel(d.label.replace(" (", "\n("), rotation=0, ha="right", va="center",
                      labelpad=12, fontsize=8.5, color=vz.INK, linespacing=1.4)
        ax.set_xlim(xlim[0] - 0.5, xlim[1] + 0.5)
        ax.set_ylim(0, ymax * 1.12)
        ax.axvline(p2.REF_YEAR, color=vz.INK, lw=0.8, ls=":", zorder=4)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="upper left", fontsize=7.6, ncol=2,
                      borderaxespad=0.3)
    axes[-1].set_xlabel("year")

    fig.suptitle(title, fontsize=12, color=vz.INK, y=1 - 0.28 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             f"JJA heatwave days, land mean over {hw.MEAN_BOX_LABEL}  ·  ≥3 days with Tmax "
             "anomaly above the calendar-day 90th percentile (Xu et al. 2026), 1991–2020 "
             f"reference of the same dataset\n{note}",
             ha="center", va="top", fontsize=8.0, color=vz.INK_SOFT, linespacing=1.4)
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


def main():
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    series = load_series()
    full = (min(int(s.year.min()) for s in series.values()),
            max(int(s.year.max()) for s in series.values()))
    title = "European summer heatwave days, season by season"
    plot(series, full, p2.FIG_DIR / "fig11_heatwave_days_timeseries.png", title,
         "dotted line: 2026  ·  ERA5 2003, 2015 and 2026 marked  ·  "
         "MPI-ESM1.2-ER omitted (no daily output)")
    plot(series, SHORT, p2.FIG_DIR / "fig11_heatwave_days_timeseries_short.png",
         f"{title}, {SHORT[0]}–{SHORT[1]}",
         "dotted line: 2026  ·  ERA5 2003, 2015 and 2026 marked; dashed: ERA5 least-squares "
         "linear trend  ·  MPI-ESM1.2-ER omitted (no daily output)", trend=True)
    for key, s in series.items():
        s = s.mean("member") if "member" in s.dims else s
        ref = s.sel(year=slice(*hw.REF_PERIOD)).mean()
        print(f"  {key:16s} {int(s.dropna('year').year.min())}-{int(s.dropna('year').year.max())}  "
              f"1991-2020 mean {float(ref):.1f} days")


if __name__ == "__main__":
    main()
