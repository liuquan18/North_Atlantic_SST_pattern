"""
Figure 18 -- return period of every observed summer in the MPI-GE climate of
its year (03_mpige_window_rp.py).

(a) ERA5 heatwave days as observed (black) and with the long-term
    trend removed (grey) -- the series and the quadratic trend of fig 17,
    taken from its own loader so the two figures cannot drift apart.
(b) the "map": for each year, the stationary GEV fitted across the 50
    MPI-GE members in the 10 years around it (500 seasons) -- its 1-in-2,
    1-in-10 and 1-in-100 year heatwave days, shaded between -- with the
    observed ERA5 series on top. Faint horizontal lines: the 2003 and 2026
    values, whose return periods across all climates are drawn in (c).
(c) the return period of each ERA5 summer under the GEV of its year (log
    axis): line with 95% bootstrap band. Dashed / dotted:
    how rare a summer of the 2026 / 2003 size would be in the MPI-GE climate
    of every year (same GEVs), labelled in the 1950s, 2003 and 2026 -- the
    cross-climate comparison.

Index and layout follow fig 17 (scripts/plotting/fig17_era5_hwd_timeseries.py):
land-mean JJA heatwave days over src.heatwave.MEAN_BOX, 180 mm wide, caption
below, white background. ERA5 is black in every row. The same figure for TXx (hottest day of land-mean
daily Tmax) is fig 19: fig19_era5_txx_return_period_mpige.py calls main("txx").
Return periods beyond RP_CAP (or above a bounded GEV's upper end) are
extrapolation far past 500 seasons and are drawn at RP_CAP.

main(index, method="rank") draws the figure without any fitted distribution
and without return periods (05_mpige_window_rank.py), in two panels: (a) the
observed ERA5 series on the 50th / 90th / 99th percentiles of the 500
MPI-GE summers, (b) the percentile rank (0-100%) of each ERA5 summer among
them, binomial 95% interval, 50% line -- figs 20 (heatwave days) and 21 (TXx).
"""
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import src.extremes as ex
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig17_era5_hwd_timeseries import TREND_DEGREE
from fig17_era5_hwd_timeseries import load as load_era5

DIR = ROOT / "data" / "heatwave_stats"
MARKED = [2003, p2.REF_YEAR]
#: ERA5 in every row; in (a) the trend-removed series is the lighter one
OBS = vz.INK
DETRENDED = vz.INK_MUTED
BAND = vz.COLORS["MPI-GE"]
#: (c): years at which the fixed-value curves are labelled; 1950s = mean RP of 1950-1959
LABEL_AT = {"1950s": (1950, 1959), "2003": (2003, 2003), "2026": (2026, 2026)}
FIXED_STYLE = {2026: (0, (4, 2)), 2003: (0, (1, 1.6))}
#: rank figures: the further fixed-value curves (summers ranked >= 90%, e.g. 1947)
OTHER_STYLE = (0, (5, 1.5, 1, 1.5))
#: (c): backing of the labels drawn over curves
WHITE = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6)
#: (c): fixed-value curves are drawn no higher than this (years)
RP_CAP = 1e4
FIG_SIZE = (7.09, 7.6)
FS = 7.5
CAPTION_SPACE = 0.17
CAPTION_TITLE = "How rare was each observed summer in the climate of its year?"
CAPTION_TITLE_RANK = "Where does each observed summer sit among the simulated summers of its climate?"
CAPTION_HWD = (
    f"(a) ERA5 June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL}: as observed "
    "(black) and as departure from a quadratic long-term trend fitted to 1940–2026 (grey), as in "
    "Fig. 17. (b) Heatwave days that MPI-ESM1.2-LR reaches once in 2, 10 and 100 summers, from a "
    "GEV fitted across its 50 members in the 10 years around each year (y−5 to y+4, 500 "
    "summers); black: ERA5; faint lines: the 2003 and 2026 values. (c) Return period of each "
    "ERA5 summer under the GEV of its year, with 95% bootstrap interval (black). Dashed / dotted: "
    "return period of a "
    "summer as large as 2026 / 2003 in the MPI-GE climate of every year. A heatwave is ≥3 days "
    "with daily maximum temperature above the calendar-day 90th percentile of 1991–2020 of the "
    "same dataset (Xu et al. 2026). MPI-GE is used as simulated: its summers vary about 1.4 "
    "times as much as ERA5's, which shortens return periods, and it warms more slowly than "
    "observed after about 2015.")
CAPTION_TXX = (
    "(a) ERA5 TXx — the hottest day of June–August in daily maximum temperature averaged over "
    f"land in {hw.MEAN_BOX_LABEL} — as anomaly from its 1991–2020 mean (black) and as departure "
    "from a quadratic long-term trend fitted to 1940–2026 (grey). (b) TXx anomaly (against "
    "MPI-GE's own 1991–2020 mean) that MPI-ESM1.2-LR reaches once in 2, 10 and 100 summers, from "
    "a GEV fitted across its 50 members in the 10 years around each year (y−5 to y+4, 500 "
    "summers); black: ERA5; faint lines: the 2003 and 2026 values. (c) Return period of each "
    "ERA5 summer under the GEV of its year, with 95% bootstrap interval (black). Dashed / dotted: "
    "return period of a "
    "TXx as high as 2026 / 2003 in the MPI-GE climate of every year; drawn at 10,000 years where "
    "the fit gives more or puts the value above its upper bound. MPI-GE TXx varies about 1.1 "
    "times as much as ERA5's and warms more slowly than observed after about 2015.")
RANK_B = ("blue: {what}: the 50th, 90th and 99th percentiles (solid) and the mean (dashed) of "
          "its 50 members in the 10 years around each year (y−5 to y+4, 500 summers); the "
          "percentiles are exceeded by half, a tenth and a hundredth of them, no fitted "
          "distribution. (b) Percentile rank of each observed ERA5 "
          "summer among those 500: the share of MPI-GE summers below it (ties counted half), "
          "with 95% binomial interval (grey); 100% is above all 500, the thin line at 50% the "
          "MPI-GE median. Marked: 2003, 2026 and every summer ranked at or above 90%. Grey "
          "lines: the percentile rank of a summer as large as each marked one in every year's "
          "window (legend). ")
CAPTION_HWD_RANK = (
    f"(a) ERA5 June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL} (black), as "
    "in Fig. 17; " + RANK_B.format(what="heatwave days of MPI-ESM1.2-LR") +
    "A heatwave is ≥3 days with daily maximum temperature above the calendar-day 90th percentile "
    "of 1991–2020 of the same dataset (Xu et al. 2026). MPI-GE is used as simulated: its summers "
    "vary about 1.4 times as much as ERA5's, which lowers the ranks of large observed "
    "summers, and it warms more slowly than observed after about 2015.")
CAPTION_TXX_RANK = (
    "(a) ERA5 TXx — the hottest day of June–August in daily maximum temperature averaged over "
    f"land in {hw.MEAN_BOX_LABEL} — as anomaly from its 1991–2020 mean (black); " +
    RANK_B.format(what="TXx anomaly of MPI-ESM1.2-LR (against its own 1991–2020 mean)") +
    "MPI-GE TXx varies about 1.1 times as much as ERA5's and warms more slowly than observed "
    "after about 2015.")

#: per method: input tables, output name suffix, legend name of the estimate, the
#: three level columns of (b) with their labels (the last one dashed)
METHODS = {
    "gev": dict(windows="mpige_window_gev", rp="era5_rp_in_mpige", out="", name="GEV",
                levels=(("rl2", "1-in-2"), ("rl10", "1-in-10"), ("rl100", "1-in-100"))),
    "rank": dict(windows="mpige_window_rank", rp="era5_rank_in_mpige", out="_rank", name="rank",
                 levels=(("q50", "50th pct"), ("q90", "90th pct"), ("q99", "99th pct"))),
}
OUT_RANK = {"hwd": "fig20_era5_return_period_mpige_rank",
            "txx": "fig21_era5_txx_return_period_mpige_rank"}

#: per index: input suffix, output name, axis label, value format, panel (a) loader, caption
SPECS = {
    "hwd": dict(suffix="", out="fig18_era5_return_period_mpige", ylabel="Heatwave days (JJA)",
                fmt=lambda v: f"{v:.1f} d", series="hwd", caption=CAPTION_HWD,
                caption_rank=CAPTION_HWD_RANK,
                tag_a="ERA5, observed and with the quadratic trend removed"),
    "txx": dict(suffix="_txx", out="fig19_era5_txx_return_period_mpige",
                ylabel="TXx anomaly (K)", fmt=lambda v: f"{v:+.1f} K", series="txx",
                caption=CAPTION_TXX, caption_rank=CAPTION_TXX_RANK,
                tag_a="ERA5 TXx, anomaly and trend removed"),
}


def load_series(kind):
    """ERA5 years, recorded series and departure from a quadratic trend (as fig 17)."""
    if kind == "hwd":
        return load_era5()
    t = pd.read_csv(DIR / "txx_table.csv")
    e = t[t.dataset == "ERA5"].sort_values("year")
    years, raw = e.year.values, e.txx_anom.values
    trend = np.polyval(np.polyfit(years.astype(float), raw, TREND_DEGREE), years)
    return years, raw, raw - trend


def padded(lo, hi, below=0.1, above=0.15):
    span = hi - lo
    return lo - below * span, hi + above * span


def label(ax, x, y, text, dy=3, **kw):
    kw.setdefault("color", vz.INK)
    ax.annotate(text, (x, y), xytext=(0, dy), textcoords="offset points", ha="center",
                va="bottom", fontsize=FS, zorder=6, **kw)


def fmt_rp(v):
    return f"{v:,.0f} yr" if v >= 10 else f"{v:.1f} yr"


def window_gevs(gev):
    return [ex.NSGEV(r.mu, 0.0, float(np.log(r.sigma)), 0.0, r.xi) for _, r in gev.iterrows()]


def panel_series(ax, years, raw, anom, spec):
    ax.plot(years, raw, color=OBS, lw=1.2, zorder=3, label="observed")
    ax.plot(years, anom, color=DETRENDED, lw=1.0, zorder=2, label="trend removed")
    ax.axhline(0, color=vz.INK_SOFT, lw=0.5, zorder=1)
    for y in MARKED:
        k = int(np.flatnonzero(years == y)[0])
        ax.plot(y, raw[k], "o", ms=3.5, color=OBS, mec="white", mew=0.6, zorder=4)
        label(ax, y, raw[k], str(y), fontweight="bold")
        # the trend-removed value's label goes below its point when the two nearly coincide
        close = abs(raw[k] - anom[k]) < 0.12 * (raw.max() - min(raw.min(), anom.min()))
        if close:
            ax.annotate(str(y), (y, anom[k]), xytext=(0, -3), textcoords="offset points",
                        ha="center", va="top", fontsize=FS, color=DETRENDED, zorder=6)
        else:
            label(ax, y, anom[k], str(y), dy=2.5, color=DETRENDED)
    ax.set_ylim(*padded(min(anom.min(), raw.min(), 0), max(raw.max(), anom.max())))
    ax.set_ylabel(spec["ylabel"], fontsize=FS)
    ax.legend(loc="upper left", bbox_to_anchor=(0.56, 1.03), frameon=False, fontsize=FS,
              handlelength=1.6, ncol=2)


def draw_levels(ax, gev, levels, styles=("-", "-", (0, (3, 1.5))), extra=()):
    """
    The window's three levels as shaded blue bands, plus any `extra` (column,
    name, linestyle) lines, labelled at the right end; labels that would touch
    are pushed apart upwards.
    """
    yr = gev.year.values
    (c1, _), (c2, _), (c3, _) = levels
    ax.fill_between(yr, gev[c1], gev[c2], color=BAND, alpha=0.30, lw=0, zorder=1)
    ax.fill_between(yr, gev[c2], gev[c3], color=BAND, alpha=0.13, lw=0, zorder=1)
    lines = [(col, name, ls) for (col, name), ls in zip(levels, styles)] + list(extra)
    for col, _, ls in lines:
        ax.plot(yr, gev[col], color=BAND, lw=0.8, ls=ls, zorder=2)
    gap = 0.09 * (gev[c3].max() - gev[c1].min())
    ends = sorted((float(gev[col].iloc[-1]), name) for col, name, _ in lines)
    y_prev = -np.inf
    for y, name in ends:
        y_prev = max(y, y_prev + gap)
        ax.text(yr[-1] + 1.2, y_prev, name, va="center", fontsize=FS - 0.5, color=BAND)


def panel_map(ax, gev, rp, values, spec, levels):
    draw_levels(ax, gev, levels)
    (c1, _), _, (c3, _) = levels
    # for ev, v in values.items():
        # ax.axhline(v, color=vz.INK_MUTED, lw=0.5, ls=FIXED_STYLE[ev], zorder=1)
    ax.plot(rp.year, rp.value, color=vz.INK, lw=1.1, zorder=3)
    for y in MARKED:
        label(ax, y, values[y], str(y), fontweight="bold")
        ax.plot(y, values[y], "o", ms=3.5, color=vz.INK, mec="white", mew=0.6, zorder=4)
    ax.set_ylabel(spec["ylabel"], fontsize=FS)
    lo, hi = min(rp.value.min(), gev[c1].min()), max(rp.value.max(), gev[c3].max())
    ax.set_ylim(0 if lo >= 0 else padded(lo, hi)[0], hi + 0.12 * (hi - min(lo, 0)))


def panel_series_levels(ax, win, years, raw, spec, levels, marked):
    """(a) for method="rank": the observed ERA5 series on the MPI-GE window percentiles."""
    draw_levels(ax, win, levels, styles=("-", "-", "-"), extra=(("mean", "mean", "--"),))
    ax.plot(years, raw, color=OBS, lw=1.2, zorder=3)
    for y in marked:
        k = int(np.flatnonzero(years == y)[0])
        ax.plot(y, raw[k], "o", ms=3.5, color=OBS, mec="white", mew=0.6, zorder=4)
        label(ax, y, raw[k], str(y), fontweight="bold")
    ax.set_ylim(*padded(min(raw.min(), win[levels[0][0]].min(), 0),
                        max(raw.max(), win[levels[-1][0]].max())))
    ax.set_ylabel(spec["ylabel"], fontsize=FS)
    handles = [
        plt.Line2D([], [], color=OBS, lw=1.2, label="ERA5"),
        plt.Rectangle((0, 0), 1, 1, color=BAND, alpha=0.3, lw=0, label="MPI-GE percentiles"),
        plt.Line2D([], [], color=BAND, lw=0.8, ls="--", label="MPI-GE mean"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.05, 1.03), frameon=False,
              fontsize=FS, handlelength=1.6, ncol=3)


def fixed_curves(win, values):
    """
    Return period of a 2003- / 2026-sized summer under every window's GEV, and
    where it is only a lower bound: beyond RP_CAP or the GEV's upper end.
    """
    fits = window_gevs(win)
    raw = {ev: np.array([f.return_period(v, 0.0) for f in fits]) for ev, v in values.items()}
    return ({ev: np.minimum(c, RP_CAP) for ev, c in raw.items()},
            {ev: c >= RP_CAP for ev, c in raw.items()})


def panel_rp(bx, gev, rp, values, fixed, bound, spec, estimate):
    yr = gev.year.values
    for ev, curve in fixed.items():
        bx.plot(yr, curve, color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[ev], zorder=2)
        for name, (a, b) in LABEL_AT.items():
            sel = (yr >= a) & (yr <= b)
            v = float(np.exp(np.mean(np.log(curve[sel]))))   # geometric mean over the decade
            x = (a + b) / 2
            if name == str(ev):          # the event's own year: that is the main line's label
                continue
            bx.plot(x, v, "o", ms=2.6, color=vz.INK_SOFT, zorder=3)
            # a decade that contains capped windows only gives a lower bound
            txt = ("≥" if bound[ev][sel].any() else "") + fmt_rp(v)
            # at 2003 the event's own bold label sits just below: go above the point
            up = name == "2003"
            bx.annotate(txt + (" (1950s)" if name == "1950s" else ""), (x, v),
                        xytext=(0, 4 if up else -4), textcoords="offset points", ha="center",
                        va="bottom" if up else "top",
                        fontsize=FS - 1, color=vz.INK_SOFT, zorder=5, bbox=WHITE)

    finite = np.isfinite(rp.rp)
    bx.fill_between(rp.year, rp.rp_lo, np.minimum(rp.rp_hi, RP_CAP), color=OBS, alpha=0.15,
                    lw=0, zorder=1)
    bx.plot(rp.year[finite], rp.rp[finite], color=OBS, lw=1.2, zorder=3)
    bx.scatter(rp.year[finite], rp.rp[finite], s=7, color=OBS, zorder=4, lw=0)
    for y in MARKED:
        r = rp[rp.year == y].iloc[0]
        txt = "∞" if r.get("above_upper_end", False) else fmt_rp(r.rp)
        label(bx, y, r.rp if np.isfinite(r.rp) else r.rp_lo, f"{y}\n{txt}", dy=4,
              fontweight="bold", linespacing=1.1, bbox=WHITE)
    bx.set_yscale("log")
    top = min(max(np.nanmax(fixed[2026]), np.nanmax(rp.rp_hi[np.isfinite(rp.rp_hi)])), RP_CAP)
    # headroom above the curves for the legend (log axis: a decade or so)
    bx.set_ylim(1, max(1e5, top * 100))
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
             if t <= top * 3]
    bx.set_yticks(ticks, [f"{t:,}" for t in ticks])
    bx.set_ylabel("Return period (years)", fontsize=FS)
    handles = [
        plt.Line2D([], [], color=OBS, lw=1.2, marker="o", ms=2.5, label=f"ERA5 summer, {estimate}"),
        plt.Rectangle((0, 0), 1, 1, color=OBS, alpha=0.15, lw=0, label="95% interval"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2026],
                   label=f"2026-like summer ({spec['fmt'](values[2026])})"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2003],
                   label=f"2003-like summer ({spec['fmt'](values[2003])})"),
    ]
    bx.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.9), ncol=2,
              frameon=False, fontsize=FS - 0.5, handlelength=2.2, columnspacing=1.2,
              borderaxespad=0.2)


def panel_pct(bx, win, rk, values, spec, marked):
    """(c) for method="rank": percentile rank of each ERA5 summer, 0-100% (linear)."""
    yr = win.year.values
    for ev in values:
        curve = win[f"pct_fixed_{ev}"].values
        bx.plot(yr, curve, color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE.get(ev, OTHER_STYLE),
                zorder=2)
        for name, (a, b) in LABEL_AT.items():
            sel = (yr >= a) & (yr <= b)
            v, x = float(curve[sel].mean()), (a + b) / 2
            # near the event's own year its bold label speaks for the curve; in the 1950s
            # all curves crowd 100%, so only the 2003 and 2026 ones are labelled there
            if abs(x - ev) < 10 or (name == "1950s" and ev not in FIXED_STYLE):
                continue
            bx.plot(x, v, "o", ms=2.6, color=vz.INK_SOFT, zorder=3)
            # both curves crowd 100%: the 2026 one is labelled above, the 2003 one below;
            # labels run rightwards, clear of the marked summers' bold labels (1947, 1994)
            # to their left: at 2003 the 2026 label steps right of the event's own label,
            # at the right end it goes beside the point, into the margin
            up = ev == 2026
            if name == "2026":
                xy, ha, va = (5, 0), "left", "center"
            else:
                xy = (16 if name == "2003" else -6, 4 if up else -4)
                ha = "left"
                va = "bottom" if up else "top"
            bx.annotate(f"{v:.1f}%" + (" (1950s)" if name == "1950s" else ""), (x, v),
                        xytext=xy, textcoords="offset points", ha=ha, va=va,
                        fontsize=FS - 1, color=vz.INK_SOFT, zorder=5, bbox=WHITE)

    bx.fill_between(rk.year, rk.pct_lo, rk.pct_hi, color=OBS, alpha=0.15, lw=0, zorder=1)
    bx.plot(rk.year, rk.pct, color=OBS, lw=1.2, zorder=3)
    bx.scatter(rk.year, rk.pct, s=7, color=OBS, zorder=4, lw=0)
    for y in marked:
        r = rk[rk.year == y].iloc[0]
        label(bx, y, r.pct, f"{y}\n{r.pct:.1f}%", dy=4, fontweight="bold", linespacing=1.1,
              bbox=WHITE)
    # the axis itself spans 0-100%; the space above it holds the legend
    bx.set_ylim(0, 128)
    bx.set_yticks(range(0, 101, 20), [f"{t}%" for t in range(0, 101, 20)])
    bx.spines["left"].set_bounds(0, 100)
    bx.axhline(50, color=vz.INK_SOFT, lw=0.5, zorder=0)
    bx.set_ylabel("Rank among MPI-GE summers", fontsize=FS)
    handles = [
        plt.Line2D([], [], color=OBS, lw=1.2, marker="o", ms=2.5, label="ERA5 summer"),
        plt.Rectangle((0, 0), 1, 1, color=OBS, alpha=0.15, lw=0, label="95% interval"),
    ] + [plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE.get(ev, OTHER_STYLE),
                    label=f"{ev}-like summer ({spec['fmt'](values[ev])})")
         for ev in sorted(values, reverse=True)]
    bx.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.02), ncol=3,
              frameon=False, fontsize=FS - 0.5, handlelength=2.2, columnspacing=1.2,
              borderaxespad=0.2)


def tag(ax, letter):
    ax.text(0.01, 0.97, f"({letter})", transform=ax.transAxes, va="top", fontsize=FS + 0.5,
            fontweight="bold")


def main(index="hwd", method="gev"):
    spec, meth = SPECS[index], METHODS[method]
    vz.use_style()
    gev = pd.read_csv(DIR / f"{meth['windows']}{spec['suffix']}.csv")
    rp = pd.read_csv(DIR / f"{meth['rp']}{spec['suffix']}.csv")
    years, raw, anom = load_series(spec["series"])
    values = {y: float(rp.loc[rp.year == y, "value"].iloc[0]) for y in MARKED}

    # the rank figure merges the series and the levels into one panel: shorter, with the
    # caption kept at the same height in inches
    height = FIG_SIZE[1] - (0 if method == "gev" else 1.3)
    scale = FIG_SIZE[1] / height
    caption_space = CAPTION_SPACE * scale
    if method == "gev":
        fig, axes = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True,
                                 gridspec_kw=dict(height_ratios=[0.85, 1, 1.25], hspace=0.14))
        sx, ax, bx = axes
        panel_series(sx, years, raw, anom, spec)
        panel_map(ax, gev, rp, values, spec, meth["levels"])
        panel_rp(bx, gev, rp, values, *fixed_curves(gev, values), spec, meth["name"])
    else:
        fig, axes = plt.subplots(2, 1, figsize=(FIG_SIZE[0], height), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.15, 1.25], hspace=0.14))
        ax, bx = axes
        # 05 writes a fixed-value curve for 2003, 2026 and every summer ranked >= 90%
        marked = sorted(int(c.rsplit("_", 1)[1]) for c in gev if c.startswith("pct_fixed_"))
        values = {y: float(rp.loc[rp.year == y, "value"].iloc[0]) for y in marked}
        panel_series_levels(ax, gev, years, raw, spec, meth["levels"], marked)
        panel_pct(bx, gev, rp, values, spec, marked)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.9, top=0.98, bottom=caption_space + 0.055 * scale)
    for a, letter in zip(axes, "abc"):
        tag(a, letter)

    for a in axes:
        a.set_facecolor("white")
        a.grid(False)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_linewidth(0.6)
            a.spines[side].set_color(vz.INK)
        a.tick_params(labelsize=FS, direction="out", length=2.5, width=0.6, pad=2,
                      colors=vz.INK)
    bx.set_xlim(gev.year.iloc[0] - 1, gev.year.iloc[-1] + 2)

    title = CAPTION_TITLE if method == "gev" else CAPTION_TITLE_RANK
    fig.text(0.012, caption_space - 0.005 * scale, title, ha="left", va="top", fontsize=FS,
             fontweight="bold", color=vz.INK)
    caption = spec["caption"] if method == "gev" else spec["caption_rank"]
    fig.text(0.012, caption_space - 0.033 * scale, textwrap.fill(caption, 122), ha="left",
             va="top", fontsize=FS, color=vz.INK, linespacing=1.45)
    out = p2.FIG_DIR / f"{spec['out'] if method == 'gev' else OUT_RANK[index]}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
