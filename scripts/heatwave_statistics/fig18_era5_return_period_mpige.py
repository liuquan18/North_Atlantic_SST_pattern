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
#: (c): backing of the labels drawn over curves
WHITE = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6)
#: (c): fixed-value curves are drawn no higher than this (years)
RP_CAP = 1e4
FIG_SIZE = (7.09, 7.6)
FS = 7.5
CAPTION_SPACE = 0.17
CAPTION_TITLE = "How rare was each observed summer in the climate of its year?"
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

#: per index: input suffix, output name, axis label, value format, panel (a) loader, caption
SPECS = {
    "hwd": dict(suffix="", out="fig18_era5_return_period_mpige", ylabel="Heatwave days (JJA)",
                fmt=lambda v: f"{v:.1f} d", series="hwd", caption=CAPTION_HWD,
                tag_a="ERA5, observed and with the quadratic trend removed"),
    "txx": dict(suffix="_txx", out="fig19_era5_txx_return_period_mpige",
                ylabel="TXx anomaly (K)", fmt=lambda v: f"{v:+.1f} K", series="txx",
                caption=CAPTION_TXX,
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


def panel_map(ax, gev, rp, values, spec):
    yr = gev.year.values
    ax.fill_between(yr, gev.rl2, gev.rl10, color=BAND, alpha=0.30, lw=0, zorder=1)
    ax.fill_between(yr, gev.rl10, gev.rl100, color=BAND, alpha=0.13, lw=0, zorder=1)
    for col, ls, name in (("rl2", "-", "1-in-2"), ("rl10", "-", "1-in-10"),
                          ("rl100", (0, (3, 1.5)), "1-in-100")):
        ax.plot(yr, gev[col], color=BAND, lw=0.8, ls=ls, zorder=2)
        ax.text(yr[-1] + 1.2, gev[col].iloc[-1], name, va="center", fontsize=FS - 0.5,
                color=BAND)
    # for ev, v in values.items():
        # ax.axhline(v, color=vz.INK_MUTED, lw=0.5, ls=FIXED_STYLE[ev], zorder=1)
    ax.plot(rp.year, rp.value, color=vz.INK, lw=1.1, zorder=3)
    for y in MARKED:
        label(ax, y, values[y], str(y), fontweight="bold")
        ax.plot(y, values[y], "o", ms=3.5, color=vz.INK, mec="white", mew=0.6, zorder=4)
    ax.set_ylabel(spec["ylabel"], fontsize=FS)
    lo, hi = min(rp.value.min(), gev.rl2.min()), max(rp.value.max(), gev.rl100.max())
    ax.set_ylim(0 if lo >= 0 else padded(lo, hi)[0], hi + 0.12 * (hi - min(lo, 0)))


def panel_rp(bx, gev, rp, values, fits, spec):
    yr = gev.year.values
    # how rare a 2003- / 2026-sized summer would be in every year's climate;
    # beyond RP_CAP (or above a bounded fit's upper end) it is drawn at RP_CAP
    raw_rp = {ev: np.array([f.return_period(v, 0.0) for f in fits]) for ev, v in values.items()}
    fixed = {ev: np.minimum(c, RP_CAP) for ev, c in raw_rp.items()}
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
            txt = ("≥" if (raw_rp[ev][sel] >= RP_CAP).any() else "") + fmt_rp(v)
            # at 2003 the event's own bold label sits just below: go above the point
            up = name == "2003"
            bx.annotate(txt + (" (1950s)" if name == "1950s" else ""), (x, v),
                        xytext=(0, 4 if up else -4), textcoords="offset points", ha="center",
                        va="bottom" if up else "top",
                        fontsize=FS - 1, color=vz.INK_SOFT, zorder=5, bbox=WHITE)

    finite = np.isfinite(rp.rp)
    bx.fill_between(rp.year, rp.rp_lo, rp.rp_hi, color=OBS, alpha=0.15, lw=0, zorder=1)
    bx.plot(rp.year[finite], rp.rp[finite], color=OBS, lw=1.2, zorder=3)
    bx.scatter(rp.year[finite], rp.rp[finite], s=7, color=OBS, zorder=4, lw=0)
    for y in MARKED:
        r = rp[rp.year == y].iloc[0]
        txt = "∞" if r.above_upper_end else fmt_rp(r.rp)
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
        plt.Line2D([], [], color=OBS, lw=1.2, marker="o", ms=2.5, label="ERA5 summer, GEV"),
        plt.Rectangle((0, 0), 1, 1, color=OBS, alpha=0.15, lw=0, label="95% interval"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2026],
                   label=f"2026-like summer ({spec['fmt'](values[2026])})"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2003],
                   label=f"2003-like summer ({spec['fmt'](values[2003])})"),
    ]
    bx.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.9), ncol=2,
              frameon=False, fontsize=FS - 0.5, handlelength=2.2, columnspacing=1.2,
              borderaxespad=0.2)


def tag(ax, letter):
    ax.text(0.01, 0.97, f"({letter})", transform=ax.transAxes, va="top", fontsize=FS + 0.5,
            fontweight="bold")


def main(index="hwd"):
    spec = SPECS[index]
    vz.use_style()
    gev = pd.read_csv(DIR / f"mpige_window_gev{spec['suffix']}.csv")
    rp = pd.read_csv(DIR / f"era5_rp_in_mpige{spec['suffix']}.csv")
    years, raw, anom = load_series(spec["series"])
    values = {y: float(rp.loc[rp.year == y, "value"].iloc[0]) for y in MARKED}
    fits = window_gevs(gev)

    fig, (sx, ax, bx) = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True,
                                     gridspec_kw=dict(height_ratios=[0.85, 1, 1.25],
                                                      hspace=0.14))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.9, top=0.98, bottom=CAPTION_SPACE + 0.055)

    panel_series(sx, years, raw, anom, spec)
    panel_map(ax, gev, rp, values, spec)
    panel_rp(bx, gev, rp, values, fits, spec)
    tag(sx, "a", )
    tag(ax, "b", )
    tag(bx, "c", )

    for a in (sx, ax, bx):
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

    fig.text(0.012, CAPTION_SPACE - 0.005, CAPTION_TITLE, ha="left", va="top", fontsize=FS,
             fontweight="bold", color=vz.INK)
    fig.text(0.012, CAPTION_SPACE - 0.033, textwrap.fill(spec["caption"], 122), ha="left",
             va="top", fontsize=FS, color=vz.INK, linespacing=1.45)
    out = p2.FIG_DIR / f"{spec['out']}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
