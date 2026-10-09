"""
Figure 18 -- return period of every observed summer in the MPI-GE climate of
its year (03_mpige_window_rp.py).

(a) ERA5 heatwave days as observed (light) and with the long-term trend
    removed (highlighted) -- the series and the quadratic trend of fig 17,
    taken from its own loader so the two figures cannot drift apart.
(b) the "map": for each year, the stationary GEV fitted across the 50
    MPI-GE members in the 10 years around it (500 seasons) -- its 1-in-2,
    1-in-10 and 1-in-100 year heatwave days, shaded between -- with the
    observed ERA5 series on top. Faint horizontal lines: the 2003 and 2026
    values, whose return periods across all climates are drawn in (c).
(c) the return period of each ERA5 summer under the GEV of its year (log
    axis): line with 95% bootstrap band; small open circles: the empirical
    estimate (k of 500 MPI-GE seasons at or above the observed value), with
    an arrow where none was and only a lower bound exists. Dashed / dotted:
    how rare a summer of the 2026 / 2003 size would be in the MPI-GE climate
    of every year (same GEVs), labelled in the 1950s, 2003 and 2026 -- the
    cross-climate comparison.

Index and layout follow fig 17 (scripts/plotting/fig17_era5_hwd_timeseries.py):
land-mean JJA heatwave days over src.heatwave.MEAN_BOX, 180 mm wide, caption
below, white background.
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
from fig17_era5_hwd_timeseries import LIGHT, TREND_DEGREE
from fig17_era5_hwd_timeseries import load as load_era5

DIR = ROOT / "data" / "heatwave_stats"
MARKED = [2003, p2.REF_YEAR]
HOT = vz.heatwave_cmap()(0.80)
BAND = vz.COLORS["MPI-GE"]
#: (c): years at which the fixed-value curves are labelled; 1950s = mean RP of 1950-1959
LABEL_AT = {"1950s": (1950, 1959), "2003": (2003, 2003), "2026": (2026, 2026)}
FIXED_STYLE = {2026: (0, (4, 2)), 2003: (0, (1, 1.6))}
#: (c): backing of the labels drawn over curves
WHITE = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6)
FIG_SIZE = (7.09, 7.6)
FS = 7.5
CAPTION_SPACE = 0.17
CAPTION_TITLE = "How rare was each observed summer in the climate of its year?"
CAPTION = (
    f"(a) ERA5 June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL}: as observed "
    "(grey) and as departure from a quadratic long-term trend fitted to 1940–2026 (red), as in "
    "Fig. 17. (b) Heatwave days that MPI-ESM1.2-LR reaches once in 2, 10 and 100 summers, from a "
    "GEV fitted across its 50 members in the 10 years around each year (y−5 to y+4, 500 "
    "summers); black: ERA5; faint lines: the 2003 and 2026 values. (c) Return period of each "
    "ERA5 summer under the GEV of its year, with 95% bootstrap interval; open circles: the plain "
    "count of MPI-GE summers at or above the observed value. Dashed / dotted: return period of a "
    "summer as large as 2026 / 2003 in the MPI-GE climate of every year. A heatwave is ≥3 days "
    "with daily maximum temperature above the calendar-day 90th percentile of 1991–2020 of the "
    "same dataset (Xu et al. 2026). MPI-GE is used as simulated: its summers vary about 1.4 "
    "times as much as ERA5's, which shortens return periods, and it warms more slowly than "
    "observed after about 2015.")


def label(ax, x, y, text, dy=3, **kw):
    kw.setdefault("color", vz.INK)
    ax.annotate(text, (x, y), xytext=(0, dy), textcoords="offset points", ha="center",
                va="bottom", fontsize=FS, zorder=6, **kw)


def fmt_rp(v):
    return f"{v:,.0f} yr" if v >= 10 else f"{v:.1f} yr"


def window_gevs(gev):
    return [ex.NSGEV(r.mu, 0.0, float(np.log(r.sigma)), 0.0, r.xi) for _, r in gev.iterrows()]


def panel_series(ax, years, raw, anom):
    ax.plot(years, raw, color=LIGHT, lw=1.0, zorder=2, label="observed")
    ax.plot(years, anom, color=HOT, lw=1.3, zorder=3, label="trend removed")
    ax.axhline(0, color=vz.INK_SOFT, lw=0.5, zorder=1)
    for y in MARKED:
        k = int(np.flatnonzero(years == y)[0])
        ax.plot(y, anom[k], "o", ms=3.5, color=HOT, mec="white", mew=0.6, zorder=4)
        label(ax, y, anom[k], str(y), fontweight="bold")
        label(ax, y, raw[k], str(y), dy=2.5, color=vz.INK_MUTED)
    ax.set_ylim(min(anom.min(), 0) * 1.2, raw.max() * 1.15)
    ax.set_ylabel("Heatwave days (JJA)", fontsize=FS)
    ax.legend(loc="upper left", bbox_to_anchor=(0.56, 1.03), frameon=False, fontsize=FS,
              handlelength=1.6, ncol=2)


def panel_map(ax, gev, rp, values):
    yr = gev.year.values
    ax.fill_between(yr, gev.rl2, gev.rl10, color=BAND, alpha=0.30, lw=0, zorder=1)
    ax.fill_between(yr, gev.rl10, gev.rl100, color=BAND, alpha=0.13, lw=0, zorder=1)
    for col, ls, name in (("rl2", "-", "1-in-2"), ("rl10", "-", "1-in-10"),
                          ("rl100", (0, (3, 1.5)), "1-in-100")):
        ax.plot(yr, gev[col], color=BAND, lw=0.8, ls=ls, zorder=2)
        ax.text(yr[-1] + 1.2, gev[col].iloc[-1], name, va="center", fontsize=FS - 0.5,
                color=BAND)
    for ev, v in values.items():
        ax.axhline(v, color=vz.INK_MUTED, lw=0.5, ls=FIXED_STYLE[ev], zorder=1)
    ax.plot(rp.year, rp.value, color=vz.INK, lw=1.1, zorder=3)
    for y in MARKED:
        label(ax, y, values[y], str(y), fontweight="bold")
        ax.plot(y, values[y], "o", ms=3.5, color=vz.INK, mec="white", mew=0.6, zorder=4)
    ax.set_ylabel("Heatwave days (JJA)", fontsize=FS)
    ax.set_ylim(0, max(rp.value.max(), gev.rl100.max()) * 1.12)


def panel_rp(bx, gev, rp, values, fits):
    yr = gev.year.values
    # how rare a 2003- / 2026-sized summer would be in every year's climate
    fixed = {ev: np.array([f.return_period(v, 0.0) for f in fits]) for ev, v in values.items()}
    for ev, curve in fixed.items():
        bx.plot(yr, curve, color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[ev], zorder=2)
        for name, (a, b) in LABEL_AT.items():
            sel = (yr >= a) & (yr <= b)
            v = float(np.exp(np.mean(np.log(curve[sel]))))   # geometric mean over the decade
            x = (a + b) / 2
            if name == str(ev):          # the event's own year: that is the main line's label
                continue
            bx.plot(x, v, "o", ms=2.6, color=vz.INK_SOFT, zorder=3)
            bx.annotate(fmt_rp(v) + (" (1950s)" if name == "1950s" else ""), (x, v),
                        xytext=(0, -4), textcoords="offset points", ha="center", va="top",
                        fontsize=FS - 1, color=vz.INK_SOFT, zorder=5, bbox=WHITE)

    finite = np.isfinite(rp.rp)
    bx.fill_between(rp.year, rp.rp_lo, rp.rp_hi, color=HOT, alpha=0.22, lw=0, zorder=1)
    bx.plot(rp.year[finite], rp.rp[finite], color=HOT, lw=1.2, zorder=3)
    bx.scatter(rp.year[finite], rp.rp[finite], s=7, color=HOT, zorder=4, lw=0)
    emp = rp[rp.k > 0]
    bx.scatter(emp.year, emp.rp_emp, s=12, facecolor="none", edgecolor=vz.INK_MUTED, lw=0.6,
               zorder=3)
    for _, r in rp[rp.k == 0].iterrows():
        bx.annotate("", (r.year, r.rp_emp_lo * 2.2), (r.year, r.rp_emp_lo),
                    arrowprops=dict(arrowstyle="-|>", color=vz.INK_MUTED, lw=0.6))
    for y in MARKED:
        r = rp[rp.year == y].iloc[0]
        txt = "∞" if r.above_upper_end else fmt_rp(r.rp)
        label(bx, y, r.rp if np.isfinite(r.rp) else r.rp_lo, f"{y}\n{txt}", dy=4,
              fontweight="bold", linespacing=1.1, bbox=WHITE)
    bx.set_yscale("log")
    top = max(np.nanmax(fixed[2026]), np.nanmax(rp.rp_hi[np.isfinite(rp.rp_hi)]))
    # headroom above the curves for the legend (log axis: a decade or so)
    bx.set_ylim(1, max(1e5, top * 30))
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000) if t <= top * 3]
    bx.set_yticks(ticks, [f"{t:,}" for t in ticks])
    bx.set_ylabel("Return period (years)", fontsize=FS)
    handles = [
        plt.Line2D([], [], color=HOT, lw=1.2, marker="o", ms=2.5, label="ERA5 summer, GEV"),
        plt.Rectangle((0, 0), 1, 1, color=HOT, alpha=0.22, lw=0, label="95% interval"),
        plt.Line2D([], [], ls="", marker="o", ms=3.5, mfc="none", mec=vz.INK_MUTED, mew=0.6,
                   label="count of 500 MPI-GE summers"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2026],
                   label=f"2026-sized summer ({values[2026]:.1f} d)"),
        plt.Line2D([], [], color=vz.INK_SOFT, lw=0.9, ls=FIXED_STYLE[2003],
                   label=f"2003-sized summer ({values[2003]:.1f} d)"),
    ]
    bx.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.9), ncol=2,
              frameon=False, fontsize=FS - 0.5, handlelength=2.2, columnspacing=1.2,
              borderaxespad=0.2)


def tag(ax, letter, text):
    ax.text(0.01, 0.97, f"({letter})", transform=ax.transAxes, va="top", fontsize=FS + 0.5,
            fontweight="bold")
    ax.text(0.045, 0.97, text, transform=ax.transAxes, va="top", fontsize=FS, color=vz.INK_SOFT)


def main():
    vz.use_style()
    gev = pd.read_csv(DIR / "mpige_window_gev.csv")
    rp = pd.read_csv(DIR / "era5_rp_in_mpige.csv")
    years, raw, anom = load_era5()
    values = {y: float(rp.loc[rp.year == y, "value"].iloc[0]) for y in MARKED}
    fits = window_gevs(gev)

    fig, (sx, ax, bx) = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True,
                                     gridspec_kw=dict(height_ratios=[0.85, 1, 1.25],
                                                      hspace=0.14))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.9, top=0.98, bottom=CAPTION_SPACE + 0.055)

    panel_series(sx, years, raw, anom)
    panel_map(ax, gev, rp, values)
    panel_rp(bx, gev, rp, values, fits)
    tag(sx, "a", f"ERA5, observed and with the quadratic trend removed")
    tag(ax, "b", "ERA5 against the MPI-GE climate of each year")
    tag(bx, "c", "how often MPI-GE reaches the observed value in that climate")

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
    fig.text(0.012, CAPTION_SPACE - 0.033, textwrap.fill(CAPTION, 122), ha="left", va="top",
             fontsize=FS, color=vz.INK, linespacing=1.45)
    out = p2.FIG_DIR / "fig18_era5_return_period_mpige.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
