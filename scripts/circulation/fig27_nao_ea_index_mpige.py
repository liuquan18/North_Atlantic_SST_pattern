"""
Figure 27 -- the summer NAO and EA indices of ERA5 against the MPI-GE summers
of their climate, as fig 22 does for heatwave days: each ERA5 summer
(1940-2026, black) on the 1st / 10th / 50th / 90th / 99th percentiles and mean
of the 500 MPI-GE member-summers in the 10 years around it (y-5 to y+4).
Indices from 01_nao_ea_eof.py: standardised PCs of each dataset's own EOFs,
ERA5 detrended, MPI-GE with the ensemble mean removed, so both are internal
variability in units of their own standard deviation. Two-sided indices, so
the lower percentiles are drawn too. Marked: the summers marked in fig 22.
The method is written out in the framed box (fig26_nao_ea_patterns.METHOD).
"""
import importlib
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "heatwave_statistics"))
eof = importlib.import_module("01_nao_ea_eof")
fig26 = importlib.import_module("fig26_nao_ea_patterns")
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
import src.pattern_2026 as p2
import src.viz2026 as vz

OUT = p2.FIG_DIR / "fig27_nao_ea_index_mpige.png"
YEARS = eof.FIT_YEARS
WIN = (5, 4)
#: (column, label) from the lowest to the highest
LEVELS = [("q1", "1st pct"), ("q10", "10th pct"), ("q50", "50th pct"),
          ("q90", "90th pct"), ("q99", "99th pct")]
CAPTION_TITLE = "The summer NAO and EA of each observed summer against the simulated summers of its climate"
CAPTION = (
    "(a) ERA5 JJA NAO index (black), 1940–2026; blue: the same index of MPI-ESM1.2-LR: the 1st, "
    "10th, 50th, 90th and 99th percentiles (solid) and the mean (dashed) of its 50 members in the "
    "10 years around each year (y−5 to y+4, 500 summers). (b) The same for the EA index. Both "
    "indices are standardised PCs of the dataset's own EOFs (see box), so they are internal "
    "variability in units of their own standard deviation; ERA5 detrended, MPI-GE with the "
    "ensemble mean removed. Marked: 2003, 2026 and every summer ranked at or above 90% of its "
    "window in heatwave days (Fig. 22).")


def window_levels(ge, years):
    """Percentiles and mean of the MPI-GE index over each year's 10-year window."""
    rows = []
    for y in years:
        x = ge.sel(year=slice(y - WIN[0], y + WIN[1])).values.ravel()
        x = x[np.isfinite(x)]
        assert len(x) == 500, f"{y}: {len(x)} summers"
        rows.append(dict(year=y, mean=x.mean(),
                         **{c: np.percentile(x, float(c[1:])) for c, _ in LEVELS}))
    return pd.DataFrame(rows)


def draw_levels(ax, win, span):
    """Bands 10th-90th (dark) and 1st-10th / 90th-99th (light), lines, end labels."""
    yr = win.year.values
    ax.fill_between(yr, win.q10, win.q90, color=fig18.BAND, alpha=0.30, lw=0, zorder=1)
    for lo, hi in (("q1", "q10"), ("q90", "q99")):
        ax.fill_between(yr, win[lo], win[hi], color=fig18.BAND, alpha=0.13, lw=0, zorder=1)
    lines = [(c, name, "-") for c, name in LEVELS] + [("mean", "mean", "--")]
    for c, _, ls in lines:
        ax.plot(yr, win[c], color=fig18.BAND, lw=0.8, ls=ls, zorder=2)
    gap = 0.06 * span
    y_prev = -np.inf
    for y, name in sorted((float(win[c].iloc[-1]), name) for c, name, _ in lines
                          if c != "q50"):   # the median sits on the mean
        y_prev = max(y, y_prev + gap)
        ax.text(yr[-1] + 1.2, y_prev, name, va="center", fontsize=fig18.FS - 0.5,
                color=fig18.BAND)


def panel(ax, win, x, marked, ylabel):
    lo, hi = fig18.padded(min(x.min(), win.q1.min()), max(x.max(), win.q99.max()))
    draw_levels(ax, win, hi - lo)
    ax.axhline(0, color=vz.INK_MUTED, lw=0.5, zorder=0)
    ax.plot(x.year, x.values, color=fig18.OBS, lw=1.2, zorder=3)
    for y in marked:
        v = float(x.sel(year=y))
        ax.plot(y, v, "o", ms=3.5, color=fig18.OBS, mec="white", mew=0.6, zorder=4)
        fig18.label(ax, y, v, str(y), fontweight="bold")
    ax.set_ylim(lo, hi)
    ax.set_ylabel(ylabel, fontsize=fig18.FS)


def rank(win_ge, ge, y, v):
    """Percentage of the window's MPI-GE summers below v."""
    x = ge.sel(year=slice(y - WIN[0], y + WIN[1])).values.ravel()
    return 100 * np.mean(x[np.isfinite(x)] < v)


def main():
    vz.use_style()
    meth = fig18.METHODS["rank"]
    win_hwd = pd.read_csv(fig18.DIR / f"{meth['windows']}.csv")
    marked = sorted(int(c.rsplit("_", 1)[1]) for c in win_hwd if c.startswith("pct_fixed_"))
    years = list(range(YEARS[0], YEARS[1] + 1))
    rows = []
    with xr.open_dataset(eof.OUT) as ds:
        for mode in ("NAO", "EA"):
            era5 = ds[f"ERA5|{mode}"].sel(year=years).load()
            ge = ds[f"MPI-GE|{mode}"].load()
            rows.append((mode, era5, ge, window_levels(ge, years)))

    text = textwrap.fill(CAPTION, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    box_in = 0.145 * fig26.method_lines(126) + 0.45
    n = len(rows)
    top = 0.985
    height = (box_in + caption_in + fig22.BOTTOM_IN + fig22.ROW_IN * n) / top
    fig, axes = plt.subplots(n, 1, figsize=(fig18.FIG_SIZE[0], height), sharex=True,
                             gridspec_kw=dict(hspace=0.14))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.88, top=top,
                        bottom=(box_in + caption_in + fig22.BOTTOM_IN) / height)
    for ax, (mode, era5, _, win) in zip(axes, rows):
        panel(ax, win, era5, marked, f"JJA {mode} index (std)")
    axes[0].legend(handles=[
        plt.Line2D([], [], color=fig18.OBS, lw=1.2, label="ERA5"),
        plt.Rectangle((0, 0), 1, 1, color=fig18.BAND, alpha=0.3, lw=0, label="MPI-GE percentiles"),
        plt.Line2D([], [], color=fig18.BAND, lw=0.8, ls="--", label="MPI-GE mean")],
        loc="upper left", bbox_to_anchor=(0.05, 1.03), frameon=False, fontsize=fig18.FS,
        handlelength=1.6, ncol=3)
    for a, letter in zip(axes, "ab"):
        fig18.tag(a, letter)
    fig18.style_axes(axes)
    axes[-1].set_xlim(YEARS[0] - 1, YEARS[1] + 2)

    y0 = (box_in + caption_in) / height
    fig.text(0.012, y0 - 0.04 / height, CAPTION_TITLE, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, y0 - fig22.CAPTION_HEAD_IN / height, text, ha="left", va="top",
             fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig26.method_box(fig, (box_in - 0.35) / height, box_in, height, width=126)
    fig.savefig(OUT, dpi=300, facecolor="white")
    fig.savefig(OUT.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {OUT}")
    for mode, era5, ge, win in rows:
        print(f" {mode}:")
        for y in marked:
            v = float(era5.sel(year=y))
            w = win.loc[win.year == y].iloc[0]
            print(f"  {y}: {v:+.2f} std, rank {rank(win, ge, y, v):.0f}% of its MPI-GE window "
                  f"(10th {w.q10:+.2f}, 90th {w.q90:+.2f})")
        print(f"  ERA5 std {float(era5.std()):.2f}; MPI-GE window 10-90th width "
              f"{(win.q90 - win.q10).mean():.2f} (mean over years)")


if __name__ == "__main__":
    main()
