"""
Figure 22 box difference -- the box-difference index MED − NAL and its two box
means, each observed summer against the MPI-GE summers of its climate.

(a) fig 22's former panel (c): ERA5 MED − NAL, Mediterranean (30–45°N, 0–20°E)
    minus North Atlantic (40–50°N, 60–20°W) JJA SST anomaly
    (data/sst_pattern_corr/box_difference.nc, 03_box_difference.py), on the
    50th / 90th / 99th percentiles and mean of the 500 MPI-GE summers in the 10
    years around each year.
(b) the same for the Mediterranean box mean alone and
(c) for the North Atlantic box mean alone ('<dataset>|med', '<dataset>|nal' of
    the same file). Neither the index nor the box means have any spatial mean
    removed: they are area means of the anomaly against each dataset's own
    1991–2020. The warming common to both boxes cancels in (a) and shows in
    (b) and (c), which share one y axis so that the two trends compare.

Marked and right-hand box plots as in fig 22.
"""
import importlib
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
import src.pattern_2026 as p2
import src.viz2026 as vz

OUT = p2.FIG_DIR / "fig22_boxdiff_mpige.png"
#: (variable suffix in box_difference.nc, y label) of each panel
ROWS = [("", f"{p2.DIPOLE_LABEL} SST anomaly (°C)"),
        ("|med", "Mediterranean SST\nanomaly (°C)"),
        ("|nal", "North Atlantic SST\nanomaly (°C)")]
CAPTION_TITLE = "The Mediterranean–North Atlantic SST contrast and its two boxes against the simulated summers of its climate"
CAPTION = (
    f"(a) ERA5 box-difference index {p2.DIPOLE_LABEL} (black): JJA SST anomaly averaged over "
    f"the Mediterranean ({p2.DIPOLE_BOXES['med']['detail']}) minus that over the North "
    f"Atlantic ({p2.DIPOLE_BOXES['nal']['detail']}), ocean cells on the common 1° grid, "
    "anomalies against 1991–2020 of the same dataset; blue: the same index in MPI-ESM1.2-LR: "
    "the 50th, 90th and 99th percentiles (solid) and the mean (dashed) of its 50 members in the "
    "10 years around each year (y−5 to y+4, 500 summers). (b) The same for the Mediterranean "
    "box mean alone and (c) for the North Atlantic box mean alone (same y axis): no spatial "
    "mean is removed, so the warming that cancels in (a) shows here. "
    "Marked: 2003, 2026 and every summer ranked at or above 90% of its window in heatwave days. "
    f"Right: the same quantities in ERA5, {fig22.BOX_YEARS[0]}–{p2.REF_YEAR}, and in the "
    "higher-resolution runs, all "
    f"historical + scenario summers {fig22.BOX_YEARS[0]}–{fig22.BOX_YEARS[1]} where available, "
    "pooled over members (ICON-EPOC: 1 run, 1991–2025; ICON-EERIE: 3, r1 1991–2030, r2 to "
    "2014, r3 to 2020; MPI-ESM1.2-ER: 3, 1991–2030; "
    "sample size in brackets; box 25th–75th percentile, line median, whiskers 5th–95th "
    "percentile, dots beyond).")


def figure(win_years, series, marked, boxes):
    """
    One panel per ROWS entry: the ERA5 `series` (DataArray over year) on the
    windows (`win_years`) of the MPI-GE values, with the column of box plots
    `boxes` (from fig22.box_samples) on the right.
    """
    vz.use_style()
    meth = fig18.METHODS["rank"]
    with xr.open_dataset(p2.dipole_file()) as ds:
        wins = [fig22.window_levels(ds["MPI-GE" + suffix].load(), win_years)
                for suffix, _ in ROWS]
    n = len(ROWS)

    text = textwrap.fill(CAPTION, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    top = 0.98
    height = (caption_in + fig22.BOTTOM_IN + fig22.ROW_IN * n) / top
    fig, axes = plt.subplots(
        n, 2, figsize=(fig18.FIG_SIZE[0], height), sharey="row",
        gridspec_kw=dict(height_ratios=[1.15] * n, hspace=0.14, wspace=0.15,
                         width_ratios=[1, 0.07 * len(fig22.BOX_KEYS) + 0.03]))
    mains, sides = axes[:, 0], axes[:, 1]
    for a in mains[1:]:
        a.sharex(mains[0])
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.985, top=top,
                        bottom=(caption_in + fig22.BOTTOM_IN) / height)

    for ax, x, win, (_, ylabel) in zip(mains, series, wins, ROWS):
        fig18.panel_series_levels(ax, win, x.year.values, x.values, dict(ylabel=ylabel),
                                  meth["levels"], [y for y in marked if y in x.year], "ERA5")
    for ax in mains[1:]:
        ax.get_legend().remove()
    for a, letter in zip(mains, "abc"):
        fig18.tag(a, letter)
    fig18.style_axes(mains)
    fig18.style_axes(sides)
    for a in mains[:-1]:
        a.tick_params(labelbottom=False)
    for i, (c, samples) in enumerate(zip(sides, boxes)):
        fig22.panel_boxes(c, samples, names=i == n - 1)

    # every row wide enough for its boxes; (b) and (c) on one common range
    lims = []
    for ax, samples in zip(mains, boxes):
        lo, hi = ax.get_ylim()
        pooled = np.concatenate([v for v in samples.values() if len(v)])
        pad = 0.05 * (hi - lo)
        lims.append((min(lo, pooled.min() - pad), max(hi, pooled.max() + pad)))
    lims[1] = lims[2] = (min(lims[1][0], lims[2][0]), max(lims[1][1], lims[2][1]))
    for ax, lim in zip(mains, lims):
        ax.set_ylim(lim)
    mains[-1].set_xlim(win_years[0] - 1, win_years[-1] + 2)

    fig.text(0.012, (caption_in - 0.04) / height, CAPTION_TITLE, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left",
             va="top", fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig.savefig(OUT, dpi=300, facecolor="white")
    fig.savefig(OUT.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {OUT}")
    for x, win, (_, ylabel) in zip(series, wins, ROWS):
        print(f" {ylabel.replace(chr(10), ' ')}:")
        for y in marked:
            if y in x.year:
                w = win.loc[win.year == y].iloc[0]
                print(f"  {y}: {float(x.sel(year=y)):+.2f}; MPI-GE window mean "
                      f"{w['mean']:+.2f}, 90th {w['q90']:+.2f}, 99th {w['q99']:+.2f}")
        a, b = win.iloc[0], win.iloc[-1]
        print(f"  MPI-GE window mean {int(a.year)} {a['mean']:+.2f} -> "
              f"{int(b.year)} {b['mean']:+.2f}")


def main():
    meth = fig18.METHODS["rank"]
    win_hwd = pd.read_csv(fig18.DIR / f"{meth['windows']}.csv")
    years, _, _ = fig18.load_series("hwd")
    marked = sorted(int(c.rsplit("_", 1)[1]) for c in win_hwd if c.startswith("pct_fixed_"))
    with xr.open_dataset(p2.dipole_file()) as ds:
        series = [ds["ERA5" + suffix].load().sel(year=slice(int(years[0]), p2.REF_YEAR))
                  .dropna("year") for suffix, _ in ROWS]
    boxes = fig22.box_samples([(p2.dipole_file(), suffix) for suffix, _ in ROWS], hwd=False)
    figure(win_hwd.year.values, series, marked, boxes)


if __name__ == "__main__":
    main()
