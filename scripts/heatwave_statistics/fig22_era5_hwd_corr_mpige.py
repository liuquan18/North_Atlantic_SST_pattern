"""
Figure 22 -- heatwave days and the SST pattern, each observed summer against
the MPI-GE summers of its climate (no fitted distribution).

(a) as fig 20 (a): ERA5 heatwave days on the 50th / 90th / 99th percentiles
    and mean of the 500 MPI-GE summers in the 10 years around each year
    (05_mpige_window_rank.py hwd).
(b) the same for the North Atlantic + Mediterranean SST pattern: the pattern
    correlation of each summer with ERA5 JJA 2026 (data/sst_pattern_corr,
    full box), ERA5 in black, the MPI-GE window percentiles and mean in blue.
    2026 itself is drawn too, at r = 1 by construction.

The box-difference index MED − NAL, formerly panel (c), and its two box means
are in fig22_boxdiff_mpige.py.

Marked: the summers marked in fig 20 (2003, 2026 and every summer ranked at
or above 90% in heatwave days).

A narrow column on the right of each panel: box plots of the same quantity in
the higher-resolution runs (BOX_KEYS), all their historical + scenario summers
of BOX_YEARS that exist, pooled over members, on the panel's y axis; and
ERA5 over the same years for comparison (2026 included, at r = 1 in the
correlation row).
MPI-ESM1.2-ER has no daily output, so its
heatwave-day box is empty.
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
w03 = importlib.import_module("03_mpige_window_rp")
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

QUANTILES = [50, 90, 99]
#: right-hand box plots: dataset key -> tick label, and the years pooled
BOX_KEYS = {"ERA5": "ERA5", "ICON-EPOC-hist": "ICON-\nEPOC", "EERIE": "ICON-\nEERIE", "MPI-ER": "MPI-\nER"}
BOX_YEARS = (1991, 2030)
OUT = p2.FIG_DIR / "fig22_era5_hwd_corr_mpige.png"
#: layout in inches: one row of panels (with its share of the gap), the strip
#: under the last row for its tick labels, and the caption's title + one line
ROW_IN, BOTTOM_IN, CAPTION_HEAD_IN, LINE_IN = 2.23, 0.42, 0.25, 0.115
CAPTION_WIDTH = 122
CAPTION_TITLE = "Heatwave days and the SST pattern of each observed summer against the simulated summers of its climate"
CAPTION = (
    f"(a) ERA5 June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL} (black), as in "
    "Fig. 20; blue: heatwave days of MPI-ESM1.2-LR: the 50th, 90th and 99th percentiles (solid) "
    "and the mean (dashed) of its 50 members in the 10 years around each year (y−5 to y+4, 500 "
    "summers). (b) The same for the SST pattern: centred pattern correlation of each JJA SST "
    "anomaly with that of ERA5 JJA 2026 over 20–60°N, 80°W–40°E (area-weighted box mean "
    "removed, anomalies against 1991–2020 of the same dataset); ERA5 in black, 2026 at r = 1 "
    "by construction. "
    "Marked: 2003, 2026 and every summer ranked at or above 90% of its window in heatwave days. "
    f"Right: the same quantities in ERA5, {BOX_YEARS[0]}–{p2.REF_YEAR}, and in the "
    "higher-resolution runs, all "
    f"historical + scenario summers {BOX_YEARS[0]}–{BOX_YEARS[1]} where available, pooled over "
    "members, keeping the summers that have both quantities (ICON-EPOC: 1 run, "
    "1991–2025; ICON-EERIE: 3, r1 1991–2030, r2 to 2014, r3 to 2020; "
    "MPI-ESM1.2-ER: 3, 1991–2030; "
    "sample size in brackets; box 25th–75th percentile, line median, whiskers 5th–95th "
    "percentile, dots beyond); "
    "MPI-ESM1.2-ER has no daily output for heatwave days.")


def window_levels(mpige, years):
    """50th / 90th / 99th percentiles and mean of the MPI-GE values in each year's window."""
    rows = []
    for year in years:
        x = mpige.sel(year=slice(year - w03.WIN_BEFORE, year + w03.WIN_AFTER)).values.ravel()
        x = x[np.isfinite(x)]
        rows.append(dict(year=year, n=len(x), mean=float(x.mean()),
                         **{f"q{q}": float(np.percentile(x, q)) for q in QUANTILES}))
    win = pd.DataFrame(rows)
    assert (win.n == 500).all(), "every window should hold 10 years x 50 members"
    return win


def _by_member_year(da):
    """A (member, year) Series of `da`; member 1 for a dataset without members."""
    dim = p2.member_dim(da)
    if dim is None:
        da = da.expand_dims(member=[1])
    else:
        da = da.rename({dim: "member"})
    return da.to_series()


def box_samples(sst, hwd=True):
    """
    Heatwave days (if `hwd`) and the SST quantities `sst`, a list of (file,
    variable suffix), e.g. (p2.dipole_file(), "|med") for the variables
    '<dataset>|med', of each BOX_KEYS dataset, BOX_YEARS, members pooled. Only
    the (member, year) summers that have all of them are kept, so every row's
    box of a dataset holds the same summers (EERIE r2 has heatwave days to 2020
    but SST only to 2014); a dataset without heatwave days (MPI-ER) keeps the
    summers that have the SST values. One {key: values} per quantity, heatwave
    days first.
    """
    a, b = BOX_YEARS
    t = pd.read_csv(fig18.DIR / "season_table.csv")
    t = t[t.year.between(a, b)]
    cols = []
    for path, suffix in sst:
        with xr.open_dataset(path) as ds:
            cols.append({k: _by_member_year(ds[k + suffix].sel(year=slice(a, b)).load())
                         for k in BOX_KEYS})
    tables = {k: pd.DataFrame({
        **({"hwd": t[t.dataset == k].set_index(["member", "year"]).hwd} if hwd else {}),
        **{i: c[k] for i, c in enumerate(cols)},
    }) for k in BOX_KEYS}
    for k, x in tables.items():
        need = list(range(len(cols))) if not hwd or x.hwd.isna().all() else list(x.columns)
        tables[k] = x.dropna(subset=need)
    return [{k: x[c].dropna().values for k, x in tables.items()}
            for c in next(iter(tables.values())).columns]


def panel_boxes(cx, samples, names=True):
    """
    Box plots of `samples` (key -> values; empty: no box) on the panel's y
    axis; if `names`, labelled with the run's name and its sample size.
    """
    for i, (key, x) in enumerate(samples.items()):
        color = vz.COLORS[key]
        if len(x) == 0:
            cx.text(i, 0.5, "n/a", transform=cx.get_xaxis_transform(), ha="center",
                    va="center", fontsize=fig18.FS - 1, color=vz.INK_SOFT)
            continue
        cx.boxplot(x, positions=[i], widths=0.6, whis=(5, 95), patch_artist=True,
                   boxprops=dict(facecolor=color, alpha=0.35, edgecolor=color, lw=0.8),
                   medianprops=dict(color=color, lw=1.2),
                   whiskerprops=dict(color=color, lw=0.8), capprops=dict(color=color, lw=0.8),
                   flierprops=dict(marker="o", ms=2, mfc=color, mec="none"))
    cx.set_xlim(-0.6, len(samples) - 0.4)
    labels = [f"{BOX_KEYS[k]}\n({len(x)})" if names else "" for k, x in samples.items()]
    cx.set_xticks(range(len(samples)), labels)
    cx.tick_params(axis="x", labelsize=fig18.FS - 1.5)
    cx.tick_params(axis="y", left=False, labelleft=False)
    cx.spines["left"].set_visible(False)


def figure(win_hwd, years, hwd, r, marked, name, title, caption, out, boxes=None):
    """
    The panels for one dataset `name`: its heatwave days `hwd` (`years`) on the
    MPI-GE windows `win_hwd` and its pattern correlation `r` (DataArray over
    year) on the same windows of the MPI-GE correlations. `boxes`, from
    box_samples(), adds the column of box plots on the right.
    """
    vz.use_style()
    spec, meth = fig18.SPECS["hwd"], fig18.METHODS["rank"]
    # (series, MPI-GE window levels, ylabel) of each SST row
    sst = []
    with xr.open_dataset(p2.corr_file(p2.PLOT_VARIANTS[0])) as ds:
        sst.append((r, window_levels(ds["MPI-GE"].load(), win_hwd.year.values),
                    "Pattern correlation with JJA 2026"))
    n = 1 + len(sst)

    text = textwrap.fill(caption, CAPTION_WIDTH)
    caption_in = CAPTION_HEAD_IN + LINE_IN * (text.count("\n") + 1)
    top = 0.98
    height = (caption_in + BOTTOM_IN + ROW_IN * n) / top
    caption_space = caption_in / height
    grid = dict(height_ratios=[1.15] * n, hspace=0.14)
    if boxes is None:
        fig, mains = plt.subplots(n, 1, figsize=(fig18.FIG_SIZE[0], height), sharex=True,
                                  gridspec_kw=grid)
        sides = [None] * n
        right = 0.9
    else:
        # the main panels' level labels sit in the gap before the box column
        fig, axes = plt.subplots(n, 2, figsize=(fig18.FIG_SIZE[0], height), sharey="row",
                                 gridspec_kw=dict(grid, width_ratios=[1, 0.07 * len(BOX_KEYS) + 0.03], wspace=0.15))
        mains, sides = axes[:, 0], axes[:, 1]
        for a, c in zip(mains[1:], sides[1:]):
            a.sharex(mains[0])
        right = 0.985
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=right, top=top,
                        bottom=caption_space + BOTTOM_IN / height)

    ax = mains[0]
    fig18.panel_series_levels(ax, win_hwd, years, hwd, spec, meth["levels"], marked, name)
    for bx, (x, win, ylabel) in zip(mains[1:], sst):
        fig18.panel_series_levels(bx, win, x.year.values, x.values, dict(ylabel=ylabel),
                                  meth["levels"], [y for y in marked if y in x.year], name)
        bx.get_legend().remove()
    for a, letter in zip(mains, "abc"):
        fig18.tag(a, letter)
    fig18.style_axes(mains)
    for a in mains[:-1]:
        a.tick_params(labelbottom=False)
    if boxes is not None:
        fig18.style_axes(sides)
        for i, (c, samples) in enumerate(zip(sides, boxes)):
            panel_boxes(c, samples, names=i == n - 1)
    mains[-1].set_xlim(win_hwd.year.iloc[0] - 1, win_hwd.year.iloc[-1] + 2)

    fig.text(0.012, (caption_in - 0.04) / height, title, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - CAPTION_HEAD_IN) / height, text, ha="left",
             va="top", fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")
    for x, win, ylabel in sst:
        print(f" {ylabel}:")
        for y in marked:
            if y in x.year:
                w = win.loc[win.year == y].iloc[0]
                print(f"  {y}: {float(x.sel(year=y)):+.2f}; MPI-GE window mean "
                      f"{w['mean']:+.2f}, 90th {w['q90']:+.2f}, 99th {w['q99']:+.2f}")


def main():
    meth = fig18.METHODS["rank"]
    win_hwd = pd.read_csv(fig18.DIR / f"{meth['windows']}.csv")
    years, raw, _ = fig18.load_series("hwd")
    marked = sorted(int(c.rsplit("_", 1)[1]) for c in win_hwd if c.startswith("pct_fixed_"))
    with xr.open_dataset(p2.corr_file(p2.PLOT_VARIANTS[0])) as ds:
        era5_r = ds["ERA5"].load()
    era5_r = era5_r.sel(year=slice(int(years[0]), p2.REF_YEAR)).dropna("year")
    figure(win_hwd, years, raw, era5_r, marked, "ERA5", CAPTION_TITLE, CAPTION, OUT,
           boxes=box_samples([(p2.corr_file(p2.PLOT_VARIANTS[0]), "")]))


if __name__ == "__main__":
    main()
