"""
Figure 23 -- as figure 22, with every higher-resolution historical + scenario
run in place of ERA5, one line per realization, over YEARS (1991-2030).

(a) heatwave days of each run on the 50th / 90th / 99th percentiles and mean
    of the 500 MPI-GE summers in the 10 years around each year (computed here
    from season_table.csv, since 05_mpige_window_rank.py stops at 2026).
(b) the pattern correlation of each run's JJA SST anomaly with ERA5 JJA 2026,
    on the same windows of the MPI-GE pattern correlations.

Runs (RUNS): ICON-EPOC (hist + SSP5-8.5, 1 run), ICON-EERIE (hist + SSP2-4.5,
r1-r3) and MPI-ESM1.2-ER (hist + SSP5-8.5, r2 and r3; r1 has no daily output
before 2014, so no heatwave reference). Colour = model, line style = member.
A run that has no heatwave days yet (MPI-ER while its daily Tmax is being
recalled from tape) is drawn in (b) only. No summers are marked.
"""
import importlib
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

YEARS = (1991, 2030)
#: dataset key -> (legend name, members drawn)
RUNS = {"ICON-EPOC-hist": ("ICON-EPOC", [1]),
        "EERIE": ("ICON-EERIE", [1, 2, 3]),
        "MPI-ER": ("MPI-ER", [2, 3])}
MEMBER_STYLE = {1: "-", 2: (0, (4, 1.5)), 3: (0, (1, 1.2))}
OUT = p2.FIG_DIR / "fig23_epoc_hwd_corr_mpige.png"
CAPTION_TITLE = "Heatwave days and the SST pattern of the higher-resolution runs against the MPI-GE summers of their climate"
CAPTION = (
    f"As Fig. 22, {YEARS[0]}–{YEARS[1]}, with the higher-resolution historical + scenario runs "
    "in place of ERA5, one line per realisation: ICON-EPOC (SSP5-8.5 from 2015, 10 km "
    "atmosphere / 5 km ocean, 1 run, to 2025), ICON-EERIE (SSP2-4.5 from 2015, 10 km / 5 km, "
    "r1 to 2050, r2 SST to 2014 and heatwave days to 2020, r3 to 2020) and MPI-ESM1.2-ER "
    "(SSP5-8.5 from 2015, ~100 km atmosphere / ~10 km ocean, r2 and r3); line style = member. "
    f"(a) June–August heatwave days averaged over land in {hw.MEAN_BOX_LABEL}, each run against "
    "its own 1991–2020 threshold; blue: heatwave days of MPI-ESM1.2-LR: the 50th, 90th and "
    "99th percentiles (solid) and the mean (dashed) of its 50 members in the 10 years around "
    "each year (y−5 to y+4, 500 summers). (b) The same for the SST pattern: centred pattern "
    "correlation of each JJA SST anomaly with that of ERA5 JJA 2026 over 20–60°N, 80°W–40°E "
    "(area-weighted box mean removed, anomalies against 1991–2020 of the same dataset).")


def run_series(da, member):
    """One member's (year, value) Series of `da` over YEARS, NaN years dropped."""
    dim = p2.member_dim(da)
    if dim is not None:
        da = da.sel({dim: member})
    return da.sel(year=slice(*YEARS)).to_series().dropna()


def load():
    """
    MPI-GE heatwave days as a (member, year) DataArray, and per panel a dict
    (key, member) -> Series of every run in RUNS that has data there.
    """
    t = pd.read_csv(fig18.DIR / "season_table.csv")
    mpige = t[t.dataset == "MPI-GE"].pivot(index="member", columns="year", values="hwd")
    mpige = xr.DataArray(mpige.values, coords=dict(member=mpige.index, year=mpige.columns),
                         dims=("member", "year"))
    hwd, corr = {}, {}
    for key, (_, members) in RUNS.items():
        for m in members:
            x = t[(t.dataset == key) & (t.member == m) & t.year.between(*YEARS)]
            if len(x):
                hwd[key, m] = x.set_index("year").hwd.sort_index()
    with xr.open_dataset(p2.corr_file(p2.PLOT_VARIANTS[0])) as ds:
        mpige_r = ds["MPI-GE"].load()
        for key, (_, members) in RUNS.items():
            for m in members:
                corr[key, m] = run_series(ds[key].load(), m)
    return mpige, mpige_r, hwd, corr


def panel(ax, win, runs, ylabel):
    """MPI-GE window levels `win` with one line per run in `runs`."""
    values = pd.concat(runs.values())
    lo, hi = fig18.padded(min(values.min(), win.q50.min()), max(values.max(), win.q99.max()))
    fig18.draw_levels(ax, win, fig18.METHODS["rank"]["levels"], styles=("-", "-", "-"),
                      extra=(("mean", "mean", "--"),), span=hi - lo)
    for (key, m), x in runs.items():
        ax.plot(x.index, x.values, color=vz.COLORS[key], ls=MEMBER_STYLE[m], lw=1.1, zorder=3)
    ax.set_ylim(lo, hi)
    ax.set_ylabel(ylabel, fontsize=fig18.FS)


def figure(rows, title, caption, out, share=()):
    """
    One panel per (MPI-GE window levels, {(key, member): Series}, ylabel) in
    `rows`, a shared legend on the first; the panels whose indices are in
    `share` get one common y range.
    """
    vz.use_style()
    text = textwrap.fill(caption, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    top, n = 0.98, len(rows)
    height = (caption_in + fig22.BOTTOM_IN + fig22.ROW_IN * n) / top
    fig, axes = plt.subplots(n, 1, figsize=(fig18.FIG_SIZE[0], height), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.15] * n, hspace=0.14))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.9, top=top,
                        bottom=(caption_in + fig22.BOTTOM_IN) / height)

    for ax, (win, runs, ylabel) in zip(axes, rows):
        panel(ax, win, runs, ylabel)
    if share:
        lims = [axes[i].get_ylim() for i in share]
        for i in share:
            axes[i].set_ylim(min(l[0] for l in lims), max(l[1] for l in lims))
    handles = [plt.Line2D([], [], color=vz.COLORS[key], ls=MEMBER_STYLE[m], lw=1.1,
                          label=name if len(members) == 1 else f"{name} r{m}")
               for key, (name, members) in RUNS.items() for m in members]
    handles += [plt.Rectangle((0, 0), 1, 1, color=fig18.BAND, alpha=0.3, lw=0,
                              label="MPI-GE percentiles"),
                plt.Line2D([], [], color=fig18.BAND, lw=0.8, ls="--", label="MPI-GE mean")]
    axes[0].legend(handles=handles, loc="upper left", bbox_to_anchor=(0.05, 1.03),
                   frameon=False, fontsize=fig18.FS - 0.5, handlelength=2.2, ncol=4)
    lo, hi = axes[0].get_ylim()
    axes[0].set_ylim(lo, hi + 0.3 * (hi - lo))    # headroom for the two-row legend
    for a, letter in zip(axes, "abcdef"):
        fig18.tag(a, letter)
    fig18.style_axes(axes)
    for a in axes[:-1]:
        a.tick_params(labelbottom=False)
    axes[-1].set_xlim(YEARS[0] - 0.5, YEARS[1] + 0.5)

    fig.text(0.012, (caption_in - 0.04) / height, title, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left",
             va="top", fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")
    for (_, runs, ylabel) in rows:
        print(f" {ylabel.replace(chr(10), ' ')}:")
        for (key, m), x in runs.items():
            print(f"  {key} r{m}: {x.index.min()}-{x.index.max()} ({len(x)} summers), "
                  f"max {x.max():+.2f} in {x.idxmax()}")


def main():
    mpige, mpige_r, hwd, corr = load()
    years = list(range(YEARS[0], YEARS[1] + 1))
    missing = [f"{name} r{m}" for k, (name, ms) in RUNS.items() for m in ms if (k, m) not in hwd]
    caption = CAPTION + (f" No heatwave days yet for {', '.join(missing)}." if missing else "")
    figure([(fig22.window_levels(mpige, years), hwd, fig18.SPECS["hwd"]["ylabel"]),
            (fig22.window_levels(mpige_r, years), corr, "Pattern correlation with JJA 2026")],
           CAPTION_TITLE, caption, OUT)


if __name__ == "__main__":
    main()
