"""
Figure 24 -- co-variability of the Mediterranean (MED) and mid-latitude North
Atlantic (NAL) JJA SST box means of fig22_boxdiff_mpige.py
(data/sst_pattern_corr/box_difference.nc, '<dataset>|med', '<dataset>|nal'),
over PERIOD (each dataset's years within it; ERA5 also over its whole record
from 1940), after removing a linear trend from both series so that the warming
common to both boxes does not count as co-variability.

(a) correlation r of MED and NAL;
(b) the standard deviations sigma_MED and sigma_NAL: r leaves out how large
    the variations are, the covariance is r * sigma_MED * sigma_NAL;
(c) the variance of the box-difference index split into its parts,
    var(MED - NAL) = sigma_MED^2 + sigma_NAL^2 - 2 cov(MED, NAL):
    stacked bar sigma_MED^2 + sigma_NAL^2, arrow -2 cov to the diamond at
    var(MED - NAL). A positive covariance narrows the index, a negative one
    widens it.

Bars: ERA5 and the higher-resolution runs of fig 23 (RUNS), one per
realization. MPI-GE: box plots of its 50 members in (a) and (b); in (c) the
member means (means add up, medians do not), the diamond's whisker the
members' 5th-95th percentile of var(MED - NAL).
"""
import importlib
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
fig23 = importlib.import_module("fig23_epoc_hwd_corr_mpige")
import src.pattern_2026 as p2
import src.viz2026 as vz

PERIOD = fig23.YEARS
#: alpha of the MED and NAL parts of a bar in (b) and (c)
ALPHA = {"med": 0.85, "nal": 0.35}
OUT = p2.FIG_DIR / "fig24_med_nal_correlation.png"
CAPTION_TITLE = "How closely, and how strongly, do Mediterranean and North Atlantic summer SSTs vary together?"
CAPTION = (
    "JJA SST anomaly averaged over the Mediterranean (MED, "
    f"{p2.DIPOLE_BOXES['med']['detail']}) and over the North Atlantic (NAL, "
    f"{p2.DIPOLE_BOXES['nal']['detail']}), ocean cells on the common 1° grid, anomalies "
    f"against 1991–2020 of the same dataset, over {PERIOD[0]}–{PERIOD[1]} or the part of it "
    "each dataset covers (years below the names), linear trend removed from both series; "
    "ERA5 also over its whole record 1940–2026 (grey). (a) Temporal correlation of MED and "
    "NAL. (b) Standard deviation of MED (dark) and NAL (light); the covariance is r·σ_MED·σ_NAL. "
    "(c) Variance of the box-difference index MED − NAL (diamond) = σ²_MED (dark) + σ²_NAL "
    "(light) − 2 cov(MED, NAL) (arrow; downward when the boxes vary together). "
    "Higher-resolution historical + scenario runs: ICON-EPOC (SSP5-8.5), ICON-EERIE r1–r3 "
    "(SSP2-4.5) and MPI-ESM1.2-ER r2 and r3 (SSP5-8.5). MPI-ESM1.2-LR (historical + "
    "SSP2-4.5): its 50 members, in (a) and (b) as box plots (box 25th–75th percentile, line "
    "median, whiskers 5th–95th percentile, dots beyond), in (c) the member means, whisker "
    "5th–95th percentile of the members' variance of MED − NAL.")


def stats(x, y, period=PERIOD):
    """
    r, sigma_MED, sigma_NAL, cov and var(MED - NAL) of two year Series over
    their common years in `period`, each linearly detrended, and the years used.
    """
    x, y = x.loc[period[0]:period[1]].dropna(), y.loc[period[0]:period[1]].dropna()
    years = x.index.intersection(y.index)
    t = years.values.astype(float)
    x = x[years].values - np.polyval(np.polyfit(t, x[years].values, 1), t)
    y = y[years].values - np.polyval(np.polyfit(t, y[years].values, 1), t)
    c = np.cov(x, y)          # ddof=1 throughout, so var(x - y) adds up exactly
    return dict(r=c[0, 1] / np.sqrt(c[0, 0] * c[1, 1]), sm=np.sqrt(c[0, 0]),
                sn=np.sqrt(c[1, 1]), cov=c[0, 1], vd=np.var(x - y, ddof=1),
                span=(int(years.min()), int(years.max())))


def series(ds, key, member=None):
    """(MED, NAL) year Series of one dataset / member."""
    out = []
    for suffix in ("|med", "|nal"):
        da = ds[key + suffix]
        dim = p2.member_dim(da)
        if dim is not None:
            da = da.sel({dim: member})
        out.append(da.to_series())
    return out


def boxes(ax, values, positions, alphas, w):
    """MPI-GE box plots of `values` (one array per position)."""
    for v, x, a in zip(values, positions, alphas):
        ax.boxplot(v, positions=[x], widths=w, whis=(5, 95), patch_artist=True,
                   boxprops=dict(facecolor=fig18.BAND, alpha=a, edgecolor=fig18.BAND, lw=0.8),
                   medianprops=dict(color=fig18.BAND, lw=1.2),
                   whiskerprops=dict(color=fig18.BAND, lw=0.8),
                   capprops=dict(color=fig18.BAND, lw=0.8),
                   flierprops=dict(marker="o", ms=2, mfc=fig18.BAND, mec="none"))


def main():
    vz.use_style()
    with xr.open_dataset(p2.dipole_file()) as ds:
        ds = ds.load()
    members = ds[p2.member_dim(ds["MPI-GE|med"])].values
    ge = [stats(*series(ds, "MPI-GE", m)) for m in members]
    ge = {k: np.array([s[k] for s in ge]) for k in ("r", "sm", "sn", "cov", "vd")}
    # (label, colour, statistics) of ERA5 and every run
    era5 = series(ds, "ERA5")
    runs = [("ERA5", vz.INK, stats(*era5)),
            ("ERA5", vz.INK_MUTED, stats(*era5, period=(None, None)))]
    for key, (name, ms) in fig23.RUNS.items():
        for m in ms:
            runs.append((name if len(ms) == 1 else f"{name} r{m}", vz.COLORS[key],
                         stats(*series(ds, key, m))))

    text = textwrap.fill(CAPTION, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    top, rows_in, ticks_in = 0.98, 3 * 2.0 + 2 * 0.25, 0.6
    height = (caption_in + ticks_in + rows_in) / top
    fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(fig18.FIG_SIZE[0], height), sharex=True,
                                     gridspec_kw=dict(hspace=0.15))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.1, right=0.985, top=top, bottom=(caption_in + ticks_in) / height)
    pos = np.arange(1, len(runs) + 1)
    w, dx = 0.55, 0.15

    # (a) correlation
    boxes(ax, [ge["r"]], [0], [0.5], w)
    for i, (_, color, s) in zip(pos, runs):
        ax.bar(i, s["r"], w, color=color, alpha=0.85, lw=0)
    ax.set_ylim(-0.6, 0.8)
    ax.set_ylabel("Correlation r", fontsize=fig18.FS)

    # (b) standard deviations
    boxes(bx, [ge["sm"], ge["sn"]], [-dx, dx], [ALPHA["med"], ALPHA["nal"]], 0.28)
    for i, (_, color, s) in zip(pos, runs):
        bx.bar(i - dx, s["sm"], 0.3, color=color, alpha=ALPHA["med"], lw=0)
        bx.bar(i + dx, s["sn"], 0.3, color=color, alpha=ALPHA["nal"], lw=0)
    bx.set_ylim(0, None)
    bx.set_ylabel("Standard deviation (°C)", fontsize=fig18.FS)
    bx.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=vz.INK_SOFT, alpha=ALPHA[k], lw=0,
                                     label=label)
                       for k, label in (("med", "MED"), ("nal", "NAL"))],
              loc="upper right", frameon=False, fontsize=fig18.FS - 0.5, ncol=2)

    # (c) var(MED - NAL) = sigma_MED^2 + sigma_NAL^2 - 2 cov
    mean = {k: float(v.mean()) for k, v in ge.items() if k != "r"}
    parts = [("MPI-GE", fig18.BAND, dict(sm2=float((ge["sm"] ** 2).mean()),
                                         sn2=float((ge["sn"] ** 2).mean()),
                                         cov=mean["cov"], vd=mean["vd"]))]
    parts += [(label, color, dict(sm2=s["sm"] ** 2, sn2=s["sn"] ** 2, cov=s["cov"], vd=s["vd"]))
              for label, color, s in runs]
    for i, (_, color, p) in zip(np.arange(len(parts)), parts):
        cx.bar(i, p["sm2"], w, color=color, alpha=ALPHA["med"], lw=0)
        cx.bar(i, p["sn2"], w, bottom=p["sm2"], color=color, alpha=ALPHA["nal"], lw=0)
        total = p["sm2"] + p["sn2"]
        if abs(p["vd"] - total) > 1e-3:
            cx.annotate("", (i + 0.2, p["vd"]), (i + 0.2, total),
                        arrowprops=dict(arrowstyle="-|>", color=vz.INK, lw=0.9,
                                        mutation_scale=7, shrinkA=0, shrinkB=0))
        cx.plot(i + 0.2, p["vd"], "D", ms=4.5, color=vz.INK, mec="white", mew=0.5, zorder=4)
    lo, hi = np.percentile(ge["vd"], [5, 95])
    cx.errorbar(0.2, mean["vd"], yerr=[[mean["vd"] - lo], [hi - mean["vd"]]], fmt="none",
                ecolor=fig18.BAND, elinewidth=0.9, capsize=2.5, zorder=3)
    cx.set_ylim(0, None)
    cx.set_ylabel("Variance (°C²)", fontsize=fig18.FS)
    cx.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=vz.INK_SOFT, alpha=ALPHA["med"], lw=0, label="σ²_MED"),
        plt.Rectangle((0, 0), 1, 1, color=vz.INK_SOFT, alpha=ALPHA["nal"], lw=0, label="σ²_NAL"),
        plt.Line2D([], [], color=vz.INK, lw=0.9, marker=r"$\downarrow$", ms=7, ls="none",
                   label="−2 cov"),
        plt.Line2D([], [], color=vz.INK, marker="D", ms=4.5, ls="none", mec="white",
                   label="var(MED − NAL)")],
        loc="upper right", frameon=False, fontsize=fig18.FS - 0.5, ncol=4)

    for a, letter in zip((ax, bx, cx), "abc"):
        a.axhline(0, color=vz.INK_SOFT, lw=0.6, zorder=0)
        fig18.tag(a, letter)
    fig18.style_axes((ax, bx, cx))
    ticks = [f"MPI-GE\n50 × {PERIOD[0]}–\n{PERIOD[1]}"]
    ticks += [f"{label.replace(' ', chr(10), 1)}\n{s['span'][0]}–\n{s['span'][1]}"
              for label, _, s in runs]
    cx.set_xticks(range(len(ticks)), ticks)
    for a in (ax, bx, cx):
        a.tick_params(axis="x", labelsize=fig18.FS - 1, length=0)
    cx.set_xlim(-0.6, len(ticks) - 0.4)
    # headroom for the legends
    for a in (bx, cx):
        a.set_ylim(0, a.get_ylim()[1] * 1.18)

    fig.text(0.012, (caption_in - 0.04) / height, CAPTION_TITLE, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left",
             va="top", fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig.savefig(OUT, dpi=300, facecolor="white")
    fig.savefig(OUT.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {OUT}")
    print("linearly detrended:        r   sd_MED  sd_NAL   cov    var(MED-NAL)")
    print(f"  MPI-GE median {'':10s} {np.median(ge['r']):+.2f}  {np.median(ge['sm']):.2f}"
          f"    {np.median(ge['sn']):.2f}   {np.median(ge['cov']):+.3f}  {np.median(ge['vd']):.3f}")
    for label, _, s in runs:
        print(f"  {label:14s} {s['span'][0]}-{s['span'][1]} {s['r']:+.2f}  {s['sm']:.2f}"
              f"    {s['sn']:.2f}   {s['cov']:+.3f}  {s['vd']:.3f}")


if __name__ == "__main__":
    main()
