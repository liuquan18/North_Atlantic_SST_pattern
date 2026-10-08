"""
Figure 8 -- does a 2026-like SST pattern come with more European heatwave days?

One scatter panel per dataset with daily Tmax (ERA5, MPI-GE, km-scale ICON
(EPOC) and EERIE; MPI-ESM1.2-ER has none). Every season (and member) in
PERIOD = 1990-2026 is one point -- the window the observations cover, so
the models are compared over the same climate and not scored on their
late-century seasons. A record that ends earlier (EPOC, until its JJA 2026
is processed; EERIE r2/r3) contributes the seasons it has. ERA5 2026 is a
point like any other and enters the fit and its r (at x = 1 by construction):

  x   pattern correlation of its JJA SST pattern with ERA5 JJA 2026
      (box mean removed, 1 deg common ocean mask; results/corr_spatial.nc)
  y   land-mean JJA heatwave days over southern Europe, 35-60N, 10W-30E
      (src.heatwave.MEAN_BOX; Xu et al. 2026 definition, 1991-2020
      reference of the same dataset)

Every run writes two versions: as is, and with the forced trend removed from
the heatwave days (file suffix _detrended). The pattern correlation is left
as it is: the box mean is removed before correlating, so the uniform warming
is already out of it. detrend(): MPI-GE minus its 50-member mean, year by
year; every other record (and EERIE member) minus its own least-squares line
over all its seasons in the panel.

With an argument `atl` or `med`, x is the pattern correlation scored over
that half of the SST box alone (results/corr_spatial_<half>.nc; file suffix
_atl / _med), and the ringed model seasons are re-selected on it.

Points are coloured by the season's mean heatwave intensity: the mean Tmax
anomaly over all European land heatwave days of that season (area-weighted
cumulative heat / area-weighted heatwave days -- the "intensity" row of
fig1_heatwaves.py, pooled over the map). Seasons without any heatwave day
have no intensity and are drawn hollow. One colour scale serves all panels;
MPI-GE's 1,800 seasons are drawn smaller. In each model panel the ringed
point is the season shown in fig1_heatwave_analogues.py (select_season:
high on both axes, or for MPI-GE furthest along the fitted line); in ERA5, JJA 2015, 2003 and 2026 are ringed -- 2026 at r = 1 by
construction, and included in the fit. Each panel carries the least-squares line of y on x with its Pearson r, two-sided p-value and slope
(heatwave days per +0.1 of pattern correlation). The
p-value treats every season as independent -- reasonable across members and
years of JJA, but it ignores any year-to-year persistence.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap, Normalize

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, land_mean_series, scored_over

KEYS = ["ERA5", "MPI-GE", "ICON-EPOC-hist", "EERIE"]
#: observed seasons marked in the ERA5 panel; 2026 is the reference season
#: (r = 1 by construction) and is part of the fit like every other season
ERA5_MARKED = [p2.REF_YEAR, 2015, 2003]
#: label offsets (points, ha) -- 2003 and 2015 sit almost on top of each other
#: how each model's season is selected (default "both"):
#:   "both"       top_on_both -- highest lower percentile rank of r and heatwave days
#:   "along_fit"  along_fit   -- highest r among seasons close to the fitted line;
#:                for MPI-GE, where "both" lands on the long hot tail of 1,800
#:                seasons rather than on the relationship the line describes
#:   "max_r_hot"  max_r_hot   -- highest r among seasons with above-median heatwave
#:                days; for EERIE, where "both" passes over its far-right,
#:                hot season (1997 r1: highest r of all, 82nd percentile in
#:                heatwave days) for a much lower-r one that is slightly hotter
SELECTION = {"MPI-GE": "along_fit", "EERIE": "max_r_hot"}
SELECTION_NOTE = {"along_fit": "highest r close to the fitted line",
                  "max_r_hot": "highest r among seasons with above-median heatwave days"}
#: "close to the fitted line": |residual| within this many residual std devs
FIT_BAND = 0.5
MARK_OFFSET = {2003: (8, 2, "left"), 2015: (0, -19, "center"), p2.REF_YEAR: (-8, 4, "right")}
PERIOD = (1990, p2.REF_YEAR)
#: intensity: the project's heatwave colours, minus the near-white end that
#: would make dots vanish on the surface
INTENSITY_CMAP = LinearSegmentedColormap.from_list(
    "intensity", vz.heatwave_cmap()(np.linspace(0.22, 1.0, 256)))


def season_table(key, corr, period=PERIOD):
    """One row per scored season: year, member, r, heatwave days, mean intensity."""
    with xr.open_dataset(hw.metrics_file(key)) as ds:
        hwd = land_mean_series(ds.hwd.load())
        # same land cells and weights in both means, so the ratio is
        # sum(w * cumulative heat) / sum(w * heatwave days)
        intensity = (land_mean_series(ds.hwcum.load()) / hwd).where(hwd > 0)
    r = corr[key]
    mdim = p2.member_dim(r)
    if mdim:
        r = r.rename({mdim: "member"})
    r, hwd, intensity = xr.align(r, hwd, intensity, join="inner")
    df = (xr.Dataset({"r": r, "hwd": hwd, "intensity": intensity}).to_dataframe()
          .reset_index().dropna(subset=["r", "hwd"]))
    if "member" not in df:
        df["member"] = np.nan
    return df.loc[df.year.between(*period), ["year", "member", "r", "hwd", "intensity"]]


def detrend(t, key):
    """
    The season table with the forced part of the heatwave days removed:
    the ensemble mean for MPI-GE, otherwise a linear fit in time over every
    season in `t`, one per member (for ERA5 that includes 2026).
    """
    t = t.copy()
    cols = ["hwd"]
    if t.member.nunique() >= p2.ENSEMBLE_MEAN_MIN_MEMBERS:
        t[cols] = t[cols] - t.groupby("year")[cols].transform("mean")
        return t
    for _, g in t.groupby(t.member.fillna(-1)):
        for c in cols:
            b = np.polyfit(g.year, g[c], 1)
            t.loc[g.index, c] = g[c] - np.polyval(b, g.year)
    return t


def top_on_both(t):
    """
    The season whose LOWER percentile rank -- of pattern correlation and of
    heatwave days, within the dataset -- is highest: near the top of both.
    """
    t = t.copy()
    pct = t[["r", "hwd"]].rank(pct=True)
    t["rank_r"], t["rank_hwd"] = pct.r, pct.hwd
    t["score"] = pct.min(axis=1)
    return t.loc[t.score.idxmax()]


def along_fit(t, fit):
    """The highest-r season among those within FIT_BAND residual std devs of the fit."""
    t = t.copy()
    resid = t.hwd - (fit.intercept + fit.slope * t.r)
    t["rank_r"] = t.r.rank(pct=True)
    t["rank_hwd"] = t.hwd.rank(pct=True)
    near = t[resid.abs() <= FIT_BAND * resid.std()]
    return near.loc[near.r.idxmax()]


def max_r_hot(t):
    """The highest-r season among those with more heatwave days than the dataset median."""
    t = t.copy()
    t["rank_r"] = t.r.rank(pct=True)
    t["rank_hwd"] = t.hwd.rank(pct=True)
    hot = t[t.hwd > t.hwd.median()]
    return hot.loc[hot.r.idxmax()]


def select_season(key, t, fit):
    """The season a model panel rings and figure 1 maps (see SELECTION)."""
    rule = SELECTION.get(key)
    if rule == "along_fit":
        return along_fit(t, fit)
    if rule == "max_r_hot":
        return max_r_hot(t)
    return top_on_both(t)


def selection_caption():
    """'MPI-GE: ...; EERIE: ...' -- the datasets that do not use the default rule."""
    return "; ".join(f"{k}: {SELECTION_NOTE[v]}" for k, v in SELECTION.items())


def scatter_intensity(ax, t, norm, s):
    """Seasons coloured by mean heatwave intensity; seasons with no heatwave hollow."""
    has = t.intensity.notna()
    ax.scatter(t.r[~has], t.hwd[~has], s=s, facecolor="none", edgecolor=vz.INK_MUTED,
               linewidths=0.7, zorder=2)
    return ax.scatter(t.r[has], t.hwd[has], s=s, c=t.intensity[has], cmap=INTENSITY_CMAP,
                      norm=norm, edgecolor=vz.SURFACE, linewidths=0.3, zorder=2.6)


def main(half=None):
    corr = xr.open_dataset(p2.corr_file(VARIANT, half)).load()
    for detrended in (False, True):
        plot(corr, half, detrended)


def plot(corr, half, detrended):
    vz.use_style()
    tables = {k: season_table(k, corr) for k in KEYS}
    # the ringed model season is always chosen on the raw values, so both
    # versions ring the season the heatwave-analogue figure 1 maps; the
    # detrended version draws it at its detrended position
    chosen = {k: select_season(k, t, stats.linregress(t.r, t.hwd))
              for k, t in tables.items() if k != "ERA5"}
    if detrended:
        tables = {k: detrend(t, k) for k, t in tables.items()}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}

    ymax = max(t.hwd.max() for t in tables.values())
    ymin = min(t.hwd.min() for t in tables.values())
    allint = pd.concat([t.intensity for t in tables.values()]).dropna()
    inorm = Normalize(*np.floor(np.percentile(allint, [2, 98]) * 2) / 2 + [0, 0.5])
    fig, axes = plt.subplots(1, len(KEYS), figsize=(3.3 * len(KEYS) + 1.2, 5.0),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.9, top=0.68, bottom=0.115, wspace=0.08)

    for ax, key in zip(axes, KEYS):
        t = tables[key].sort_values("year")
        n = len(t)
        sc = scatter_intensity(ax, t, inorm, s=7 if n > 1000 else 24)
        if key == "ERA5":
            full = t.set_index("year")
            marks = [(full.loc[y], y, f"JJA {y}") for y in ERA5_MARKED]
        else:
            c = chosen[key]
            same = (t.year == c.year) & ((t.member == c.member) | (t.member.isna() & (c.member != c.member)))
            sel = t[same].iloc[0]
            lab = f"JJA {int(sel.year)}" + (f" r{int(sel.member)}" if sel.member == sel.member else "")
            marks = [(sel, None, lab)]
        for row, y, lab in marks:
            ax.scatter(row.r, row.hwd, s=90, facecolor="none", edgecolor=vz.INK,
                       linewidths=1.5, zorder=4)
            dx, dy, ha = MARK_OFFSET.get(y, (8, 2, "left"))
            if ha == "left" and row.r > 0.75:     # near the right edge: label on the left
                dx, ha = -abs(dx), "right"
            ax.annotate(lab, (row.r, row.hwd), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, va="bottom", fontsize=7.5, color=vz.INK, zorder=5,
                        bbox=dict(facecolor=vz.SURFACE, alpha=0.8, edgecolor="none",
                                  boxstyle="round,pad=0.15"))
        f = fits[key]
        xs = np.array([t.r.min(), t.r.max()])
        ax.plot(xs, f.intercept + f.slope * xs, color=vz.INK, lw=1.6, zorder=3)

        d = p2.DATASETS[key]
        ax.set_title(f"{d.label}\n{d.resolution}", fontsize=8.8, color=vz.INK, linespacing=1.3)
        p_txt = "p < 0.001" if f.pvalue < 0.001 else f"p = {f.pvalue:.3f}"
        ax.text(0.03, 0.97, f"r = {f.rvalue:+.2f}  ({p_txt})\n"
                            f"slope = {f.slope / 10:+.2f} days per +0.1",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.35,
                zorder=6, bbox=dict(facecolor=vz.SURFACE, alpha=0.85, edgecolor="none",
                                    boxstyle="round,pad=0.25"))
        ax.grid(False)
        ax.set_xlabel("SST pattern correlation with ERA5 JJA 2026", fontsize=8.2)
        if detrended:
            ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
        ax.tick_params(labelsize=7.5)

    cax = fig.add_axes([0.915, 0.115, 0.012, 0.565])
    cb = fig.colorbar(sc, cax=cax, extend="both")
    cb.set_label("mean heatwave intensity\n(Tmax anomaly on heatwave days, °C)",
                 fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    axes[0].set_ylabel("S. European land-mean JJA heatwave days"
                       + ("\nanomaly from forced part" if detrended else ""), fontsize=8.5)
    # -0.75 holds the whole-box and Atlantic correlations; the Mediterranean
    # ones reach further out, so widen to the data where needed
    rmin = min(t.r.min() for t in tables.values())
    axes[0].set_xlim(min(-0.75, rmin - 0.05), 1.08)
    # headroom: the stats text sits above the data
    if detrended:
        axes[0].set_ylim(ymin - 1, ymax + 0.3 * (ymax - ymin))
    else:
        axes[0].set_ylim(-0.5, ymax * 1.2)


    fig.suptitle("Do 2026-like summer SST patterns come with more European heatwave days? "
                 f"JJA {PERIOD[0]}–{PERIOD[1]}{scored_over(half)}"
                 + ("  ·  trend removed from heatwave days" if detrended else ""),
                 fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.48, 0.92,
             "one point per JJA season (and member), coloured by mean heatwave intensity "
             "(mean Tmax anomaly over the season's European land heatwave days)\n"
             f"x: pattern correlation over {p2.half_label(half)}  ·  y: land-mean heatwave "
             f"days, {hw.MEAN_BOX_LABEL} (Xu et al. 2026)  ·  both against 1991–2020 of the same dataset\n"
             + ("heatwave days with the forced part removed: MPI-GE minus its 50-member mean; ERA5, EPOC and each "
                f"EERIE member minus a linear fit over each record's seasons in {PERIOD[0]}–{PERIOD[1]}  ·  "
                "pattern correlation as is (box mean already removed)\n"
                if detrended else "")
             + f"ringed: ERA5 {p2.REF_YEAR} (reference, r = 1, included in the fit), 2015, 2003; models: "
             + "the season shown in the heatwave-analogue figure 1"
             + (", selected on raw heatwave days" if detrended else "")
             + f" ({selection_caption()})"
             + "\nline: least-squares fit  ·  r = Pearson correlation of the linear relationship, "
             "two-sided p-value  ·  slope: heatwave days per +0.1 of pattern correlation",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.4)

    out = p2.FIG_DIR / (f"fig8_heatwave_days_vs_pattern_corr{p2.half_suffix(half)}"
                        + ("_detrended" if detrended else "") + ".png")
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    for k, f in fits.items():
        print(f"  {k:16s} n={len(tables[k]):5d}  r = {f.rvalue:+.3f}  p = {f.pvalue:.2g}  "
              f"slope = {f.slope / 10:+.2f} days per +0.1  intercept = {f.intercept:.2f}")


if __name__ == "__main__":
    main(p2.half_from_argv(sys.argv))
