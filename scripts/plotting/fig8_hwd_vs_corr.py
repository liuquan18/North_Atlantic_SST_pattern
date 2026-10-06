"""
Figure 8 -- does a 2026-like SST pattern come with more European heatwave days?

One scatter panel per dataset with daily Tmax (ERA5, MPI-GE, km-scale ICON
(EPOC) and EERIE; MPI-ESM1.2-ER has none). Every season (and member) in
PERIOD = 1990-2025 is one point -- the window the observations and EPOC
cover, so the models are compared over the same climate and not scored on
their late-century seasons:

  x   pattern correlation of its JJA SST pattern with ERA5 JJA 2026
      (box mean removed, 1 deg common ocean mask; results/corr_spatial.nc)
  y   land-mean JJA heatwave days over the whole European map
      (Xu et al. 2026 definition, 1991-2020 reference of the same dataset)

Points are coloured by the season's mean heatwave intensity: the mean Tmax
anomaly over all European land heatwave days of that season (area-weighted
cumulative heat / area-weighted heatwave days -- the "intensity" row of
fig1_heatwaves.py, pooled over the map). Seasons without any heatwave day
have no intensity and are drawn hollow. One colour scale serves all panels;
MPI-GE's 1,800 seasons are drawn smaller. In each model panel the ringed
point is the season shown in fig1_heatwave_analogues.py (select_season:
high on both axes, or for MPI-GE furthest along the fitted line); in ERA5, JJA 2015, 2003 and 2026 are ringed -- 2026 drawn as the
reference, outside PERIOD, at r = 1 by construction and left out of the fit. Each panel carries the least-squares line of y on x with its Pearson r, two-sided p-value and slope
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
from fig1_patterns_simulations import VARIANT, land_mean_series

KEYS = ["ERA5", "MPI-GE", "ICON-EPOC-hist", "EERIE"]
#: observed seasons marked in the ERA5 panel; 2026 is outside PERIOD and is
#: drawn as the reference (r = 1 by construction), not used in the fit
ERA5_MARKED = [p2.REF_YEAR, 2015, 2003]
#: label offsets (points, ha) -- 2003 and 2015 sit almost on top of each other
#: how each model's season is selected (default "both"):
#:   "both"       top_on_both -- highest lower percentile rank of r and heatwave days
#:   "along_fit"  along_fit   -- highest r among seasons close to the fitted line;
#:                for MPI-GE, where "both" lands on the long hot tail of 1,800
#:                seasons rather than on the relationship the line describes
SELECTION = {"MPI-GE": "along_fit"}
#: "close to the fitted line": |residual| within this many residual std devs
FIT_BAND = 0.5
MARK_OFFSET = {2003: (8, 2, "left"), 2015: (0, -19, "center"), p2.REF_YEAR: (-8, 4, "right")}
PERIOD = (1990, 2025)
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


def select_season(key, t, fit):
    """The season a model panel rings and figure 1 maps (see SELECTION)."""
    return along_fit(t, fit) if SELECTION.get(key) == "along_fit" else top_on_both(t)


def scatter_intensity(ax, t, norm, s):
    """Seasons coloured by mean heatwave intensity; seasons with no heatwave hollow."""
    has = t.intensity.notna()
    ax.scatter(t.r[~has], t.hwd[~has], s=s, facecolor="none", edgecolor=vz.INK_MUTED,
               linewidths=0.7, zorder=2)
    return ax.scatter(t.r[has], t.hwd[has], s=s, c=t.intensity[has], cmap=INTENSITY_CMAP,
                      norm=norm, edgecolor=vz.SURFACE, linewidths=0.3, zorder=2.6)


def main():
    vz.use_style()
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    tables = {k: season_table(k, corr) for k in KEYS}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}

    ymax = max(t.hwd.max() for t in tables.values())
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
            full = season_table(key, corr, (PERIOD[0], p2.REF_YEAR)).set_index("year")
            ref = full.loc[[p2.REF_YEAR]].reset_index()
            scatter_intensity(ax, ref, inorm, s=24)
            marks = [(full.loc[y], y, f"JJA {y}") for y in ERA5_MARKED]
        else:
            sel = select_season(key, t, fits[key])
            lab = f"JJA {int(sel.year)}" + (f" r{int(sel.member)}" if sel.member == sel.member else "")
            marks = [(sel, None, lab)]
        for row, y, lab in marks:
            ax.scatter(row.r, row.hwd, s=90, facecolor="none", edgecolor=vz.INK,
                       linewidths=1.5, zorder=4)
            dx, dy, ha = MARK_OFFSET.get(y, (8, 2, "left"))
            ax.annotate(lab, (row.r, row.hwd), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, va="bottom", fontsize=7.5, color=vz.INK, zorder=5)
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
        ax.tick_params(labelsize=7.5)

    cax = fig.add_axes([0.915, 0.115, 0.012, 0.565])
    cb = fig.colorbar(sc, cax=cax, extend="both")
    cb.set_label("mean heatwave intensity\n(Tmax anomaly on heatwave days, °C)",
                 fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    axes[0].set_ylabel("European land-mean JJA heatwave days", fontsize=8.5)
    axes[0].set_xlim(-0.75, 1.08)
    axes[0].set_ylim(-0.5, ymax * 1.2)   # headroom: the stats text sits above the data


    fig.suptitle("Do 2026-like summer SST patterns come with more European heatwave days? "
                 f"JJA {PERIOD[0]}–{PERIOD[1]}",
                 fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.48, 0.92,
             "one point per JJA season (and member), coloured by mean heatwave intensity "
             "(mean Tmax anomaly over the season's European land heatwave days)\n"
             "x: pattern correlation, box mean removed, 30–60°N, 80°W–40°E  ·  y: land-mean heatwave "
             "days, 35–70°N, 12°W–42°E (Xu et al. 2026)  ·  both against 1991–2020 of the same dataset\n"
             f"ringed: ERA5 {p2.REF_YEAR} (reference, r = 1, not in the fit), 2015, 2003; models: the season "
             "shown in the heatwave-analogue figure 1 (MPI-GE: highest r close to the fitted line)\nline: least-squares fit  ·  r = Pearson correlation of the linear relationship, "
             "two-sided p-value  ·  slope: heatwave days per +0.1 of pattern correlation",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.4)

    out = p2.FIG_DIR / "fig8_heatwave_days_vs_pattern_corr.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    for k, f in fits.items():
        print(f"  {k:16s} n={len(tables[k]):5d}  r = {f.rvalue:+.3f}  p = {f.pvalue:.2g}  "
              f"slope = {f.slope / 10:+.2f} days per +0.1  intercept = {f.intercept:.2f}")


if __name__ == "__main__":
    main()
