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

The large ensemble (DENSITY_KEYS) is drawn as its joint probability density
instead: 1,800 overlapping dots hide where the seasons actually are. Filled
contours enclose the most probable 25/50/75/90/95 % of seasons (highest-density
regions of a Gaussian KDE); only the seasons outside the outermost contour
are drawn as dots.

Points and the density wear each dataset's colour from figure 2. Each panel carries the
least-squares line of y on x with its Pearson r and two-sided p-value. The
p-value treats every season as independent -- reasonable across members and
years of JJA, but it ignores any year-to-year persistence.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, land_mean_series

KEYS = ["ERA5", "MPI-GE", "ICON-EPOC-hist", "EERIE"]
PERIOD = (1990, 2025)
DENSITY_KEYS = {"MPI-GE"}
#: probability mass enclosed by each density contour, outermost first
HDR_MASS = [0.95, 0.90, 0.75, 0.50, 0.25]
#: probability shades, light (outer) -> dark (inner), in MPI-GE's blue
DENSITY_CMAP = LinearSegmentedColormap.from_list("density", ["#d6e6f8", vz.COLORS["MPI-GE"], "#0f3a6e"])


def season_table(key, corr):
    """One row per scored season: year, member, r, heatwave days."""
    with xr.open_dataset(hw.metrics_file(key)) as ds:
        hwd = land_mean_series(ds.hwd.load())
    r = corr[key]
    mdim = p2.member_dim(r)
    if mdim:
        r = r.rename({mdim: "member"})
    r, hwd = xr.align(r, hwd, join="inner")
    df = xr.Dataset({"r": r, "hwd": hwd}).to_dataframe().reset_index().dropna(subset=["r", "hwd"])
    if "member" not in df:
        df["member"] = np.nan
    return df.loc[df.year.between(*PERIOD), ["year", "member", "r", "hwd"]]


def draw_density(ax, x, y, cmap):
    """
    Joint probability density of (x, y) as nested highest-density regions.

    Each contour is the density level above which HDR_MASS of the seasons lie,
    so the bands read directly as "the most likely p % of seasons".
    Returns the outlying seasons (outside the outermost contour).
    """
    kde = stats.gaussian_kde(np.vstack([x, y]))
    at_points = kde(np.vstack([x, y]))
    levels = [np.quantile(at_points, 1 - m) for m in HDR_MASS]
    gx, gy = np.meshgrid(np.linspace(x.min() - 0.15, x.max() + 0.15, 200),
                         np.linspace(max(y.min() - 3, -0.5), y.max() + 3, 200))
    dens = kde(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    shades = cmap(np.linspace(0.15, 0.85, len(levels)))
    ax.contourf(gx, gy, dens, levels=levels + [dens.max() * 1.01], colors=shades, zorder=2)
    ax.contour(gx, gy, dens, levels=levels, colors=[vz.SURFACE], linewidths=0.6, zorder=2.5)
    outside = at_points < levels[0]
    return outside, levels, shades


def main():
    vz.use_style()
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    tables = {k: season_table(k, corr) for k in KEYS}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}

    ymax = max(t.hwd.max() for t in tables.values())
    fig, axes = plt.subplots(1, len(KEYS), figsize=(3.3 * len(KEYS) + 1.2, 4.8),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.70, bottom=0.12, wspace=0.08)

    for ax, key in zip(axes, KEYS):
        t = tables[key].sort_values("year")
        n = len(t)
        if key in DENSITY_KEYS:
            outside, _, shades = draw_density(ax, t.r.values, t.hwd.values, DENSITY_CMAP)
            ax.scatter(t.r[outside], t.hwd[outside], s=5, color=vz.INK_MUTED, linewidths=0,
                       zorder=2)
            handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in shades[::-1]]
            ax.legend(handles, [f"{m:.0%}" for m in HDR_MASS[::-1]], title="most likely",
                      loc="upper right", fontsize=6.8, title_fontsize=6.8, frameon=False,
                      handlelength=1.2, handleheight=0.9, borderaxespad=0.3)
        else:
            ax.scatter(t.r, t.hwd, s=22, color=vz.COLORS[key], alpha=0.8, linewidths=0,
                       zorder=2)
        f = fits[key]
        xs = np.array([t.r.min(), t.r.max()])
        ax.plot(xs, f.intercept + f.slope * xs, color=vz.INK, lw=1.6, zorder=3)

        d = p2.DATASETS[key]
        ax.set_title(f"{d.label}\n{d.resolution}", fontsize=8.8, color=vz.INK, linespacing=1.3)
        p_txt = "p < 0.001" if f.pvalue < 0.001 else f"p = {f.pvalue:.3f}"
        ax.text(0.03, 0.97, f"r = {f.rvalue:+.2f}  ({p_txt})",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.35)
        ax.grid(True, lw=0.5)
        ax.set_xlabel("SST pattern correlation with ERA5 JJA 2026", fontsize=8.2)
        ax.tick_params(labelsize=7.5)

    axes[0].set_ylabel("European land-mean JJA heatwave days", fontsize=8.5)
    axes[0].set_xlim(-0.75, 0.85)
    axes[0].set_ylim(-0.5, ymax * 1.06)


    fig.suptitle("Do 2026-like summer SST patterns come with more European heatwave days? "
                 f"JJA {PERIOD[0]}–{PERIOD[1]}",
                 fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.5, 0.92,
             "one point per JJA season (and member)  ·  MPI-GE: joint probability "
             "density, contours enclose the most likely 25–95 % of seasons, dots outside 95 %\n"
             "x: pattern correlation, box mean removed, 30–60°N, 80°W–40°E  ·  y: land-mean heatwave "
             "days, 35–70°N, 12°W–42°E (Xu et al. 2026)  ·  both against 1991–2020 of the same dataset\n"
             "line: least-squares fit  ·  r = Pearson correlation of the linear relationship, "
             "two-sided p-value",
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
