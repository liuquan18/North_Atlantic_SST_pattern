"""
Figure 3 (historical variant) -- the same split violins as figure 3, but with
the single-model ensembles compared over the observed era instead of over the
end of the century, and ERA5 itself as the first row.

  ERA5, MPI-GE,    upper: 1990-2025, lower: pre-1980, i.e. as much of
  MPI-ER           1850-1979 as the record covers (ERA5 1940-1979,
                   MPI-ER 1950-1979)
  EPOC, EERIE      upper: forced run over 1990-2025, lower: whole control run
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import gaussian_kde

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz
HALF_WIDTH = 0.40   # max extent of each half-violin, in row units
RECENT = (1990, 2025)
EARLY = (1850, 1979)


def _years(da):
    """Years with at least one finite value."""
    finite = np.isfinite(da)
    mdim = p2.member_dim(da)
    if mdim:
        finite = finite.any(mdim)
    return da.year.values[finite.values]


def _vals(da):
    v = da.values.ravel()
    return v[np.isfinite(v)]


def _period(da, y0, y1):
    """Values and an "a-b" label over the years of [y0, y1] the record covers."""
    yrs = [y for y in _years(da) if y0 <= y <= y1 and y != p2.REF_YEAR]
    return _vals(da.sel(year=yrs)), f"{yrs[0]}–{yrs[-1]}"


def half_violin(ax, vals, y, color, upper):
    """One half of a horizontal violin plus its summary marks."""
    grid = np.linspace(vals.min(), vals.max(), 200)
    dens = gaussian_kde(vals)(grid)
    dens = dens / dens.max() * HALF_WIDTH
    sign = 1 if upper else -1
    ax.fill_between(grid, y, y + sign * dens, facecolor=color,
                    alpha=0.38 if upper else 0.14, lw=0, zorder=2)
    ax.plot(grid, y + sign * dens, color=color, lw=1.0, zorder=3,
            ls="-" if upper else (0, (4, 2)))

    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    yb = y + sign * 0.07
    ax.plot([q1, q3], [yb, yb], color=color, lw=2.6, solid_capstyle="butt", zorder=5)
    ax.plot([med], [yb], marker="|", ms=7, color=vz.SURFACE, mew=1.6, zorder=6)


def build_hist_rows(corr):
    """(label, color, (upper_vals, upper_label), (lower_vals, lower_label))"""
    rows = []
    for key in ["ERA5", "MPI-GE", "MPI-ER"]:
        if key in corr:
            rows.append((p2.DATASETS[key].label, vz.COLORS[key],
                         _period(corr[key], *RECENT), _period(corr[key], *EARLY)))
    pairs = [("ICON-EPOC-hist", "ICON-EPOC-ctrl", "km-scale ICON (EPOC)",
              "hist+ssp585", "control"),
             ("EERIE", "EERIE-ctrl", "EERIE ICON-ESM-ER",
              "hist+ssp245", "control")]
    for fk, ck, label, fl, cl in pairs:
        if fk in corr and ck in corr:
            f, c = corr[fk], corr[ck]
            # forced runs over the same 1990-2025 window as the rows above
            # (EPOC covers exactly that; EERIE runs 1950-2050)
            fv, flab = _period(f, *RECENT)
            rows.append((label, vz.COLORS[fk],
                         (fv, f"{fl} {flab}"),
                         (_vals(c), f"{cl} {_years(c)[0]}\u2013{_years(c)[-1]}")))
    return rows


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    corr = {v: xr.open_dataset(p2.RESULT_DIR / f"corr_{v}.nc") for v in p2.VARIANTS}

    variants = p2.PLOT_VARIANTS
    n_rows = len(build_hist_rows(corr[variants[0]]))
    fig_h = 1.15 * n_rows + 1.65
    fig, axes = plt.subplots(1, len(variants),
                             figsize=(5.6 * len(variants) + 3.4, fig_h),
                             sharey=True, squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(wspace=0.06, top=1 - 1.05 / fig_h, bottom=0.55 / fig_h,
                        left=0.205, right=0.985)

    for j, variant in enumerate(variants):
        ax = axes[j]
        rows = build_hist_rows(corr[variant])
        era5 = best[(best.variant == variant) & (best.dataset == "ERA5")].iloc[0]
        obs_ref, obs_year = float(era5.best_r), int(era5.best_year)

        ax.axvline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
        ax.grid(axis="x", zorder=0)

        for i, (label, color, (up, up_lab), (lo, lo_lab)) in enumerate(rows):
            y = n_rows - 1 - i
            ax.axhline(y, color=vz.INK_MUTED, lw=0.4, xmin=0.02, xmax=0.98, zorder=1)
            half_violin(ax, up, y, color, upper=True)
            half_violin(ax, lo, y, color, upper=False)
            ax.text(-0.88, y + 0.24, f"{up_lab}  (n={up.size:,})", fontsize=7.2,
                    color=vz.INK_SOFT, va="bottom", ha="left")
            ax.text(-0.88, y - 0.24, f"{lo_lab}  (n={lo.size:,})", fontsize=7.2,
                    color=vz.INK_SOFT, va="top", ha="left")
            if label == p2.DATASETS["ERA5"].label:
                # the best observed season is the upper half's dot
                yb = y + 0.07
                ax.plot([obs_ref], [yb], marker="o", ms=4.5, color=vz.INK, zorder=6)
                ax.annotate(f"JJA {obs_year}", (obs_ref, yb), xytext=(0, 7),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=7.6, color=vz.INK)

        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([r[0].replace(" (", "\n(") for r in reversed(rows)],
                           fontsize=8.3)
        if len(variants) > 1:
            ax.set_title(p2.VARIANTS[variant]["title"], fontsize=9.5, color=vz.INK, pad=7)
        ax.set_xlabel("pattern correlation with observed JJA 2026")
        ax.set_xlim(-0.9, 0.9)
        ax.set_ylim(-0.55, n_rows - 0.25)

    fig.suptitle("Global warming influences 2026-like summer SST pattern?",
                 fontsize=12, color=vz.INK, y=1 - 0.26 / fig_h)
    fig.text(0.5, 1 - 0.60 / fig_h,
             f"upper half: {RECENT[0]}–{RECENT[1]} / forced run  ·  "
             f"lower half: before {EARLY[1] + 1} / control  ·  "
             f"{p2.VARIANTS[variants[0]]['long']}  ·  "
             "bar = IQR, tick = median",
             ha="center", fontsize=8.2, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig3_hist_correlation_distribution.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
