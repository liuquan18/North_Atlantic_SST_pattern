"""
Figure 3 -- not just "can a model make this pattern once", but how often, and
whether that changes with the forcing.

Figure 1 shows each simulation's single best season, which by construction is
its most flattering one. This figure shows the whole distribution behind that
number as split violins: each row is one model, the upper half one forcing
state and the lower half another, so a shift of the distribution with warming
reads directly as an asymmetry between the two halves.

  MPI-GE, MPI-ER   upper: last 30 years of the record, lower: first 30 years
  EPOC, EERIE      upper: forced run (transient GHG / hist+ssp245),
                   lower: constant-forcing control

The benchmark line is taken from the observations themselves: the highest
correlation any *other* observed summer since 1940 reaches against JJA 2026.
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

#: years per half for the single-run-type rows
WINDOW = 30

HALF_WIDTH = 0.40   # max extent of each half-violin, in row units


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


def _window(da, first):
    yrs = _years(da)
    sel = yrs[:WINDOW] if first else yrs[-WINDOW:]
    return _vals(da.sel(year=sel)), f"{sel[0]}–{sel[-1]}"


def build_rows(corr):
    """(label, color, (upper_vals, upper_label), (lower_vals, lower_label))"""
    rows = []
    for key in ["MPI-GE", "MPI-ER"]:
        if key in corr:
            rows.append((p2.DATASETS[key].label, vz.COLORS[key],
                         _window(corr[key], False), _window(corr[key], True)))
    pairs = [("ICON-EPOC-hist", "ICON-EPOC-ctrl", "km-scale ICON (EPOC)",
              "transient GHG", "control"),
             ("EERIE", "EERIE-ctrl", "EERIE ICON-ESM-ER",
              "hist+ssp245", "control")]
    for fk, ck, label, fl, cl in pairs:
        if fk in corr and ck in corr:
            f, c = corr[fk], corr[ck]
            rows.append((label, vz.COLORS[fk],
                         (_vals(f), f"{fl} {_years(f)[0]}–{_years(f)[-1]}"),
                         (_vals(c), f"{cl} {_years(c)[0]}–{_years(c)[-1]}")))
    return rows


def half_violin(ax, vals, y, color, upper):
    """Draw one half of a horizontal violin plus its summary marks."""
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
    ax.plot([vals.max()], [yb], marker="o", ms=4.5,
            mfc=color if upper else vz.SURFACE, mec=color, mew=1.1, zorder=6)


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    corr = {v: xr.open_dataset(p2.RESULT_DIR / f"corr_{v}.nc") for v in p2.VARIANTS}

    variants = p2.PLOT_VARIANTS
    n_rows = len(build_rows(corr[variants[0]]))
    # margins in absolute inches, so the chrome does not grow with the number
    # of rows the way a fractional top/bottom would
    fig_h = 1.15 * n_rows + 1.65
    fig, axes = plt.subplots(1, len(variants),
                             figsize=(5.6 * len(variants) + 3.4, fig_h),
                             sharey=True, squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(wspace=0.06, top=1 - 1.05 / fig_h, bottom=0.55 / fig_h,
                        left=0.205, right=0.985)

    for j, variant in enumerate(variants):
        ax = axes[j]
        rows = build_rows(corr[variant])
        era5 = best[(best.variant == variant) & (best.dataset == "ERA5")].iloc[0]
        obs_ref, obs_year = float(era5.best_r), int(era5.best_year)

        ax.axvline(obs_ref, color=vz.INK, lw=1.0, ls=":", zorder=4)
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

        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([r[0].replace(" (", "\n(") for r in reversed(rows)],
                           fontsize=8.3)
        if len(variants) > 1:
            ax.set_title(p2.VARIANTS[variant]["title"], fontsize=9.5, color=vz.INK, pad=7)
        ax.set_xlabel("pattern correlation with observed JJA 2026")
        ax.set_xlim(-0.9, 0.9)
        ax.set_ylim(-0.55, n_rows - 0.25)
        # left of the benchmark line, in the strip above the top violin
        ax.text(obs_ref, n_rows - 0.30, f"closest observed\nanalogue: JJA {obs_year} ",
                fontsize=7.4, color=vz.INK, va="top", ha="right", linespacing=1.3)

    fig.suptitle("Global warming influences 2026-like summer SST pattern?",
                 fontsize=12, color=vz.INK, y=1 - 0.26 / fig_h)
    fig.text(0.5, 1 - 0.60 / fig_h,
             f"upper half (warmer): last {WINDOW} years / forced run  ·  "
             f"lower half: first {WINDOW} years / control  ·  "
             f"{p2.VARIANTS[variants[0]]['long']}  ·  "
             "bar = IQR, tick = median, dot = best season",
             ha="center", fontsize=8.2, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig3_correlation_distribution.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
