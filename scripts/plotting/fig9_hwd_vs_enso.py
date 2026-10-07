"""
Figure 9 -- do European heatwave days follow ENSO? The scatter of
fig8_hwd_vs_corr.py with the SST pattern correlation on x replaced by the
season's ENSO state:

  x   relative Nino-3.4 of the same JJA: Nino-3.4 (5S-5N, 170W-120W) minus
      20S-20N mean SST anomaly, 1991-2020 reference of the same dataset
      (08_nino34.sh). Subtracting the tropical mean removes the forced
      warming, which heatwave days share, so the fit is not a trend artefact.
  y   land-mean JJA heatwave days over the whole European map, as in fig 8

Same datasets, PERIOD = 1990-2025 window, intensity colouring, least-squares
fit and statistics as fig 8; the slope is in heatwave days per +1 K. In ERA5,
JJA 2026, 2015 and 2003 are ringed, 2026 drawn outside PERIOD and left out of
the fit, as in fig 8.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import Normalize
from scipy import stats

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig8_hwd_vs_corr import (ERA5_MARKED, KEYS, MARK_OFFSET, PERIOD,
                              scatter_intensity, season_table)

ENSO_FILE = p2.DATA_DIR / "nino34_jja.nc"
XLABEL = "JJA relative Niño-3.4 (K)"


def load_enso():
    """Relative Nino-3.4 per dataset, in the layout season_table expects of `corr`."""
    ds = xr.open_dataset(ENSO_FILE).load()
    return ds[[v for v in ds.data_vars if not v.endswith(("_n34", "_trop"))]]


def fit_text(ax, f, n=None):
    p_txt = "p < 0.001" if f.pvalue < 0.001 else f"p = {f.pvalue:.3f}"
    txt = f"r = {f.rvalue:+.2f}  ({p_txt})\nslope = {f.slope:+.2f} days per +1 K"
    if n is not None:
        txt += f"\nn = {n}"
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=7.8,
            color=vz.INK_SOFT, linespacing=1.35, zorder=6,
            bbox=dict(facecolor=vz.SURFACE, alpha=0.85, edgecolor="none",
                      boxstyle="round,pad=0.25"))


def draw_fit(ax, t, f):
    xs = np.array([t.r.min(), t.r.max()])
    ax.plot(xs, f.intercept + f.slope * xs, color=vz.INK, lw=1.6, zorder=3)
    ax.axvline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)


def colorbar(fig, sc, rect):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(sc, cax=cax, extend="both")
    cb.set_label("mean heatwave intensity\n(Tmax anomaly on heatwave days, °C)",
                 fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)


def main():
    vz.use_style()
    enso = load_enso()
    tables = {k: season_table(k, enso) for k in KEYS}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}

    ymax = max(t.hwd.max() for t in tables.values())
    xabs = max(np.abs(t.r).max() for t in tables.values())
    allint = pd.concat([t.intensity for t in tables.values()]).dropna()
    inorm = Normalize(*np.floor(np.percentile(allint, [2, 98]) * 2) / 2 + [0, 0.5])
    fig, axes = plt.subplots(1, len(KEYS), figsize=(3.3 * len(KEYS) + 1.2, 5.0),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.9, top=0.7, bottom=0.115, wspace=0.08)

    for ax, key in zip(axes, KEYS):
        t = tables[key].sort_values("year")
        sc = scatter_intensity(ax, t, inorm, s=7 if len(t) > 1000 else 24)
        if key == "ERA5":
            full = season_table(key, enso, (PERIOD[0], p2.REF_YEAR)).set_index("year")
            scatter_intensity(ax, full.loc[[p2.REF_YEAR]].reset_index(), inorm, s=24)
            for y in ERA5_MARKED:
                row = full.loc[y]
                ax.scatter(row.r, row.hwd, s=90, facecolor="none", edgecolor=vz.INK,
                           linewidths=1.5, zorder=4)
                dx, dy, ha = MARK_OFFSET.get(y, (8, 2, "left"))
                ax.annotate(f"JJA {y}", (row.r, row.hwd), xytext=(dx, dy),
                            textcoords="offset points", ha=ha, va="bottom", fontsize=7.5,
                            color=vz.INK, zorder=5)
        draw_fit(ax, t, fits[key])
        d = p2.DATASETS[key]
        ax.set_title(f"{d.label}\n{d.resolution}", fontsize=8.8, color=vz.INK, linespacing=1.3)
        fit_text(ax, fits[key], len(t))
        ax.grid(False)
        ax.set_xlabel(XLABEL, fontsize=8.2)
        ax.tick_params(labelsize=7.5)

    colorbar(fig, sc, [0.915, 0.115, 0.012, 0.585])
    axes[0].set_ylabel("European land-mean JJA heatwave days", fontsize=8.5)
    axes[0].set_xlim(-1.05 * xabs, 1.05 * xabs)
    axes[0].set_ylim(-0.5, ymax * 1.2)   # headroom: the stats text sits above the data

    fig.suptitle("Do European summer heatwave days follow ENSO? "
                 f"JJA {PERIOD[0]}–{PERIOD[1]}", fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.48, 0.92,
             "one point per JJA season (and member), coloured by mean heatwave intensity "
             "(mean Tmax anomaly over the season's European land heatwave days)\n"
             "x: relative Niño-3.4 of the same JJA = Niño-3.4 (5°S–5°N, 170°W–120°W) minus "
             "20°S–20°N mean SST anomaly  ·  y: land-mean heatwave days, 35–70°N, 12°W–42°E "
             "(Xu et al. 2026)\nboth against 1991–2020 of the same dataset  ·  "
             f"ringed: ERA5 {p2.REF_YEAR} (outside the window, not in the fit), 2015, 2003  ·  "
             "line: least-squares fit, Pearson r, two-sided p, slope per +1 K",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.4)

    out = p2.FIG_DIR / "fig9_heatwave_days_vs_enso.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    for k, f in fits.items():
        print(f"  {k:16s} n={len(tables[k]):5d}  r = {f.rvalue:+.3f}  p = {f.pvalue:.2g}  "
              f"slope = {f.slope:+.2f} days per +1 K")


if __name__ == "__main__":
    main()
