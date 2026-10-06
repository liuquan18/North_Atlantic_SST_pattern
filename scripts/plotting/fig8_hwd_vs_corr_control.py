"""
Figure 8 (control runs) -- the scatter of fig8_hwd_vs_corr.py for the three
km-scale runs under constant forcing:

  km-scale ICON EPOC control   (epoc2_010, constant 1990 forcing, 1990-2024)
  EERIE ICON-ESM-ER control    (eerie-control-1950, constant 1950, 1950-2050)
  ICON Sapphire sap0006        (constant 1950 forcing, 1961-2049 used)

Every season of each run is one point (no 1990-2025 window: without forcing
there is no trend to keep out). Same axes, colouring, fit and statistics as
fig8_hwd_vs_corr.py.

Caveats, both from what is on disk:
  * EPOC control: its daily-max files are on tape only, so its heatwaves are
    detected in daily MEAN temperature. The caption reports how the EPOC
    transient run's heatwave days compare between Tmean and Tmax detection.
    Its intensities are Tmean anomalies, smaller than the Tmax ones, so its
    dots sit pale on the shared colour scale.
  * sap0006: Tmax is the maximum of the 8 three-hourly instantaneous values.
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
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT, land_mean_series
from fig8_hwd_vs_corr import scatter_intensity, season_table

KEYS = ["ICON-EPOC-ctrl", "EERIE-ctrl", "ICON-sap-ctrl"]
ALL_YEARS = (0, 9999)
SUBTITLE = {
    "ICON-EPOC-ctrl": "constant 1990 forcing  ·  heatwaves from daily MEAN T",
    "EERIE-ctrl": "constant 1950 forcing",
    "ICON-sap-ctrl": "constant 1950 forcing  ·  Tmax from 3-hourly values",
}


def tmean_check():
    """EPOC transient: land-mean heatwave days from Tmax vs from Tmean detection."""
    out = {}
    for key in ("ICON-EPOC-hist", "ICON-EPOC-hist-tmean"):
        with xr.open_dataset(hw.metrics_file(key)) as ds:
            out[key] = land_mean_series(ds.hwd.load()).to_series()
    both = pd.concat(out, axis=1).dropna()
    r = np.corrcoef(both.iloc[:, 0], both.iloc[:, 1])[0, 1]
    return r, both.iloc[:, 0].mean(), both.iloc[:, 1].mean(), len(both)


def main():
    vz.use_style()
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    tables = {k: season_table(k, corr, ALL_YEARS) for k in KEYS}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}
    r_tm, d_tmax, d_tmean, n_tm = tmean_check()

    ymax = max(t.hwd.max() for t in tables.values())
    # colour scale from the two Tmax-based runs, so the Tmean one does not stretch it
    allint = pd.concat([tables[k].intensity for k in KEYS if k != "ICON-EPOC-ctrl"]).dropna()
    inorm = Normalize(*np.floor(np.percentile(allint, [2, 98]) * 2) / 2 + [0, 0.5])

    fig, axes = plt.subplots(1, len(KEYS), figsize=(3.6 * len(KEYS) + 1.2, 5.0),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.88, top=0.68, bottom=0.115, wspace=0.08)

    for ax, key in zip(axes, KEYS):
        t = tables[key].sort_values("year")
        sc = scatter_intensity(ax, t, inorm, s=24)
        f = fits[key]
        xs = np.array([t.r.min(), t.r.max()])
        ax.plot(xs, f.intercept + f.slope * xs, color=vz.INK, lw=1.6, zorder=3)

        d = p2.DATASETS[key]
        y0, y1 = int(t.year.min()), int(t.year.max())
        ax.set_title(f"{d.label}\n{d.resolution}  ·  JJA {y0}–{y1}\n{SUBTITLE[key]}",
                     fontsize=8.6, color=vz.INK, linespacing=1.3)
        p_txt = "p < 0.001" if f.pvalue < 0.001 else f"p = {f.pvalue:.3f}"
        ax.text(0.03, 0.97, f"r = {f.rvalue:+.2f}  ({p_txt})\n"
                            f"slope = {f.slope / 10:+.2f} days per +0.1\nn = {len(t)}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.8,
                color=vz.INK_SOFT, linespacing=1.35, zorder=6,
                bbox=dict(facecolor=vz.SURFACE, alpha=0.85, edgecolor="none",
                          boxstyle="round,pad=0.25"))
        ax.grid(False)
        ax.set_xlabel("SST pattern correlation with ERA5 JJA 2026", fontsize=8.2)
        ax.tick_params(labelsize=7.5)

    cax = fig.add_axes([0.895, 0.115, 0.012, 0.565])
    cb = fig.colorbar(sc, cax=cax, extend="both")
    cb.set_label("mean heatwave intensity\n(Tmax anomaly on heatwave days, °C)",
                 fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    axes[0].set_ylabel("European land-mean JJA heatwave days", fontsize=8.5)
    axes[0].set_xlim(-0.75, 0.85)
    axes[0].set_ylim(-0.5, ymax * 1.25)   # headroom: the stats text sits above the data

    fig.suptitle("Do 2026-like summer SST patterns come with more European heatwave days? "
                 "Control runs", fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.48, 0.92,
             "one point per JJA season, coloured by mean heatwave intensity  ·  x: pattern "
             "correlation with ERA5 JJA 2026, box mean removed, 30–60°N, 80°W–40°E\n"
             "y: land-mean heatwave days, 35–70°N, 12°W–42°E (Xu et al. 2026), against model "
             "years 1991–2020 of the same run  ·  line: least-squares fit, Pearson r, slope "
             "per +0.1 of pattern correlation\n"
             f"EPOC control: daily Tmax only on tape, heatwaves from daily mean T -- in the EPOC "
             f"transient run the two give r = {r_tm:+.3f} year to year ({n_tm} seasons), "
             f"on average {d_tmean:.1f} (Tmean) vs {d_tmax:.1f} (Tmax) days; its colours are Tmean intensities",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.4)

    out = p2.FIG_DIR / "fig8_heatwave_days_vs_pattern_corr_control.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    print(f"  EPOC transient, Tmax vs Tmean heatwave days: r = {r_tm:+.3f}, "
          f"means {d_tmax:.2f} vs {d_tmean:.2f} days, n = {n_tm}")
    for k, f in fits.items():
        t = tables[k]
        print(f"  {k:16s} n={len(t):4d} ({int(t.year.min())}-{int(t.year.max())})  "
              f"r = {f.rvalue:+.3f}  p = {f.pvalue:.2g}  slope = {f.slope / 10:+.2f} days per +0.1")


if __name__ == "__main__":
    main()
