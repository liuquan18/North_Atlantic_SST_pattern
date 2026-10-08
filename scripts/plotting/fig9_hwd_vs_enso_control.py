"""
Figure 9 (control runs) -- the heatwave-days vs ENSO scatter of
fig9_hwd_vs_enso.py for the three km-scale runs under constant forcing
(EPOC control, EERIE control, ICON Sapphire sap0006), every season of each
run, as in fig8_hwd_vs_corr_control.py -- with the same caveats: EPOC
control heatwaves come from daily MEAN temperature (its Tmax is on tape
only), sap0006 Tmax is the maximum of 3-hourly instantaneous values.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy import stats

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig8_hwd_vs_corr import scatter_intensity, season_table
from fig8_hwd_vs_corr_control import ALL_YEARS, KEYS, SUBTITLE, tmean_check
from fig9_hwd_vs_enso import XLABEL, colorbar, draw_fit, fit_text, load_enso


def main():
    vz.use_style()
    enso = load_enso()
    tables = {k: season_table(k, enso, ALL_YEARS) for k in KEYS}
    fits = {k: stats.linregress(t.r, t.hwd) for k, t in tables.items()}
    r_tm, d_tmax, d_tmean, n_tm = tmean_check()

    ymax = max(t.hwd.max() for t in tables.values())
    xabs = max(np.abs(t.r).max() for t in tables.values())
    # colour scale from the two Tmax-based runs, so the Tmean one does not stretch it
    allint = pd.concat([tables[k].intensity for k in KEYS if k != "ICON-EPOC-ctrl"]).dropna()
    inorm = Normalize(*np.floor(np.percentile(allint, [2, 98]) * 2) / 2 + [0, 0.5])

    fig, axes = plt.subplots(1, len(KEYS), figsize=(3.6 * len(KEYS) + 1.2, 5.0),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.88, top=0.68, bottom=0.115, wspace=0.08)

    for ax, key in zip(axes, KEYS):
        t = tables[key].sort_values("year")
        sc = scatter_intensity(ax, t, inorm, s=24)
        draw_fit(ax, t, fits[key])
        d = p2.DATASETS[key]
        y0, y1 = int(t.year.min()), int(t.year.max())
        ax.set_title(f"{d.label}\n{d.resolution}  ·  JJA {y0}–{y1}\n{SUBTITLE[key]}",
                     fontsize=8.6, color=vz.INK, linespacing=1.3)
        fit_text(ax, fits[key], len(t))
        ax.grid(False)
        ax.set_xlabel(XLABEL, fontsize=8.2)
        ax.tick_params(labelsize=7.5)

    colorbar(fig, sc, [0.895, 0.115, 0.012, 0.565])
    axes[0].set_ylabel("S. European land-mean JJA heatwave days", fontsize=8.5)
    axes[0].set_xlim(-1.05 * xabs, 1.05 * xabs)
    axes[0].set_ylim(-0.5, ymax * 1.25)   # headroom: the stats text sits above the data

    fig.suptitle("Do European summer heatwave days follow ENSO? Control runs",
                 fontsize=12, color=vz.INK, y=0.985)
    fig.text(0.48, 0.92,
             "one point per JJA season, coloured by mean heatwave intensity  ·  x: relative "
             "Niño-3.4 of the same JJA = Niño-3.4 (5°S–5°N, 170°W–120°W) minus 20°S–20°N mean "
             f"SST anomaly\ny: land-mean heatwave days, {hw.MEAN_BOX_LABEL} (Xu et al. 2026), "
             "against model years 1991–2020 of the same run  ·  line: least-squares fit, "
             "Pearson r, slope per +1 K\n"
             f"EPOC control: daily Tmax only on tape, heatwaves from daily mean T -- in the EPOC "
             f"transient run the two give r = {r_tm:+.3f} year to year ({n_tm} seasons), "
             f"on average {d_tmean:.1f} (Tmean) vs {d_tmax:.1f} (Tmax) days; its colours are "
             "Tmean intensities",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.4)

    out = p2.FIG_DIR / "fig9_heatwave_days_vs_enso_control.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    for k, f in fits.items():
        t = tables[k]
        print(f"  {k:16s} n={len(t):4d} ({int(t.year.min())}-{int(t.year.max())})  "
              f"r = {f.rvalue:+.3f}  p = {f.pvalue:.2g}  slope = {f.slope:+.2f} days per +1 K")


if __name__ == "__main__":
    main()
