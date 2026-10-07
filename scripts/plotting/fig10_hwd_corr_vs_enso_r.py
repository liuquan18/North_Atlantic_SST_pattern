"""
Figure 10 -- which predictor do European heatwave days follow more closely,
the 2026-like North Atlantic SST pattern or ENSO? One point per dataset:

  x   Pearson r of land-mean JJA heatwave days on the SST pattern
      correlation with ERA5 JJA 2026 (the fit of fig8_hwd_vs_corr.py)
  y   Pearson r of the same heatwave days on relative Nino-3.4
      (the fit of fig9_hwd_vs_enso.py)

All datasets with daily temperature -- ERA5, the transient runs (MPI-GE,
EPOC, EERIE) and the constant-forcing controls (EPOC, EERIE, sap0006) --
over the same PERIOD = 1990-2025 as fig 8 (for the controls these are model
years, which carry no forcing meaning; EERIE control and sap0006 simply
contribute those 36 of their seasons). MPI-ER has no daily Tmax.

Transient runs are filled, controls hollow in the same hue. MPI-GE (diamond)
is the pooled r over all 50 members x 36 seasons; its 50 single-member r
pairs are the small dots -- the spread a single 36-season run can show from
internal variability alone. The dashed lines mark |r| at p = 0.05 for
n = 36 (two-sided), the sample of a single run. EPOC control heatwaves come
from daily MEAN temperature (its Tmax is on tape only).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy import stats

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns_simulations import VARIANT
from fig8_hwd_vs_corr import PERIOD, season_table
from fig9_hwd_vs_enso import load_enso

#: (key, legend label, transient?) -- control hollow, transient filled
KEYS = [
    ("ERA5", "ERA5 (observed)", True),
    ("MPI-GE", "MPI-ESM1-2-LR GE (pooled, 50 members)", True),
    ("ICON-EPOC-hist", "EPOC km-scale ICON, hist+ssp585", True),
    ("ICON-EPOC-ctrl", "EPOC km-scale ICON, control (Tmean)", False),
    ("EERIE", "EERIE ICON-ESM-ER, hist+ssp245 (3 members)", True),
    ("EERIE-ctrl", "EERIE ICON-ESM-ER, control", False),
    ("ICON-sap-ctrl", "ICON Sapphire sap0006, control", False),
]
#: sap0006 has no project colour; MPI-ER's violet is free here (no Tmax)
COLORS = {**vz.COLORS, "ICON-sap-ctrl": vz.COLORS["MPI-ER"]}


def paired_table(key, corr, enso):
    """Seasons in PERIOD with both predictors: year, member, r_pat, r_enso, hwd."""
    a = season_table(key, corr).rename(columns={"r": "r_pat"})
    b = season_table(key, enso).rename(columns={"r": "r_enso"})
    return a.merge(b[["year", "member", "r_enso"]], on=["year", "member"],
                   how="inner").drop(columns="intensity")


def pearson(t):
    return (stats.pearsonr(t.r_pat, t.hwd).statistic,
            stats.pearsonr(t.r_enso, t.hwd).statistic)


def main():
    vz.use_style()
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    enso = load_enso()
    tables = {k: paired_table(k, corr, enso) for k, _, _ in KEYS}
    points = {k: pearson(t) for k, t in tables.items()}
    members = np.array([pearson(g) for _, g in tables["MPI-GE"].groupby("member")])
    n1 = len(np.arange(PERIOD[0], PERIOD[1] + 1))
    tcrit = stats.t.ppf(0.975, n1 - 2)
    rcrit = tcrit / np.sqrt(n1 - 2 + tcrit ** 2)

    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.83, bottom=0.27)

    lim = 0.75
    ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
    ax.axvline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color=vz.INK_MUTED, lw=0.6, ls=":", zorder=1)
    ax.text(0.62, 0.66, "y = x", fontsize=7, color=vz.INK_MUTED, rotation=45,
            ha="center", va="center")
    # for v in (-rcrit, rcrit):
    #     ax.axhline(v, color=vz.INK_MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
    #     ax.axvline(v, color=vz.INK_MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
    # ax.text(rcrit + 0.01, -lim + 0.02, f"p = 0.05, n = {n1}", fontsize=7,
    #         color=vz.INK_MUTED, rotation=90, ha="left", va="bottom")

    ax.scatter(members[:, 0], members[:, 1], s=12, color=COLORS["MPI-GE"], alpha=0.35,
               edgecolor="none", zorder=2, label="MPI-ESM1-2-LR GE, single members")
    for key, label, transient in KEYS:
        x, y = points[key]
        c = COLORS[key]
        marker = {"ERA5": "*", "MPI-GE": "D"}.get(key, "o")
        size = {"ERA5": 190, "MPI-GE": 60}.get(key, 85)
        ax.scatter(x, y, s=size, marker=marker,
                   facecolor=c if transient else "none", edgecolor=c,
                   linewidths=0.8 if transient else 2.0, zorder=4,
                   label=f"{label}  (n = {len(tables[key])})")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.tick_params(labelsize=8)
    ax.set_xlabel("r (heatwave days, ATL-MED SST pattern correlation)", fontsize=9)
    ax.set_ylabel("r (heatwave days, relative Niño-3.4)", fontsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, fontsize=7.3,
              frameon=False, handletextpad=0.4, columnspacing=1.2)

    fig.suptitle("European heatwave days: North Atlantic SST pattern vs ENSO, "
                 f"JJA {PERIOD[0]}–{PERIOD[1]}", fontsize=11.5, color=vz.INK, y=0.975)
    fig.text(0.54, 0.925,
             "each point: Pearson r over the seasons (and members) of one dataset  ·  "
             "x from figure 8, y from figure 9\n"
             "heatwave days: land mean, 35–70°N, 12°W–42°E (Xu et al. 2026)  ·  ENSO: "
             "Niño-3.4 minus 20°S–20°N SST anomaly\n"
             "filled: transient forcing  ·  hollow: constant-forcing control "
             "(model years 1990–2025)  ·  EPOC control from daily mean T",
             ha="center", va="top", fontsize=7.4, color=vz.INK_SOFT, linespacing=1.45)

    out = p2.FIG_DIR / "fig10_hwd_r_pattern_vs_enso.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")
    print(f"  |r| at p = 0.05 for n = {n1}: {rcrit:.3f}")
    for key, _, _ in KEYS:
        x, y = points[key]
        print(f"  {key:16s} n={len(tables[key]):5d}  r_pattern = {x:+.3f}  r_enso = {y:+.3f}")
    print(f"  MPI-GE members: r_pattern {members[:, 0].mean():+.3f} "
          f"[{members[:, 0].min():+.3f}, {members[:, 0].max():+.3f}], "
          f"r_enso {members[:, 1].mean():+.3f} "
          f"[{members[:, 1].min():+.3f}, {members[:, 1].max():+.3f}]")


if __name__ == "__main__":
    main()
