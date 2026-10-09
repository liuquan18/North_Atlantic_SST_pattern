"""
Figure 15 -- return periods of European summer heatwaves: observed and simulated.

Rows: land-mean JJA heatwave days and cumulative heat (src.heatwave.MEAN_BOX).
Horizontal lines in every panel: the ERA5 values of 2003 (dashed) and 2026
(solid).

Left, ERA5: return-level curves of the shift+scale nonstationary GEV fitted
to 1950-2026 (02_statistics.py) in the climate of the 1950s, of 2003 and of
2026, with 95% bootstrap bands; the shift fit (constant scale) in the 2026
climate dashed. Dots: every observed season moved into the 2026 climate
(same percentile under the fit), at its empirical return period -- the
WWA-style check of the fit; 2003 and 2026 ringed.

Right, simulations: MPI-GE seasons of all 50 members in 1850-1900, at ERA5's
2003 warming level, around 2026 (calendar +-10 yr) and at ERA5's 2026
warming level, as empirical return periods (no fit); the km-scale
transient runs' nonstationary GEV in their own 2026 climate (EPOC: 2025, its
last season); the control runs' stationary GEV (dashed). Each simulation is
counted against its own 1991-2020 threshold.
"""
import importlib
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.extremes as ex
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

st = importlib.import_module("02_statistics")

DIR = ROOT / "data" / "heatwave_stats"
METRICS = {"hwd": ("heatwave days", "days"), "hwcum": ("cumulative heat", "°C·day")}
EVENTS = (2003, 2026)
MAIN, SHIFT = "shift+scale 1950-2026", "shift 1950-2026"
RP = np.geomspace(1.05, 3000, 200)
COLORS = {**vz.COLORS, "ICON-sap-ctrl": vz.COLORS["MPI-ER"]}
LABEL = {"ICON-EPOC-hist": "EPOC hist (2025)", "EERIE": "EERIE hist (2026)",
         "EERIE-ctrl": "EERIE control", "ICON-sap-ctrl": "Sapphire control",
         "ICON-EPOC-ctrl": "EPOC control (Tmean)"}
#: ERA5 climates: line style, grey level
CLIM_STYLE = {"1950s": (vz.CONTEXT, (0, (1, 1.5))), "2003": (vz.INK_SOFT, (0, (4, 2))),
              "2026": (vz.INK, "-")}
#: MPI-GE windows: lightness ramp of its blue
MPIGE_ALPHA = [0.3, 0.5, 0.75, 1.0]


def gev(row):
    return ex.NSGEV(**{k: float(row[k]) for k in ("mu0", "mu1", "s0", "s1", "xi")})


def event_lines(ax, e_val, m, right=False):
    for ev, ls in ((2026, "-"), (2003, (0, (4, 2)))):
        v = e_val[(m, ev)]
        ax.axhline(v, color=vz.INK, lw=0.8, ls=ls, zorder=1)
        ax.text(RP[-1] if right else RP[0], v, f"ERA5 {ev} ", ha="right" if right else "left",
                va="bottom", fontsize=7, color=vz.INK)


def era5_panel(ax, m, seasons, pars, boot, clim, e_val):
    fit = gev(pars[(pars.metric == m) & (pars.variant == MAIN)].iloc[0])
    bt = [gev(r) for _, r in boot[(boot.metric == m) & (boot.variant == MAIN)].iterrows()]
    for c, T in clim.items():
        color, ls = CLIM_STYLE[c]
        if bt:
            lv = np.array([b.return_level(RP, T) for b in bt])
            ax.fill_between(RP, *np.nanpercentile(lv, [2.5, 97.5], axis=0), color=color,
                            alpha=0.13, lw=0)
        ax.plot(RP, fit.return_level(RP, T), color=color, ls=ls, lw=1.5,
                label=f"{c} climate")
    shift = gev(pars[(pars.metric == m) & (pars.variant == SHIFT)].iloc[0])
    ax.plot(RP, shift.return_level(RP, clim["2026"]), color=vz.INK, lw=0.9, ls=(0, (6, 2, 1, 2)),
            label="2026 climate, shift fit")
    # observations moved into the 2026 climate at their own percentile
    moved = fit.return_level(1 / np.clip(fit.sf(seasons[m].values, seasons["T"].values),
                                         1e-9, None), clim["2026"])
    vals, rps = ex.plotting_positions(moved)
    ax.scatter(rps, vals, s=10, color=vz.INK_MUTED, zorder=3, lw=0)
    for ev in EVENTS:
        mv = float(moved[seasons.year.values == ev][0])
        rp = rps[np.argmin(np.abs(vals - mv))]
        ax.scatter(rp, mv, s=48, facecolor="none", edgecolor=vz.INK, lw=1.1, zorder=4)
        ax.annotate(f"{ev}", (rp, mv), xytext=(5, -9), textcoords="offset points", fontsize=7,
                    color=vz.INK)
    event_lines(ax, e_val, m, right=True)
    return float(np.nanmax(moved[np.isin(seasons.year.values, EVENTS)]))


def model_panel(ax, m, t, mpars, mrp, e_clim, e_val):
    g = t[t.dataset == "MPI-GE"]
    windows = [("1850–1900", g[g.year.between(*st.PREIND)])]
    for ev in EVENTS:
        windows.append((f"at ERA5 {ev} warming", g[(g["T"] - e_clim[str(ev)]).abs() <= st.WL_TOL]))
        if ev == 2026:
            windows.insert(-1, (f"2016–2036", g[g.year.between(2016, 2036)]))
    for (lab, w), a in zip(windows, MPIGE_ALPHA):
        vals, rps = ex.plotting_positions(w[m].values)
        ax.plot(rps, vals, color=COLORS["MPI-GE"], alpha=a, lw=1.6,
                label=f"MPI-GE {lab} (n={len(w)})")
    for key in ("ICON-EPOC-hist", "EERIE"):
        row = mpars[(mpars.dataset == key) & (mpars.metric == m)].iloc[0]
        T = mrp[(mrp.dataset == key) & (mrp.climate == "own 2026") & (mrp.metric == m)]["T"].iloc[0]
        ax.plot(RP, gev(row).return_level(RP, T), color=COLORS[key], lw=1.5, label=LABEL[key])
    for key in st.CONTROLS:
        row = mpars[(mpars.dataset == key) & (mpars.metric == m)].iloc[0]
        ax.plot(RP, gev(row).return_level(RP, 0.0), color=COLORS[key], lw=1.3, ls=(0, (4, 2)),
                label=LABEL[key])
    event_lines(ax, e_val, m, right=True)


def main():
    vz.use_style()
    t = st.load()
    seasons = pd.read_csv(DIR / "era5_seasons.csv")
    pars = pd.read_csv(DIR / "era5_gev_params.csv")
    boot = pd.read_csv(DIR / "era5_gev_boot.csv")
    mpars = pd.read_csv(DIR / "model_gev_params.csv")
    mrp = pd.read_csv(DIR / "model_event_rp.csv")
    e = st.era5(t)
    clim = st.climates(e)
    clim = {k: clim[k] for k in ("1950s", "2003", "2026")}
    e_val = {(m, ev): float(e.loc[e.year == ev, m].iloc[0]) for m in METRICS for ev in EVENTS}

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6), sharex=True,
                             gridspec_kw=dict(hspace=0.22, wspace=0.16))
    for row, (m, (name, unit)) in enumerate(METRICS.items()):
        moved_max = era5_panel(axes[row, 0], m, seasons, pars, boot, clim, e_val)
        model_panel(axes[row, 1], m, t, mpars, mrp, clim, e_val)
        # room for 2003 moved into the 2026 climate
        top = max(e_val[(m, 2026)] * 1.6, moved_max * 1.12)
        for col, ax in enumerate(axes[row]):
            ax.set_xscale("log")
            ax.set_xlim(RP[0], RP[-1])
            ax.set_ylim(0, top)
            ax.grid(True, which="major")
            ax.set_ylabel(f"{name} ({unit})" if col == 0 else "")
            ax.set_title(f"({'abcd'[2 * row + col]}) "
                         + ("ERA5, nonstationary GEV 1950–2026" if col == 0 else "simulations")
                         + f" — {name}", loc="left")
        axes[row, 1].sharey(axes[row, 0])
    for ax in axes[1]:
        ax.set_xlabel("return period (years)")
        ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000],
                      ["1", "2", "5", "10", "20", "50", "100", "200", "500", "1000"])
    axes[0, 0].legend(loc="upper left", fontsize=7.4)
    axes[0, 1].legend(loc="upper left", fontsize=7.2, ncol=1, bbox_to_anchor=(1.02, 1.0))

    fig.suptitle("How often should a summer like 2003 or 2026 occur? Return periods in observed "
                 "and simulated climates", fontsize=11.5, color=vz.INK, x=0.01, ha="left", y=0.995)
    fig.text(0.01, 0.0,
             f"land mean {hw.MEAN_BOX_LABEL}, Xu et al. (2026) heatwaves, each dataset's own 1991–2020 "
             "threshold  ·  ERA5: shift+scale GEV (location and log-scale linear in the smoothed JJA "
             "global-mean SST anomaly), 95% bootstrap bands; dots: observed seasons moved to the 2026 "
             "climate at their percentile\nMPI-GE: empirical, all members; 'at ERA5 … warming': seasons "
             f"whose forced global-mean SST is within ±{st.WL_TOL} K of ERA5's  ·  km-scale transient: "
             "nonstationary GEV at their own 2026  ·  dashed: control runs (constant forcing), stationary GEV",
             fontsize=7.4, color=vz.INK_SOFT, va="top")
    out = p2.FIG_DIR / "fig15_heatwave_return_levels.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
