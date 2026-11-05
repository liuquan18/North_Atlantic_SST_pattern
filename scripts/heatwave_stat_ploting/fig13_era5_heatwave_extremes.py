"""
Figure 13 -- ERA5 European summer heatwaves 1940-2026, and how extreme 2003
and 2026 were for their time under different ways of removing the warming.

Rows: land-mean JJA heatwave days and cumulative heat (src.heatwave.MEAN_BOX,
southern Europe). Left: every season, with the nonstationary GEV fitted to
1950-2026 (02_statistics.py) drawn year by year -- shift+scale fit (location
and log-scale follow the smoothed global-mean SST anomaly): median and
1-in-10 / 1-in-100 year levels; shift fit (constant scale): median, dashed.
1940-1949 are grey and not fitted. Right: the robustness matrix -- rank among
the fitted seasons, z (normal quantile of the season's percentile in its own
climate) and return period in its own climate of 2003 and 2026 under every
variant in era5_robustness.csv, coloured by z. The linear-in-time rows are
the conventional detrending (constant spread).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
import src.extremes as ex
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

DIR = ROOT / "data" / "heatwave_stats"
METRICS = {"hwd": ("heatwave days", "days"), "hwcum": ("cumulative heat", "°C·day")}
EVENTS = (2003, 2026)
MAIN = "shift+scale 1950-2026"
SHIFT = "shift 1950-2026"
ROW_ORDER = [
    "linear in time 1991-2026", "linear in time 1950-2026",
    "shift 1991-2026", "shift 1950-2026", "shift 1950-2026, event left out", "shift 1940-2026",
    "shift+scale 1991-2026", "shift+scale 1950-2026", "shift+scale 1950-2026, event left out",
    "shift+scale 1940-2026", "shift+scale 1950-2026, 10-yr covariate",
    "shift+scale 1950-2026, 20-yr covariate",
]


def fmt_rp(rp):
    if not np.isfinite(rp):
        return "—"
    if rp > 1e4:
        return ">10⁴ yr"
    return f"{rp:,.0f} yr" if rp >= 10 else f"{rp:.1f} yr"


def params(df, metric, variant):
    r = df[(df.metric == metric) & (df.variant == variant)].iloc[0]
    return ex.NSGEV(**{k: float(r[k]) for k in ("mu0", "mu1", "s0", "s1", "xi")})


def record_panel(ax, m, name, unit, seasons, allyrs, pars):
    fit, shift = params(pars, m, MAIN), params(pars, m, SHIFT)
    yr, T = seasons.year.values, seasons["T"].values
    ax.fill_between(yr, fit.return_level(10, T), fit.return_level(100, T), color=vz.CONTEXT,
                    alpha=0.28, lw=0, label="1-in-10 to 1-in-100 yr", zorder=1)
    ax.plot(yr, fit.return_level(2, T), color=vz.INK_SOFT, lw=1.3, label="median", zorder=2)
    ax.plot(yr, shift.return_level(2, T), color=vz.INK_SOFT, lw=1.0, ls=(0, (4, 2)),
            label="median, shift fit", zorder=2)
    early = allyrs[allyrs.year < yr.min()]
    ax.scatter(early.year, early[m], s=13, color=vz.CONTEXT, zorder=3,
               label=f"{int(early.year.min())}–{int(early.year.max())}, not fitted")
    rest = seasons[~seasons.year.isin(EVENTS)]
    ax.scatter(rest.year, rest[m], s=15, color=vz.INK_MUTED, zorder=3)
    for ev in EVENTS:
        r = seasons[seasons.year == ev].iloc[0]
        ax.scatter(ev, r[m], s=44, color=vz.INK, zorder=4)
        ax.annotate(str(ev), (ev, r[m]), xytext=(-7 if ev == 2026 else 7, 0),
                    textcoords="offset points", ha="right" if ev == 2026 else "left",
                    va="center", fontsize=8, color=vz.INK, fontweight="bold")
    # 2003 as a reference level: the recent summers that matched or beat it
    v03 = float(seasons.loc[seasons.year == 2003, m].iloc[0])
    ax.axhline(v03, color=vz.INK, lw=0.6, ls=(0, (1, 2)), zorder=0)
    n_since = int((seasons[(seasons.year > 2003) & (seasons.year < 2026)][m] >= v03).sum())
    ax.text(1941, v03, f"2003 level · reached again in {n_since} of 2004–2025",
            va="bottom", fontsize=7, color=vz.INK_SOFT)
    ax.set_xlim(1938, 2028.5)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(f"{name} ({unit})")
    ax.grid(axis="y")


def robustness_panel(ax, rob, pars, m, norm, cmap):
    r = rob[rob.metric == m].set_index(["variant", "event"])
    rows = [v for v in ROW_ORDER if (v, 2003) in r.index]
    z = np.array([[r.loc[(v, ev), "z"] for ev in EVENTS] for v in rows])
    ax.imshow(z, cmap=cmap, norm=norm, aspect="auto")
    for i, v in enumerate(rows):
        for j, ev in enumerate(EVENTS):
            c = r.loc[(v, ev)]
            txt = f"#{int(c['rank'])}   z {c.z:.1f}"
            if np.isfinite(c.rp):
                txt += f"\n{fmt_rp(c.rp)}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.8, linespacing=1.15,
                    color=vz.SURFACE if c.z > norm.vmin + 0.62 * (norm.vmax - norm.vmin) else vz.INK)
    labels = []
    aic = pars[pars.metric == m].set_index("variant").aic
    for v in rows:
        n = int(r.loc[(v, 2003), "n"])
        extra = f", AIC {aic[v]:.0f}" if v in aic.index else ""
        short = (v.replace("shift+scale 1950-2026, ", "shift+scale 1950-2026 · ")
                  .replace("-yr covariate", "-yr T").replace("event left out", "w/o event"))
        labels.append(f"{short}  (n={n}{extra})")
    ax.set_yticks(range(len(rows)), labels, fontsize=7)
    ax.set_xticks(range(len(EVENTS)), [f"JJA {e}" for e in EVENTS], fontsize=8)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.axhline(1.5, color=vz.SURFACE, lw=2.5)
    ax.axhline(5.5, color=vz.SURFACE, lw=2.5)


def main():
    vz.use_style()
    seasons = pd.read_csv(DIR / "era5_seasons.csv")
    pars = pd.read_csv(DIR / "era5_gev_params.csv")
    rob = pd.read_csv(DIR / "era5_robustness.csv")
    table = pd.read_csv(DIR / "season_table.csv")
    allyrs = table[table.dataset == "ERA5"].sort_values("year")
    cmap = vz.heatwave_cmap()
    norm = Normalize(0, max(3.0, np.ceil(rob.z.max() * 2) / 2))

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.6),
                             gridspec_kw=dict(width_ratios=[2.0, 1], hspace=0.42, wspace=0.7))
    for row, (m, (name, unit)) in enumerate(METRICS.items()):
        record_panel(axes[row, 0], m, name, unit, seasons, allyrs, pars)
        axes[row, 0].set_title(f"({'ab'[row]}) JJA {name}, ERA5, land mean {hw.MEAN_BOX_LABEL}",
                               loc="left")
        robustness_panel(axes[row, 1], rob, pars, m, norm, cmap)
        axes[row, 1].set_title(f"({'cd'[row]}) rank · z · return period in own climate",
                               loc="left", pad=20, fontsize=8.8)
    axes[0, 0].legend(loc="upper left", fontsize=7.4, ncol=2, handlelength=2.2)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    # own axes right of the matrices, so it cannot crowd panel (d)'s title
    box = axes[1, 1].get_position()
    cax = fig.add_axes([box.x1 + 0.015, box.y0 + 0.05, 0.01, axes[0, 1].get_position().y1 - box.y0 - 0.1])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("z (σ)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.suptitle("European summer heatwaves in ERA5: how extreme were 2003 and 2026 for their time?",
                 fontsize=11.5, color=vz.INK, x=0.01, ha="left", y=0.995)
    fig.text(0.01, -0.02,
             "Xu et al. (2026) heatwaves, 1991–2020 threshold  ·  GEV covariate: 4-yr trailing mean of the "
             "JJA global-mean SST anomaly  ·  shift: location linear in it; shift+scale: log-scale too  ·  "
             "linear in time: residual of a straight line in years, constant spread\n"
             "z = normal quantile of the season's percentile in its own climate (for the linear rows: "
             "residual / residual std); rank among the fitted seasons; return period in the season's own "
             "climate",
             fontsize=7.4, color=vz.INK_SOFT, va="top")
    out = p2.FIG_DIR / "fig13_era5_heatwave_extremes.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
