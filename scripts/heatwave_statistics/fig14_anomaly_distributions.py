"""
Figure 14 -- can the simulations produce a summer as anomalous as 2003 or 2026?

Distributions of the season-to-season anomaly from each dataset's forced
mean (02_statistics.py part B: MPI-GE minus its 50-member mean; single
transient runs minus a least-squares fit on the smoothed global-mean SST;
control runs minus a linear drift), for land-mean JJA heatwave days (top)
and cumulative heat (bottom). Left: in days / degC*day; right: in units of
each sample's own standard deviation, which takes out the differences in
variability between datasets (ERA5's anomalies there are divided by ERA5's
local spread, the std within +-15 years of the event). Horizontal lines:
ERA5's 2026 (solid) and 2003 (dashed) anomaly. Above each distribution: its
size and standard deviation -- the models' spread against ERA5's. Under each distribution: the share of its seasons at or
above each line. MPI-GE is shown in 1850-1900, pooled over 1991-2026 and in
the seasons at ERA5's 2003 and 2026 warming level (forced global-mean SST within 0.05 K),
because its spread grows with warming. Controls are hatched (unforced).
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
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

st = importlib.import_module("02_statistics")

DIR = ROOT / "data" / "heatwave_stats"
METRICS = {"hwd": ("heatwave days", "days"), "hwcum": ("cumulative heat", "°C·day")}
COLORS = {**vz.COLORS, "ICON-sap-ctrl": vz.COLORS["MPI-ER"], "ERA5": vz.INK_MUTED}
SHORT = {
    "ERA5": "ERA5", "MPI-GE": "MPI-GE", "ICON-EPOC-hist": "EPOC\nhist", "EERIE": "EERIE\nhist",
    "EERIE-ctrl": "EERIE\nctrl", "ICON-sap-ctrl": "Sapphire\nctrl", "ICON-EPOC-ctrl": "EPOC ctrl\n(Tmean)",
}


def sample_label(key, label):
    if key == "MPI-GE" and "warming level" in label:
        return f"MPI-GE\nat ERA5\n{label.split()[1]}\nwarming"
    if key == "ERA5":
        return f"ERA5\n{label}"
    if key == "MPI-GE":
        return f"MPI-GE\n{label}"
    return SHORT[key]


def main():
    vz.use_style()
    t = st.load()
    a = pd.read_csv(DIR / "season_anomalies.csv")
    exc = pd.read_csv(DIR / "anomaly_exceedance.csv")
    e_clim = st.climates(st.era5(t))
    samp = st.samples(a, e_clim)
    e = a[a.dataset == "ERA5"].set_index("year")

    fig, axes = plt.subplots(2, 2, figsize=(18.5, 9.6),
                             gridspec_kw=dict(hspace=0.62, wspace=0.14))
    for row, (m, (name, unit)) in enumerate(METRICS.items()):
        e_row = exc[(exc.metric == m) & (exc.dataset == "ERA5")].iloc[0]
        sd_e = {ev: e_row[f"{ev}_era5_local_std"] for ev in (2003, 2026)}
        for col, how in enumerate(("abs", "std")):
            ax = axes[row, col]
            for i, ((key, label), g) in enumerate(samp.items()):
                s = g[f"{m}_anom"].dropna().values
                sd = s.std(ddof=1)
                v = s if how == "abs" else s / sd
                c = COLORS[key]
                parts = ax.violinplot(v, positions=[i], widths=0.8, showextrema=False)
                for b in parts["bodies"]:
                    b.set_facecolor(c)
                    b.set_edgecolor(c)
                    b.set_alpha(0.32)
                    if "ctrl" in key:
                        b.set_hatch("////")
                        b.set_linestyle((0, (4, 2)))
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                ax.vlines(i, q1, q3, color=c, lw=3.2, zorder=3)
                ax.scatter(i, med, color=vz.SURFACE, s=9, zorder=4)
                if len(v) <= 200:
                    jit = np.random.default_rng(i).uniform(-0.13, 0.13, len(v))
                    ax.scatter(i + jit, v, s=4, color=c, alpha=0.75, lw=0, zorder=2)
                r = exc[(exc.metric == m) & (exc.dataset == key) & (exc["sample"] == label)].iloc[0]
                txt = "\n".join(f"’{str(ev)[2:]}: {100 * r[f'{ev}_{how}_p']:.1f}%" for ev in (2003, 2026))
                ax.text(i, -0.02, txt, transform=ax.get_xaxis_transform(), ha="center", va="top",
                        fontsize=6.4, color=vz.INK_SOFT)
                ax.text(i, 1.0, f"n={len(v)}\nσ={sd:.1f}" if how == "abs" else f"n={len(v)}",
                        transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                        fontsize=6.2, color=vz.INK_MUTED, linespacing=1.1)
            ys = {ev: e.loc[ev, f"{m}_anom"] / (1 if how == "abs" else sd_e[ev]) for ev in (2003, 2026)}
            for ev, ls in ((2026, "-"), (2003, (0, (4, 2)))):
                y = ys[ev]
                ax.axhline(y, color=vz.INK, ls=ls, lw=1.1, zorder=1)
                # the label of the higher line above it, of the lower one below
                up = y >= ys[2003 if ev == 2026 else 2026]
                ax.text(len(samp) - 0.45, y, f" ERA5 {ev} ({y:.1f})", va="bottom" if up else "top",
                        fontsize=7.2, color=vz.INK)
            ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=0)
            ax.set_xticks(range(len(samp)),
                          [sample_label(k, lab) for (k, lab) in samp], fontsize=7)
            ax.tick_params(axis="x", length=0, pad=24)
            ax.set_xlim(-0.6, len(samp) + 0.5)
            ax.set_ylabel(f"anomaly ({unit})" if how == "abs" else "anomaly / sample std (σ)")
            ax.grid(axis="y")
            ax.set_title(f"({'abcd'[2 * row + col]}) {name}, "
                         + ("absolute" if how == "abs" else "standardised"), loc="left", pad=20)

    fig.suptitle("Season-to-season anomalies of European summer heatwaves, forced signal removed",
                 fontsize=11.5, color=vz.INK, x=0.01, ha="left", y=1.0)
    fig.text(0.01, -0.03,
             f"land mean {hw.MEAN_BOX_LABEL}, Xu et al. (2026) heatwaves with each dataset's own 1991–2020 "
             "threshold  ·  forced part: MPI-GE 50-member mean; ERA5 (1950–2026), EPOC and EERIE members: "
             "regression on smoothed global-mean SST; controls: linear drift  ·  standardised: ERA5 "
             "divided by its std within ±15 yr of the event\n"
             "percentages: share of the sample's seasons at or above ERA5's 2003 / 2026 anomaly  ·  "
             "box: interquartile range, dot: median  ·  hatched: control runs (constant forcing)  ·  "
             "EPOC control from daily mean temperature",
             fontsize=7.4, color=vz.INK_SOFT, va="top")
    out = p2.FIG_DIR / "fig14_heatwave_anomaly_distributions.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
