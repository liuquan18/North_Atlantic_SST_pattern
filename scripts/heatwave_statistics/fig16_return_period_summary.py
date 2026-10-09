"""
Figure 16 -- every return-period estimate of the 2003 and 2026 summers in one place.

Left: the ERA5 2003 value, right: the ERA5 2026 value (land-mean JJA
heatwave days filled, cumulative heat hollow), with 95% intervals
(bootstrap for the GEV fits, Clopper-Pearson for empirical counts).
"own year" = the climate of the event's own year: for ERA5 its smoothed
global-mean SST, for MPI-GE the seasons at ERA5's warming of that year or
within +-10 calendar years, for the km-scale transient runs their own
forced state of that year. Rows are grouped by source; controls are
unforced, so they give the return period without any warming (EPOC control
from daily mean temperature). An arrow marks an interval that runs past the
axis, an open-ended interval (no simulated season reached the value) or a
value beyond the fitted upper bound of the GEV (return period infinite).
All values from 02_statistics.py (era5_event_rp.csv, model_event_rp.csv).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

DIR = ROOT / "data" / "heatwave_stats"
EVENTS = (2003, 2026)
METRICS = {"hwd": ("heatwave days", 0.0), "hwcum": ("cumulative heat", 0.22)}
XMIN, XMAX = 1, 1e5
COLORS = {**vz.COLORS, "ICON-sap-ctrl": vz.COLORS["MPI-ER"]}
MAIN, SHIFT = "shift+scale 1950-2026", "shift 1950-2026"


def rows_for(ev):
    """(group, label, source, selector) per figure row, for event `ev`."""
    own = str(ev)
    return [
        ("ERA5", "shift+scale GEV · own year", "era5", dict(variant=MAIN, climate=own)),
        ("ERA5", "shift GEV · own year", "era5", dict(variant=SHIFT, climate=own)),
        ("ERA5", "shift+scale · event left out", "era5",
         dict(variant=f"{MAIN}, event left out", climate=own)),
        ("ERA5", "shift+scale · 20-yr covariate", "era5",
         dict(variant="shift+scale 1950-2026, 20-yr covariate", climate=own)),
        ("ERA5", "shift+scale · 1991–2026 fit", "era5",
         dict(variant="shift+scale 1991-2026", climate=own)),
        ("ERA5", "shift+scale · 2026 climate", "era5", dict(variant=MAIN, climate="2026")),
        ("ERA5", "shift+scale · 1950s climate", "era5", dict(variant=MAIN, climate="1950s")),
        ("ERA5", "shift · 1950s climate", "era5", dict(variant=SHIFT, climate="1950s")),
        ("MPI-GE", "at ERA5 warming of own year", "model",
         dict(dataset="MPI-GE", method="empirical", climate=f"ERA5 {ev} warming level")),
        ("MPI-GE", "own year ±10", "model",
         dict(dataset="MPI-GE", method="empirical", climate=f"{ev}±10 (calendar)")),
        ("MPI-GE", "at ERA5 2026 warming", "model",
         dict(dataset="MPI-GE", method="empirical", climate="ERA5 2026 warming level")),
        ("MPI-GE", "1850–1900", "model",
         dict(dataset="MPI-GE", method="empirical", climate="1850-1900")),
        ("MPI-GE", "NS-GEV · 1850–1900", "model",
         dict(dataset="MPI-GE", method="NS-GEV", climate="1850-1900")),
        ("ICON-EPOC-hist", "EPOC hist · NS-GEV · own year", "model",
         dict(dataset="ICON-EPOC-hist", method="NS-GEV", climate=f"own {ev}")),
        ("EERIE", "EERIE hist · NS-GEV · own year", "model",
         dict(dataset="EERIE", method="NS-GEV", climate=f"own {ev}")),
        ("EERIE-ctrl", "EERIE control · GEV", "model",
         dict(dataset="EERIE-ctrl", method="GEV", climate="control")),
        ("ICON-sap-ctrl", "Sapphire control · GEV", "model",
         dict(dataset="ICON-sap-ctrl", method="GEV", climate="control")),
        ("ICON-EPOC-ctrl", "EPOC control (Tmean) · GEV", "model",
         dict(dataset="ICON-EPOC-ctrl", method="GEV", climate="control")),
    ]


def pick(df, metric, ev, sel):
    q = (df.metric == metric) & (df.event == ev)
    if "quantity" in df:
        q &= df.quantity == "rp"
    for k, v in sel.items():
        q &= df[k] == v
    r = df[q]
    return None if r.empty else r.iloc[0]


def draw(ax, ev, era5, model):
    rows = rows_for(ev)
    y = np.arange(len(rows))[::-1].astype(float)
    groups = [r[0] for r in rows]
    for i, (group, label, src, sel) in enumerate(rows):
        color = vz.INK if group == "ERA5" else COLORS[group]
        for m, (_, dy) in METRICS.items():
            r = pick(era5 if src == "era5" else model, m, ev, sel)
            if r is None:
                continue
            yy = y[i] - dy
            est, lo, hi = r.est, r.lo, r.hi
            lo_d, hi_d = np.clip([lo, hi], XMIN, XMAX)
            ax.plot([lo_d, hi_d], [yy, yy], color=color, lw=1.3, solid_capstyle="butt")
            if not np.isfinite(hi) or hi >= XMAX:
                ax.annotate("", (XMAX * 1.25, yy), (XMAX * 0.92, yy),
                            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1))
            filled = m == "hwd"
            if np.isfinite(est) and est <= XMAX:
                ax.scatter(est, yy, s=28, zorder=3, color=color if filled else vz.SURFACE,
                           edgecolor=color, lw=1.1)
            elif not np.isfinite(est):
                ax.text(XMAX * 1.4, yy, "∞", va="center", fontsize=7.5, color=color)
            else:
                # finite but off the axis: the arrow carries it, say how far
                ax.text(XMAX * 1.4, yy, f"{est:.0e}".replace("e+0", "e"), va="center",
                        fontsize=6.5, color=color)
    # group separators and labels
    for i in range(1, len(rows)):
        if groups[i] != groups[i - 1] and groups[i - 1] in ("ERA5", "MPI-GE", "EERIE"):
            ax.axhline((y[i] + y[i - 1]) / 2, color=vz.CONTEXT, lw=0.6)
    ax.set_yticks(y, [r[1] for r in rows], fontsize=7.4)
    ax.set_xscale("log")
    ax.set_xlim(XMIN, XMAX * 2.2)
    ax.set_ylim(y.min() - 0.7, y.max() + 0.6)
    ax.grid(axis="x", which="major")
    ax.tick_params(axis="y", length=0)
    return y, rows


def main():
    vz.use_style()
    era5 = pd.read_csv(DIR / "era5_event_rp.csv")
    model = pd.read_csv(DIR / "model_event_rp.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 7.8), sharey=False,
                             gridspec_kw=dict(wspace=0.04))
    for k, (ax, ev) in enumerate(zip(axes, EVENTS)):
        y, rows = draw(ax, ev, era5, model)
        vals = {m: float(era5[(era5.metric == m) & (era5.event == ev)].value.iloc[0])
                for m in METRICS}
        ax.set_title(f"({'ab'[k]}) ERA5 JJA {ev}: {vals['hwd']:.1f} heatwave days, "
                     f"{vals['hwcum']:.0f} °C·day", loc="left")
        ax.set_xlabel("return period (years)")
        ax.set_xticks([1, 10, 100, 1000, 10000, 100000], ["1", "10", "100", "1,000", "10⁴", "10⁵"])
        if k == 1:
            ax.set_yticklabels([])
    # group labels on the left
    groups = {"ERA5": "ERA5", "MPI-GE": "MPI-GE (50 members)", "ICON-EPOC-hist": "km-scale",
              "EERIE-ctrl": "km-scale controls"}
    rows = rows_for(EVENTS[0])
    yy = np.arange(len(rows))[::-1]
    for g, lab in groups.items():
        idx = [i for i, r in enumerate(rows) if r[0] == g]
        axes[0].text(-0.5, yy[idx[0]] + 0.45, lab, transform=axes[0].get_yaxis_transform(),
                     fontsize=8, fontweight="bold", color=vz.INK, ha="left")
    h = [plt.Line2D([], [], marker="o", ls="", color=vz.INK_SOFT, label="heatwave days"),
         plt.Line2D([], [], marker="o", ls="", markerfacecolor=vz.SURFACE, color=vz.INK_SOFT,
                    label="cumulative heat")]
    axes[0].legend(handles=h, loc="lower right", fontsize=7.6)

    fig.suptitle("Return period of the observed 2003 and 2026 summers, by source and climate",
                 fontsize=11.5, color=vz.INK, x=0.01, ha="left", y=0.995)
    fig.text(0.01, 0.0,
             f"land mean {hw.MEAN_BOX_LABEL}, Xu et al. (2026) heatwaves, each dataset's own 1991–2020 "
             "threshold  ·  95% intervals: bootstrap (GEV), Clopper–Pearson (empirical)  ·  "
             "arrow: interval open or beyond the axis; ∞: value above the fitted GEV's upper bound\n"
             "own year: ERA5 at its smoothed global-mean SST of that year; MPI-GE at the same warming or "
             "±10 calendar years; km-scale transient runs at their own state of that year (EPOC 2026: "
             "its 2025)",
             fontsize=7.4, color=vz.INK_SOFT, va="top")
    out = p2.FIG_DIR / "fig16_return_period_summary.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
