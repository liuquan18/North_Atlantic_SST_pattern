"""
Figure 1 (heatwaves) -- the European summer heatwave of each season in
figure 1, described by several quantities instead of heatwave days alone.

Columns: the seasons of figure 1 (ERA5 JJA 2026, ERA5's closest earlier
summer, and the best SST-analogue season of MPI-GE, km-scale ICON (EPOC) and
EERIE). MPI-ESM1.2-ER is left out: with no daily Tmax every one of its panels
would be empty.

Rows (Xu et al. 2026 definition, 1991-2020 reference; see
scripts/european_heatwave/):

  days         JJA days that belong to a heatwave            -- how much
  events       heatwaves with at least one day in JJA        -- how often
  longest      full length of the longest of those events    -- how persistent
  intensity    mean Tmax anomaly over the JJA heatwave days  -- how hot, typically
  peak         Tmax anomaly of the hottest heatwave day      -- how hot, at worst
  cumulative   sum of those anomalies (days x intensity)     -- total extra heat
  excess heat  sum of the anomalies above the 90th-percentile threshold
               (Perkins-Kirkpatrick & Lewis 2020)            -- heat beyond "normal hot"
  onset        first heatwave day of the May-Sep window      -- how early

Cumulative heat counts each heatwave day's anomaly from the climatological
mean (Xu et al. 2026); excess heat counts only the part above the threshold,
so it is not inflated by the ~5 degC every heatwave day must clear anyway.

Intensity, peak and onset are undefined where there was no heatwave; those
land cells are grey, not blank.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.pattern_2026 as p2
import src.viz2026 as vz
from fig1_patterns import heatwave_season, land_mean, season

COLUMNS = [("ERA5", "ref"), ("ERA5", "best"), ("MPI-GE", "best"),
           ("ICON-EPOC-hist", "best"), ("EERIE", "best")]

# first day of each month on the 365-day calendar, for the onset colour bar
MONTH_DOY = {"1 May": 121, "1 Jun": 152, "1 Jul": 182, "1 Aug": 213, "1 Sep": 244, "30 Sep": 273}


def quantity_rows():
    """(name, label, units, function Dataset -> field, levels, cmap, extend, land_under, fmt)."""
    amp = vz.heatwave_cmap()
    return [
        ("days", "heatwave days", "days",
         lambda s: s.hwd, np.arange(0, 46, 5), amp, "max", False, "{:.1f} days"),
        ("events", "number of heatwaves", "events",
         lambda s: s.hwn, np.arange(-0.5, 7.5, 1), amp, "max", False, "{:.1f}"),
        ("longest", "longest heatwave", "days",
         lambda s: s.hwmax, np.array([0, 3, 5, 7, 10, 14, 21, 28, 35]), amp, "max", False,
         "{:.1f} days"),
        ("intensity", "mean heatwave intensity", "Tmax anomaly (°C)",
         lambda s: (s.hwcum / s.hwd).where(s.hwd > 0), np.arange(3, 11.1, 1), amp, "both", True,
         "{:.1f} °C"),
        ("peak", "peak heatwave intensity", "Tmax anomaly (°C)",
         lambda s: s.hwpeak, np.arange(4, 16.1, 1.5), amp, "both", True, "{:.1f} °C"),
        ("cumulative", "cumulative heat", "°C · day",
         lambda s: s.hwcum, np.array([0, 10, 25, 50, 75, 100, 150, 200, 250, 300]), amp, "max",
         False, "{:.0f} °C·day"),
        ("excess", "excess heat\n(above threshold)", "°C · day",
         lambda s: s.hwexcess, np.array([0, 5, 10, 20, 30, 40, 60, 80, 100]), amp, "max",
         False, "{:.0f} °C·day"),
        ("onset", "heatwave onset", "first heatwave day",
         lambda s: s.onset, np.arange(121, 275, 7), plt.get_cmap("viridis"), "neither", True,
         None),
    ]


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(p2.RESULT_DIR / "best_analogues.csv")

    cols = []
    for key, which in COLUMNS:
        year, mem, r = season(key, which, best)
        hws = heatwave_season(key, year, mem)
        if hws is None:
            print(f"!! {key} JJA {year}: no heatwave metrics, column skipped")
            continue
        cols.append(dict(key=key, which=which, year=year, member=mem, r=r, hw=hws))

    rows = quantity_rows()
    n, m = len(cols), len(rows)
    panel_w = 2.3
    panel_h = panel_w / vz.panel_aspect(vz.EUROPE_EXTENT)
    fig_w = panel_w * n + 2.3
    fig_h = panel_h * m + 1.9

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(m, n + 1, figure=fig, width_ratios=[1] * n + [0.04],
                  hspace=0.10, wspace=0.06,
                  left=0.085, right=0.93, top=1 - 1.35 / fig_h, bottom=0.45 / fig_h)

    for i, (name, label, units, fn, levels, cmap, extend, under, fmt) in enumerate(rows):
        norm = BoundaryNorm(levels, cmap.N, extend=extend)
        im = None
        for j, c in enumerate(cols):
            ax = fig.add_subplot(gs[i, j], projection=vz.europe_projection())
            field = fn(c["hw"])
            im = vz.draw_land_map(ax, field, cmap, norm, land_under=under,
                                  labels_bottom=(i == m - 1), labels_left=(j == 0)) or im

            if i == 0:
                d = p2.DATASETS[c["key"]]
                if c["key"] == "ERA5":
                    sub = "reference season" if c["which"] == "ref" else "closest earlier summer"
                else:
                    sub = d.resolution
                when = f"JJA {c['year']}" + (f" · r{c['member']}" if c["member"] else "")
                rtxt = "" if c["r"] is None else f"  (SST r = {c['r']:+.2f})"
                ax.set_title(f"{d.label}\n{sub}\n{when}{rtxt}", pad=5, fontsize=8.6,
                             color=vz.INK, linespacing=1.3)

            if fmt is not None:
                vz.panel_tag(ax, "land mean " + fmt.format(land_mean(field)), loc="upper right")
            if j == 0:
                ax.text(-0.24, 0.5, label, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=9.2, color=vz.INK)

        cax = fig.add_subplot(gs[i, -1])
        if name == "onset":
            cb = fig.colorbar(im, cax=cax, ticks=list(MONTH_DOY.values()))
            cb.ax.set_yticklabels(list(MONTH_DOY.keys()))
        elif name == "events":
            cb = fig.colorbar(im, cax=cax, extend=extend, ticks=np.arange(0, 8))
        else:
            cb = fig.colorbar(im, cax=cax, extend=extend, ticks=levels)
        cb.set_label(units, fontsize=7.8, color=vz.INK_SOFT)
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_linewidth(0)

    fig.suptitle("European summer heatwaves in JJA 2026 and in the closest SST-analogue "
                 "season of each dataset", fontsize=12, color=vz.INK, y=1 - 0.30 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             "heatwave: ≥3 consecutive days with Tmax anomaly above the calendar-day 90th "
             "percentile (15-day window), detected May–Sep, land only (Xu et al. 2026)  ·  "
             "1991–2020 reference of the same dataset  ·  grey land: no heatwave  ·  "
             "MPI-ESM1.2-ER omitted (no daily output)",
             ha="center", fontsize=8, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig1_heatwaves_2026.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
