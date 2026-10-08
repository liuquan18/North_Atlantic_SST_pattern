"""
Figure 4 -- what the reference pattern actually is, and how unusual 2026 was.

Upper block: June, July and August 2026 separately plus the JJA mean, all with
the box mean removed. A seasonal-mean pattern is only worth correlating against
if it persists through the summer rather than being carried by one month, and
this block is what lets the reader check that.

Lower panel: the area-averaged JJA anomaly over the box for every year since
1940, against the global-ocean mean. The gap between the two curves is the part
of 2026 that is regional rather than global warming -- which is exactly what the
"global mean removed" variant of the pattern isolates.
"""
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

MONTHS = {6: "June 2026", 7: "July 2026", 8: "August 2026"}
VMAX_CAP = 2.0


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)

    monthly = p2._single_var(xr.open_dataset(
        p2.DATA_DIR / f"era5_{p2.REF_YEAR}_monthly_anom_na025.nc")).squeeze(drop=True)
    jja = p2.make_pattern(p2.load_anomaly("ERA5"), p2.load_gmsst("ERA5"),
                          "spatial").sel(year=p2.REF_YEAR)
    box = p2._single_var(xr.open_dataset(p2.RESULT_DIR / "era5_box_mean_anom.nc"))
    gm = p2.load_gmsst("ERA5")

    panels = [(p2.remove_spatial_mean(
        monthly.sel(time=monthly["time.month"] == m).squeeze(drop=True)), MONTHS[m])
        for m in (6, 7, 8)]
    panels.append((jja, "JJA 2026 mean  —  the reference pattern"))

    vmax = min(VMAX_CAP,
               max(float(np.nanpercentile(np.abs(f.values), 98)) for f, _ in panels))
    vmax = np.round(vmax * 2) / 2
    levels = np.linspace(-vmax, vmax, 21)
    cmap = vz.anomaly_cmap()
    proj = vz.map_projection()

    panel_w = 5.0
    panel_h = panel_w / vz.panel_aspect()
    ts_h = 2.7
    fig_w = 2 * panel_w + 1.3
    fig_h = 2 * panel_h + ts_h + 1.95

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(3, 3, figure=fig,
                  height_ratios=[panel_h, panel_h, ts_h],
                  width_ratios=[1, 1, 0.035],
                  hspace=0.30, wspace=0.06,
                  left=0.062, right=0.95,
                  top=1 - 1.25 / fig_h, bottom=0.62 / fig_h)

    for k, (field, label) in enumerate(panels):
        i, j = divmod(k, 2)
        ax = fig.add_subplot(gs[i, j], projection=proj)
        im = vz.draw_map(ax, field, levels, cmap,
                         labels_bottom=(i == 1), labels_left=(j == 0))
        ax.set_title(label, fontsize=9, color=vz.INK, pad=4,
                     fontweight="bold" if k == 3 else "normal")

    cax = fig.add_subplot(gs[0:2, 2])
    _pos = cax.get_position()          # a full-height bar next to two short maps
    _h = _pos.height * 0.62            # reads as a slab; shorten and centre it
    cax.set_position([_pos.x0, _pos.y0 + (_pos.height - _h) / 2, _pos.width, _h])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.set_label("SST anomaly, box mean removed (°C)", fontsize=8.2, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0)

    # --- amplitude in context ----------------------------------------------
    ax = fig.add_subplot(gs[2, :])
    ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
    ax.grid(axis="y", zorder=0)
    ax.plot(gm.year, gm, color=vz.COLORS["MPI-GE"], lw=1.5, zorder=2,
            label="global ocean mean")
    ax.plot(box.year, box, color=vz.INK, lw=1.5, zorder=3,
            label=f"North Atlantic + Mediterranean box mean ({p2.REGION_LABEL})")

    v = float(box.sel(year=p2.REF_YEAR))
    summary = json.loads((p2.RESULT_DIR / "summary.json").read_text())
    ax.plot([p2.REF_YEAR], [v], marker="o", ms=6, mfc=vz.INK, mec=vz.SURFACE,
            mew=1.3, zorder=5, clip_on=False)
    ax.annotate(f"2026: {v:+.2f} °C\n(rank {summary['era5_box_mean_rank']} of "
                f"{summary['era5_n_years']} since 1940)",
                (p2.REF_YEAR, v), textcoords="offset points", xytext=(10, 0),
                ha="left", va="center", fontsize=8.2, color=vz.INK,
                linespacing=1.35, zorder=6)
    ax.set_xlabel("year")
    ax.set_ylabel("JJA SST anomaly (°C)")
    # leave room to the right of 2026 for its label
    ax.set_xlim(int(box.year.min()), int(box.year.max()) + 12)
    ax.legend(loc="upper left", fontsize=8.3)

    fig.suptitle("The observed reference: JJA 2026 sea surface temperature over the "
                 "North Atlantic and Mediterranean",
                 fontsize=12, color=vz.INK, y=1 - 0.28 / fig_h)
    fig.text(0.5, 1 - 0.62 / fig_h,
             "ERA5, anomalies against the 1991–2020 JJA climatology  ·  "
             "July and August 2026 assembled from hourly ERA5T",
             ha="center", fontsize=8.3, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig4_era5_2026_context.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}  (colour scale ±{vmax} °C)")


if __name__ == "__main__":
    main()
