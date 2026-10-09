"""
Figure 5 -- what the km-scale runs actually resolve.

The pattern correlations in figures 1-3 are computed on a 1 deg grid, which is
deliberately blind to the difference between a 5 km ocean and a 1 deg one. This
figure shows what that common grid throws away, in the two places where ocean
resolution is known to matter for North Atlantic + Mediterranean SST: the Gulf
Stream / North Atlantic Current separation, and the Mediterranean basin.

Every dataset is drawn on the same 0.05 deg canvas, so the blockiness of a panel
is that model's own resolution rather than a plotting choice. The point of the
EERIE column is the one that is easy to miss: ICON-ESM-ER runs a 5 km ocean, but
the CMORized archive used here is stored at 0.25 deg, so its panel shows what is
*available*, not what the model resolved.

Colour ramp: cmocean "thermal", the oceanographic standard for absolute SST. It
is multi-hue but strictly monotonic in lightness -- the property that the
no-rainbow rule exists to protect -- so fronts stay readable in greyscale and
for colour-vision-deficient readers.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.gridspec import GridSpec

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

ZOOMS = [
    ("gulfstream", "Gulf Stream extension", [-75, -35, 33, 48]),
    ("medsea", "Mediterranean", [-2, 30, 33, 45]),
]
COLS = [
    ("ERA5", "ERA5", "0.25° reanalysis"),
    ("MPI-GE", "MPI-ESM1-2-LR", "~1° ocean"),
    ("EERIE", "EERIE ICON-ESM-ER", "5 km ocean, archived at 0.25°"),
    ("ICON-EPOC", "km-scale ICON (EPOC)", "5 km ocean, native"),
]


def load_zoom(name, zoom):
    f = p2.DATA_DIR / f"zoom_{zoom}_{name}.nc"
    if not f.exists():
        return None
    da = p2._single_var(xr.open_dataset(f)).squeeze(drop=True)
    if float(da.max()) > 100:      # ERA5 arrives in kelvin
        da = da - 273.15
    return da


def main():
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)

    try:
        import cmocean
        cmap = cmocean.cm.thermal
    except ImportError:
        cmap = plt.get_cmap("magma")

    fields = {(c[0], z[0]): load_zoom(c[0], z[0]) for c in COLS for z in ZOOMS}
    missing = [k for k, v in fields.items() if v is None]
    if missing:
        raise SystemExit(f"missing zoom files for {missing}; run scripts/pre_process/06_kmscale_zoom.sh first")

    panel_w = 3.0
    fig_w = panel_w * len(COLS) + 1.8
    heights = [panel_w / vz.panel_aspect(z[2]) for z in ZOOMS]
    fig_h = sum(heights) + 2.3

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(len(ZOOMS), len(COLS) + 1, figure=fig,
                  width_ratios=[1] * len(COLS) + [0.04],
                  height_ratios=heights, hspace=0.30, wspace=0.06,
                  left=0.088, right=0.95,
                  top=1 - 1.45 / fig_h, bottom=0.45 / fig_h)

    proj = vz.map_projection()
    ims = []
    for i, (zkey, zlabel, extent) in enumerate(ZOOMS):
        vals = np.concatenate([fields[(c[0], zkey)].values.ravel() for c in COLS])
        vals = vals[np.isfinite(vals)]
        lo, hi = np.percentile(vals, [1, 99])
        levels = np.linspace(np.floor(lo), np.ceil(hi), 25)

        for j, (ckey, clabel, cres) in enumerate(COLS):
            ax = fig.add_subplot(gs[i, j], projection=proj)
            im = vz.draw_map(ax, fields[(ckey, zkey)], levels, cmap, extent=extent,
                             labels_bottom=True, labels_left=(j == 0),
                             lon_step=10, lat_step=5)
            if i == 0:
                ax.set_title(f"{clabel}\n{cres}", fontsize=8.8, color=vz.INK, pad=5)
            if j == 0:
                ax.text(-0.115, 0.5, zlabel, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=8.6, color=vz.INK)

        cax = fig.add_subplot(gs[i, len(COLS)])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7)
        cb.outline.set_linewidth(0)
        cb.set_label("°C", fontsize=8, color=vz.INK_SOFT)
        ims.append(im)

    fig.suptitle("What the km-scale ocean resolves that the pattern correlation cannot see",
                 fontsize=12, color=vz.INK, y=1 - 0.30 / fig_h)
    fig.text(0.5, 1 - 0.72 / fig_h,
             "JJA 1991–2020 mean sea surface temperature, every dataset conservatively remapped onto "
             "the same 0.05° canvas  ·  blockiness is the model's own resolution",
             ha="center", fontsize=8.3, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig5_kmscale_resolution.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
