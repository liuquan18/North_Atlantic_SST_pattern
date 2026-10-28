"""
Figure 25 -- the SST and Z300 pattern of the summers that sit in the top
decile of the MPI-GE summers of their climate, in MPI-GE and in each of the
higher-resolution runs.

A summer is selected when its index reaches the 90th percentile of the 500
MPI-GE summers in the 10 years around it (y-5 to y+4, as the q90 line of
fig 22/23), over YEARS (1991-2030). Two indices, one figure each (argument
`hwd` or `boxdiff`, both without):

  hwd      June-August heatwave days over land in hw.MEAN_BOX
           (data/heatwave_stats/season_table.csv; fig23_epoc_hwd_corr_mpige.py)
  boxdiff  the box-difference index MED - NAL
           (data/sst_pattern_corr/box_difference.nc; fig23_boxdiff_mpige.py)

Panel (a) is the composite of the MPI-GE summers so selected (~10% of its
2000), the others one composite per run of fig 23's RUNS, over that run's
selected summers. Shading: JJA SST anomaly with its area-weighted box mean
removed (the pattern of fig 1), averaged over the summers; contours: JJA 300
hPa geopotential height anomaly, its box mean removed summer by summer, then
averaged (data/sst_anomaly/<dataset>_jja_z300_anom_1deg.nc,
pre_process/10_model_zg.sh). Runs without upper-air output over these years
(EERIE r2 and r3, MPI-ER r3) get the SST alone. The tag is the centred pattern
correlation of the SST composite with ERA5 JJA 2026 over the box.
"""
import importlib
import sys
import textwrap
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
fig23 = importlib.import_module("fig23_epoc_hwd_corr_mpige")
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

VARIANT = p2.PLOT_VARIANTS[0]
YEARS = fig23.YEARS
QUANTILE = "q90"
ZLEVEL = 300
#: contour interval of the Z300 composites (gpm); no zero line
ZG_STEP = 10
VMAX_CAP = 1.5
NCOL = 2
#: SST and Z300 anomaly file stem of each dataset
ZG_STEM = {"MPI-GE": "mpige", "ICON-EPOC-hist": "epoc2_020",
           "EERIE": "eerie_hist_ssp245", "MPI-ER": "mpier"}
INDICES = {
    "hwd": dict(label=f"number of heatwave days over land in {hw.MEAN_BOX_LABEL} reaches",
                short="heatwave days", boxes=[hw.MEAN_BOX]),
    "boxdiff": dict(label=f"box-difference index {p2.DIPOLE_LABEL} reaches",
                    short=p2.DIPOLE_LABEL, boxes=list(p2.DIPOLE_BOXES.values())),
}


def selection(index):
    """
    The selected summers: a list of (key, member, [years]) with MPI-GE first
    (its years per member), then every RUNS member that has the index, and the
    MPI-GE window levels.
    """
    years = list(range(YEARS[0], YEARS[1] + 1))
    if index == "hwd":
        mpige, _, runs, _ = fig23.load()
    else:
        with xr.open_dataset(p2.dipole_file()) as ds:
            mpige = ds["MPI-GE"].load().rename({"member_MPI-GE": "member"})
            runs = {(key, m): fig23.run_series(ds[key].load(), m)
                    for key, (_, members) in fig23.RUNS.items() for m in members}
    q = fig22.window_levels(mpige, years).set_index("year")[QUANTILE]
    hit = mpige.sel(year=years) >= xr.DataArray(q.values, coords=dict(year=years))
    out = [("MPI-GE", m, [int(y) for y in hit.year[hit.sel(member=m).values].values])
           for m in hit.member.values]
    for (key, m), x in runs.items():
        x = x[x.index.isin(years)]
        out.append((key, int(m), [int(y) for y in x.index if x[y] >= q[y]]))
    return out, q


def sst_pattern(key, member, years):
    """Box-mean-removed JJA SST anomaly of each selected summer (dim `year`)."""
    da = p2.load_anomaly(key)
    if "member" in da.dims:
        da = da.sel(member=member)
    da = da.sel(year=years).load()
    return p2.make_pattern(da, None, VARIANT)


def zg_pattern(key, member, years):
    """
    Z300 anomaly of each selected summer minus its area-weighted mean over the
    analysis box (as fig1_patterns_era5.zg_field), or None without upper-air
    output for that member and those years.
    """
    f = p2.DATA_DIR / f"{ZG_STEM[key]}_jja_z{ZLEVEL}_anom_1deg.nc"
    if not f.exists():
        return None
    with xr.open_dataset(f) as ds:
        z = ds[f"z{ZLEVEL}"]
        if member not in z.member or not set(years) <= set(z.year.values):
            return None
        z = z.sel(member=member, year=years).load()
    if z.isnull().all(("lat", "lon")).any():
        return None
    box = z.sel(lon=slice(*p2.REGION["lon"]), lat=slice(*p2.REGION["lat"]))
    w = np.cos(np.deg2rad(box.lat)) * xr.ones_like(box.lon)
    return z - (box * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))


def composites(sel):
    """One panel dict per run, MPI-GE first; MPI-GE pools all its members."""
    ref = sst_pattern("ERA5", None, [p2.REF_YEAR]).squeeze("year")
    panels = []
    ge = [(m, ys) for key, m, ys in sel if key == "MPI-GE" and ys]
    runs = [("MPI-GE", None, ge)] + [(key, m, [(m, ys)] if ys else [])
                                    for key, m, ys in sel if key != "MPI-GE"]
    for key, m, parts in runs:
        n = sum(len(ys) for _, ys in parts)
        p = dict(key=key, member=m, n=n, years=[y for _, ys in parts for y in ys])
        if n:
            ssts = [sst_pattern(key, mm, ys) for mm, ys in parts]
            zgs = [zg_pattern(key, mm, ys) for mm, ys in parts]
            p["sst"] = xr.concat(ssts, "year").mean("year")
            p["zg"] = None if any(z is None for z in zgs) else xr.concat(zgs, "year").mean("year")
            p["r"] = float(p2.pattern_corr(p2.coarsen_to_1deg(p["sst"]),
                                           p2.coarsen_to_1deg(ref), centered=True))
        panels.append(p)
    return panels


def title(p):
    """Panel title: model, member, resolution."""
    if p["key"] == "MPI-GE":
        return f"MPI-ESM1.2-LR (MPI-GE), 50 members · {p2.DATASETS['MPI-GE'].resolution}"
    name = fig23.RUNS[p["key"]][0]
    return f"{name} r{p['member']} · {p2.DATASETS[p['key']].resolution}"


def year_tag(p):
    """Number of composited summers, and which ones for a single run."""
    if p["key"] == "MPI-GE":
        return f"n = {p['n']} of {50 * (YEARS[1] - YEARS[0] + 1)} summers"
    if not p["n"]:
        return ""
    return f"n = {p['n']}: " + ", ".join(map(str, p["years"]))


def figure(index, panels, out):
    vz.use_style()
    meta = INDICES[index]
    have = [p for p in panels if p["n"]]
    vmax = min(VMAX_CAP, max(float(np.nanpercentile(np.abs(p["sst"].values), 98)) for p in have))
    vmax = max(0.5, np.round(vmax * 2) / 2)
    levels = np.linspace(-vmax, vmax, 21)
    cmap = vz.anomaly_cmap()

    no_zg = [f"{fig23.RUNS[p['key']][0]} r{p['member']}" for p in have if p["zg"] is None]
    caption = (
        f"Summers ({YEARS[0]}–{YEARS[1]}) whose {meta['label']} the 90th percentile of "
        "the 500 MPI-ESM1.2-LR summers in the 10 years around them (y−5 to y+4; the q90 line of "
        "Fig. 23). (a) Composite of the MPI-GE summers so selected, all 50 members (historical + "
        "SSP2-4.5); the other panels: the selected summers of each higher-resolution historical + "
        "scenario run, one panel per realisation, years listed (ICON-EPOC: SSP5-8.5 from 2015, "
        "to 2025; ICON-EERIE: SSP2-4.5 from 2015, r1 to 2030, r2 SST to 2014, r3 to 2020; "
        "MPI-ESM1.2-ER: SSP5-8.5 from 2015). Shading: JJA SST anomaly against 1991–2020 of the "
        f"same dataset, area-weighted mean over {p2.REGION_LABEL} removed, averaged over the "
        f"summers. Contours: JJA {ZLEVEL} hPa geopotential height anomaly, treated the same way "
        f"(box mean removed summer by summer), every {ZG_STEP} gpm, negative dashed, zero "
        "omitted" + (f"; none where the run has no upper-air output for these years "
                     f"({', '.join(no_zg)})" if no_zg else "")
        + f". Dashed box{'es' if len(meta['boxes']) > 1 else ''}: where "
        f"the index is taken. r: centred pattern correlation of the SST composite with ERA5 JJA "
        f"{p2.REF_YEAR} over the box (1° grid).")
    missing = [f"{fig23.RUNS[k][0]} r{m}" for k, (_, ms) in fig23.RUNS.items() for m in ms
               if not any(p["key"] == k and p["member"] == m for p in panels)]
    if missing:
        caption += f" Not shown, no {meta['short']}: {', '.join(missing)}."
    text = textwrap.fill(caption, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)

    nrow = int(np.ceil(len(panels) / NCOL))
    width = fig18.FIG_SIZE[0]
    panel_w = (width - 0.75) / NCOL
    panel_h = panel_w / vz.panel_aspect()
    row_in = panel_h + 0.3
    cbar_in = 0.85
    top = 0.995
    height = (caption_in + cbar_in + nrow * row_in + 0.1) / top
    fig = plt.figure(figsize=(width, height))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(nrow, NCOL, left=0.06, right=0.985, top=top,
                          bottom=(caption_in + cbar_in) / height, hspace=0.3, wspace=0.08)

    im = None
    for i, (p, letter) in enumerate(zip(panels, "abcdefghij")):
        r, c = divmod(i, NCOL)
        ax = fig.add_subplot(gs[r, c], projection=vz.map_projection())
        last_row = r == nrow - 1 or i + NCOL >= len(panels)
        if p["n"]:
            im = vz.draw_map(ax, p["sst"], levels, cmap, labels_bottom=last_row,
                             labels_left=c == 0)
            if p["zg"] is not None:
                z = p["zg"]
                zmax = ZG_STEP * np.ceil(float(np.abs(z).max()) / ZG_STEP)
                zl = [v for v in np.arange(-zmax, zmax + 1, ZG_STEP) if v != 0]
                cs = ax.contour(z.lon, z.lat, z.values, levels=zl, colors=vz.INK,
                                linewidths=0.8, linestyles=["--" if v < 0 else "-" for v in zl],
                                transform=ccrs.PlateCarree(), zorder=5)
                ax.clabel(cs, fmt="%d", fontsize=5.5, inline_spacing=1)
            vz.panel_tag(ax, f"r = {p['r']:+.2f}", loc="upper right")
            vz.panel_tag(ax, year_tag(p), loc="lower left")
        else:
            vz.draw_map(ax, xr.full_like(panels[0]["sst"], np.nan), levels, cmap,
                        labels_bottom=last_row, labels_left=c == 0)
            ax.text(0.5, 0.5, f"no summer at or above\nthe MPI-GE 90th percentile",
                    transform=ax.transAxes, ha="center", va="center", fontsize=fig18.FS,
                    color=vz.INK_SOFT, zorder=6,
                    bbox=dict(facecolor="#ffffff", alpha=0.9, edgecolor="none",
                              boxstyle="round,pad=0.4"))
        for b in meta["boxes"]:
            vz.draw_box(ax, b, lw=0.9, color=vz.INK_SOFT)
        ax.set_title(f"({letter}) {title(p)}", loc="left", fontsize=fig18.FS + 0.5,
                     color=vz.INK, pad=3)

    cax = fig.add_axes([0.3, (caption_in + 0.42) / height, 0.4, 0.1 / height])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both",
                      ticks=np.arange(-vmax, vmax + 0.01, 0.5))
    cb.set_label("SST anomaly, box mean removed (°C)", fontsize=fig18.FS, color=vz.INK_SOFT)
    cb.ax.tick_params(labelsize=fig18.FS - 0.5)

    head = ("SST and Z300 pattern of the summers in the top decile of MPI-GE: "
            + ("heatwave days" if index == "hwd" else f"the {p2.DIPOLE_LABEL} SST contrast"))
    fig.text(0.012, (caption_in - 0.04) / height, head, ha="left", va="top",
             fontsize=fig18.FS, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left",
             va="top", fontsize=fig18.FS, color=vz.INK, linespacing=1.45)
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}  (SST scale ±{vmax} °C)")


def main(index):
    sel, _ = selection(index)
    panels = composites(sel)
    figure(index, panels, p2.FIG_DIR / f"fig25_composite_maps_{index}.png")
    for p in panels:
        name = "MPI-GE" if p["key"] == "MPI-GE" else f"{p['key']} r{p['member']}"
        extra = "" if not p["n"] else (f"  r = {p['r']:+.2f}  Z300 "
                                       f"{'yes' if p['zg'] is not None else 'no'}")
        print(f"  {name:20s} n = {p['n']:3d}{extra}  {p['years'] if p['key'] != 'MPI-GE' else ''}")


if __name__ == "__main__":
    for index in [a for a in sys.argv[1:] if a in INDICES] or list(INDICES):
        main(index)
