"""
Figure 28 -- heatwave days, the Mediterranean-North Atlantic SST contrast and
a blocking index over the heatwave region: the three pairwise linear
regressions, one scatter panel per pair and dataset, as fig 8 does for heatwave
days against the SST pattern correlation.

Three JJA quantities per season (and member):

  hwd   land-mean heatwave days over src.heatwave.MEAN_BOX (35-60N, 10W-30E;
        data/heatwave_stats/season_table.csv, as fig 8's y axis)
  bd    the box-difference index MED - NAL (data/sst_pattern_corr/box_difference.nc)
  slp   blocking index: JJA mean sea-level pressure averaged (cos-lat weighted,
        land and sea) over the same MEAN_BOX, anomaly against 1991-2020 of the
        same dataset and member (data/circulation/<ds>_jja_slp_2p5deg.nc,
        pre_process/11_slp.sh). Persistent anticyclones over Europe raise the
        seasonal-mean SLP there; it is a proxy for blocking, not a blocking
        detection (no reversal or persistence criterion).
  z300  the same with the JJA 300 hPa geopotential height anomaly (in dam) in
        place of SLP (data/sst_anomaly/<ds>_jja_z300_anom_1deg.nc,
        pre_process/09_era5_zg.sh and 10_model_zg.sh). The upper-level ridge
        of a blocked summer, which the surface heat low over hot land hides
        from SLP; it also grows with the warmth of the column below.

The blocking index is chosen by the argument `slp` (default) or `z300`
(CIRC); the outputs are named fig28_hwd_sst_<index>_regression*. Z300 exists
for the models from 1991 only, and for ICON-EERIE r1 only.

Pairs (pairs(), y regressed on x): hwd on the index, bd on the index, hwd on bd.
Datasets (KEYS): ERA5, MPI-GE, ICON-EPOC (hist + SSP5-8.5) and ICON-EERIE
(hist + SSP2-4.5, r1-r3 pooled) -- those with all three quantities;
MPI-ESM1.2-ER has no heatwave days. Seasons of PERIOD (1990-2026) as in fig 8;
ERA5 is drawn over its whole record 1940-2026 (black line) with the PERIOD
fit light. Only seasons that have all three quantities enter.

Two versions (file suffix _detrended): as is, and with the forced part of all
three quantities removed as fig 8's detrend(): MPI-GE minus its 50-member mean
year by year, every other record (and EERIE member) minus its own linear trend
over the seasons used. The warming raises heatwave days, and possibly SLP and
the box difference, in step; the detrended version keeps only the year-to-year
co-variability.

fig28_..._bars / _bars_detrended: the slopes and r of the three pairs per
dataset, as is and detrended, one figure each on the same axes; MPI-GE's whisker is the 5th-95th percentile of the
same statistic fitted to each member alone over PERIOD (37 seasons, the
length of the observed record), i.e. how far one realisation of that length
can stray.
"""
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
import src.heatwave as hw
import src.pattern_2026 as p2
import src.viz2026 as vz

KEYS = ["ERA5", "MPI-GE", "ICON-EPOC-hist", "EERIE"]
NAMES = {"ERA5": "ERA5", "MPI-GE": "MPI-GE", "ICON-EPOC-hist": "ICON-EPOC",
         "EERIE": "ICON-EERIE r1–r3"}
SLP_STEM = {"ERA5": "era5", "MPI-GE": "mpige", "ICON-EPOC-hist": "epoc2_020",
            "EERIE": "eerie_hist_ssp245"}
PERIOD = (1990, p2.REF_YEAR)
#: MPI-GE's late century (SSP2-4.5), a second MPI-GE bar in the bar figures;
#: anomalies and heatwave thresholds stay those of 1991-2020
LATE = (2061, 2100)
LATE_LABEL = f"MPI-GE {LATE[0]}–{LATE[1]}"
ERA5_MARKED = [1947, 2003, p2.REF_YEAR]
#: vertical label offset (points) of the marked seasons; 1947 and 2003 lie close
MARK_DY = {1947: -8, 2003: 6, p2.REF_YEAR: 6}
LABELS = {"hwd": "Heatwave days (JJA)", "bd": f"{p2.DIPOLE_LABEL} SST anomaly (°C)"}
DIR = p2.PROJECT_ROOT / "data" / "circulation"
ZG_STEM = {"ERA5": "era5", "MPI-GE": "mpige", "ICON-EPOC-hist": "epoc2_020",
           "EERIE": "eerie_hist_ssp245"}


def slp_index(key):
    """(member, year) Series of the JJA SLP anomaly averaged over hw.MEAN_BOX."""
    with xr.open_dataset(DIR / f"{SLP_STEM[key]}_jja_slp_2p5deg.nc") as ds:
        slp = ds.slp.load()
    box = slp.sel(lat=slice(*hw.MEAN_BOX["lat"]), lon=slice(*hw.MEAN_BOX["lon"]))
    w = np.cos(np.deg2rad(box.lat)) * xr.ones_like(box.lon)
    x = (box * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))
    x = x.where(slp.notnull().any(("lat", "lon")))
    x = x - x.sel(year=slice(p2.CLIM_PERIOD[0], p2.CLIM_PERIOD[1])).mean("year")
    if "member" not in x.dims:
        x = x.expand_dims(member=[1])
    return x.to_series()


def z300_index(key):
    """(member, year) Series of the JJA Z300 anomaly (dam) averaged over hw.MEAN_BOX."""
    with xr.open_dataset(p2.DATA_DIR / f"{ZG_STEM[key]}_jja_z300_anom_1deg.nc") as ds:
        z = ds.z300.load()
    if key == "ERA5":
        z = z.squeeze("plev", drop=True)
        z = z.assign_coords(year=z.time.dt.year).swap_dims(time="year").drop_vars("time")
        z = z.expand_dims(member=[1])
    box = z.sel(lat=slice(*hw.MEAN_BOX["lat"]), lon=slice(*hw.MEAN_BOX["lon"]))
    w = np.cos(np.deg2rad(box.lat)) * xr.ones_like(box.lon)
    x = (box * w).sum(("lat", "lon")) / w.sum(("lat", "lon")) / 10
    return x.where(z.notnull().any(("lat", "lon"))).to_series()


#: the two blocking indices: short name, unit, axis label, loader, caption
#: text, EERIE members that have it
CIRC = {
    "slp": dict(name="SLP", unit="hPa", label="SLP anomaly over the\nheatwave box (hPa)",
                load=slp_index, eerie="r1–r3",
                desc=("SLP: JJA mean sea-level pressure averaged over the same "
                      f"{hw.MEAN_BOX_LABEL} box, land and sea, a proxy for blocking over "
                      "Europe (high = anticyclonic).")),
    "z300": dict(name="Z300", unit="dam", label="Z300 anomaly over the\nheatwave box (dam)",
                 load=z300_index, eerie="r1",
                 desc=("Z300: JJA 300 hPa geopotential height anomaly averaged over the same "
                       f"{hw.MEAN_BOX_LABEL} box, in dam, a proxy for the upper-level ridge of "
                       "blocking (it also rises with the warmth of the column below). Models "
                       "from 1991; ICON-EERIE r1 only, the other members have no upper-air "
                       "output.")),
}


def pairs(var):
    """(x, y, slope unit, short label) of each regression for blocking index `var`."""
    c = CIRC[var]
    return [("blk", "hwd", f"days per {c['unit']}", f"Heatwave days on {c['name']}"),
            ("blk", "bd", f"°C per {c['unit']}", f"{p2.DIPOLE_LABEL} on {c['name']}"),
            ("bd", "hwd", "days per °C", f"Heatwave days on {p2.DIPOLE_LABEL}")]


def names(var):
    """Panel / bar names of the datasets; ICON-EERIE named by the members used."""
    return dict(NAMES, EERIE=f"ICON-EERIE {CIRC[var]['eerie']}")


def box_difference(key):
    """(member, year) Series of the box-difference index."""
    with xr.open_dataset(p2.dipole_file()) as ds:
        da = ds[key].load()
    dim = p2.member_dim(da)
    da = da.expand_dims(member=[1]) if dim is None else da.rename({dim: "member"})
    return da.to_series()


def season_table(var):
    """One row per dataset, member and year with all three quantities."""
    t = pd.read_csv(fig18.DIR / "season_table.csv")
    rows = []
    for key in KEYS:
        h = t[t.dataset == key].set_index(["member", "year"]).hwd
        d = pd.DataFrame({"hwd": h, "bd": box_difference(key), "blk": CIRC[var]["load"](key)})
        d = d.dropna().reset_index()
        # the models over PERIOD only, so that the detrending sees the seasons
        # of the panel (as fig 8); ERA5 keeps its whole record
        if key != "ERA5":
            keep = d.year.between(*PERIOD)
            if key == "MPI-GE":
                keep |= d.year.between(*LATE)
            d = d[keep]
        d.insert(0, "dataset", key)
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def detrend(t):
    """All three quantities with their forced part removed (as fig 8's detrend)."""
    t = t.copy()
    cols = ["hwd", "bd", "blk"]
    for key, g in t.groupby("dataset"):
        if g.member.nunique() >= p2.ENSEMBLE_MEAN_MIN_MEMBERS:
            t.loc[g.index, cols] = g[cols] - g.groupby("year")[cols].transform("mean")
            continue
        for _, gm in g.groupby("member"):
            for c in cols:
                b = np.polyfit(gm.year, gm[c], 1)
                t.loc[gm.index, c] = gm[c] - np.polyval(b, gm.year)
    return t


def fit(x, y):
    """(slope, r, p) of the least-squares line of y on x."""
    res = stats.linregress(x, y)
    return res.slope, res.rvalue, res.pvalue


def subsets(t, var):
    """(label, key, rows) of every dataset fit: ERA5 over its whole record and PERIOD."""
    out = []
    for key in KEYS:
        g = t[t.dataset == key]
        if key == "ERA5":
            out.append((f"ERA5 {g.year.min()}–{g.year.max()}", key, g))
        out.append((names(var)[key] if key != "ERA5" else f"ERA5 {PERIOD[0]}–{PERIOD[1]}", key,
                    g[g.year.between(*PERIOD)]))
        if key == "MPI-GE":
            out.append((LATE_LABEL, key, g[g.year.between(*LATE)]))
    return out


def scatter(t, var, detrended, out):
    vz.use_style()
    prs = pairs(var)
    labels = dict(LABELS, blk=CIRC[var]["label"])
    n_row, n_col = len(prs), len(KEYS)
    fs = fig18.FS
    caption = (
        "One point per JJA season (and member). Heatwave days: land mean over "
        f"{hw.MEAN_BOX_LABEL} (as Fig. 8). {p2.DIPOLE_LABEL}: Mediterranean "
        f"({p2.DIPOLE_BOXES['med']['detail']}) minus North Atlantic "
        f"({p2.DIPOLE_BOXES['nal']['detail']}) JJA SST anomaly. {CIRC[var]['desc']} All "
        "against 1991–2020 of the same dataset. "
        f"Seasons {PERIOD[0]}–{PERIOD[1]} (ICON-EPOC to 2025, ICON-EERIE r1 to 2026"
        + (" and r2 to 2014, r3 to 2020, its members pooled" if var == "slp" else "") + "); ERA5 all seasons 1940–2026 (black line, first "
        f"statistics) and {PERIOD[0]}–{PERIOD[1]} (grey line, second); marked "
        + ", ".join(map(str, ERA5_MARKED)) + ". Lines: least-squares fit of y on x; r: Pearson "
        "correlation, two-sided p treating every season as independent."
        + (" Detrended: MPI-GE minus its 50-member mean year by year, every other record (and "
           "ICON-EERIE member) minus its own linear trend, for all three quantities."
           if detrended else ""))
    text = textwrap.fill(caption, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    panel_in = 1.55
    top_in = 0.3
    height = top_in + n_row * (panel_in + 0.42) + 0.15 + caption_in
    fig, axes = plt.subplots(n_row, n_col, figsize=(fig18.FIG_SIZE[0], height), sharey="row")
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.1, right=0.985, top=1 - top_in / height,
                        bottom=(caption_in + 0.42) / height, hspace=0.42, wspace=0.08)
    for i, (xv, yv, unit, _) in enumerate(prs):
        for j, key in enumerate(KEYS):
            ax = axes[i, j]
            g = t[t.dataset == key]
            if key != "ERA5":
                g = g[g.year.between(*PERIOD)]
            color = vz.COLORS[key]
            big = len(g) > 300
            ax.scatter(g[xv], g[yv], s=3 if big else 9, color=color, alpha=0.35 if big else 0.75,
                       lw=0, zorder=2)
            lines = [(g, color, 1.1)]
            if key == "ERA5":
                lines.append((g[g.year.between(*PERIOD)], vz.INK_MUTED, 0.9))
                for y in ERA5_MARKED:
                    s = g[g.year == y]
                    if len(s):
                        ax.scatter(s[xv], s[yv], s=26, facecolor="none", edgecolor=vz.INK,
                                   lw=0.9, zorder=4)
                        # left of points far to the right, so the label stays in its panel
                        left = s[xv].iloc[0] > np.percentile(g[xv], 75)
                        ax.annotate(str(y), (s[xv].iloc[0], s[yv].iloc[0]),
                                    xytext=(-5 if left else 5, MARK_DY.get(y, 3)),
                                    ha="right" if left else "left", va="center",
                                    textcoords="offset points", fontsize=fs - 1.5, color=vz.INK)
            else:
                lines = [(g[g.year.between(*PERIOD)], vz.INK, 1.1)]
            tags = []
            for gg, c, lw in lines:
                slope, r, p = fit(gg[xv], gg[yv])
                xs = np.linspace(gg[xv].min(), gg[xv].max(), 2)
                ax.plot(xs, gg[yv].mean() + slope * (xs - gg[xv].mean()), color=c, lw=lw, zorder=3)
                ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
                tags.append(f"r = {r:+.2f} ({ptxt})\n{slope:+.2f} {unit}")
            ax.text(0.16 if j == 0 else 0.03, 0.97, "\n".join(tags), transform=ax.transAxes,
                    va="top", ha="left",
                    fontsize=fs - 2, color=vz.INK, linespacing=1.25,
                    bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1))
            if i == 0:
                ax.set_title(names(var)[key], fontsize=fs, color=vz.INK, pad=3)
            ax.set_xlabel(labels[xv], fontsize=fs - 1)
            if j == 0:
                ax.set_ylabel(labels[yv], fontsize=fs - 1)
            ax.tick_params(labelsize=fs - 1.5)
            ax.axhline(0, color=vz.INK_MUTED, lw=0.4, zorder=0)
            ax.axvline(0, color=vz.INK_MUTED, lw=0.4, zorder=0)
        fig18.tag(axes[i, 0], "abc"[i])
        # head room for the statistics above the cloud (the row shares its y axis)
        lo, hi = axes[i, 0].get_ylim()
        axes[i, 0].set_ylim(lo, hi + 0.45 * (hi - lo))
    fig18.style_axes(axes.ravel())
    title = ("Heatwave days, the Mediterranean–North Atlantic SST contrast and "
             f"{CIRC[var]['name']} over the heatwave region" + (", detrended" if detrended else ""))
    fig.text(0.012, (caption_in - 0.04) / height, title, ha="left", va="top", fontsize=fs,
             fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left", va="top",
             fontsize=fs, color=vz.INK, linespacing=1.45)
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")


#: short tick labels of the bar figure
TICK = {f"ERA5 1940–{p2.REF_YEAR}": "ERA5 1940–", f"ERA5 {PERIOD[0]}–{PERIOD[1]}": f"ERA5 {PERIOD[0]}–",
        "MPI-GE": f"MPI-GE {PERIOD[0]}–", LATE_LABEL: f"MPI-GE {LATE[0]}–{LATE[1]}",
        "ICON-EPOC": "ICON-EPOC", "ICON-EERIE r1–r3": "ICON-EERIE", "ICON-EERIE r1": "ICON-EERIE r1"}


def bar_color(label, key):
    """
    (colour, alpha): ERA5 over its whole record black, over PERIOD grey; models
    their colour, MPI-GE's late century lighter.
    """
    if key == "ERA5" and not label.startswith("ERA5 1940"):
        return vz.INK_MUTED, 0.85
    return vz.COLORS[key], 0.4 if label == LATE_LABEL else 0.85


def member_spread(t, var, years):
    """5th / 95th percentile of each pair's slope and r over MPI-GE members, `years`."""
    g = t[(t.dataset == "MPI-GE") & t.year.between(*years)]
    out = {}
    for xv, yv, _, _ in pairs(var):
        f = np.array([fit(m[xv], m[yv]) for _, m in g.groupby("member")])
        out[xv, yv] = {"slope": np.percentile(f[:, 0], [5, 95]), "r": np.percentile(f[:, 1], [5, 95])}
    return out


def results(t_raw, t_det, var):
    """Slope, r and p of every pair and dataset fit, as is and detrended; MPI-GE member spread."""
    rows = []
    for version, t in (("as is", t_raw), ("detrended", t_det)):
        for label, key, g in subsets(t, var):
            for xv, yv, _, _ in pairs(var):
                slope, r, p = fit(g[xv], g[yv])
                rows.append(dict(version=version, label=label, key=key, pair=(xv, yv),
                                 slope=slope, r=r, p=p, n=len(g)))
    # version -> MPI-GE bar label -> pair -> statistic -> (5th, 95th percentile)
    spread = {version: {names(var)["MPI-GE"]: member_spread(t, var, PERIOD),
                        LATE_LABEL: member_spread(t, var, LATE)}
              for version, t in (("as is", t_raw), ("detrended", t_det))}
    return pd.DataFrame(rows), spread


def slope_limits(res, spread, var):
    """Per pair, one slope axis range for both versions (bars and MPI-GE whiskers)."""
    lims = {}
    for xv, yv, _, _ in pairs(var):
        v = list(res[res.pair == (xv, yv)].slope) + [0.0]
        v += [b for s in spread.values() for m in s.values() for b in m[xv, yv]["slope"]]
        lo, hi = min(v), max(v)
        pad = 0.08 * (hi - lo)
        lims[xv, yv] = (lo - pad if lo < 0 else 0, hi + 0.15 * (hi - lo) if hi > 0 else 0)
    # one symmetric r range for all pairs, in steps of 0.25, at least +-0.75
    r = list(res.r.abs()) + [abs(b) for s in spread.values() for m in s.values()
                             for v in m.values() for b in v["r"]]
    lims["r"] = max(0.75, np.ceil(max(r) / 0.25) * 0.25)
    return lims


def bars(res, spread, var, version, lims, out):
    """Slope (top) and r (bottom) of the three pairs per dataset, one version."""
    vz.use_style()
    prs = pairs(var)
    fs = fig18.FS
    res = res[res.version == version]
    labels = list(dict.fromkeys(res.label))
    detrended = version == "detrended"

    caption = (
        "Slope (top, in the units of y per unit of x) and Pearson r (bottom) of the three "
        f"regressions of Fig. 28{' (detrended)' if detrended else ''} per dataset"
        + (", with the forced part of all three quantities removed (MPI-GE minus its ensemble "
           "mean year by year, every other record and ICON-EERIE member minus its own linear "
           "trend)" if detrended else ", as is (no trend removed)")
        + ". Dot at the end of a bar: p < 0.05 (seasons treated as independent). "
        "ERA5 over its whole record black, over the common period grey. MPI-GE: "
        f"all 50 members pooled over {PERIOD[0]}–{PERIOD[1]}; whisker the 5th–95th percentile of "
        f"the same statistic fitted to each member alone over {PERIOD[0]}–{PERIOD[1]}, the "
        "range one realisation as long as the observed record can show. Light blue: MPI-GE in "
        f"{LATE[0]}–{LATE[1]} (SSP2-4.5), all 50 members pooled, whisker from single members over "
        "the same 40 years; its heatwave days still count against the 1991–2020 threshold and its "
        "other quantities against 1991–2020"
        + (". ICON-EERIE members pooled. " if var == "slp" else ". ")
        + f"{CIRC[var]['desc']} Heatwave days over land in the same box; {p2.DIPOLE_LABEL} the "
        "SST box difference. "
        "Same axes as the " + ("raw" if detrended else "detrended") + " version of this figure.")
    text = textwrap.fill(caption, fig22.CAPTION_WIDTH)
    caption_in = fig22.CAPTION_HEAD_IN + fig22.LINE_IN * (text.count("\n") + 1)
    panel_in, top_in, tick_in = 1.75, 0.3, 0.85
    height = top_in + 2 * panel_in + 0.35 + tick_in + caption_in
    fig, axes = plt.subplots(2, len(prs), figsize=(fig18.FIG_SIZE[0], height), sharex=True)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.08, right=0.99, top=1 - top_in / height,
                        bottom=(caption_in + tick_in) / height, hspace=0.18, wspace=0.32)
    x = np.arange(len(labels))
    w = 0.62
    for j, (xv, yv, unit, name) in enumerate(prs):
        sub = res[res.pair == (xv, yv)].set_index("label")
        for i, stat in enumerate(("slope", "r")):
            ax = axes[i, j]
            for xi, lab in enumerate(labels):
                row = sub.loc[lab]
                color, alpha = bar_color(lab, row.key)
                ax.bar(xi, row[stat], w, color=color, lw=0, alpha=alpha)
                if row.key == "MPI-GE":
                    lo, hi = spread[version][lab][xv, yv][stat]
                    ax.errorbar(xi, row[stat], yerr=[[row[stat] - lo], [hi - row[stat]]],
                                fmt="none", ecolor=vz.INK, elinewidth=0.7, capsize=2)
                if row.p < 0.05:
                    ax.plot(xi, row[stat], "o", ms=2.6, color=vz.INK, mec="white", mew=0.4,
                            zorder=4)
            ax.axhline(0, color=vz.INK_SOFT, lw=0.6)
            if i == 0:
                ax.set_title(name, fontsize=fs, color=vz.INK, pad=3)
                ax.set_ylabel(f"Slope ({unit})", fontsize=fs - 0.5)
                ax.set_ylim(*lims[xv, yv])
            else:
                ax.set_ylabel("Correlation r", fontsize=fs - 0.5)
                ax.set_ylim(-lims["r"], lims["r"])
            ax.tick_params(labelsize=fs - 1.5)
            fig18.tag(ax, "abcdef"[i * len(prs) + j])
    for ax in axes[-1]:
        ax.set_xticks(x, [TICK.get(l, l) for l in labels], fontsize=fs - 2, rotation=35,
                      ha="right", rotation_mode="anchor")
    fig18.style_axes(axes.ravel())
    fig.text(0.012, (caption_in - 0.04) / height,
             "Regression slopes and correlations between heatwave days, the SST contrast and "
             + CIRC[var]["name"]
             + (", detrended" if detrended else ""),
             ha="left", va="top", fontsize=fs, fontweight="bold", color=vz.INK)
    fig.text(0.012, (caption_in - fig22.CAPTION_HEAD_IN) / height, text, ha="left", va="top",
             fontsize=fs, color=vz.INK, linespacing=1.45)
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"wrote {out}")


def main(var="slp"):
    t = season_table(var)
    print(t.groupby("dataset").agg(n=("year", "size"), y0=("year", "min"), y1=("year", "max"),
                                   members=("member", "nunique")))
    t_det = detrend(t)
    base = p2.FIG_DIR / f"fig28_hwd_sst_{var}_regression"
    scatter(t, var, False, base.with_suffix(".png"))
    scatter(t_det, var, True, base.parent / (base.name + "_detrended.png"))
    res, spread = results(t, t_det, var)
    lims = slope_limits(res, spread, var)
    bars(res, spread, var, "as is", lims, base.parent / (base.name + "_bars.png"))
    bars(res, spread, var, "detrended", lims, base.parent / (base.name + "_bars_detrended.png"))
    pd.set_option("display.width", 160)
    res["pair"] = res.pair.map(lambda p: f"{p[1]}~{p[0]}")
    print(res.pivot_table(index=["pair", "label"], columns="version",
                          values=["slope", "r", "p"]).round(3))


if __name__ == "__main__":
    main(*[a for a in sys.argv[1:] if a in CIRC][:1])
