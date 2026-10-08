"""
Figure 12 -- how much of the resemblance to JJA 2026 comes from the Atlantic
and how much from the Mediterranean? Season-by-season pattern correlation
with ERA5 JJA 2026, scored three ways:

  box   the whole analysis box, 20-60N, 80W-40E (results/corr_spatial.nc)
  atl   the Atlantic half alone, west of 0E       (corr_spatial_atl.nc)
  med   the Mediterranean half alone, east of 0E  (corr_spatial_med.nc)

One row per dataset, each on its own time axis so the short km-scale records
are not slivers on MPI-GE's 1850-2100. Each row reports the Pearson
correlation over time of the Atlantic and of the Mediterranean series with
the whole-box series: which half the whole-box score follows.

Ensembles: the r values are pooled over every member and season; the lines
show a single member -- the one holding the dataset's whole-box best
analogue -- with that member's own r alongside. ERA5's 2026 is included
(r = 1 in all three by construction).
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

VARIANT = p2.PLOT_VARIANTS[0]
KEYS = ["ERA5", "MPI-GE", "MPI-ER", "ICON-EPOC-hist", "EERIE",
        "ICON-EPOC-ctrl", "EERIE-ctrl", "ICON-sap-ctrl"]
REGIONS = {
    None: dict(label=f"whole box ({p2.REGION_LABEL})", color=vz.INK, lw=1.5),
    "atl": dict(label=f"Atlantic ({p2.SPLIT_REGIONS['atl']['detail']})",
                color=vz.COLORS["MPI-GE"], lw=1.1),
    "med": dict(label=f"Mediterranean ({p2.SPLIT_REGIONS['med']['detail']})",
                color=vz.COLORS["ICON-EPOC-hist"], lw=1.1),
}


def load():
    """{half: Dataset of correlation series}, member dims renamed to 'member'."""
    out = {}
    for half in REGIONS:
        ds = xr.open_dataset(p2.corr_file(VARIANT, half)).load()
        out[half] = {k: (ds[k].rename({p2.member_dim(ds[k]): "member"})
                         if p2.member_dim(ds[k]) else ds[k])
                     for k in KEYS if k in ds}
    return out


def pearson(a, b):
    """Pearson r over every (member,) season where both are finite."""
    a, b = xr.align(a, b, join="inner")
    x, y = a.values.ravel(), b.values.ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[ok], y[ok])[0, 1]), int(ok.sum())


def main():
    vz.use_style()
    corr = load()
    best = pd.read_csv(p2.best_file())
    keys = [k for k in KEYS if k in corr[None]]

    fig_h = 1.35 * len(keys) + 1.9
    fig, axes = plt.subplots(len(keys), 1, figsize=(10.5, fig_h), squeeze=False)
    axes = axes[:, 0]
    fig.subplots_adjust(hspace=0.42, top=1 - 1.25 / fig_h, bottom=0.55 / fig_h,
                        left=0.17, right=0.80)

    rows = []
    for ax, key in zip(axes, keys):
        series = {h: corr[h][key] for h in REGIONS}
        r_atl, n = pearson(series["atl"], series[None])
        r_med, _ = pearson(series["med"], series[None])

        shown, member_note = series, ""
        if "member" in series[None].dims:
            m = best[(best.variant == VARIANT) & (best.dataset == key)].best_member.iloc[0]
            m = int(m)
            shown = {h: s.sel(member=m) for h, s in series.items()}
            ra_m, _ = pearson(shown["atl"], shown[None])
            rm_m, _ = pearson(shown["med"], shown[None])
            nm = len(vz.valid_members(series[None]))
            member_note = (f"\npooled over {nm} members; lines: member r{m}\n"
                           f"(member r{m} alone: Atl {ra_m:+.2f}, Med {rm_m:+.2f})")

        ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
        ax.grid(axis="y", zorder=0)
        for h, st in REGIONS.items():
            s = shown[h].dropna("year")
            ax.plot(s.year, s, color=st["color"], lw=st["lw"], zorder=3 if h is None else 2)
        if key == "ERA5":
            ax.axvline(p2.REF_YEAR, color=vz.INK, lw=0.8, ls=":", zorder=4)

        yrs = series[None].dropna("year", how="all").year
        ax.set_xlim(int(yrs.min()) - 0.5, int(yrs.max()) + 0.5)
        ax.set_ylim(-1.0, 1.05)
        ax.set_yticks([-0.5, 0, 0.5, 1])
        ax.tick_params(labelsize=7.5)
        d = p2.DATASETS[key]
        ax.set_ylabel(d.label.replace(" (", "\n("), rotation=0, ha="right", va="center",
                      labelpad=12, fontsize=8.3, color=vz.INK, linespacing=1.35)
        ax.text(1.015, 0.5,
                f"r(Atl, box) = {r_atl:+.2f}\nr(Med, box) = {r_med:+.2f}\n"
                f"n = {n}{member_note}",
                transform=ax.transAxes, ha="left", va="center", fontsize=7.6,
                color=vz.INK_SOFT, linespacing=1.35)
        rows.append(dict(dataset=key, n=n, r_atl_box=round(r_atl, 3), r_med_box=round(r_med, 3)))
    axes[-1].set_xlabel("year")

    handles = [Line2D([], [], color=st["color"], lw=st["lw"] + 0.3, label=st["label"])
               for st in REGIONS.values()]
    fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 1 - 0.78 / fig_h))
    fig.suptitle("Which half carries the resemblance to JJA 2026? Pattern correlation over the "
                 "whole box, the Atlantic and the Mediterranean",
                 fontsize=11.5, color=vz.INK, y=1 - 0.22 / fig_h)
    fig.text(0.5, 1 - 0.5 / fig_h,
             f"centred pattern correlation with ERA5 JJA {p2.REF_YEAR}, 1° common ocean grid, each "
             "half scored on its own  ·  right: Pearson r over time of each half's series with the "
             "whole-box series  ·  dotted: 2026",
             ha="center", fontsize=7.9, color=vz.INK_SOFT)

    out = p2.FIG_DIR / "fig12_corr_atl_med_timeseries.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    pd.DataFrame(rows).to_csv(p2.RESULT_DIR / "corr_halves_vs_box.csv", index=False)
    print(f"wrote {out}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
