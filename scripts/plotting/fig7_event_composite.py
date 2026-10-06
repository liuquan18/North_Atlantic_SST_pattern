"""
Figure 7 -- is there decadal persistence in the resemblance to JJA 2026?

Figure 2 hints at slow swings in the pattern correlation, but the ensemble
envelope hides any single member's history. This figure asks the question
directly with a superposed-epoch composite: find every season where a record
starts to look like 2026 (r reaches a per-record threshold from below), align
those seasons at lag 0, and look at the 20 seasons either side. If the
resemblance is carried by slow ocean variability, the composite should rise
ahead of lag 0 and decay over several years after it; if it is weather-like,
it should spike at lag 0 and nowhere else.

Two choices keep the composite honest:

  * The forced part of each series is removed first (pattern_2026.
    internal_component). The r >= 0.5 seasons cluster after 2000 in every
    forced record, so a composite of raw r would slope upward around lag 0 from
    the trend alone and look like persistence.
  * Each panel carries a null band: the 10-90% range of the composite mean of
    the same number of onsets drawn at random from the same record. A
    composite inside the band at lag k is what chance alignment would give.

The internal component is smoothed with a 3-yr running mean (SMOOTH) before it is
composited, so the composite shows the slow part of each event's history rather
than year-to-year scatter. Events are still picked on the annual values. The
smoothing spreads the event season itself over lags -1..+1 (a box of 1/3 its
height, drawn dashed), so only a composite that stands above that box near lag
0, or outside the null band beyond |lag| 1, says anything about persistence.

The threshold is set per record (THRESHOLD): one cut for all would leave the
single runs with one to three events, too few to composite. Each panel states
its threshold and event count; with a handful of events the null band is wide
and the composite is suggestive at best.
The last panel gives the same answer from every season rather than the
selected ones: the lag autocorrelation of the internal component.
"""
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2
import src.viz2026 as vz

VARIANT = p2.PLOT_VARIANTS[0]
#: panels, row by row: observations and the forced runs on top; below them the
#: autocorrelation summary, the two controls (each under its forced run) and
#: MPI-ER under MPI-GE. "ACF" marks the autocorrelation panel.
LAYOUT = ["ERA5", "ICON-EPOC-hist", "EERIE", "MPI-GE",
          "ACF", "ICON-EPOC-ctrl", "EERIE-ctrl", "MPI-ER"]
#: event threshold on the annual correlation, per record
THRESHOLD = {"MPI-GE": 0.6, "MPI-ER": 0.5, "ERA5": 0.4,
             "ICON-EPOC-hist": 0.4, "ICON-EPOC-ctrl": 0.4,
             "EERIE": 0.4, "EERIE-ctrl": 0.4}
LAG = 20
#: lags shown; the composite and its null band are computed out to LAG
XLIM = 15
#: quantiles of the random-onset composite drawn as the null band
NULL_RANGE = (0.1, 0.9)
ACF_LAG = 15
#: centred running mean applied to the internal component before compositing;
#: events are still picked on the annual values, which is what "a season that
#: resembles 2026" means
SMOOTH = 3


def member_series(da: xr.DataArray) -> list[np.ndarray]:
    """One 1-D array per member (or a single one), on the full year axis."""
    mdim = p2.member_dim(da)
    if mdim is None:
        return [da.values]
    return [da.sel({mdim: m}).values for m in vz.valid_members(da, mdim)]


def pooled_acf(series: list[np.ndarray], max_lag: int) -> np.ndarray:
    """Autocorrelation pooled over members: all (x_t, x_t+k) pairs of all members together."""
    out = np.empty(max_lag + 1)
    for k in range(max_lag + 1):
        a = np.concatenate([s[: len(s) - k] for s in series])
        b = np.concatenate([s[k:] for s in series])
        ok = np.isfinite(a) & np.isfinite(b)
        out[k] = np.corrcoef(a[ok], b[ok])[0, 1]
    return out


def composite(key: str, r: xr.DataArray) -> dict:
    if key == "ERA5":
        # 2026 against itself is 1 by construction -- the reference, not an event
        r = r.where(r["year"] != p2.REF_YEAR)
    anom = p2.internal_component(r, control=key.endswith("-ctrl"))
    raw, internal = member_series(r), member_series(anom)

    threshold = THRESHOLD[key]
    rule = f"r ≥ {threshold:.2f}"
    onsets = [p2.event_onsets(x, threshold) for x in raw]

    # annual lag-0 height, for the "spike alone" reference: a single-season
    # spike of height h, smoothed, is a box of h/SMOOTH over |lag| <= SMOOTH//2
    spike = float(np.nanmean(np.concatenate(
        [x[o] for x, o in zip(internal, onsets) if len(o)])))
    smoothed = [p2.running_mean(x, SMOOTH) for x in internal]
    windows = np.concatenate([p2.epoch_windows(x, o, LAG) for x, o in zip(smoothed, onsets)])
    n = len(windows)
    band = p2.random_onset_band(smoothed, n, LAG, quantiles=NULL_RANGE)
    count = np.isfinite(windows).sum(0)
    with np.errstate(invalid="ignore"):
        mean = np.where(count >= min(3, n), np.nanmean(windows, 0), np.nan)
    n_seasons = int(sum(np.isfinite(x).sum() for x in internal))
    return dict(windows=windows, mean=mean, band=band, n=n, rule=rule, spike=spike,
                acf=pooled_acf(internal, ACF_LAG), n_seasons=n_seasons,
                n_series=len(internal))


def draw_composite(ax, key: str, c: dict):
    """One composite panel: null band, single events, composite mean, spike-alone reference."""
    lags = np.arange(-LAG, LAG + 1)
    color = vz.COLORS[key]
    ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
    ax.axvline(0, color=vz.INK_MUTED, lw=0.6, ls=":", zorder=1)
    ax.grid(axis="y", zorder=0)
    ax.fill_between(lags, c["band"][0], c["band"][1], color=vz.CONTEXT, alpha=0.35,
                    lw=0, zorder=2, label=f"random onsets, {NULL_RANGE[0]:.0%}–{NULL_RANGE[1]:.0%}")

    # many events: faint, so the mean stays readable; few: each one matters
    alpha = 0.16 if c["n"] > 15 else 0.45
    for w in c["windows"]:
        ax.plot(lags, w, color=color, lw=0.6, alpha=alpha, zorder=3)
    ax.plot([], [], color=color, lw=0.6, alpha=0.6, label="single events")
    ax.plot(lags, c["mean"], color=color, lw=2.0, ls=vz.LINESTYLES[key], zorder=5,
            label="composite mean")
    # what the composite would be if the event season alone stood out
    box = np.where(np.abs(lags) <= SMOOTH // 2, c["spike"] / SMOOTH, 0.0)
    ax.plot(lags, box, color=vz.INK, lw=0.9, ls=(0, (3, 2)), zorder=4,
            label=f"event season alone, {SMOOTH}-yr smoothed")

    # lags beyond the smoothing's reach of the event season where the
    # composite leaves the null band -- ~20% of lags by chance alone
    out = (np.abs(lags) > SMOOTH // 2) & ((c["mean"] > c["band"][1]) | (c["mean"] < c["band"][0]))
    ax.plot(lags[out], c["mean"][out], ls="none", marker="o", ms=3.6, mfc=color,
            mec=vz.SURFACE, mew=0.8, zorder=6)

    ax.set_title(p2.DATASETS[key].label, fontsize=9, color=vz.INK, pad=5, loc="left")
    vz.panel_tag(ax, f"{c['n']} onsets of {c['rule']}", loc="upper left")
    ax.set_xlim(-XLIM, XLIM)
    ax.set_ylim(-0.32, 0.42)
    ax.set_xticks(np.arange(-XLIM, XLIM + 1, 5))
    ax.set_xlabel("years from onset")


def main():
    # windows past a record's end are all-NaN at some lags; that is expected
    warnings.filterwarnings("ignore", "Mean of empty slice", RuntimeWarning)
    vz.use_style()
    p2.FIG_DIR.mkdir(parents=True, exist_ok=True)
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc")
    keys = [k for k in LAYOUT if k in corr.data_vars]
    res = {k: composite(k, corr[k]) for k in keys}

    ncol = 4
    nrow = int(np.ceil(len(LAYOUT) / ncol))
    fig_h = 3.05 * nrow + 1.25
    fig, axes = plt.subplots(nrow, ncol, figsize=(15.6, fig_h), squeeze=False)
    fig.subplots_adjust(hspace=0.42, wspace=0.17, top=1 - 1.1 / fig_h,
                        bottom=0.55 / fig_h, left=0.05, right=0.99)
    flat = axes.ravel()

    for key in keys:
        ax = flat[LAYOUT.index(key)]
        draw_composite(ax, key, res[key])
        # first column only: further in, the label would run into the
        # neighbouring panel, and every composite shares the same axis
        if LAYOUT.index(key) % ncol == 0:
            ax.set_ylabel(f"r, forced part removed,\n{SMOOTH}-yr running mean")

    flat[LAYOUT.index(keys[0])].legend(loc="lower left", fontsize=7.2, borderaxespad=0.2)

    # lag autocorrelation from every season, not just the events
    ax = flat[LAYOUT.index("ACF")]
    k = np.arange(ACF_LAG + 1)
    # the white-noise 95% range for the shortest and an ~ERA5-length record;
    # MPI-GE pools 50 members and its range is narrower than the line width
    for n_yr, a in ((35, 0.18), (86, 0.32)):
        ax.fill_between(k, -1.96 / np.sqrt(n_yr), 1.96 / np.sqrt(n_yr),
                        color=vz.CONTEXT, alpha=a, lw=0, zorder=1)
    ax.text(ACF_LAG, 1.96 / np.sqrt(35) + 0.02, "white noise, 95%: 35 yr",
            fontsize=6.8, color=vz.INK_SOFT, ha="right", va="bottom")
    ax.text(ACF_LAG, 1.96 / np.sqrt(86) - 0.02, "86 yr", fontsize=6.8,
            color=vz.INK_SOFT, ha="right", va="top")
    ax.axhline(0, color=vz.INK_MUTED, lw=0.6, zorder=1)
    for key in keys:
        ax.plot(k[1:], res[key]["acf"][1:], color=vz.COLORS[key], ls=vz.LINESTYLES[key],
                lw=1.9 if key == "MPI-GE" else 1.1, marker="o", ms=2.6, zorder=3,
                label=key)
    ax.set_xlim(0.5, ACF_LAG + 0.5)
    ax.set_ylim(-0.6, 0.6)
    ax.set_xticks([1, 5, 10, 15])
    ax.set_xlabel("lag (years)")
    ax.set_ylabel("autocorrelation of r,\nforced part removed")
    ax.set_title("Lag autocorrelation, all annual seasons", fontsize=9, color=vz.INK, pad=5,
                 loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=6.8, ncol=4,
              borderaxespad=0, handlelength=1.8, columnspacing=1.0)
    for i, name in enumerate(LAYOUT):
        if name != "ACF" and name not in keys:
            flat[i].set_visible(False)

    fig.suptitle("Do seasons that resemble JJA 2026 come in decadal spells?",
                 fontsize=12, color=vz.INK, y=1 - 0.25 / fig_h)
    fig.text(0.5, 1 - 0.5 / fig_h,
             f"pattern correlation with ERA5 JJA 2026 (box mean removed), aligned on the season it "
             f"first reaches the threshold  ·  {SMOOTH}-yr running mean before compositing\n"
             "forced part removed: 50-member mean (MPI-GE), "
             "record mean (controls), linear/quadratic trend (other runs)  ·  "
             f"dots: outside the random-onset band, |lag| > {SMOOTH // 2}",
             ha="center", va="top", fontsize=7.8, color=vz.INK_SOFT, linespacing=1.5)

    for key in keys:
        c = res[key]
        print(f"{key:15s} n={c['n']:3d} ({c['rule']})  acf1-3={np.round(c['acf'][1:4], 2)}")

    out = p2.FIG_DIR / "fig7_event_composite.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
