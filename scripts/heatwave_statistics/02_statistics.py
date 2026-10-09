"""
How extreme were the European summers of 2003 and 2026, in the observations
and in the simulations? Three questions, one output table each (all in
data/heatwave_stats/):

A. ERA5 alone -- era5_seasons.csv, era5_robustness.csv, era5_gev_params.csv,
   era5_gev_boot.csv, era5_event_rp.csv
   A nonstationary GEV (src.extremes) is fitted with the smoothed global-mean
   SST anomaly T as covariate, in two forms: "shift" (location linear in T,
   constant scale -- the usual choice for temperature) and "shift+scale"
   (log-scale linear in T too). The second matters here: with the threshold
   fixed at 1991-2020, heatwave days spread out as the climate warms (ERA5
   detrended std 1.1 days in 1950-89, 3.0 in 2007-26), and how extreme 2026
   was depends on whether its larger spread is allowed for. Every season
   gets its return period IN ITS OWN CLIMATE and the matching normal z.
   VARIANTS also vary the fit period and leave the event out of the fit;
   era5_robustness.csv collects the rank, z and own-climate return period of
   2003 and 2026 under each, plus the conventional linear-in-time
   detrending. For each variant, the return period of the 2003 / 2026 value
   in the climate of 2026, of 2003 and of the 1950s (mean T 1950-1959), the
   probability ratios, and how much weaker an equally rare 1950s season was.
   95% intervals from bootstrap refits (N_BOOT for the two main variants).

B. Internal variability -- anomaly_exceedance.csv
   Anomaly = season minus forced mean. Forced mean: MPI-GE 50-member mean
   year by year; single transient runs (ERA5, EPOC, EERIE members) least
   squares on T; control runs a linear drift in time. How often does each
   simulation produce an anomaly at least as large as ERA5's 2003 / 2026
   (absolute, in days or degC*day), and how often one as large in units of its
   own standard deviation (standardised; ERA5's anomaly is divided by its
   LOCAL spread, the residual std within +-LOCAL_HALF years of the event,
   not the 1950-2026 pooled one that the calm early decades shrink)? MPI-GE
   is also scored in 1850-1900 and in the seasons whose forced T matches
   ERA5's 2003 / 2026 warming, since its spread grows with T.

C. Return periods of the observed values in the simulated climates --
   model_event_rp.csv, model_gev_params.csv
   * MPI-GE, empirical (no fit): seasons of all members in 1850-1900, around
     2003 and 2026 (calendar +-10 years), and at ERA5's 2003 / 2026 warming
     level (forced T within +-WL_TOL of ERA5's smoothed T).
   * MPI-GE and the km-scale transient runs: the same NS-GEV as ERA5,
     evaluated at the dataset's own T of 2003 and 2026 and at ERA5's.
   * control runs: stationary GEV (drift removed, mean kept) and empirical
     counts -- the return period in an unforced climate.
   Each simulation uses its own 1991-2020 heatwave threshold, so the
   observed value means "this many days above the local 90th percentile" in
   every dataset; the reference-period means are listed for comparison.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
import src.extremes as ex

DIR = ROOT / "data" / "heatwave_stats"
METRICS = ["hwd", "hwcum"]
EVENTS = [2003, 2026]
ERA5_FIT = (1950, 2026)
#: ERA5 fit variants: (name, period, scale_trend, leave the event out,
#: covariate smoothing in years). The first two are the main ones (full
#: bootstrap, kept for the figures). The longer smoothings test whether the
#: 2023-26 global SST surge -- partly El Nino, i.e. internal -- counting as
#: "the 2026 climate" is what makes 2026 look less rare.
VARIANTS = [
    ("shift+scale 1950-2026", (1950, 2026), True, False, 4),
    ("shift 1950-2026", (1950, 2026), False, False, 4),
    ("shift+scale 1950-2026, event left out", (1950, 2026), True, True, 4),
    ("shift 1950-2026, event left out", (1950, 2026), False, True, 4),
    ("shift+scale 1950-2026, 10-yr covariate", (1950, 2026), True, False, 10),
    ("shift+scale 1950-2026, 20-yr covariate", (1950, 2026), True, False, 20),
    ("shift+scale 1940-2026", (1940, 2026), True, False, 4),
    ("shift 1940-2026", (1940, 2026), False, False, 4),
    ("shift+scale 1991-2026", (1991, 2026), True, False, 4),
    ("shift 1991-2026", (1991, 2026), False, False, 4),
]
MAIN_VARIANTS = 2
LOCAL_HALF = 15      # years, ERA5 local spread for the standardised anomaly
FIFTIES = (1950, 1959)
#: bootstrap sizes; HWS_N_BOOT scales all of them down for a quick test run
N_BOOT = int(os.environ.get("HWS_N_BOOT", 1000))
N_BOOT_SENS = max(N_BOOT * 3 // 10, 10)
WL_TOL = 0.05        # K, warming-level matching of MPI-GE seasons
CAL_HALF = 10        # years, calendar windows around 2003 / 2026
PREIND = (1850, 1900)
TRANSIENT_FIT = {"MPI-GE": (1850, 2100), "ICON-EPOC-hist": (1990, 2026), "EERIE": (1950, 2050)}
CONTROLS = ["EERIE-ctrl", "ICON-sap-ctrl", "ICON-EPOC-ctrl"]
#: bootstrap processes: a whole compute node under SLURM, a few on a login node
JOBS = min(64, int(os.environ.get("SLURM_CPUS_ON_NODE", 8)))


def load():
    return pd.read_csv(DIR / "season_table.csv")


def era5(t, period=ERA5_FIT, window=ex.SMOOTH_YEARS):
    """ERA5 seasons of `period`; T re-smoothed over `window` years if not the default."""
    e = t[t.dataset == "ERA5"].sort_values("year").reset_index(drop=True)
    if window != ex.SMOOTH_YEARS:
        e["T"] = ex.smooth_covariate(e.gmsst.values, e.year.values, window)
    return e[e.year.between(*period)].reset_index(drop=True)


def climates(e):
    """Named covariate values on ERA5's T axis."""
    T = e.set_index("year")["T"]
    return {"2026": float(T.loc[2026]), "2003": float(T.loc[2003]),
            "1950s": float(T.loc[FIFTIES[0]:FIFTIES[1]].mean())}


def z_of(fit, x, T):
    """Normal-equivalent of a value's percentile in its own climate."""
    return stats.norm.ppf(np.clip(fit.cdf(x, T), 1e-12, 1 - 1e-12))


# --------------------------------------------------------------------------
# A. ERA5
# --------------------------------------------------------------------------
def diagnostics(fit, x, T, scale_trend):
    """
    AIC, a KS test of the fitted percentiles against uniform, and the spread
    of the normal-equivalent z in the first and second half of the record
    (with Levene's p): a fit that misses a change in spread leaves z narrow
    early and wide late.
    """
    k = 5 if scale_trend else 4
    nll = ex._nll([fit.mu0, fit.mu1, fit.s0, fit.s1, fit.xi], x, T, True)
    u = fit.cdf(x, T)
    z = stats.norm.ppf(np.clip(u, 1e-9, 1 - 1e-9))
    h = len(z) // 2
    return dict(aic=2 * nll + 2 * k, ks_p=stats.kstest(u, "uniform").pvalue,
                z_std_early=z[:h].std(ddof=1), z_std_late=z[h:].std(ddof=1),
                levene_p=stats.levene(z[:h], z[h:]).pvalue)


def own_climate_scores(fit, x, T):
    """Return period in its own climate, normal z and rank of every season."""
    rp = fit.return_period(x, T)
    z = z_of(fit, x, T)
    return rp, z, pd.Series(z).rank(ascending=False).astype(int).values


def linear_ranks(e, m, period):
    """Rank and z of each season after removing a linear trend in time (constant spread)."""
    w = e[e.year.between(*period)]
    r = w[m] - np.polyval(np.polyfit(w.year, w[m], 1), w.year)
    z = r / r.std(ddof=1)
    return dict(zip(w.year, z)), dict(zip(w.year, z.rank(ascending=False).astype(int)))


def era5_part(t):
    e_all = era5(t, (1900, 2100))
    clim = climates(e_all)
    e = era5(t)
    seasons = e[["year", "T"] + METRICS].copy()
    params, events, boots, robust = [], [], [], []
    for m in METRICS:
        # conventional detrending, for comparison
        for period in ((1950, 2026), (1991, 2026)):
            z, rank = linear_ranks(e_all, m, period)
            tag = "lin" if period[0] == 1950 else "lin91"
            seasons[f"{m}_z_{tag}"] = seasons.year.map(z)
            seasons[f"{m}_rank_{tag}"] = seasons.year.map(rank)
            for ev in EVENTS:
                robust.append(dict(metric=m, variant=f"linear in time {period[0]}-{period[1]}",
                                   event=ev, n=len(z), rank=rank[ev], z=z[ev], rp=np.nan,
                                   xi=np.nan, s1=np.nan))

        for vi, (name, period, scale_trend, leave_out, window) in enumerate(VARIANTS):
            n_boot = N_BOOT if vi < MAIN_VARIANTS else N_BOOT_SENS
            w = era5(t, period, window)
            # climates on this variant's own covariate
            vclim = climates(era5(t, (1900, 2100), window))
            for ev in (EVENTS if leave_out else [None]):
                fw = w if ev is None else w[w.year != ev]
                x, T = fw[m].values, fw["T"].values
                fit = ex.fit_nsgev(x, T, scale_trend=scale_trend)
                boot = ex.bootstrap_nsgev(x, T, n=n_boot, seed=1 + vi, start=fit,
                                          scale_trend=scale_trend, jobs=JOBS)
                # scores of every season of the period (the left-out one included)
                rp, z, rank = own_climate_scores(fit, w[m].values, w["T"].values)
                sc = pd.DataFrame(dict(year=w.year.values, rp=rp, z=z, rank=rank)).set_index("year")
                for e_ev in (EVENTS if ev is None else [ev]):
                    xe = float(w.loc[w.year == e_ev, m].iloc[0])
                    robust.append(dict(metric=m, variant=name, event=e_ev, n=len(x),
                                       rank=int(sc.loc[e_ev, "rank"]), z=sc.loc[e_ev, "z"],
                                       rp=sc.loc[e_ev, "rp"], xi=fit.xi, s1=fit.s1))
                    events += event_rows(m, e_ev, xe, fit, boot, vclim, name)
                if ev is None:
                    row = dict(metric=m, variant=name, n=len(x), **fit.as_dict(),
                               **diagnostics(fit, x, T, scale_trend))
                    for k in fit.as_dict():
                        row[f"{k}_lo"], row[f"{k}_hi"] = ex.ci([getattr(b, k) for b in boot])
                    params.append(row)
                if vi < MAIN_VARIANTS:
                    tag = "" if vi == 0 else "_shift"
                    seasons[f"{m}_median{tag}"] = fit.return_level(2, seasons["T"])
                    seasons[f"{m}_rp_own{tag}"] = seasons.year.map(sc.rp)
                    seasons[f"{m}_z{tag}"] = seasons.year.map(sc.z)
                    seasons[f"{m}_rank{tag}"] = seasons.year.map(sc["rank"])
                    boots += [dict(metric=m, variant=name, i=i, **b.as_dict())
                              for i, b in enumerate(boot)]
            print(f"  ERA5 {m:5s} {name:40s} " + "  ".join(
                f"{r['event']}: #{r['rank']} z={r['z']:.2f} RP={r['rp']:.0f}"
                for r in robust[-2:]))
    seasons.to_csv(DIR / "era5_seasons.csv", index=False, float_format="%.5g")
    pd.DataFrame(robust).to_csv(DIR / "era5_robustness.csv", index=False, float_format="%.5g")
    pd.DataFrame(params).to_csv(DIR / "era5_gev_params.csv", index=False, float_format="%.5g")
    pd.DataFrame(boots).to_csv(DIR / "era5_gev_boot.csv", index=False, float_format="%.6g")
    ev = pd.DataFrame(events)
    ev.to_csv(DIR / "era5_event_rp.csv", index=False, float_format="%.5g")
    return seasons, ev


def event_rows(metric, event, xe, fit, boot, clim, variant):
    """Return period of value xe in each named climate, PRs, 1950s-equivalent value."""
    out = []
    rp = {c: fit.return_period(xe, T) for c, T in clim.items()}
    rpb = {c: np.array([b.return_period(xe, T) for b in boot]) for c, T in clim.items()}
    for c, T in clim.items():
        lo, hi = ex.ci_log(rpb[c])
        out.append(dict(metric=metric, event=event, value=xe, variant=variant, quantity="rp",
                        climate=c, T=T, est=rp[c], lo=lo, hi=hi))
    for a, b in (("2026", "1950s"), ("2003", "1950s"), ("2026", "2003")):
        # probability ratio p_a / p_b
        pr = float(ex.probability_ratio(1 / rp[a], 1 / rp[b]))
        lo, hi = ex.ci_log(ex.probability_ratio(1 / rpb[a], 1 / rpb[b]))
        out.append(dict(metric=metric, event=event, value=xe, variant=variant, quantity="pr",
                        climate=f"{a} vs {b}", T=np.nan, est=pr, lo=lo, hi=hi))
    # the 1950s season that is as rare as the event was in its own year
    T_ev = clim["2026"] if event == 2026 else clim["2003"]
    rp_own = fit.return_period(xe, T_ev)
    dx = xe - fit.return_level(rp_own, clim["1950s"])
    dxb = [xe - b.return_level(b.return_period(xe, T_ev), clim["1950s"]) for b in boot]
    lo, hi = ex.ci(dxb)
    out.append(dict(metric=metric, event=event, value=xe, variant=variant,
                    quantity="change vs 1950s", climate=f"{event}", T=T_ev, est=dx, lo=lo, hi=hi))
    return out


# --------------------------------------------------------------------------
# B. anomalies from the forced mean
# --------------------------------------------------------------------------
def anomalies(t):
    """Season table plus `<m>_anom` columns and the sample each season is scored in."""
    rows = []
    for key, g in t.groupby("dataset", sort=False):
        g = g.copy()
        if key == "MPI-GE":
            for m in METRICS:
                g[f"{m}_anom"] = g[m] - g.groupby("year")[m].transform("mean")
        elif key in CONTROLS:
            for m in METRICS:
                b = np.polyfit(g.year, g[m], 1)
                g[f"{m}_anom"] = g[m] - np.polyval(b, g.year)
        else:
            if key == "ERA5":
                g = g[g.year.between(*ERA5_FIT)].copy()
            for m in METRICS:
                for _, gm in g.groupby("member"):
                    b = np.polyfit(gm["T"], gm[m], 1)
                    g.loc[gm.index, f"{m}_anom"] = gm[m] - np.polyval(b, gm["T"])
        rows.append(g)
    return pd.concat(rows)


def samples(a, e_clim):
    """(dataset, sample label) -> rows of `a` making up that sample."""
    out = {}
    for key, g in a.groupby("dataset", sort=False):
        if key == "ERA5":
            out[(key, f"{ERA5_FIT[0]}-{ERA5_FIT[1]}")] = g
        elif key == "MPI-GE":
            out[(key, "1850-1900")] = g[g.year.between(*PREIND)]
            out[(key, "1991-2026")] = g[g.year.between(1991, 2026)]
            for ev in EVENTS:
                Tt = e_clim[str(ev)]
                out[(key, f"ERA5 {ev} warming level")] = g[(g["T"] - Tt).abs() <= WL_TOL]
        elif key in CONTROLS:
            out[(key, f"all {g.year.nunique()} years")] = g
        else:
            out[(key, "1990-2026")] = g[g.year.between(1990, 2026)]
    return out


def local_std(e, m, year):
    """Std of ERA5's anomalies within +-LOCAL_HALF years of `year` (e indexed by year)."""
    w = e.loc[year - LOCAL_HALF:year + LOCAL_HALF, f"{m}_anom"]
    return float(w.std(ddof=1))


def anomaly_part(t):
    a = anomalies(t)
    e = a[a.dataset == "ERA5"].set_index("year")
    e_clim = climates(era5(t))
    rows = []
    for (key, label), g in samples(a, e_clim).items():
        for m in METRICS:
            s = g[f"{m}_anom"].dropna().values
            row = dict(metric=m, dataset=key, kind=g.kind.iloc[0], sample=label, n=len(s),
                       years=f"{g.year.min()}-{g.year.max()}", std=s.std(ddof=1),
                       skew=stats.skew(s), max=s.max(), q99=np.quantile(s, 0.99))
            for ev in EVENTS:
                ae = float(e.loc[ev, f"{m}_anom"])
                sd_e = local_std(e, m, ev)
                row[f"{ev}_era5_local_std"] = sd_e
                for how, thr in (("abs", ae), ("std", ae / sd_e * s.std(ddof=1))):
                    r = ex.empirical_return_period(s, thr)
                    row[f"{ev}_{how}_thr"] = thr
                    row[f"{ev}_{how}_k"] = r["k"]
                    row[f"{ev}_{how}_p"] = r["p"]
                    row[f"{ev}_{how}_rp"] = r["rp"]
                    row[f"{ev}_{how}_rp_lo"] = r["rp_lo"]
                    row[f"{ev}_{how}_rp_hi"] = r["rp_hi"]
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(DIR / "anomaly_exceedance.csv", index=False, float_format="%.5g")
    a.to_csv(DIR / "season_anomalies.csv", index=False, float_format="%.5g")
    return df


# --------------------------------------------------------------------------
# C. return periods of the observed values in the simulations
# --------------------------------------------------------------------------
def model_part(t):
    e = era5(t)
    e_val = {(m, ev): float(e.loc[e.year == ev, m].iloc[0]) for m in METRICS for ev in EVENTS}
    e_clim = climates(e)
    rows, params = [], []

    def add(dataset, method, climate, m, ev, est, lo, hi, n, T=np.nan, k=np.nan):
        rows.append(dict(metric=m, event=ev, value=e_val[(m, ev)], dataset=dataset,
                         method=method, climate=climate, T=T, n=n, k=k, est=est, lo=lo, hi=hi))

    # MPI-GE empirical
    g = t[t.dataset == "MPI-GE"]
    windows = {"1850-1900": g[g.year.between(*PREIND)]}
    for ev in EVENTS:
        windows[f"{ev}±{CAL_HALF} (calendar)"] = g[g.year.between(ev - CAL_HALF, ev + CAL_HALF)]
        Tt = e_clim[str(ev)]
        wl = g[(g["T"] - Tt).abs() <= WL_TOL]
        windows[f"ERA5 {ev} warming level"] = wl
        print(f"  MPI-GE at ERA5 {ev} warming level (T = {Tt:+.2f} K): "
              f"years {wl.year.min()}-{wl.year.max()}, {len(wl)} seasons")
    for name, w in windows.items():
        for m in METRICS:
            for ev in EVENTS:
                r = ex.empirical_return_period(w[m], e_val[(m, ev)])
                add("MPI-GE", "empirical", name, m, ev, r["rp"], r["rp_lo"], r["rp_hi"],
                    r["n"], T=w["T"].mean(), k=r["k"])

    # NS-GEV for the transient simulations
    for key, period in TRANSIENT_FIT.items():
        g = t[(t.dataset == key) & t.year.between(*period)].dropna(subset=["T"])
        Tyear = g.groupby("year")["T"].first()
        own = {ev: float(Tyear.loc[min(ev, Tyear.index.max())]) for ev in EVENTS}
        n_boot = N_BOOT // 2 if key == "MPI-GE" else N_BOOT
        for m in METRICS:
            fit = ex.fit_nsgev(g[m].values, g["T"].values)
            boot = ex.bootstrap_nsgev(g[m].values, g["T"].values, n=n_boot, seed=4, jobs=JOBS)
            params.append(dict(dataset=key, metric=m, n=len(g), T_min=g["T"].min(),
                               T_max=g["T"].max(), **fit.as_dict()))
            clim = {**{f"own {ev}": own[ev] for ev in EVENTS},
                    **{f"ERA5 {ev} warming level": e_clim[str(ev)] for ev in EVENTS}}
            if key == "MPI-GE":
                clim["1850-1900"] = float(g[g.year.between(*PREIND)]["T"].mean())
            for ev in EVENTS:
                for c, T in clim.items():
                    x = e_val[(m, ev)]
                    lo, hi = ex.ci_log([b.return_period(x, T) for b in boot])
                    add(key, "NS-GEV", c, m, ev, fit.return_period(x, T), lo, hi, len(g), T=T)
            print(f"  {key} {m}: NS-GEV {fit.as_dict()}  T {g['T'].min():+.2f}..{g['T'].max():+.2f}")

    # control runs: stationary, drift removed (mean kept), plus empirical counts
    for key in CONTROLS:
        g = t[t.dataset == key].sort_values("year")
        for m in METRICS:
            x = g[m].values
            x = x - np.polyval(np.polyfit(g.year, x, 1), g.year) + x.mean()
            fit = ex.fit_nsgev(x)
            boot = ex.bootstrap_nsgev(x, n=N_BOOT, seed=5, jobs=JOBS)
            params.append(dict(dataset=key, metric=m, n=len(x), **fit.as_dict()))
            for ev in EVENTS:
                xe = e_val[(m, ev)]
                lo, hi = ex.ci_log([b.return_period(xe, 0.0) for b in boot])
                add(key, "GEV", "control", m, ev, fit.return_period(xe, 0.0), lo, hi, len(x))
                r = ex.empirical_return_period(x, xe)
                add(key, "empirical", "control", m, ev, r["rp"], r["rp_lo"], r["rp_hi"],
                    r["n"], k=r["k"])

    df = pd.DataFrame(rows)
    df.to_csv(DIR / "model_event_rp.csv", index=False, float_format="%.5g")
    pd.DataFrame(params).to_csv(DIR / "model_gev_params.csv", index=False, float_format="%.5g")

    # reference-period climatology of each dataset, for reading the absolute values
    ref = (t[t.year.between(1991, 2020)].groupby("dataset", sort=False)[METRICS]
           .agg(["mean", "std"]))
    ref.columns = [f"{a}_{b}" for a, b in ref.columns]
    ref.to_csv(DIR / "reference_1991_2020.csv", float_format="%.4g")
    print(ref.round(2).to_string())
    return df


def main():
    t = load()
    print("== A. ERA5 nonstationary GEV")
    era5_part(t)
    print("== B. anomalies from the forced mean")
    anomaly_part(t)
    print("== C. observed values in the simulated climates")
    model_part(t)
    print(f"wrote tables to {DIR}")


if __name__ == "__main__":
    main()
