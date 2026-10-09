"""
Return period of every ERA5 summer in the MPI-GE climate of its year.

A GEV describes a stationary sample, so it is fitted where stationarity
holds: across the 50 MPI-GE members inside a 10-year window, 500 seasons of
one background climate. "1-in-100 years" then means 1 in 100 seasons of that
climate. Rolling the window through time gives a stationary GEV -- a
return-level curve -- for every background climate, and each ERA5 season is
scored against the window centred on its own calendar year.

  window   MPI-GE years y-5 .. y+4 for ERA5 year y (10 x 50 = 500 seasons);
           the trend of the ensemble mean inside the window is reported as a
           check on the no-trend assumption, not removed
  fit      stationary GEV, src.extremes.fit_nsgev (|xi| <= 0.5)
  RP       1 / P(X >= ERA5 value) under that fit; 95% interval from N_BOOT
           bootstrap refits of the 500 seasons; infinite if the value lies
           above a bounded (xi < 0) fit's upper end
  check    empirical k of 500 with a Clopper-Pearson interval -- no fit

Index: land-mean JJA heatwave days over src.heatwave.MEAN_BOX (as fig 17),
each dataset against its own 1991-2020 threshold. MPI-GE is used as is: it
is ~1.4x as variable as ERA5 (1991-2020 std 3.7 vs 2.7 days), which makes
ERA5 values look more frequent (shorter return periods); and it warms more
slowly than observed after ~2015, so recent ERA5 years are scored against a
cooler model climate.

Input : data/heatwave_stats/season_table.csv (01_season_table.py)
Output: data/heatwave_stats/mpige_window_gev.csv   one row per window
        data/heatwave_stats/era5_rp_in_mpige.csv   one row per ERA5 season
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
import src.extremes as ex

DIR = ROOT / "data" / "heatwave_stats"
METRIC = "hwd"
#: window offsets around the ERA5 year: y-5 .. y+4
WIN_BEFORE, WIN_AFTER = 5, 4
RETURN_PERIODS = [2, 5, 10, 20, 50, 100, 500]
N_BOOT = int(os.environ.get("HWS_N_BOOT", 1000))
JOBS = min(64, int(os.environ.get("SLURM_CPUS_ON_NODE", 8)))


def window(mpige, year):
    return mpige[mpige.year.between(year - WIN_BEFORE, year + WIN_AFTER)]


def main():
    t = pd.read_csv(DIR / "season_table.csv")
    mpige = t[t.dataset == "MPI-GE"]
    era5 = t[t.dataset == "ERA5"].sort_values("year").set_index("year")[METRIC]

    gev_rows, rp_rows = [], []
    for year, xe in era5.items():
        w = window(mpige, year)
        x = w[METRIC].values
        fit = ex.fit_nsgev(x)
        boot = ex.bootstrap_nsgev(x, n=N_BOOT, seed=year, start=fit, jobs=JOBS)
        ens = w.groupby("year")[METRIC].mean()
        trend = np.polyfit(ens.index, ens.values, 1)[0] * (WIN_BEFORE + WIN_AFTER)

        gev_rows.append(dict(
            year=year, first=year - WIN_BEFORE, last=year + WIN_AFTER, n=len(x),
            mu=fit.mu0, sigma=float(fit.scale(0)), xi=fit.xi,
            mean=x.mean(), std=x.std(ddof=1), max=x.max(), trend_in_window=trend,
            **{f"rl{rp}": float(fit.return_level(rp, 0)) for rp in RETURN_PERIODS},
            **{f"rl{rp}_lo": ex.ci([b.return_level(rp, 0) for b in boot])[0]
               for rp in RETURN_PERIODS},
            **{f"rl{rp}_hi": ex.ci([b.return_level(rp, 0) for b in boot])[1]
               for rp in RETURN_PERIODS}))

        rp = float(fit.return_period(xe, 0))
        lo, hi = ex.ci_log([b.return_period(xe, 0) for b in boot])
        emp = ex.empirical_return_period(x, xe)
        rp_rows.append(dict(year=year, value=xe, rp=rp, rp_lo=lo, rp_hi=hi,
                            above_upper_end=bool(np.isinf(rp)), k=emp["k"], n=emp["n"],
                            rp_emp=emp["rp"], rp_emp_lo=emp["rp_lo"], rp_emp_hi=emp["rp_hi"]))

    gev = pd.DataFrame(gev_rows)
    rps = pd.DataFrame(rp_rows)
    gev.to_csv(DIR / "mpige_window_gev.csv", index=False, float_format="%.5g")
    rps.to_csv(DIR / "era5_rp_in_mpige.csv", index=False, float_format="%.5g")

    assert (gev.n == 500).all(), "every window should hold 10 years x 50 members"
    print(f"  {len(gev)} windows of {gev.n.iloc[0]} seasons, {N_BOOT} bootstrap refits each")
    print(f"  ensemble-mean trend inside the window: max {gev.trend_in_window.abs().max():.2f} "
          f"days (window std {gev['std'].min():.1f}-{gev['std'].max():.1f})")
    top = rps.replace(np.inf, 1e9).sort_values("rp", ascending=False).head(8)
    for _, r in top.iterrows():
        rp_txt = "inf" if r.above_upper_end else f"{r.rp:.0f}"
        print(f"  {int(r.year)}: {r.value:5.1f} days  RP {rp_txt:>6s} yr "
              f"({r.rp_lo:.0f}-{r.rp_hi:.0f})  empirical {int(r.k)}/{int(r.n)}")
    print(f"wrote {DIR / 'mpige_window_gev.csv'}, {DIR / 'era5_rp_in_mpige.csv'}")


if __name__ == "__main__":
    main()
