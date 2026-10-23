"""
As 03_mpige_window_rp.py, but without any fitted distribution and without
return periods: where each ERA5 summer sits, from small to large, among the
500 MPI-GE summers of its window (50 members x the years y-5 .. y+4).

  below    MPI-GE summers in the window below the ERA5 value
  equal    MPI-GE summers equal to it (heatwave days has many ties, mostly 0)
  pct      percentile rank 100 * (below + equal / 2) / n: 0% under every
           MPI-GE summer, 100% above all of them
  95% CI   Clopper-Pearson: lower bound from `below`, upper from
           `below + equal`, so ties widen the interval instead of hiding in it
  levels   the window's 50th / 90th / 99th percentiles (q50 / q90 / q99),
           exceeded by half, a tenth and a hundredth of its summers
  fixed    the percentile rank of the 2003 and 2026 ERA5 values in every
           window (the cross-climate curves of the figures)

Index: argument `hwd` (default) or `txx`, as for 03.
Output: data/heatwave_stats/mpige_window_rank<suffix>.csv  one row per window
        data/heatwave_stats/era5_rank_in_mpige<suffix>.csv one row per ERA5 season
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

w03 = importlib.import_module("03_mpige_window_rp")

DIR = ROOT / "data" / "heatwave_stats"
QUANTILES = [50, 90, 99]
FIXED_YEARS = [2003, 2026]


def percentile_rank(sample, value, level=0.95):
    """Percentile rank (%) of `value` among `sample`, ties counted half, with an interval."""
    s = np.asarray(sample, float)
    s = s[np.isfinite(s)]
    n, below, equal = len(s), int(np.sum(s < value)), int(np.sum(s == value))
    a = 1 - level
    lo = 0.0 if below == 0 else stats.beta.ppf(a / 2, below, n - below + 1)
    k = below + equal
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return dict(below=below, equal=equal, n=n, pct=100 * (below + equal / 2) / n,
                pct_lo=100 * float(lo), pct_hi=100 * float(hi))


def main(index="hwd"):
    table, col, suffix = w03.METRICS[index]
    t = pd.read_csv(DIR / table)
    mpige = t[t.dataset == "MPI-GE"]
    era5 = t[t.dataset == "ERA5"].sort_values("year").set_index("year")[col]
    fixed_vals = {y: float(era5.loc[y]) for y in FIXED_YEARS}

    win_rows, pct_rows = [], []
    for year, xe in era5.items():
        x = w03.window(mpige, year)[col].values
        win_rows.append(dict(
            year=year, n=len(x),
            **{f"q{q}": float(np.percentile(x, q)) for q in QUANTILES},
            **{f"pct_fixed_{y}": percentile_rank(x, v)["pct"] for y, v in fixed_vals.items()}))
        pct_rows.append(dict(year=year, value=xe, **percentile_rank(x, xe)))

    win = pd.DataFrame(win_rows)
    pcts = pd.DataFrame(pct_rows)
    assert (win.n == 500).all(), "every window should hold 10 years x 50 members"
    win.to_csv(DIR / f"mpige_window_rank{suffix}.csv", index=False, float_format="%.5g")
    pcts.to_csv(DIR / f"era5_rank_in_mpige{suffix}.csv", index=False, float_format="%.5g")
    print(f"  [{index}] {len(win)} windows of 500 summers; ERA5 seasons above the whole "
          f"ensemble: {list(pcts.year[pcts.below == pcts.n])}")
    for _, r in pcts.sort_values("pct", ascending=False).head(6).iterrows():
        print(f"  {int(r.year)}: {r.value:6.2f}  below {int(r.below):3d}  equal {int(r.equal):3d}"
              f"  {r.pct:5.1f}% ({r.pct_lo:.1f}-{r.pct_hi:.1f})")
    for y in FIXED_YEARS:
        print(f"  {y} value in the windows of 1950s / 2003 / 2026: " + " / ".join(
            f"{win.loc[win.year.between(a, b), f'pct_fixed_{y}'].mean():.1f}%"
            for a, b in ((1950, 1959), (2003, 2003), (2026, 2026))))
    print(f"wrote {DIR / f'mpige_window_rank{suffix}.csv'}, {DIR / f'era5_rank_in_mpige{suffix}.csv'}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
