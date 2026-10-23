"""
As 03_mpige_window_rp.py, but without any fitted distribution: the return
period of each ERA5 summer comes purely from its rank among the 500 MPI-GE
summers of its window (50 members x the years y-5 .. y+4).

  k        MPI-GE summers in the window at or above the ERA5 value
  rank     k + 1: the ERA5 summer's place, from the top, among the 501
           values (the 500 plus itself)
  RP       (n + 1) / rank = 501 / (k + 1) years -- the Weibull plotting
           position; finite even for k = 0, where it is 501 years, "beyond
           the whole ensemble" (a lower bound in effect)
  95% CI   from the Clopper-Pearson interval of k / n
  levels   the 1-in-2 / 10 / 100 year heatwave days or TXx of a window are its
           50th / 90th / 99th percentiles
  fixed    the same rank-based return period of the 2003 and 2026 ERA5
           values in every window (the cross-climate curves of the figures)

Index: argument `hwd` (default) or `txx`, as for 03.
Output: data/heatwave_stats/mpige_window_rank<suffix>.csv  one row per window
        data/heatwave_stats/era5_rank_in_mpige<suffix>.csv one row per ERA5 season
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

import src.extremes as ex

w03 = importlib.import_module("03_mpige_window_rp")

DIR = ROOT / "data" / "heatwave_stats"
RETURN_PERIODS = [2, 10, 100]
FIXED_YEARS = [2003, 2026]


def rank_rp(sample, value):
    """Rank-based return period of `value` among `sample`, with a 95% interval."""
    k, n, _, p_lo, p_hi = ex.exceedance(sample, value)
    rp = (n + 1) / (k + 1)
    return dict(k=k, n=n, rp=rp, rp_lo=1 / p_hi, rp_hi=np.inf if p_lo == 0 else 1 / p_lo)


def main(index="hwd"):
    table, col, suffix = w03.METRICS[index]
    t = pd.read_csv(DIR / table)
    mpige = t[t.dataset == "MPI-GE"]
    era5 = t[t.dataset == "ERA5"].sort_values("year").set_index("year")[col]
    fixed_vals = {y: float(era5.loc[y]) for y in FIXED_YEARS}

    win_rows, rp_rows = [], []
    for year, xe in era5.items():
        x = w03.window(mpige, year)[col].values
        win_rows.append(dict(
            year=year, n=len(x),
            **{f"rl{rp}": float(np.quantile(x, 1 - 1 / rp)) for rp in RETURN_PERIODS},
            **{f"rp_fixed_{y}": rank_rp(x, v)["rp"] for y, v in fixed_vals.items()},
            **{f"k_fixed_{y}": rank_rp(x, v)["k"] for y, v in fixed_vals.items()}))
        rp_rows.append(dict(year=year, value=xe, **rank_rp(x, xe)))

    win = pd.DataFrame(win_rows)
    rps = pd.DataFrame(rp_rows)
    assert (win.n == 500).all(), "every window should hold 10 years x 50 members"
    win.to_csv(DIR / f"mpige_window_rank{suffix}.csv", index=False, float_format="%.5g")
    rps.to_csv(DIR / f"era5_rank_in_mpige{suffix}.csv", index=False, float_format="%.5g")
    print(f"  [{index}] {len(win)} windows of 500 summers; ERA5 seasons beyond the whole "
          f"ensemble (k = 0): {list(rps.year[rps.k == 0])}")
    for _, r in rps.sort_values("rp", ascending=False).head(6).iterrows():
        print(f"  {int(r.year)}: {r.value:6.2f}  k = {int(r.k):3d}  RP {r.rp:6.1f} yr "
              f"({r.rp_lo:.0f}-{r.rp_hi:.0f})")
    for y in FIXED_YEARS:
        print(f"  {y} value in the windows of 1950s / 2003 / 2026: " + " / ".join(
            f"{win.loc[win.year.between(a, b), f'rp_fixed_{y}'].pipe(lambda s: np.exp(np.log(s).mean())):.0f} yr"
            for a, b in ((1950, 1959), (2003, 2003), (2026, 2026))))
    print(f"wrote {DIR / f'mpige_window_rank{suffix}.csv'}, {DIR / f'era5_rank_in_mpige{suffix}.csv'}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
