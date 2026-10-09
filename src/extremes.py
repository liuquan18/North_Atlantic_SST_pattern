"""
Extreme-value tools for the heatwave statistics (scripts/heatwave_statistics/).

The seasonal heatwave indices (land-mean JJA heatwave days, cumulative heat)
are right-skewed and their spread grows as the climate warms -- with the
threshold fixed at 1991-2020, MPI-GE's heatwave-day standard deviation rises
from 2.4 days in 1850-1900 to 5.4 days in 2016-2036. A fit that only shifts
the distribution would therefore overstate how rare a recent season is. The
nonstationary GEV here lets location AND scale follow a warming covariate
(smoothed global-mean SST anomaly, the "GMST" of attribution studies):

    mu(T)    = mu0 + mu1 * T
    sigma(T) = exp(s0 + s1 * T)
    xi       constant

Return periods are 1 / P(X >= x | T). GEV shape is bounded to |xi| <= 0.5:
with 36-90 seasons an unconstrained xi wanders to values no physical
heatwave index has.

Empirical exceedance probabilities (large ensemble, control runs) carry a
Clopper-Pearson interval, so "no season exceeded" is reported as a lower
bound on the return period rather than infinity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

XI_BOUND = 0.5
SMOOTH_YEARS = 4   # trailing running mean of the covariate, as in WWA-type studies


# --------------------------------------------------------------------------
# covariate
# --------------------------------------------------------------------------
def smooth_covariate(values, years, window: int = SMOOTH_YEARS) -> np.ndarray:
    """Trailing `window`-year running mean (fewer years at the start of a record)."""
    s = pd.Series(np.asarray(values, float), index=np.asarray(years))
    return s.sort_index().rolling(window, min_periods=1).mean().reindex(s.index).values


# --------------------------------------------------------------------------
# nonstationary GEV
# --------------------------------------------------------------------------
@dataclass
class NSGEV:
    mu0: float
    mu1: float
    s0: float
    s1: float
    xi: float

    def loc(self, t):
        return self.mu0 + self.mu1 * np.asarray(t, float)

    def scale(self, t):
        return np.exp(self.s0 + self.s1 * np.asarray(t, float))

    def _dist(self, t):
        # scipy's genextreme uses c = -xi
        return stats.genextreme(-self.xi, loc=self.loc(t), scale=self.scale(t))

    def sf(self, x, t):
        return self._dist(t).sf(x)

    def cdf(self, x, t):
        return self._dist(t).cdf(x)

    def return_period(self, x, t):
        with np.errstate(divide="ignore"):
            return 1.0 / self.sf(x, t)

    def return_level(self, period, t):
        return self._dist(t).isf(1.0 / np.asarray(period, float))

    def as_dict(self):
        return dict(mu0=self.mu0, mu1=self.mu1, s0=self.s0, s1=self.s1, xi=self.xi)


def gev_logpdf(x, mu, sigma, xi):
    """GEV log-density in the climate convention (xi > 0: heavy upper tail).

    Plain numpy -- scipy's frozen-distribution call costs more than the
    arithmetic for the ~10^4-season MPI-GE fits. -inf outside the support.
    """
    z = (x - mu) / sigma
    if abs(xi) < 1e-6:
        return -np.log(sigma) - z - np.exp(-z)
    tz = 1.0 + xi * z
    with np.errstate(invalid="ignore", divide="ignore"):
        out = -np.log(sigma) - (1.0 + 1.0 / xi) * np.log(tz) - tz ** (-1.0 / xi)
    return np.where(tz > 0, out, -np.inf)


def _nll(theta, x, t, scale_trend):
    mu0, mu1, s0, s1, xi = theta if scale_trend else (*theta[:3], 0.0, theta[3])
    if abs(xi) > XI_BOUND:
        return np.inf
    ll = gev_logpdf(x, mu0 + mu1 * t, np.exp(s0 + s1 * t), xi)
    if not np.all(np.isfinite(ll)):
        return np.inf
    return -ll.sum()


def fit_nsgev(x, t=None, *, scale_trend: bool = True, start: NSGEV | None = None) -> NSGEV:
    """
    Maximum-likelihood fit. `t=None` (or scale_trend with a constant t) gives a
    stationary GEV. Several starting shapes are tried and the best likelihood
    wins; `start` (e.g. the full-sample fit, for a bootstrap refit) is tried
    first and alone if it converges.
    """
    x = np.asarray(x, float)
    t = np.zeros_like(x) if t is None else np.asarray(t, float)
    stationary = np.ptp(t) == 0
    if stationary:
        scale_trend = False
    # start: linear regression for the location, residual spread for the scale
    if stationary:
        a, b = 0.0, x.mean()
    else:
        a, b = np.polyfit(t, x, 1)
    resid = x - (a * t + b)
    s_start = np.log(max(resid.std() * np.sqrt(6) / np.pi, 1e-3))
    starts = []
    if start is not None:
        starts.append([start.mu0, start.mu1, start.s0, start.s1, start.xi] if scale_trend
                      else [start.mu0, start.mu1, start.s0, start.xi])
    for xi0 in (-0.2, 0.0, 0.15):
        mu_start = b - 0.5772 * np.exp(s_start) if xi0 == 0 else b
        starts.append([mu_start, a, s_start, 0.0, xi0] if scale_trend
                      else [mu_start, a, s_start, xi0])
    best = None
    for k, th0 in enumerate(starts):
        if stationary:
            th0[1] = 0.0
        # a start that leaves seasons outside the GEV support has an infinite
        # likelihood everywhere around it: the simplex would only wander
        if not np.isfinite(_nll(th0, x, t, scale_trend)):
            continue
        # the simplex meets the inf barrier at |xi| > XI_BOUND: harmless inf - inf
        with np.errstate(invalid="ignore"):
            r = optimize.minimize(_nll, th0, args=(x, t, scale_trend), method="Nelder-Mead",
                                  options=dict(maxiter=5000, maxfev=5000, xatol=1e-7,
                                               fatol=1e-9))
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
            best = r
        if start is not None and k == 0 and r.success and np.isfinite(r.fun):
            break
    if best is None:
        raise RuntimeError("GEV fit failed for every starting point")
    th = best.x
    if not scale_trend:
        th = [th[0], th[1], th[2], 0.0, th[3]]
    if stationary:
        th[1] = 0.0
    return NSGEV(*map(float, th))


def _refit(args):
    x, t, scale_trend, start = args
    try:
        return fit_nsgev(x, t, scale_trend=scale_trend, start=start)
    except RuntimeError:
        return None


def bootstrap_nsgev(x, t=None, *, n: int = 1000, seed: int = 0, scale_trend: bool = True,
                    start: NSGEV | None = None, jobs: int = 1) -> list[NSGEV]:
    """
    Refits on seasons resampled with replacement (pairs of index and covariate),
    warm-started from `start` (default: the full-sample fit), over `jobs`
    processes. The resamples depend only on `seed`, not on `jobs`.
    """
    x = np.asarray(x, float)
    t = np.zeros_like(x) if t is None else np.asarray(t, float)
    if start is None:
        start = fit_nsgev(x, t, scale_trend=scale_trend)
    rng = np.random.default_rng(seed)
    idx = [rng.integers(0, len(x), len(x)) for _ in range(n)]
    tasks = [(x[i], t[i], scale_trend, start) for i in idx]
    if jobs > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(jobs) as pool:
            fits = list(pool.map(_refit, tasks, chunksize=max(1, n // (4 * jobs))))
    else:
        fits = [_refit(a) for a in tasks]
    return [f for f in fits if f is not None]


def ci(values, level: float = 0.95):
    """Central interval of a bootstrap sample, ignoring non-finite values."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)] if np.isfinite(v).any() else v
    q = (1 - level) / 2
    return float(np.nanquantile(v, q)), float(np.nanquantile(v, 1 - q))


def ci_log(values, level: float = 0.95, cap: float = 1e6):
    """
    As ci(), for return periods and ratios: infinities are capped at `cap`,
    undefined values (NaN, e.g. a 0/0 probability ratio) are dropped.
    """
    v = np.asarray(values, float)
    v = np.clip(v[~np.isnan(v)], None, cap)
    if len(v) == 0:
        return np.nan, np.nan
    q = (1 - level) / 2
    return float(np.quantile(v, q)), float(np.quantile(v, 1 - q))


def probability_ratio(p_a, p_b):
    """p_a / p_b elementwise; inf where only p_b is 0, NaN where both are."""
    p_a, p_b = np.asarray(p_a, float), np.asarray(p_b, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((p_a == 0) & (p_b == 0), np.nan, p_a / p_b)


# --------------------------------------------------------------------------
# empirical exceedance
# --------------------------------------------------------------------------
def exceedance(sample, threshold, level: float = 0.95):
    """
    P(X >= threshold) in a sample, with a Clopper-Pearson interval.
    Returns (k, n, p, p_lo, p_hi).
    """
    s = np.asarray(sample, float)
    s = s[np.isfinite(s)]
    n, k = len(s), int(np.sum(s >= threshold))
    a = 1 - level
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return k, n, k / n, float(lo), float(hi)


def empirical_return_period(sample, threshold, level: float = 0.95):
    """
    Return period 1/p with its interval from exceedance(); with no exceedance
    the best estimate is NaN and the interval is (1/p_hi, inf).
    """
    k, n, p, lo, hi = exceedance(sample, threshold, level)
    rp = np.inf if k == 0 else 1 / p
    return dict(k=k, n=n, p=p, rp=rp, rp_lo=1 / hi, rp_hi=np.inf if lo == 0 else 1 / lo)


def plotting_positions(x):
    """Sorted values and their empirical return periods, (n + 1) / rank (Weibull)."""
    s = np.sort(np.asarray(x, float)[np.isfinite(x)])[::-1]
    rp = (len(s) + 1) / np.arange(1, len(s) + 1)
    return s, rp
