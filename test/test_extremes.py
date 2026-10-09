"""Tests for src/extremes.py."""
import numpy as np
import pytest
from scipy import stats

import src.extremes as ex


def test_smooth_covariate_trailing():
    years = np.arange(2000, 2006)
    v = np.array([0, 4, 8, 12, 16, 20], float)
    s = ex.smooth_covariate(v, years, window=4)
    np.testing.assert_allclose(s, [0, 2, 4, 6, 10, 14])


def test_smooth_covariate_keeps_input_order():
    years = np.array([2002, 2000, 2001])
    s = ex.smooth_covariate([2.0, 0.0, 1.0], years, window=2)
    np.testing.assert_allclose(s, [1.5, 0.0, 0.5])


def test_stationary_fit_recovers_parameters():
    rng = np.random.default_rng(1)
    x = stats.genextreme(-0.1, loc=5, scale=2).rvs(4000, random_state=rng)
    f = ex.fit_nsgev(x)
    assert f.mu1 == 0 and f.s1 == 0
    assert f.loc(0) == pytest.approx(5, abs=0.1)
    assert f.scale(0) == pytest.approx(2, rel=0.05)
    assert f.xi == pytest.approx(0.1, abs=0.04)


def test_nonstationary_fit_recovers_trends():
    rng = np.random.default_rng(2)
    t = np.repeat(np.linspace(-0.5, 1.0, 40), 100)
    mu, sigma = 4 + 6 * t, np.exp(0.5 + 0.4 * t)
    x = stats.genextreme(0.1, loc=mu, scale=sigma).rvs(random_state=rng)
    f = ex.fit_nsgev(x, t)
    assert f.mu1 == pytest.approx(6, abs=0.3)
    assert f.s1 == pytest.approx(0.4, abs=0.08)
    assert f.xi == pytest.approx(-0.1, abs=0.04)


def test_return_period_and_level_are_inverse():
    f = ex.NSGEV(mu0=5, mu1=3, s0=0.5, s1=0.2, xi=0.05)
    for t in (0.0, 0.5):
        x = f.return_level(50, t)
        assert f.return_period(x, t) == pytest.approx(50)
    # warming makes a fixed level more frequent
    assert f.return_period(15, 0.5) < f.return_period(15, 0.0)


def test_shape_is_bounded():
    rng = np.random.default_rng(3)
    x = rng.standard_cauchy(60) ** 2
    assert abs(ex.fit_nsgev(x).xi) <= ex.XI_BOUND


def test_exceedance_interval():
    k, n, p, lo, hi = ex.exceedance(np.arange(100), 90)
    assert (k, n) == (10, 100) and p == pytest.approx(0.1)
    assert lo < 0.1 < hi


def test_empirical_return_period_without_exceedance_is_a_bound():
    r = ex.empirical_return_period(np.zeros(500), 1.0)
    assert r["k"] == 0 and np.isinf(r["rp"]) and np.isinf(r["rp_hi"])
    # one-sided 97.5% bound for k=0: 1 - 0.025**(1/n), rp_lo ~ n / 3.7
    assert 100 < r["rp_lo"] < 200


def test_plotting_positions():
    s, rp = ex.plotting_positions([3.0, 1.0, 2.0, np.nan])
    np.testing.assert_allclose(s, [3, 2, 1])
    np.testing.assert_allclose(rp, [4, 2, 4 / 3])


def test_numpy_logpdf_matches_scipy():
    x = np.linspace(-2, 12, 50)
    for xi in (-0.3, 0.0, 0.2):
        ref = stats.genextreme.logpdf(x, -xi, loc=3, scale=1.5)
        np.testing.assert_allclose(ex.gev_logpdf(x, 3, 1.5, xi), ref, rtol=1e-6, atol=1e-9)


def test_bootstrap_parallel_matches_serial():
    rng = np.random.default_rng(4)
    t = np.linspace(-0.4, 0.6, 80)
    x = stats.genextreme(0.1, loc=5 + 4 * t, scale=1.5).rvs(random_state=rng)
    a = ex.bootstrap_nsgev(x, t, n=12, seed=7)
    b = ex.bootstrap_nsgev(x, t, n=12, seed=7, jobs=3)
    assert len(a) == len(b) == 12
    for fa, fb in zip(a, b):
        assert fa.mu1 == pytest.approx(fb.mu1)


def test_probability_ratio_edge_cases():
    pr = ex.probability_ratio([0.1, 0.1, 0.0], [0.01, 0.0, 0.0])
    assert pr[0] == pytest.approx(10) and np.isinf(pr[1]) and np.isnan(pr[2])
    lo, hi = ex.ci_log([1.0, np.inf, np.nan, 3.0], cap=100)
    assert hi == pytest.approx(100, rel=0.1) and lo >= 1.0
