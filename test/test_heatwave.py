"""Tests for the heatwave definition in src/heatwave.py."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.heatwave as hw


def test_noleap_doy_ignores_feb29():
    t = pd.to_datetime(["2003-05-01", "2004-05-01", "2004-02-29", "2026-09-30"])
    assert list(hw.noleap_doy(t)) == [121, 121, 60, 273]


def test_run_length():
    m = np.array([0, 1, 1, 0, 1, 1, 1, 1, 0, 1], bool)[:, None]
    assert list(hw.run_length(m)[:, 0]) == [0, 2, 2, 0, 4, 4, 4, 4, 0, 1]


def _year(start="2003-04-01", end="2003-10-31"):
    t = pd.date_range(start, end, freq="D")
    return t, hw.noleap_doy(t)


def test_three_day_minimum_and_jja_counting():
    t, doy = _year()
    a = np.zeros((len(t), 4))
    thr = np.ones_like(a)
    ix = {d: i for i, d in enumerate(t.strftime("%m-%d"))}

    # point 0: 2-day spell in July -> not a heatwave
    a[ix["07-10"]:ix["07-12"], 0] = 2
    # point 1: 5-day spell in July, anomaly 3 -> 5 days, 1 event, cumulative 15
    a[ix["07-10"]:ix["07-15"], 1] = 3
    # point 2: event 29 May - 3 Jun: 3 JJA days, 1 event, full length 6
    a[ix["05-29"]:ix["06-04"], 2] = 2
    # point 3: exceedance in April only (outside detection window) -> nothing
    a[ix["04-10"]:ix["04-20"], 3] = 5

    m = hw.season_metrics(a, thr, doy)
    assert list(m["hwd"]) == [0, 5, 3, 0]
    assert list(m["hwn"]) == [0, 1, 1, 0]
    assert np.allclose(m["hwcum"], [0, 15, 6, 0])
    assert list(m["hwmax"]) == [0, 5, 6, 0]
    # threshold is 1: excess = (3-1)*5 at point 1, (2-1)*3 JJA days at point 2
    assert np.allclose(m["hwexcess"], [0, 10, 3, 0])
    assert np.allclose(m["hwpeak"][1:3], [3, 2])
    assert np.isnan(m["hwpeak"][[0, 3]]).all()
    assert m["onset"][1] == hw.noleap_doy(pd.to_datetime(["2003-07-10"]))[0]
    assert np.isnan(m["onset"][0])


def test_event_starting_on_june_first_counted_once():
    t, doy = _year()
    a = np.zeros((len(t), 1))
    i = list(t.strftime("%m-%d")).index("06-01")
    a[i:i + 4, 0] = 2
    m = hw.season_metrics(a, np.ones_like(a), doy)
    assert m["hwn"][0] == 1 and m["hwd"][0] == 4


def test_event_cut_at_window_end():
    """Exceedances running past 30 Sep only count up to the window edge."""
    t, doy = _year()
    a = np.zeros((len(t), 1))
    i = list(t.strftime("%m-%d")).index("09-29")
    a[i:i + 10, 0] = 2          # 29 Sep .. 8 Oct: only 2 days inside the window
    m = hw.season_metrics(a, np.ones_like(a), doy)
    assert np.isnan(m["onset"][0])


def test_threshold_is_90th_percentile_of_15_day_pool():
    years = range(1991, 2021)
    t = pd.DatetimeIndex(np.concatenate(
        [pd.date_range(f"{y}-04-01", f"{y}-10-31") for y in years]))
    doy = hw.noleap_doy(t)
    rng = np.random.default_rng(0)
    tmax = rng.normal(20, 3, size=(len(t), 2))
    ref = np.ones(len(t), bool)

    clim = hw.daily_climatology(tmax, doy, ref)
    anom = tmax - clim[doy]
    thr = hw.percentile_threshold(anom, doy, ref, [182])

    pool = anom[np.abs(doy - 182) <= 7]
    assert pool.shape[0] == 30 * 15
    assert np.allclose(thr[182], np.percentile(pool, 90, axis=0))
    assert np.isnan(thr[181]).all()
