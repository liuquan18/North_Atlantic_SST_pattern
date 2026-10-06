"""Tests for the 2026 pattern-correlation machinery in src/pattern_2026.py."""
import sys

import numpy as np
import pytest
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2


@pytest.fixture
def grid():
    """The common 0.25 deg analysis grid: 30-60N, 80W-40E."""
    lat = np.arange(30.125, 60, 0.25)
    lon = np.arange(-79.875, 40, 0.25)
    return lat, lon


@pytest.fixture
def ones(grid):
    lat, lon = grid
    return xr.DataArray(np.ones((len(lat), len(lon))),
                        coords=[("lat", lat), ("lon", lon)])


def test_coarsen_preserves_constant(ones):
    c = p2.coarsen_to_1deg(ones)
    assert c.shape == (30, 120)
    assert np.allclose(c.values, 1.0)


def test_coarsen_is_area_weighted(ones):
    """A field equal to cos(lat) must coarsen to the cos-weighted block mean."""
    field = ones * np.cos(np.deg2rad(ones.lat))
    coarse = p2.coarsen_to_1deg(field)

    block = field.isel(lat=slice(0, 4), lon=slice(0, 4))
    w = np.cos(np.deg2rad(block.lat)) * xr.ones_like(block)
    expected = float((block * w).sum() / w.sum())
    assert float(coarse.values[0, 0]) == pytest.approx(expected, rel=1e-9)


def test_coarsen_drops_cells_below_ocean_fraction(ones):
    """A 1 deg cell that is mostly land must not enter the correlation."""
    mostly_land = ones.where(~((ones.lat < 31.0) & (ones.lon < -79.0)))
    assert np.isnan(p2.coarsen_to_1deg(mostly_land).values[0, 0])

    # 12 of 16 sub-cells valid is above the 0.5 threshold and must survive
    partial = ones.copy(deep=True)
    partial[0:1, 0:4] = np.nan
    assert np.isfinite(p2.coarsen_to_1deg(partial).values[0, 0])


@pytest.fixture
def noise(grid):
    lat, lon = grid
    rng = np.random.RandomState(0)
    return xr.DataArray(rng.randn(len(lat), len(lon)),
                        coords=[("lat", lat), ("lon", lon)])


def test_centered_correlation_bounds(noise):
    assert float(p2.pattern_corr(noise, noise, centered=True)) == pytest.approx(1.0)
    assert float(p2.pattern_corr(noise, -noise, centered=True)) == pytest.approx(-1.0)


def test_centered_correlation_ignores_offsets(noise):
    """The 'spatial mean removed' metric must not care about a uniform offset."""
    assert float(p2.pattern_corr(noise + 5, noise, centered=True)) == pytest.approx(1.0)


def test_uncentered_correlation_sees_offsets(noise):
    """The 'global mean removed' metric keeps basin-mean warmth as signal."""
    shifted = float(p2.pattern_corr(noise + 5, noise, centered=False))
    assert abs(shifted - 1.0) > 0.1


def test_correlation_broadcasts_over_member_and_year(noise):
    """One call must score the whole ensemble-by-year matrix."""
    stack = noise.expand_dims(year=[2000, 2001, 2002]).expand_dims(member=[1, 2])
    r = p2.pattern_corr(stack, noise, centered=True)
    assert set(r.dims) == {"member", "year"}
    assert np.allclose(r.values, 1.0)


def test_remove_spatial_mean_zeroes_the_area_mean(noise):
    out = p2.remove_spatial_mean(noise)
    w = np.cos(np.deg2rad(out.lat)) * xr.ones_like(out.lon)
    assert float((out * w).sum() / w.sum()) == pytest.approx(0.0, abs=1e-12)


def test_common_ocean_mask_is_an_intersection(ones):
    a = ones.where(ones.lon > -70)
    b = ones.where(ones.lat < 55)
    mask = p2.common_ocean_mask({"a": a, "b": b})
    assert bool(mask.sel(lon=0.125, lat=40.125, method="nearest"))
    assert not bool(mask.sel(lon=-75.125, lat=40.125, method="nearest"))
    assert not bool(mask.sel(lon=0.125, lat=58.125, method="nearest"))


# --- event composites --------------------------------------------------------
def test_event_onsets_counts_a_spell_once():
    x = np.array([0.1, 0.6, 0.7, 0.2, 0.5, np.nan, 0.9, 0.1])
    # 0.6 begins a two-season spell; 0.5 is a fresh onset; 0.9 follows a gap,
    # so whether it is an onset is unknown and it is skipped
    assert list(p2.event_onsets(x, 0.5)) == [1, 4]


def test_event_onsets_skips_a_spell_open_at_the_start():
    assert list(p2.event_onsets(np.array([0.8, 0.9, 0.1, 0.6]), 0.5)) == [3]


def test_epoch_windows_pad_past_the_record():
    w = p2.epoch_windows(np.arange(5.0), np.array([0, 4]), 2)
    assert np.array_equal(w[0], [np.nan, np.nan, 0, 1, 2], equal_nan=True)
    assert np.array_equal(w[1], [2, 3, 4, np.nan, np.nan], equal_nan=True)


def test_internal_component_large_ensemble_removes_member_mean():
    years = np.arange(1900, 2000)
    rng = np.random.default_rng(1)
    r = xr.DataArray(rng.normal(size=(12, 100)) + 0.01 * (years - 1900),
                     coords=[("member_X", np.arange(12)), ("year", years)])
    a = p2.internal_component(r)
    assert np.allclose(a.mean("member_X"), 0)


def test_internal_component_single_run_removes_trend_not_oscillation():
    years = np.arange(1940, 2026)
    cycle = 0.2 * np.sin(2 * np.pi * (years - 1940) / 65)
    r = xr.DataArray(0.004 * (years - 1940) + cycle, coords=[("year", years)])
    a = p2.internal_component(r)
    # a linear fit (86-yr record) leaves most of a 65-yr cycle in place
    assert np.corrcoef(a, cycle)[0, 1] > 0.9
    assert abs(float(a.mean())) < 1e-10


def test_random_onset_band_brackets_zero_for_noise():
    rng = np.random.default_rng(2)
    series = [rng.normal(size=200) for _ in range(5)]
    lo, hi = p2.random_onset_band(series, 30, 5, n_draws=300)
    assert np.all(lo < 0) and np.all(hi > 0)


def test_largest_events_takes_the_top_and_declusters():
    a = np.array([0.1, 0.9, 0.8, 0.2, 0.1, 0.1, 0.7])
    b = np.array([0.85, np.nan, 0.3])
    # 0.8 at a[2] sits next to a[1]=0.9 and is skipped
    out = p2.largest_events([a, b], 3, min_separation=3)
    assert [list(o) for o in out] == [[1, 6], [0]]
