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
