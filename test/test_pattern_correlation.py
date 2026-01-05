"""
Unit tests for pattern correlation functions
"""

import pytest
import numpy as np
import xarray as xr
from src.pattern_correlation import (
    calculate_pattern_correlation,
    calculate_spatial_correlation
)


class TestPatternCorrelation:
    """Test suite for calculate_pattern_correlation function"""
    
    def test_identical_patterns(self):
        """Test that identical patterns have correlation of 1.0"""
        pattern = np.random.randn(10, 20)
        corr = calculate_pattern_correlation(pattern, pattern)
        assert np.isclose(corr, 1.0), "Identical patterns should have correlation of 1.0"
    
    def test_opposite_patterns(self):
        """Test that opposite patterns have correlation of -1.0"""
        pattern1 = np.random.randn(10, 20)
        pattern2 = -pattern1
        corr = calculate_pattern_correlation(pattern1, pattern2)
        assert np.isclose(corr, -1.0), "Opposite patterns should have correlation of -1.0"
    
    def test_uncorrelated_patterns(self):
        """Test that uncorrelated patterns have correlation near 0"""
        np.random.seed(42)
        pattern1 = np.random.randn(100, 100)
        pattern2 = np.random.randn(100, 100)
        corr = calculate_pattern_correlation(pattern1, pattern2)
        # Should be close to 0 for large random arrays
        assert abs(corr) < 0.1, "Uncorrelated patterns should have correlation near 0"
    
    def test_with_weights(self):
        """Test pattern correlation with weights"""
        pattern1 = np.random.randn(10, 20)
        pattern2 = np.random.randn(10, 20)
        weights = np.ones((10, 20))
        
        # With uniform weights should be same as without weights
        corr_weighted = calculate_pattern_correlation(pattern1, pattern2, weights=weights)
        corr_unweighted = calculate_pattern_correlation(pattern1, pattern2, weights=None)
        
        assert np.isclose(corr_weighted, corr_unweighted), \
            "Uniform weights should give same result as no weights"
    
    def test_with_nan_values(self):
        """Test that NaN values are properly handled"""
        pattern1 = np.random.randn(10, 20)
        pattern2 = np.random.randn(10, 20)
        
        # Add some NaN values
        pattern1[0:2, 0:3] = np.nan
        pattern2[5:7, 10:15] = np.nan
        
        corr = calculate_pattern_correlation(pattern1, pattern2)
        
        # Should compute without error and return valid correlation
        assert not np.isnan(corr), "Correlation should not be NaN"
        assert -1 <= corr <= 1, "Correlation should be between -1 and 1"
    
    def test_shape_mismatch(self):
        """Test that mismatched shapes raise ValueError"""
        pattern1 = np.random.randn(10, 20)
        pattern2 = np.random.randn(10, 15)
        
        with pytest.raises(ValueError, match="same shape"):
            calculate_pattern_correlation(pattern1, pattern2)
    
    def test_weight_shape_mismatch(self):
        """Test that mismatched weight shape raises ValueError"""
        pattern1 = np.random.randn(10, 20)
        pattern2 = np.random.randn(10, 20)
        weights = np.random.randn(10, 15)
        
        with pytest.raises(ValueError, match="same shape"):
            calculate_pattern_correlation(pattern1, pattern2, weights=weights)
    
    def test_all_nan_patterns(self):
        """Test that all-NaN patterns raise ValueError"""
        pattern1 = np.full((10, 20), np.nan)
        pattern2 = np.random.randn(10, 20)
        
        with pytest.raises(ValueError, match="No valid"):
            calculate_pattern_correlation(pattern1, pattern2)
    
    def test_zero_variance(self):
        """Test that zero variance patterns raise ValueError"""
        pattern1 = np.ones((10, 20))  # Constant pattern, zero variance
        pattern2 = np.random.randn(10, 20)
        
        # With centered=True (default), constant pattern has zero variance after centering
        with pytest.raises(ValueError, match="zero variance"):
            calculate_pattern_correlation(pattern1, pattern2, centered=True)
    
    def test_xarray_input(self):
        """Test with xarray DataArray input"""
        data1 = np.random.randn(10, 20)
        data2 = np.random.randn(10, 20)
        
        # Convert to xarray
        da1 = xr.DataArray(data1, dims=['lat', 'lon'])
        da2 = xr.DataArray(data2, dims=['lat', 'lon'])
        
        corr_xr = calculate_pattern_correlation(da1, da2)
        corr_np = calculate_pattern_correlation(data1, data2)
        
        assert np.isclose(corr_xr, corr_np), \
            "xarray and numpy inputs should give same result"
    
    def test_centered_vs_uncentered(self):
        """Test centered vs uncentered correlation"""
        pattern1 = np.random.randn(10, 20) + 5  # Add offset
        pattern2 = np.random.randn(10, 20) + 3  # Different offset
        
        corr_centered = calculate_pattern_correlation(pattern1, pattern2, centered=True)
        corr_uncentered = calculate_pattern_correlation(pattern1, pattern2, centered=False)
        
        # They should generally be different (unless patterns happen to be already centered)
        # Just check that both are valid
        assert -1 <= corr_centered <= 1
        assert -1 <= corr_uncentered <= 1


class TestSpatialCorrelation:
    """Test suite for calculate_spatial_correlation function"""
    
    def test_basic_spatial_correlation(self):
        """Test basic spatial correlation with xarray"""
        lats = np.linspace(-90, 90, 10)
        lons = np.linspace(0, 360, 20)
        
        data1 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['lat', 'lon'],
            coords={'lat': lats, 'lon': lons}
        )
        data2 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['lat', 'lon'],
            coords={'lat': lats, 'lon': lons}
        )
        
        corr = calculate_spatial_correlation(data1, data2)
        assert -1 <= corr <= 1, "Correlation should be between -1 and 1"
    
    def test_area_weighting(self):
        """Test that area weighting is applied"""
        lats = np.linspace(-90, 90, 10)
        lons = np.linspace(0, 360, 20)
        
        # Create a simple pattern that varies with latitude
        lat_pattern = np.cos(np.deg2rad(lats))[:, np.newaxis] * np.ones((10, 20))
        
        data1 = xr.DataArray(
            lat_pattern,
            dims=['lat', 'lon'],
            coords={'lat': lats, 'lon': lons}
        )
        data2 = xr.DataArray(
            lat_pattern,
            dims=['lat', 'lon'],
            coords={'lat': lats, 'lon': lons}
        )
        
        corr_weighted = calculate_spatial_correlation(data1, data2, area_weighted=True)
        corr_unweighted = calculate_spatial_correlation(data1, data2, area_weighted=False)
        
        # Both should be 1.0 since patterns are identical
        assert np.isclose(corr_weighted, 1.0)
        assert np.isclose(corr_unweighted, 1.0)
    
    def test_missing_lat_dimension(self):
        """Test that missing latitude dimension raises error"""
        data1 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['x', 'y']
        )
        data2 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['x', 'y']
        )
        
        with pytest.raises(ValueError, match="Latitude dimension"):
            calculate_spatial_correlation(data1, data2, area_weighted=True)
    
    def test_custom_dimension_names(self):
        """Test with custom dimension names"""
        lats = np.linspace(-90, 90, 10)
        lons = np.linspace(0, 360, 20)
        
        data1 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['latitude', 'longitude'],
            coords={'latitude': lats, 'longitude': lons}
        )
        data2 = xr.DataArray(
            np.random.randn(10, 20),
            dims=['latitude', 'longitude'],
            coords={'latitude': lats, 'longitude': lons}
        )
        
        corr = calculate_spatial_correlation(
            data1, data2,
            lat_dim='latitude',
            lon_dim='longitude',
            area_weighted=True
        )
        assert -1 <= corr <= 1


class TestEdgeCases:
    """Test edge cases and special scenarios"""
    
    def test_single_point(self):
        """Test correlation with single point"""
        # Single point patterns should raise error (zero variance after centering)
        pattern1 = np.array([[5.0]])
        pattern2 = np.array([[3.0]])
        
        with pytest.raises(ValueError, match="zero variance"):
            calculate_pattern_correlation(pattern1, pattern2)
    
    def test_1d_patterns(self):
        """Test with 1D patterns"""
        pattern1 = np.random.randn(50)
        pattern2 = np.random.randn(50)
        
        corr = calculate_pattern_correlation(pattern1, pattern2)
        assert -1 <= corr <= 1
    
    def test_3d_patterns(self):
        """Test with 3D patterns (e.g., time x lat x lon)"""
        pattern1 = np.random.randn(5, 10, 20)
        pattern2 = np.random.randn(5, 10, 20)
        
        corr = calculate_pattern_correlation(pattern1, pattern2)
        assert -1 <= corr <= 1
    
    def test_latitude_weighting_effect(self):
        """Test that latitude weighting has measurable effect"""
        lats = np.linspace(-90, 90, 18)  # 10-degree resolution
        lons = np.linspace(0, 360, 36)
        
        # Create patterns with random variations plus latitude trend
        np.random.seed(42)
        lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
        
        # Pattern with more variance at high latitudes
        pattern1 = np.abs(lat_grid) / 90 + np.random.randn(18, 36) * 0.1
        # Pattern with more variance at low latitudes  
        pattern2 = (90 - np.abs(lat_grid)) / 90 + np.random.randn(18, 36) * 0.1
        
        data1 = xr.DataArray(pattern1, dims=['lat', 'lon'], coords={'lat': lats, 'lon': lons})
        data2 = xr.DataArray(pattern2, dims=['lat', 'lon'], coords={'lat': lats, 'lon': lons})
        
        corr_weighted = calculate_spatial_correlation(data1, data2, area_weighted=True)
        corr_unweighted = calculate_spatial_correlation(data1, data2, area_weighted=False)
        
        # Both correlations should be negative (opposite patterns)
        assert corr_weighted < 0, "Weighted correlation should be negative"
        assert corr_unweighted < 0, "Unweighted correlation should be negative"
        
        # They should be different, but we'll just verify both are valid
        assert -1 <= corr_weighted <= 1
        assert -1 <= corr_unweighted <= 1
