"""
Pattern Correlation Functions

This module provides functions for calculating pattern correlation between
SST (Sea Surface Temperature) patterns.
"""

import numpy as np
import xarray as xr
from typing import Union, Tuple


def calculate_pattern_correlation(
    pattern1: Union[np.ndarray, xr.DataArray],
    pattern2: Union[np.ndarray, xr.DataArray],
    weights: Union[np.ndarray, xr.DataArray, None] = None,
    centered: bool = True
) -> float:
    """
    Calculate the pattern correlation between two spatial patterns.
    
    The pattern correlation is a measure of similarity between two spatial patterns.
    It is calculated as the Pearson correlation coefficient between the two patterns,
    optionally weighted by grid cell area or other weights.
    
    Parameters
    ----------
    pattern1 : np.ndarray or xr.DataArray
        First spatial pattern (e.g., SST anomalies)
    pattern2 : np.ndarray or xr.DataArray
        Second spatial pattern (e.g., SST anomalies from model)
    weights : np.ndarray, xr.DataArray, or None, optional
        Weights for each grid point (e.g., cosine of latitude for area weighting).
        If None, all points are equally weighted. Default is None.
    centered : bool, optional
        If True, center the patterns by removing their weighted means before 
        calculating correlation. Default is True.
    
    Returns
    -------
    float
        Pattern correlation coefficient (between -1 and 1)
    
    Examples
    --------
    >>> import numpy as np
    >>> pattern1 = np.random.randn(10, 20)
    >>> pattern2 = np.random.randn(10, 20)
    >>> corr = calculate_pattern_correlation(pattern1, pattern2)
    >>> print(f"Pattern correlation: {corr:.3f}")
    
    >>> # With latitude weights
    >>> lats = np.linspace(-90, 90, 10)
    >>> weights = np.cos(np.deg2rad(lats))[:, np.newaxis]
    >>> corr = calculate_pattern_correlation(pattern1, pattern2, weights=weights)
    
    Notes
    -----
    - NaN values are automatically masked out in the calculation
    - Patterns should have the same shape
    - For SST patterns, it's common to use cosine of latitude as weights
      to account for grid cell area differences
    """
    # Convert xarray DataArrays to numpy arrays if needed
    if isinstance(pattern1, xr.DataArray):
        pattern1 = pattern1.values
    if isinstance(pattern2, xr.DataArray):
        pattern2 = pattern2.values
    if isinstance(weights, xr.DataArray):
        weights = weights.values
    
    # Ensure patterns are numpy arrays
    pattern1 = np.asarray(pattern1)
    pattern2 = np.asarray(pattern2)
    
    # Check that patterns have the same shape
    if pattern1.shape != pattern2.shape:
        raise ValueError(
            f"Patterns must have the same shape. "
            f"Got pattern1: {pattern1.shape}, pattern2: {pattern2.shape}"
        )
    
    # Flatten patterns for easier computation
    pattern1_flat = pattern1.flatten()
    pattern2_flat = pattern2.flatten()
    
    # Create mask for valid (non-NaN) values
    valid_mask = ~(np.isnan(pattern1_flat) | np.isnan(pattern2_flat))
    
    if not np.any(valid_mask):
        raise ValueError("No valid (non-NaN) data points found in both patterns")
    
    # Apply mask
    pattern1_valid = pattern1_flat[valid_mask]
    pattern2_valid = pattern2_flat[valid_mask]
    
    # Handle weights
    if weights is None:
        weights_valid = np.ones_like(pattern1_valid)
    else:
        weights = np.asarray(weights)
        if weights.shape != pattern1.shape:
            raise ValueError(
                f"Weights must have the same shape as patterns. "
                f"Got weights: {weights.shape}, patterns: {pattern1.shape}"
            )
        weights_flat = weights.flatten()
        weights_valid = weights_flat[valid_mask]
    
    # Normalize weights
    weights_valid = weights_valid / np.sum(weights_valid)
    
    # Center the patterns if requested
    if centered:
        pattern1_mean = np.sum(pattern1_valid * weights_valid)
        pattern2_mean = np.sum(pattern2_valid * weights_valid)
        pattern1_centered = pattern1_valid - pattern1_mean
        pattern2_centered = pattern2_valid - pattern2_mean
    else:
        pattern1_centered = pattern1_valid
        pattern2_centered = pattern2_valid
    
    # Calculate weighted covariance
    covariance = np.sum(weights_valid * pattern1_centered * pattern2_centered)
    
    # Calculate weighted standard deviations
    std1 = np.sqrt(np.sum(weights_valid * pattern1_centered**2))
    std2 = np.sqrt(np.sum(weights_valid * pattern2_centered**2))
    
    # Calculate correlation
    if std1 == 0 or std2 == 0 or np.abs(std1) < 1e-10 or np.abs(std2) < 1e-10:
        raise ValueError("One or both patterns have zero variance")
    
    correlation = covariance / (std1 * std2)
    
    return correlation


def calculate_spatial_correlation(
    data1: xr.DataArray,
    data2: xr.DataArray,
    lat_dim: str = 'lat',
    lon_dim: str = 'lon',
    area_weighted: bool = True
) -> float:
    """
    Calculate spatial pattern correlation between two xarray DataArrays.
    
    This is a convenience wrapper around calculate_pattern_correlation
    specifically for xarray DataArrays with latitude/longitude dimensions.
    
    Parameters
    ----------
    data1 : xr.DataArray
        First spatial field
    data2 : xr.DataArray
        Second spatial field
    lat_dim : str, optional
        Name of the latitude dimension. Default is 'lat'.
    lon_dim : str, optional
        Name of the longitude dimension. Default is 'lon'.
    area_weighted : bool, optional
        If True, weight by cosine of latitude. Default is True.
    
    Returns
    -------
    float
        Pattern correlation coefficient
    
    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> data1 = xr.DataArray(
    ...     np.random.randn(10, 20),
    ...     dims=['lat', 'lon'],
    ...     coords={'lat': np.linspace(-90, 90, 10), 'lon': np.linspace(0, 360, 20)}
    ... )
    >>> data2 = xr.DataArray(
    ...     np.random.randn(10, 20),
    ...     dims=['lat', 'lon'],
    ...     coords={'lat': np.linspace(-90, 90, 10), 'lon': np.linspace(0, 360, 20)}
    ... )
    >>> corr = calculate_spatial_correlation(data1, data2)
    """
    # Calculate area weights if requested
    if area_weighted:
        if lat_dim not in data1.dims:
            raise ValueError(f"Latitude dimension '{lat_dim}' not found in data1")
        
        # Get latitude values
        lats = data1[lat_dim].values
        
        # Calculate cosine weights
        weights = np.cos(np.deg2rad(lats))
        
        # Broadcast weights to match data shape
        # Create a weight array that matches the data dimensions
        weight_array = xr.DataArray(
            weights,
            dims=[lat_dim],
            coords={lat_dim: data1[lat_dim]}
        )
        # Broadcast to match full data shape
        weight_array = weight_array * xr.ones_like(data1)
    else:
        weight_array = None
    
    # Calculate correlation
    correlation = calculate_pattern_correlation(data1, data2, weights=weight_array)
    
    return correlation
