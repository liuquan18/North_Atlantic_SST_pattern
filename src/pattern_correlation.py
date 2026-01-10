"""
Pattern Correlation Functions

This module provides functions for calculating pattern correlation between
SST (Sea Surface Temperature) patterns.
"""

import numpy as np
import xarray as xr
from typing import Union


# %%
def calculate_spatial_correlation(field1, field2):
    """
    Calculate spatial correlation between two fields.

    Parameters:
    -----------
    field1 : xarray.DataArray
        First spatial field (2D or 3D with time)
    field2 : xarray.DataArray
        Second spatial field (2D or 3D with time), must have same spatial dimensions as field1

    Returns:
    --------
    float or xarray.DataArray
        Spatial correlation coefficient. If inputs have time dimension, returns time series of correlations.
    """
    # Flatten spatial dimensions
    # Get the spatial dimension names (exclude 'time' if present)

    weights = np.cos(np.deg2rad(field1.lat))
    weights.name = "weights"
    weights = xr.DataArray(weights, dims=["lat"], coords={"lat": field1.lat})
    # Broadcast weights to match field1's spatial dimensions
    weights = weights * xr.ones_like(field1.isel(time=0) if 'time' in field1.dims else field1)


    spatial_dims = [dim for dim in field1.dims if dim != "time"]

    # Stack spatial dimensions into a single dimension
    field1_flat = field1.stack(space=spatial_dims)
    field2_flat = field2.stack(space=spatial_dims)
    weights = weights.stack(space=spatial_dims)

    # Remove NaN values (ocean mask)
    valid_mask = ~(field1_flat.isnull() | field2_flat.isnull())
    field1_valid = field1_flat.where(valid_mask)
    field2_valid = field2_flat.where(valid_mask)

    # Calculate correlation along spatial dimension
    # Using xarray's correlation method
    correlation = xr.corr(field1_valid, field2_valid, dim="space", weights=weights)

    return correlation
