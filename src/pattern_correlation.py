

import numpy as np
import xarray as xr


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

    weights = np.cos(np.deg2rad(field2.lat))
    weights.name = "weights"
    weights = xr.DataArray(weights, dims=["lat"], coords={"lat": field1.lat})
    # Broadcast weights to match field1's spatial dimensions
    weights = weights * xr.ones_like(
        field2.isel(time=0) if "time" in field2.dims else field2
    )

    spatial_dims = [dim for dim in field2.dims if dim != "time"]

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


def spatial_corr_model_era(model_data, era_data):
    """
    Calculate spatial correlation between model data (with time) and ERA5 data (single time).

    Parameters:
    -----------
    model_data : xarray.DataArray
        Model data with dimensions (time, lat, lon) or (time, lat, lon, ...)
    era_data : xarray.DataArray
        ERA5 reference data with dimensions (lat, lon) - single time snapshot

    Returns:
    --------
    xarray.DataArray
        Time series of spatial correlation coefficients between each model time step
        and the ERA5 reference pattern
    """

    # Option 1: Simple loop over time (most reliable)
    correlations = []
    for t in range(len(model_data.time)):
        model_slice = model_data.isel(time=t)
        corr_t = calculate_spatial_correlation(era_data, model_slice)
        correlations.append(corr_t)

    # Combine into DataArray with time coordinate
    result = xr.DataArray(correlations, dims=["time"], coords={"time": model_data.time})

    return result


# %%
