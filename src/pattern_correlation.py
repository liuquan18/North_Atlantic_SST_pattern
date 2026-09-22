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


def calculate_spatial_correlation_keepareamean(
    field1,
    field2,
    *,
    lat_dim="lat",
    time_dim="time",
    area_weighted=True,
):
    """
    Calculate spatial pattern correlation without removing spatial mean.

    This computes an uncentered (cosine-like) pattern correlation:
    sum(w * x * y) / sqrt(sum(w * x^2) * sum(w * y^2))

    Parameters
    ----------
    field1 : xarray.DataArray
        First spatial field.
    field2 : xarray.DataArray
        Second spatial field.
    lat_dim : str, optional
        Latitude dimension name, by default "lat".
    time_dim : str, optional
        Time dimension name to exclude from spatial flattening, by default "time".
    area_weighted : bool, optional
        If True, apply cosine(latitude) area weighting, by default True.

    Returns
    -------
    xarray.DataArray
        Uncentered spatial pattern correlation. If input has a time dimension,
        returns a time series.
    """
    if not isinstance(field1, xr.DataArray) or not isinstance(field2, xr.DataArray):
        raise TypeError("field1 and field2 must be xarray.DataArray objects")

    spatial_dims = [dim for dim in field2.dims if dim != time_dim]
    if len(spatial_dims) == 0:
        raise ValueError("No spatial dimensions found to compute pattern correlation")

    # Fields remapped onto "the same" grid by independent cdo pipelines can
    # carry lat/lon coordinate values that agree to ~1e-6 but not bit-for-bit
    # (join="exact" would reject that as a mismatch). Snap field1's spatial
    # coordinates onto field2's once they're confirmed to agree within
    # tolerance, then align on those now-identical labels.
    for dim in spatial_dims:
        if dim in field1.coords and dim in field2.coords:
            if field1.sizes.get(dim) != field2.sizes.get(dim) or not np.allclose(
                field1[dim].values, field2[dim].values, rtol=0, atol=1e-4
            ):
                raise ValueError(
                    f"Coordinate '{dim}' differs beyond tolerance between field1 and field2"
                )
            field1 = field1.assign_coords({dim: field2[dim].values})

    field1, field2 = xr.align(field1, field2, join="exact")

    field1_flat = field1.stack(space=spatial_dims)
    field2_flat = field2.stack(space=spatial_dims)

    valid_mask = ~(field1_flat.isnull() | field2_flat.isnull())
    field1_valid = field1_flat.where(valid_mask)
    field2_valid = field2_flat.where(valid_mask)

    if area_weighted:
        if lat_dim not in spatial_dims:
            raise ValueError(
                f"Latitude dimension '{lat_dim}' not found in spatial dims {spatial_dims}"
            )

        lat_weights = np.cos(np.deg2rad(field2[lat_dim]))
        lat_weights.name = "weights"

        spatial_template = (
            field2.isel({time_dim: 0}) if time_dim in field2.dims else field2
        )
        weights = lat_weights * xr.ones_like(spatial_template)
    else:
        spatial_template = (
            field2.isel({time_dim: 0}) if time_dim in field2.dims else field2
        )
        weights = xr.ones_like(spatial_template)

    weights_flat = weights.stack(space=spatial_dims).where(valid_mask)

    numerator = (weights_flat * field1_valid * field2_valid).sum(
        dim="space", skipna=True
    )
    denominator = np.sqrt(
        (weights_flat * field1_valid**2).sum(dim="space", skipna=True)
        * (weights_flat * field2_valid**2).sum(dim="space", skipna=True)
    )

    correlation = numerator / denominator
    return correlation.where(denominator > 0)


def spatial_corr_model_era(model_data, era_data, rm_spatial_mean=True):
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
    # xr.corr broadcasts over any dims not in `dim=`, so the whole time
    # series is computed in one call -- no need to loop and correlate
    # year by year (verified to give identical results, ~50-100x fewer
    # dask task-graph evaluations for a typical ensemble member).
    if rm_spatial_mean:
        return calculate_spatial_correlation(model_data, era_data)
    else:
        return calculate_spatial_correlation_keepareamean(model_data, era_data)


# %%
