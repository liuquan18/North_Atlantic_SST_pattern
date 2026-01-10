# %%
from pdb import main
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# %%
# Data directory
data_dir = Path(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/MPI_GE_CMIP6/pattern_corr"
)

# Output directory
output_dir = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern/figures")
output_dir.mkdir(parents=True, exist_ok=True)


# %%
def read_ensemble_data(ens_num):
    """
    Read all pattern correlation files for a single ensemble member.

    Parameters:
    -----------
    ens_num : int
        Ensemble member number (1-50)

    Returns:
    --------
    xarray.DataArray
        Time series of pattern correlation values
    """
    ens_dir = data_dir / f"r{ens_num}i1p1f1"

    # Get all NetCDF files sorted by filename
    nc_files = sorted(ens_dir.glob("*.nc"))

    # Open all files and concatenate along time dimension
    ds_list = []
    for nc_file in nc_files:
        ds = xr.open_dataset(nc_file)
        # Pattern correlation should be stored as a variable
        # Assuming the variable name is 'ts' or similar
        ds_list.append(ds)

    # Concatenate along time dimension
    ds_combined = xr.concat(ds_list, dim="time")

    # Extract the correlation values (should be a single value per time step)
    # The variable name might vary, so we'll get the first data variable
    var_name = list(ds_combined.data_vars)[-1]
    corr_values = ds_combined[var_name]

    # Squeeze out spatial dimensions (should be 1x1 after fldcor)
    corr_values = corr_values.squeeze()
    corr_values["ens"] = ens_num

    return corr_values


# %%

# Storage for all data
all_times = []
all_ens_nums = []
all_corr_values = []

# Loop over all ensemble members
for ens_num in range(1, 51):
    print(f"Processing ensemble {ens_num}/50...")

    corr_data = read_ensemble_data(ens_num)
    all_corr_values.append(corr_data)

corr_data = xr.concat(all_corr_values, dim="ens")
# %%

corr_data.mean(dim="ens").plot(figsize=(12, 6))
# %%
obs = xr.open_dataset(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_jja.nc"
).var34
# %%
obs = obs.squeeze()
# %%
obs.plot()
# %%
model = xr.open_dataset(
    "/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r1i1p1f1/ts_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc"
).ts
# %%
model = model.isel(time=0).squeeze()
# %%
model.plot()
# %%
obs.mean()
# %%
model.mean()
# %%
obs_anom = obs - obs.mean()
# %%
model_anom = model - model.mean()
# %%
obs_anom.plot()
# %%
model_anom.plot()


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


# %%
# Calculate spatial correlation between obs and model (raw fields)
print("Calculating spatial correlation between obs and model (raw fields)...")
corr_raw = calculate_spatial_correlation(obs, model)
print(f"Spatial correlation (raw): {corr_raw.values:.4f}")

# %%
# Calculate spatial correlation between obs_anom and model_anom (anomaly fields)
print(
    "\nCalculating spatial correlation between obs_anom and model_anom (anomaly fields)..."
)
corr_anom = calculate_spatial_correlation(obs_anom, model_anom)
print(f"Spatial correlation (anomaly): {corr_anom.values:.4f}")

# %%
