# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# %%
data_2023_14 = xr.open_dataset(
    "/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r15i1p1f1/ts_Amon_MPI-ESM1-2-LR_ssp245_r15i1p1f1_gn_201501-203412.nc"
)
# %%
data_2023_14 = data_2023_14.sel(time="2030")
# %%
data_2023_14_mean = data_2023_14.mean(dim="time").ts
# %%
# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Plot 1: Original data
data_2023_14_mean.plot(
    ax=ax1,
    levels=np.arange(-1.8, 2.0, 0.1),
    cmap="RdBu_r",
    cbar_kwargs={"label": "SST Anomaly (K)"},
)
ax1.set_title("Original Data")

# Plot 2: Spatial mean removed
spatial_mean = data_2023_14_mean.mean(dim=["lat", "lon"])
spatial_anomaly = data_2023_14_mean - spatial_mean

spatial_anomaly.plot(
    ax=ax2,
    levels=np.arange(-1.8, 2.0, 0.1),
    cmap="RdBu_r",
    cbar_kwargs={"label": "SST Spatial Anomaly (K)"},
)
ax2.set_title(f"Spatial Mean Removed (mean={spatial_mean.values:.3f} K)")

plt.tight_layout()
# %%
corr = xr.open_dataset(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr_JJA/mpi_era5_patterncorr_ens15.nc"
)
# %%
corr.sel(time=2030, method="nearest").__xarray_dataarray_variable__.values
# %%
