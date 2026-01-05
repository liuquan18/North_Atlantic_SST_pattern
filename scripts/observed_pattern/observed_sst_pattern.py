# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# %%
# Define paths
data_dir = "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/"
files = {
    "JJA": "anom_2023_jja_natl.nc",
    "June": "anom_2023_jun_natl.nc",
    "July": "anom_2023_jul_natl.nc",
    "August": "anom_2023_aug_natl.nc",
}

# %%
# Setup plot
fig, axes = plt.subplots(
    2, 2, figsize=(12, 10), subplot_kw={"projection": ccrs.PlateCarree()}
)
axes = axes.flatten()

# Plot each period
for i, (period, filename) in enumerate(files.items()):
    ax = axes[i]
    filepath = data_dir + filename

    try:
        ds = xr.open_dataset(filepath)
        # Find the variable name (assuming it's the first data variable)
        var_name = list(ds.data_vars)[-1]
        data = ds[var_name]

        # Handle time dimension if it exists (squeeze it out)
        if "time" in data.dims:
            data = data.squeeze()

        # Plot
        im = ax.contourf(
            data.lon,
            data.lat,
            data,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            levels=np.linspace(-2, 2, 21),
            extend="both",
        )

        ax.add_feature(cfeature.COASTLINE)
        ax.gridlines(
            draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
        )
        ax.set_title(f"SST Anomaly 2023 - {period}")

        # Set extent to North Atlantic (approximate based on user's box)
        # User used sellonlatbox,280,360,0,70. 280 is -80.
        ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

    except FileNotFoundError:
        print(f"File not found: {filepath}")
        ax.text(
            0.5, 0.5, "File not found", ha="center", va="center", transform=ax.transAxes
        )

# Add colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax, label="SST Anomaly (K)")

plt.suptitle(
    "North Atlantic SST Anomalies (2023) vs 1991-2020 Climatology", fontsize=16
)
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/doc/observed_sst_anomaly_2023.png",
    bbox_inches="tight",
)
plt.show()

# %%
