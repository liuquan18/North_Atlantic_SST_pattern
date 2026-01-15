# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# %%
# Define paths
data_dir = "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/"
monthly_file = "anom_2023_06_2023_08.nc"
jja_file = "anom_2023_jja.nc"

# %%
# Setup plot
fig, axes = plt.subplots(
    2, 2, figsize=(12, 10), subplot_kw={"projection": ccrs.PlateCarree()}
)
axes = axes.flatten()

# Read monthly data (Jun, Jul, Aug)
try:
    ds_monthly = xr.open_dataset(data_dir + monthly_file)
    var_name = list(ds_monthly.data_vars)[-1]
    data_monthly = ds_monthly[var_name]

    # Plot each month (Jun, Jul, Aug) with spatial mean removed
    month_names = ["June", "July", "August"]
    for i in range(3):
        ax = axes[i]
        data = data_monthly.isel(time=i)

        # Remove spatial mean
        spatial_mean = data.mean(dim=["lat", "lon"])
        spatial_anomaly = data - spatial_mean

        im = spatial_anomaly.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            vmin=-2,
            vmax=2,
            add_colorbar=False,
            extend="both",
        )

        ax.add_feature(cfeature.COASTLINE)
        ax.gridlines(
            draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
        )
        ax.set_title(
            f"SST Spatial Anomaly 2023 - {month_names[i]}\n(Spatial mean removed: {spatial_mean.values:.3f} K)"
        )
        ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

except FileNotFoundError:
    print(f"File not found: {data_dir + monthly_file}")
    for i in range(3):
        axes[i].text(
            0.5,
            0.5,
            "File not found",
            ha="center",
            va="center",
            transform=axes[i].transAxes,
        )

# Plot JJA seasonal mean with spatial mean removed
try:
    ax = axes[3]
    ds_jja = xr.open_dataset(data_dir + jja_file)
    var_name = list(ds_jja.data_vars)[-1]
    data_jja = ds_jja[var_name]

    # Handle time dimension if it exists
    if "time" in data_jja.dims:
        data_jja = data_jja.squeeze()

    # Remove spatial mean
    spatial_mean_jja = data_jja.mean(dim=["lat", "lon"])
    spatial_anomaly_jja = data_jja - spatial_mean_jja

    im = spatial_anomaly_jja.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=-2,
        vmax=2,
        add_colorbar=False,
        extend="both",
    )

    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(
        draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
    )
    ax.set_title(
        f"SST Spatial Anomaly 2023 - JJA (Seasonal Mean)\n(Spatial mean removed: {spatial_mean_jja.values:.3f} K)"
    )
    ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

except FileNotFoundError:
    print(f"File not found: {data_dir + jja_file}")
    axes[3].text(
        0.5,
        0.5,
        "File not found",
        ha="center",
        va="center",
        transform=axes[3].transAxes,
    )

# Add colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax, label="SST Spatial Anomaly (K)")

plt.suptitle(
    "North Atlantic SST Spatial Anomalies (2023) - Spatial Mean Removed", fontsize=16
)
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/doc/observed_sst_spatial_anomaly_2023.png",
    bbox_inches="tight",
)
plt.show()

# %%
