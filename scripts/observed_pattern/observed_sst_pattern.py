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

    # Plot each month (Jun, Jul, Aug)
    month_names = ["June", "July", "August"]
    for i in range(3):
        ax = axes[i]
        data = data_monthly.isel(time=i)

        im = data.plot(
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
        ax.set_title(f"SST Anomaly 2023 - {month_names[i]}")
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

# Plot JJA seasonal mean
try:
    ax = axes[3]
    ds_jja = xr.open_dataset(data_dir + jja_file)
    var_name = list(ds_jja.data_vars)[-1]
    data_jja = ds_jja[var_name]

    # Handle time dimension if it exists
    if "time" in data_jja.dims:
        data_jja = data_jja.squeeze()

    im = data_jja.plot(
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
    ax.set_title(f"SST Anomaly 2023 - JJA (Seasonal Mean)")
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
# %%
# function to do spatial normalize
def spatial_normalize(data):
    spatial_mean = data.mean(dim=["lat", "lon"])
    spatial_anomaly = data - spatial_mean
    spatial_std = spatial_anomaly.std(dim=["lat", "lon"])
    normalized_data = spatial_anomaly / spatial_std
    return normalized_data


# %%
# Plot spatially normalized maps
fig2, axes2 = plt.subplots(
    2, 2, figsize=(12, 10), subplot_kw={"projection": ccrs.PlateCarree()}
)
axes2 = axes2.flatten()

# Read and normalize monthly data (Jun, Jul, Aug)
try:
    ds_monthly = xr.open_dataset(data_dir + monthly_file)
    var_name = list(ds_monthly.data_vars)[-1]
    data_monthly = ds_monthly[var_name]

    # Plot each month (Jun, Jul, Aug) normalized
    month_names = ["June", "July", "August"]
    for i in range(3):
        ax = axes2[i]
        data = data_monthly.isel(time=i)
        normalized_data = spatial_normalize(data)

        im2 = normalized_data.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            vmin=-3,
            vmax=3,
            add_colorbar=False,
            extend="both",
        )

        ax.add_feature(cfeature.COASTLINE)
        ax.gridlines(
            draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
        )
        ax.set_title(f"Normalized SST Anomaly 2023 - {month_names[i]}")
        ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

except FileNotFoundError:
    print(f"File not found: {data_dir + monthly_file}")
    for i in range(3):
        axes2[i].text(
            0.5,
            0.5,
            "File not found",
            ha="center",
            va="center",
            transform=axes2[i].transAxes,
        )

# Plot JJA seasonal mean normalized
try:
    ax = axes2[3]
    ds_jja = xr.open_dataset(data_dir + jja_file)
    var_name = list(ds_jja.data_vars)[-1]
    data_jja = ds_jja[var_name]

    # Handle time dimension if it exists
    if "time" in data_jja.dims:
        data_jja = data_jja.squeeze()

    normalized_jja = spatial_normalize(data_jja)

    im2 = normalized_jja.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=-3,
        vmax=3,
        add_colorbar=False,
        extend="both",
    )

    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(
        draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
    )
    ax.set_title(f"Normalized SST Anomaly 2023 - JJA (Seasonal Mean)")
    ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

except FileNotFoundError:
    print(f"File not found: {data_dir + jja_file}")
    axes2[3].text(
        0.5,
        0.5,
        "File not found",
        ha="center",
        va="center",
        transform=axes2[3].transAxes,
    )

# Add colorbar
cbar_ax2 = fig2.add_axes([0.92, 0.15, 0.02, 0.7])
fig2.colorbar(im2, cax=cbar_ax2, label="Normalized SST Anomaly (σ)")

plt.suptitle("Spatially Normalized North Atlantic SST Anomalies (2023)", fontsize=16)
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/doc/observed_sst_normalized_2023.png",
    bbox_inches="tight",
)
plt.show()

# %%
