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
    ax.set_title("SST Anomaly 2023 - JJA (Seasonal Mean)")
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

    return spatial_anomaly


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

    normalized_jja = data_jja

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
    ax.set_title("Normalized SST Anomaly 2023 - JJA (Seasonal Mean)")
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

plt.suptitle("remove spatial mean North Atlantic SST Anomalies (2023)", fontsize=16)
# plt.savefig(
#     "/work/mh0033/m300883/North_Atlantic_SST_pattern/doc/observed_sst_normalized_2023.png",
#     bbox_inches="tight",
# )
plt.show()

# %%
########
# JJA-mean only, original, remove the patial mean, remove the global mean
########
JJA_original = xr.open_dataset(data_dir + "anom_2023_jja.nc").var34.squeeze()
JJA_rm_spam = spatial_normalize(JJA_original)
# %%
JJA_rm_glbm = xr.open_dataset(data_dir + "anom_2023_jja_rm_glbm.nc").var34.squeeze()

# %%
fig3, axes3 = plt.subplots(
    1, 3, figsize=(16, 5), subplot_kw={"projection": ccrs.PlateCarree()}
)

fields = [JJA_original, JJA_rm_spam, JJA_rm_glbm]
titles = ["JJA Original", "JJA Remove Spatial Mean", "JJA Remove Global Mean"]

levels = np.arange(-2, 2.1, 0.2)

for ax, field, title in zip(axes3, fields, titles):
    im3 = field.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        levels=levels,
        add_colorbar=False,
        extend="both",
    )
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(
        draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0
    )
    ax.set_title(title)
    ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())

fig3.subplots_adjust(right=0.90)
cbar_ax3 = fig3.add_axes([0.92, 0.18, 0.015, 0.64])
fig3.colorbar(im3, cax=cbar_ax3, label="SST Anomaly (K)")
plt.suptitle("JJA Pattern Comparison", fontsize=16)
plt.show()

# %%
ens14_2023_jja = xr.open_dataset(
    "/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r14i1p1f1/ts_Amon_MPI-ESM1-2-LR_ssp245_r14i1p1f1_gn_201501-203412.nc"
).ts
# %%
ens14_2023_jja = ens14_2023_jja.groupby("time.year").mean("time").sel(year=2023)
# %%
ens14_2023_jja_rm_spam = spatial_normalize(ens14_2023_jja)

# %%
ens14_2023_jja_rm_glbm = xr.open_dataset(
    "/scratch/m/m300883/nalt/MPI_GE_CMIP6/removed_glbm/r14i1p1f1/ts_Amon_MPI-ESM1-2-LR_ssp245_r14i1p1f1_gn_201501-203412.nc"
).ts
# %%
ens14_2023_jja_rm_glbm = (
    ens14_2023_jja_rm_glbm.groupby("time.year").mean("time").sel(year=2023)
)

# %%
fig4, axes4 = plt.subplots(
    2, 3, figsize=(18, 9), subplot_kw={"projection": ccrs.PlateCarree()}
)

compare_fields = [
    [JJA_original, JJA_rm_spam, JJA_rm_glbm],
    [ens14_2023_jja, ens14_2023_jja_rm_spam, ens14_2023_jja_rm_glbm],
]
col_titles = ["Original", "Remove Spatial Mean", "Remove Global Mean"]
row_labels = ["ERA5 (2023 JJA)", "Ensemble 14 (2023 JJA)"]

levels_compare = np.arange(-2, 2.1, 0.2)

for r in range(2):
    for c in range(3):
        ax = axes4[r, c]
        im4 = compare_fields[r][c].plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            levels=levels_compare,
            add_colorbar=False,
            extend="both",
        )
        ax.add_feature(cfeature.COASTLINE)
        ax.gridlines(
            draw_labels=False, dms=True, x_inline=False, y_inline=False, linewidth=0
        )
        ax.set_extent([-80, 0, 0, 70], crs=ccrs.PlateCarree())
        if r == 0:
            ax.set_title(col_titles[c], fontsize=12)
        if c == 0:
            ax.text(
                -0.15,
                0.5,
                row_labels[r],
                transform=ax.transAxes,
                rotation=90,
                va="center",
                ha="center",
                fontsize=12,
            )

fig4.subplots_adjust(right=0.90, wspace=0.05, hspace=0.10)
cbar_ax4 = fig4.add_axes([0.92, 0.20, 0.015, 0.60])
fig4.colorbar(im4, cax=cbar_ax4, label="SST Anomaly (K)")
plt.suptitle("ERA5 vs Ensemble 14 (2023 JJA)", fontsize=16)
plt.show()

# %%
