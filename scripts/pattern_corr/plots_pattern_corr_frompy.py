# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import glob
import cartopy.crs as ccrs

# %%
base_dir = "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr/"

# %%
corrs = xr.open_mfdataset(base_dir + "*.nc", combine="nested", concat_dim="ens")
# %%
corrs = corrs.__xarray_dataarray_variable__
# %%
ds_year = corrs.mean(dim="ens").groupby("time.year").mean(dim="time")
# %%
fig, ax = plt.subplots(figsize=(8, 6))
ds_year.plot(ax=ax)
ax.set_title("Pattern Correlation average over all ensemble members")
# %%
# Prepare data for combined plot
# Add decade grouping
decades = (corrs.time.dt.year // 10) * 10

# Count records above 0.4 for each decade
ds_above_threshold = corrs > 0.4
ds_decade_count = ds_above_threshold.groupby(decades).sum(dim=["time", "ens"])

# Prepare data for violin plot - flatten time and ens dimensions for each decade
data_by_decade = []
decade_labels = []

for decade in np.unique(decades.values):
    decade_mask = decades == decade
    decade_data = corrs.isel(time=decade_mask).values.flatten()
    data_by_decade.append(decade_data)
    decade_labels.append(str(decade))

# %%
# Create figure with two rows
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7))

# First row: Violin plot
violin_parts = ax1.violinplot(
    data_by_decade,
    positions=range(len(decade_labels)),
    showmeans=True,
    showmedians=True,
)
ax1.set_xticks(range(len(decade_labels)))
ax1.set_xticklabels(decade_labels, rotation=45)
ax1.set_xlabel("Decade", fontsize=12)
ax1.set_ylabel("Pattern Correlation", fontsize=12)
ax1.set_title(
    "Distribution of Pattern Correlation",
    fontsize=13,
)
ax1.grid(True, alpha=0.3, axis="y")

# Second row: Count above 0.4
ds_decade_count.plot(ax=ax2, marker="o", linewidth=2, markersize=8)
ax2.set_xlabel("Decade", fontsize=12)
ax2.set_ylabel("Count of Records Above 0.4", fontsize=12)
ax2.set_title(
    "Count of Pattern Correlation Records Above 0.4 by Decade",
    fontsize=13,
)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
# %%
files_first = glob.glob(
    "/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r*i1p1f1/*_199001-200912.nc"
)
# %%
files_first = sorted(files_first, key=lambda x: int(x.split("/r")[1].split("i")[0]))
# %%
anom_first = xr.open_mfdataset(files_first, combine="nested", concat_dim="ens")
# %%
anom_first["ens"] = np.arange(1, 51)
# %%
# Find events in the first decade where pattern correlation > 0.4
first_decade = np.min(decades.values)
first_decade_mask = decades == first_decade

# Get pattern correlation values for first decade
corr_first = corrs.sel(time=first_decade_mask)

# Find indices where pattern correlation > 0.4
high_corr_mask = corr_first > 0.4

# Get the time and ensemble indices where condition is met
high_corr_times = []
high_corr_ens = []

for ens_idx in range(corr_first.shape[0]):  # iterate over ensemble dimension
    for time_idx in range(corr_first.shape[1]):  # iterate over time dimension
        if high_corr_mask.values[ens_idx, time_idx]:
            high_corr_times.append(corr_first.time.values[time_idx])
            high_corr_ens.append(ens_idx + 1)  # ens is 1-indexed

print(
    f"Found {len(high_corr_times)} events in first decade ({first_decade}s) with pattern correlation > 0.4"
)

# %%
# Select corresponding anomaly data for these events
selected_events = []

for time_val, ens_val in zip(high_corr_times, high_corr_ens):
    event_data = anom_first.sel(time=time_val, ens=ens_val, method="nearest")
    selected_events.append(event_data)

# Stack and average all selected events
events_stacked = xr.concat(selected_events, dim="event")
avg_pattern = events_stacked.mean(dim="event")

# %%
# Plot the spatial map of the average pattern
fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": ccrs.PlateCarree()})

avg_pattern.ts.plot(
    ax=ax,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    center=0,
    cbar_kwargs={"label": "SST Anomaly (K)", "shrink": 0.8},
)

ax.coastlines()
ax.gridlines(draw_labels=True, alpha=0.3)
ax.set_title(
    f"Average SST Anomaly Pattern for Events with Pattern Correlation > 0.4\nFirst Decade ({first_decade}s), N={len(high_corr_times)} events",
    fontsize=14,
)

plt.tight_layout()
# %%
