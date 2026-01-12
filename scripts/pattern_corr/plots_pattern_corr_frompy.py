# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import glob
import cartopy.crs as ccrs

# %%
# Read pattern correlation data with proper ensemble indexing
base_dir = "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr/"
corr_files = sorted(glob.glob(base_dir + "mpi_era5_patterncorr_ens*.nc"))
# Extract ensemble numbers from filenames
corr_ens_nums = [int(f.split("ens")[1].split(".nc")[0]) for f in corr_files]

corrs = xr.open_mfdataset(corr_files, combine="nested", concat_dim="ens")
corrs = corrs.__xarray_dataarray_variable__
# Assign proper ensemble coordinate
corrs["ens"] = corr_ens_nums

# %%
# Read anomaly data with proper ensemble indexing

anom_dir = "/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/"
anom_data = []
for ens_num in corr_ens_nums:
    ens_dir = f"{anom_dir}r{ens_num}i1p1f1/"
    ens_files = sorted(glob.glob(ens_dir + "*.nc"))
    ens_ds = xr.open_mfdataset(ens_files, combine="by_coords")
    ens_ds["ens"] = ens_num
    anom_data.append(ens_ds)
# Combine all ensembles into a single dataset
anom_data = xr.concat(anom_data, dim="ens")

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
ds_above_threshold = corrs >= 0.5
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
# add horizontal line at 0.5
ax1.axhline(0.5, color="red", linestyle="dotted", linewidth=1)

# Second row: Count above 0.4
ds_decade_count.plot(ax=ax2, marker="o", linewidth=2, markersize=8)
ax2.set_xlabel("Decade", fontsize=12)
ax2.set_ylabel("Count of Records Above 0.5", fontsize=12)
ax2.set_title(
    "Count of Pattern Correlation Records Above 0.5 by Decade",
    fontsize=13,
)
ax2.grid(True, alpha=0.3)
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/figures/violin_count_pattern_corr.png",
    dpi=300,
)

plt.tight_layout()
# %%
# Find events where pattern correlation >= 0.5 for all decades between 1990s to 2040s
threshold = 0.5
decades = (corrs.time.dt.year // 10) * 10

# Filter for decades between 1990 and 2040
decade_mask = (decades >= 1990) & (decades <= 2040)
corr_filtered = corrs.where(decade_mask, drop=True)

# Find where correlation >= threshold using index-based method
high_corr_mask = corr_filtered >= threshold

# Get indices using np.where
ens_indices, time_indices = np.where(high_corr_mask.values)

print(
    f"Found {len(ens_indices)} events in decades 1990s-2040s with pattern correlation >= {threshold}"
)

# %%
# Select anomaly data using indices (no for loops)
# Get the actual time and ens values for the found indices
selected_times = corr_filtered.time.values[time_indices]
selected_ens = corr_filtered.ens.values[ens_indices]
selected_corr_values = corr_filtered.values[ens_indices, time_indices]

# Stack selections to create a multi-index selection
# Use xr.concat to gather all events efficiently
selected_events = []
for i in range(len(ens_indices)):
    event = anom_data.sel(time=selected_times[i], ens=selected_ens[i])
    selected_events.append(event)

events_stacked = xr.concat(selected_events, dim="event")

# Compute average pattern
avg_pattern = events_stacked.mean(dim="event")

# Find strongest pattern
max_corr_idx = np.argmax(selected_corr_values)
max_corr_value = selected_corr_values[max_corr_idx]
strongest_pattern = selected_events[max_corr_idx]
strongest_time = selected_times[max_corr_idx]
strongest_ens = selected_ens[max_corr_idx]

print(f"Strongest correlation: {max_corr_value:.4f}")
print(
    f"Strongest event: Ensemble {strongest_ens}, Time {np.datetime_as_string(strongest_time, unit='M')}"
)

# %%
# Plot the spatial maps: average and strongest pattern
fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(16, 6), subplot_kw={"projection": ccrs.PlateCarree()}
)

# First column: Average pattern
avg_pattern.ts.plot(
    ax=ax1,
    transform=ccrs.PlateCarree(),
    levels=np.linspace(-2, 2, 21),
    cmap="RdBu_r",
    center=0,
    cbar_kwargs={"label": "SST Anomaly (K)", "shrink": 0.8},
    extend="both",
)
ax1.coastlines()
ax1.gridlines(draw_labels=True, alpha=0.3)
ax1.set_title(
    f"Average Pattern\n(N={len(ens_indices)} events, 1990s-2040s)",
    fontsize=13,
)

# Second column: Strongest pattern
strongest_pattern.ts.plot(
    ax=ax2,
    transform=ccrs.PlateCarree(),
    levels=np.linspace(-2, 2, 21),
    cmap="RdBu_r",
    center=0,
    cbar_kwargs={"label": "SST Anomaly (K)", "shrink": 0.8},
    extend="both",
)
ax2.coastlines()
ax2.gridlines(draw_labels=True, alpha=0.3)
ax2.set_title(
    f"Strongest Pattern (Corr={max_corr_value:.3f})\n(Ens {strongest_ens}, Time: {np.datetime_as_string(strongest_time, unit='M')})",
    fontsize=13,
)

plt.suptitle(
    f"SST Anomaly Patterns for Pattern Correlation >= {threshold} (1990s-2040s)",
    fontsize=15,
    y=1.02,
)
plt.tight_layout()
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/figures/first_decade_high_corr_patterns.png",
    dpi=300,
)
# %%
# Plot all individual events
n_events = len(selected_events)
ncols = min(4, n_events)  # Maximum 4 columns
nrows = int(np.ceil(n_events / ncols))

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(5 * ncols, 4 * nrows),
    subplot_kw={"projection": ccrs.PlateCarree()},
)

# Flatten axes array for easy iteration
if n_events == 1:
    axes = [axes]
elif nrows == 1:
    axes = axes
else:
    axes = axes.flatten()

# Plot each event
for i in range(n_events):
    ax = axes[i]
    im = selected_events[i].ts.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=np.linspace(-2, 2, 21),
        cmap="RdBu_r",
        center=0,
        add_colorbar=False,
        extend="both",
    )
    ax.coastlines()
    ax.gridlines(draw_labels=False, alpha=0.3)
    ax.set_title(
        f"Event {i+1}: Corr={selected_corr_values[i]:.3f}\nEns {selected_ens[i]}, {np.datetime_as_string(selected_times[i], unit='M')}",
        fontsize=11,
    )

# Hide any unused subplots
for i in range(n_events, len(axes)):
    axes[i].set_visible(False)

# Add a single colorbar for all subplots
cbar = fig.colorbar(
    im,
    ax=axes[10:],
    label="SST Anomaly (K)",
    shrink=0.8,
    pad=0.02,
    orientation="horizontal",
)

plt.suptitle(
    f"All {n_events} Events with Pattern Correlation >= {threshold} (1990s-2040s)",
    fontsize=15,
    y=0.995,
)
plt.tight_layout()
plt.savefig(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/figures/all_high_corr_events.png",
    dpi=300,
    bbox_inches="tight",
)
# %%
