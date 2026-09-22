"""
Updated figures for the tos-based pattern correlation (see
doc/pattern_corr_JJA_tos.md for the write-up this supports).

Mirrors the two-methodology structure of doc/pattern_corr_JJA.pptx:
  1. remove spatial mean before correlating
  2. remove global mean before correlating (extended back to 1850)
but now against Omon/tos (actual SST) instead of Amon/ts + land mask.
"""

import glob
import warnings

import cartopy.crs as ccrs
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

PROJECT_ROOT = "/work/mh0033/m300883/North_Atlantic_SST_pattern"
FIG_DIR = f"{PROJECT_ROOT}/figures"
THRESHOLD = 0.4

VARIANTS = {
    "method1": {
        "label": "Method 1: spatial mean removed",
        "corr_dir": f"{PROJECT_ROOT}/data/pattern_corr_JJA_tos",
        "anom_dir": "/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/anomaly",
        "suffix": "tos",
    },
    "method2": {
        "label": "Method 2: global mean removed",
        "corr_dir": f"{PROJECT_ROOT}/data/pattern_corr_JJA_tos_noglm",
        "anom_dir": "/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/removed_glbm",
        "suffix": "tos_noglm",
    },
}


def load_corr(corr_dir):
    files = sorted(glob.glob(f"{corr_dir}/mpi_era5_patterncorr_ens*.nc"))
    ens_nums = [int(f.split("ens")[1].split(".nc")[0]) for f in files]
    ds = xr.open_mfdataset(files, combine="nested", concat_dim="ens")
    corrs = ds[list(ds.data_vars)[-1]]
    corrs["ens"] = ens_nums
    return corrs.load()


def load_anomaly_field(anom_dir, ens_num, year):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        ds = xr.open_mfdataset(
            f"{anom_dir}/r{ens_num}i1p1f1/*.nc", combine="by_coords", data_vars="minimal"
        ).tos
    ds = ds.groupby("time.year").mean("time")
    return ds.sel(year=year).load()


def plot_avg_timeseries(corrs_by_method):
    fig, ax = plt.subplots(figsize=(9, 5))
    for key, corrs in corrs_by_method.items():
        corrs.mean(dim="ens").plot(ax=ax, label=VARIANTS[key]["label"])
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("MPI-GE tos pattern correlation with ERA5 2023 JJA SST\n(ensemble mean, 50 members)")
    ax.set_ylabel("Pattern correlation")
    ax.set_xlabel("Year")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/pattern_corr_avg_ens_JJA_tos.png", dpi=300)
    plt.close(fig)


def plot_violin_and_count(key, corrs):
    variant = VARIANTS[key]
    decades = (corrs.time // 10) * 10
    above = corrs >= THRESHOLD
    decade_count = above.groupby(decades).sum(dim=["time", "ens"])

    data_by_decade, labels = [], []
    for decade in np.unique(decades.values):
        data_by_decade.append(corrs.isel(time=(decades == decade).values).values.flatten())
        labels.append(str(decade))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))
    ax1.violinplot(data_by_decade, positions=range(len(labels)), showmeans=True, showmedians=True)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=45)
    ax1.set_ylabel("Pattern correlation")
    ax1.set_title(f"Distribution of pattern correlation -- {variant['label']} (tos)")
    ax1.axhline(THRESHOLD, color="red", linestyle="dotted", linewidth=1)
    ax1.grid(True, alpha=0.3, axis="y")

    decade_count.plot(ax=ax2, marker="o", linewidth=2, markersize=7)
    ax2.set_xlabel("Decade")
    ax2.set_ylabel(f"Count of records >= {THRESHOLD}")
    ax2.set_title(f"Count of pattern-correlation records >= {THRESHOLD} by decade")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/violin_count_pattern_corr_JJA_{variant['suffix']}.png", dpi=300)
    plt.close(fig)


def plot_high_corr_events(key, corrs):
    variant = VARIANTS[key]
    high = corrs >= THRESHOLD
    ens_idx, time_idx = np.where(high.values)
    n_events = len(ens_idx)
    print(f"  {variant['label']}: {n_events} events with corr >= {THRESHOLD} (1850-2100)")
    if n_events == 0:
        return

    years = corrs.time.values[time_idx]
    ens_nums = corrs.ens.values[ens_idx]
    corr_vals = corrs.values[ens_idx, time_idx]

    order = np.argsort(corr_vals)[::-1]
    years, ens_nums, corr_vals = years[order], ens_nums[order], corr_vals[order]

    events = [
        load_anomaly_field(variant["anom_dir"], int(e), int(y))
        for e, y in zip(ens_nums, years)
    ]
    events_stacked = xr.concat(events, dim="event")
    avg_pattern = events_stacked.mean(dim="event")
    strongest = events[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, field, title in [
        (ax1, avg_pattern, f"Average pattern\n(N={n_events} events)"),
        (ax2, strongest, f"Strongest match\n(corr={corr_vals[0]:.3f}, ens {ens_nums[0]}, {years[0]})"),
    ]:
        field.plot(
            ax=ax, transform=ccrs.PlateCarree(), levels=np.linspace(-2, 2, 21),
            cmap="RdBu_r", center=0, cbar_kwargs={"label": "tos anomaly (K)", "shrink": 0.8}, extend="both",
        )
        ax.coastlines()
        ax.gridlines(draw_labels=True, alpha=0.3)
        ax.set_title(title)

    fig.suptitle(f"MPI-GE tos anomaly patterns, corr >= {THRESHOLD} -- {variant['label']}")
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/high_corr_patterns_JJA_{variant['suffix']}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    corrs_by_method = {key: load_corr(v["corr_dir"]) for key, v in VARIANTS.items()}

    print("Plotting ensemble-mean time series...")
    plot_avg_timeseries(corrs_by_method)

    for key, corrs in corrs_by_method.items():
        print(f"Plotting violin/count -- {VARIANTS[key]['label']}...")
        plot_violin_and_count(key, corrs)

    for key, corrs in corrs_by_method.items():
        print(f"Plotting high-correlation composite patterns -- {VARIANTS[key]['label']}...")
        plot_high_corr_events(key, corrs)

    print("Done.")
