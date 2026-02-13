# %%
import xarray as xr
import numpy as np
import warnings

# %%
import src.pattern_correlation as pc

spatial_corr_model_era = pc.spatial_corr_model_era
# %%
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
normalize = False

# Debug: Verify MPI is working
print(f"MPI initialized: Rank {rank} of {size} processes", flush=True)
if rank == 0:
    print(
        f"Total ensembles to process: 50, distributed across {size} ranks", flush=True
    )
# %%
#era5_data = "anom_2023_jja.nc"
era5_data = "anom_2023_jja_rm_glbm.nc"

# %%
ERA5_data = xr.open_dataset(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/" + era5_data
).var34.squeeze()


# %%
def read_ensemble_data(ens_num):
    # Suppress FutureWarning about compat parameter
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        model_data = xr.open_mfdataset(
            # f"/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r{ens_num}i1p1f1/*.nc",
            f"/scratch/m/m300883/nalt/MPI_GE_CMIP6/removed_glbm/r{ens_num}i1p1f1/*.nc",
            combine="by_coords",
            data_vars="minimal",
        ).ts
    model_data = model_data.groupby("time.year").mean("time")
    # rename year to time
    model_data = model_data.rename({"year": "time"})
    return model_data

#%%
# function to do spatial normalize
def spatial_normalize(data):
    spatial_mean = data.mean(dim=["lat", "lon"])
    spatial_anomaly = data - spatial_mean
    spatial_std = spatial_anomaly.std(dim=["lat", "lon"])
    normalized_data = spatial_anomaly / spatial_std
    return normalized_data


# %%

all_ens = np.arange(1, 51)
ens_for_rank = np.array_split(all_ens, size)[rank]


# %%
for i, ens_num in enumerate(ens_for_rank):
    print(
        f"Rank {rank} processing ensemble {ens_num} ({i+1}/{len(ens_for_rank)})...",
        flush=True,
    )
    model_data = read_ensemble_data(ens_num)
    if normalize:
        model_data = spatial_normalize(model_data)
        ERA5_data = spatial_normalize(ERA5_data)
    corr = spatial_corr_model_era(model_data, ERA5_data, rm_spatial_mean=False)

    # Save the correlation results for this ensemble member
    output_file = f"/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr_JJA_noglm/mpi_era5_patterncorr_ens{ens_num:02d}.nc"
    corr.to_netcdf(output_file)
    print(f"Rank {rank} saved ensemble {ens_num} to {output_file}", flush=True)
# %%
