# %%
import xarray as xr
import numpy as np

# %%
import src.pattern_correlation as pc

spatial_corr_model_era = pc.spatial_corr_model_era
# %%
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# %%
era5_data = "anom_2023_jja.nc"
# %%
ERA5_data = xr.open_dataset(
    "/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/" + era5_data
).var34.squeeze()


# %%
def read_ensemble_data(ens_num):
    model_data = xr.open_mfdataset(
        f"/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r{ens_num}i1p1f1/*.nc",
        combine="by_coords",
        data_vars="minimal",
    ).ts
    return model_data


# %%

all_ens = np.arange(1, 51)
ens_for_rank = np.array_split(all_ens, size)[rank]


# %%
for i, ens_num in enumerate(ens_for_rank):
    print(f"Rank {rank} processing ensemble {ens_num} ({i+1}/{len(ens_for_rank)})...")
    model_data = read_ensemble_data(ens_num)
    corr = spatial_corr_model_era(model_data, ERA5_data)

    # Save the correlation results for this ensemble member
    output_file = f"/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr/mpi_era5_patterncorr_ens{ens_num:02d}.nc"
    corr.to_netcdf(output_file)
    print(f"Rank {rank} saved ensemble {ens_num} to {output_file}")
# %%
