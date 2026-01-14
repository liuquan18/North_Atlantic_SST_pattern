#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
data_2023_14 =xr.open_dataset("/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r14i1p1f1/ts_Amon_MPI-ESM1-2-LR_ssp245_r14i1p1f1_gn_201501-203412.nc")
# %%
data_2023_14 = data_2023_14.sel(time = '2023')
# %%
data_2023_14_mean = data_2023_14.mean(dim='time').ts
# %%
data_2023_14_mean.plot(levels = np.arange(-1.8,2.,0.1), cmap = 'RdBu_r')
# %%
corr = xr.open_dataset("/work/mh0033/m300883/North_Atlantic_SST_pattern/data/pattern_corr/mpi_era5_patterncorr_ens14.nc")
# %%
corr.sel(time = '2023').mean(dim = 'time').__xarray_dataarray_variable__.values
# %%
