#!/bin/bash

module load cdo/2.5.0-gcc-11.2.0
module load parallel

ens=1

hist_dir="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r${ens}i1p1f1/Amon/ts/gn/v*/"
ssp245_dir="/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r${ens}i1p1f1/Amon/ts/gn/v*/"

clim_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6/climatology/r${ens}i1p1f1/
ano_dir=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/MPI_GE_CMIP6/anomaly/r${ens}i1p1f1/

mkdir -p $clim_dir
mkdir -p $ano_dir

# Define mask file (Land Area Fraction)
# sftlf: Percentage of the grid cell occupied by land (including lakes)
mask_file="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/fx/sftlf/gn/v20190710/sftlf_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc"
ocean_mask="${clim_dir}/ocean_mask_natl.nc"

# Create regional ocean mask if it doesn't exist
# 1. Select region
# 2. Select ocean (sftlf < 1%) -> gives 1 for ocean, 0 for land
# 3. Set 0 (land) to missing value
if [ ! -f "$ocean_mask" ]; then
    echo "Creating ocean mask..."
    cdo -setctomiss,0 -ltc,1 -sellonlatbox,280,360,0,70 "$mask_file" "$ocean_mask"
fi

# to compute climatology for 1991-2020, find the files from both the historical and ssp245 directories
hist_clim_files=$(ls ${hist_dir}/ts_Amon_MPI-ESM1-2-LR_historical_r${ens}i1p1f1_gn_199*.nc ${hist_dir}/ts_Amon_MPI-ESM1-2-LR_historical_r${ens}i1p1f1_gn_20*.nc)
ssp_clim_files=$(ls ${ssp245_dir}/ts_Amon_MPI-ESM1-2-LR_ssp245_r${ens}i1p1f1_gn_2015*.nc)
clim_files="${hist_clim_files} ${ssp_clim_files}"

echo "=========================================="
echo "Calculating Climatologies (1991-2020) for MPI-ESM1-2-LR r${ens}i1p1f1"
echo "=========================================="

# 1. June Climatology
echo "Calculating June Climatology..."
# Apply ocean mask using ifthen
cdo -r -f nc -ifthen "$ocean_mask" -timmean -sellonlatbox,280,360,0,70 -selmon,6 -mergetime $clim_files "${clim_dir}/clim_jun.nc"