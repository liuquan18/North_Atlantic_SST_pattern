#!/bin/bash
module load cdo/2.5.0-gcc-11.2.0

# Define directories
era_dir="/pool/data/ERA5/E5/sf/an/1M/034"
output_dir="/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5"
inter_dir="/scratch/m/m300883/nalt/"
clim_dir="${output_dir}/climatology"
anom_dir="${output_dir}/anomaly"

# Create output directories
mkdir -p "$clim_dir"
mkdir -p "$anom_dir"

# Load CDO
module load cdo/2.5.0-gcc-11.2.0

# Target grid for remapping
target_grid="/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r1i1p1f1/ts_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc"

# Define the climatology period
start_year=1991
end_year=2020

# Construct the list of files for the climatology period
# We assume the file naming convention: E5sf00_1M_YYYY_034.grb
clim_files=""
for year in $(seq $start_year $end_year); do
    clim_files="${clim_files} ${era_dir}/E5sf00_1M_${year}_034.grb"
done

echo "=========================================="
echo "Calculating Climatologies (1991-2020)"
echo "=========================================="


# climatology calculations
cdo -r -f nc -remapbil,$target_grid -setgridtype,regular -ymonmean -selmon,6,7,8 -mergetime $clim_files ${clim_dir}/clim_nalt_1991_06-2020_08.nc

echo "=========================================="
echo "Calculating Anomalies for 2023"
echo "=========================================="
year_anom=2023
file_2023="${era_dir}/E5sf00_1M_${year_anom}_034.grb"

# Check if 2023 file exists
if [ ! -f "$file_2023" ]; then
    echo "Error: File for 2023 not found: $file_2023"
    exit 1
fi

echo "Calculating 2023 Anomaly..."

cdo -r -f nc -ymonsub -selmon,6,7,8 -remapbil,$target_grid -setgridtype,regular "$file_2023" "${clim_dir}/clim_nalt_1991_06-2020_08.nc" "${anom_dir}/anom_2023_06_2023_08.nc"


# 4. JJA Anomaly
echo "Calculating JJA Anomaly for 2023..."
cdo -r -f nc -timmean "${anom_dir}/anom_2023_06_2023_08.nc" "${anom_dir}/anom_2023_jja.nc"

# remove global mean and select region
echo "Removing global mean and selecting North Atlantic region..."
cdo -r -f nc -sellonlatbox,280,360,0,70 -sub "${anom_dir}/anom_2023_jja.nc" -enlarge,"${anom_dir}/anom_2023_jja.nc" -fldmean "${anom_dir}/anom_2023_jja.nc" "${anom_dir}/anom_2023_jja_rm_glbm.nc"