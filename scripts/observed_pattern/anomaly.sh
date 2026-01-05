#!/bin/bash

# Define directories
era_dir="/pool/data/ERA5/E5/sf/an/1M/034"
output_dir="/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5"
clim_dir="${output_dir}/climatology"
anom_dir="${output_dir}/anomaly"

# Create output directories
mkdir -p "$clim_dir"
mkdir -p "$anom_dir"

# Load CDO
module load cdo/2.5.0-gcc-11.2.0

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

# 1. June Climatology
echo "Calculating June Climatology..."
cdo -O -f nc timmean -selmon,6 -mergetime $clim_files "${clim_dir}/clim_jun.nc"

# 2. July Climatology
echo "Calculating July Climatology..."
cdo -O -f nc timmean -selmon,7 -mergetime $clim_files "${clim_dir}/clim_jul.nc"

# 3. August Climatology
echo "Calculating August Climatology..."
cdo -O -f nc timmean -selmon,8 -mergetime $clim_files "${clim_dir}/clim_aug.nc"

# 4. JJA Climatology
echo "Calculating JJA Climatology..."
cdo -O -f nc timmean -selmon,6,7,8 -mergetime $clim_files "${clim_dir}/clim_jja.nc"


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

# 1. June Anomaly
echo "Calculating June 2023 Anomaly..."
# Select June 2023, subtract June Climatology
cdo -O -f nc sub -selmon,6 "$file_2023" "${clim_dir}/clim_jun.nc" "${anom_dir}/anom_2023_jun.nc"

# 2. July Anomaly
echo "Calculating July 2023 Anomaly..."
cdo -O -f nc sub -selmon,7 "$file_2023" "${clim_dir}/clim_jul.nc" "${anom_dir}/anom_2023_jul.nc"

# 3. August Anomaly
echo "Calculating August 2023 Anomaly..."
cdo -O -f nc sub -selmon,8 "$file_2023" "${clim_dir}/clim_aug.nc" "${anom_dir}/anom_2023_aug.nc"

# 4. JJA Anomaly
echo "Calculating JJA 2023 Anomaly..."
# Calculate JJA mean for 2023 first, then subtract JJA Climatology
cdo -O -f nc sub -timmean -selmon,6,7,8 "$file_2023" "${clim_dir}/clim_jja.nc" "${anom_dir}/anom_2023_jja.nc"

echo "Done."

