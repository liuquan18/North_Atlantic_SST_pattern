#!/bin/bash
module load cdo/2.5.0-gcc-11.2.0
module load parallel


Junfile=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_jun.nc
Julfile=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_jul.nc
Augfile=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_aug.nc

JJAfile=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_jja.nc

# calculate the gridcorr between all three months and JJA
echo "Calculating grid correlation between June and JJA..."

out_dir=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/self_pattern_corr/
mkdir -p $out_dir

cdo -fldcor $JJAfile $Junfile ${out_dir}/gridcorr_ERA5_jun_jja.nc

cdo -fldcor $JJAfile $Julfile ${out_dir}/gridcorr_ERA5_jul_jja.nc

cdo -fldcor $JJAfile $Augfile ${out_dir}/gridcorr_ERA5_aug_jja.nc