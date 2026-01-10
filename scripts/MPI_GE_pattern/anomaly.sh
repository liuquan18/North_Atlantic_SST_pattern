#!/bin/bash

module load cdo/2.5.0-gcc-11.2.0
module load parallel

# Function to process a single ensemble member
process_ensemble() {
    ens=$1

    hist_dir="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r${ens}i1p1f1/Amon/ts/gn/v*/"
    ssp245_dir="/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r${ens}i1p1f1/Amon/ts/gn/v*/"

    clim_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6/climatology/r${ens}i1p1f1/
    ano_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r${ens}i1p1f1/

    mkdir -p $clim_dir
    mkdir -p $ano_dir

    # Define mask file (Land Area Fraction)
    # sftlf: Percentage of the grid cell occupied by land (including lakes)
    mask_file="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/fx/sftlf/gn/v20190710/sftlf_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc"
    ocean_mask="/scratch/m/m300883/nalt/MPI_GE_CMIP6/climatology/ocean_mask_natl.nc"

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

    # 1. JJA climatology
    echo "Calculating JJA Climatology..."
    cdo -r -f nc -ifthen "$ocean_mask" -selmon,6/8 -ymonmean -sellonlatbox,280,360,0,70 -selyear,1991/2020 -mergetime $clim_files "${clim_dir}/clim_r${ens}i1p1f1_JJA.nc"

    # 2. calculate anomaly for each year in 1991-2100
    # Get all historical and SSP245 files
    ssp_all_files=$(ls ${ssp245_dir}/ts_Amon_MPI-ESM1-2-LR_ssp245_r${ens}i1p1f1_gn_*.nc)
    all_files="${hist_clim_files} ${ssp_all_files}"

    echo "Calculating JJA Anomalies for each year (1991-2100)..."
    echo "Number of files to process: $(echo $all_files | wc -w)"

    # Function to calculate anomaly for a single file
    calculate_anomaly() {
        infile=$1
        ens=$2
        ano_dir=$3
        clim_dir=$4
        ocean_mask=$5
        
        outfile="${ano_dir}$(basename $infile)"
        cdo -r -f nc -ifthen "$ocean_mask" -ymonsub -selmon,6/8 -sellonlatbox,280,360,0,70 $infile "${clim_dir}/clim_r${ens}i1p1f1_JJA.nc" "$outfile"
    }

    # Export the function for parallel
    export -f calculate_anomaly

    # Use parallel to process all files (7 parallel jobs per ensemble)
    parallel -j 7 calculate_anomaly {} $ens $ano_dir $clim_dir $ocean_mask ::: $all_files
}

# Export the function for parallel
export -f process_ensemble

# Process all ensemble members in parallel (1 to 50)
echo "=========================================="
echo "Processing ensemble members 1 to 50"
echo "=========================================="
parallel -j 1 process_ensemble ::: {1..50}
