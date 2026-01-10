#!/bin/bash

module load cdo/2.5.0-gcc-11.2.0
module load parallel

era_file=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/ERA5/anomaly/anom_2023_jja.nc

# Function to process a single ensemble member
process_ensemble() {
    ens=$1
    
    corr_dir=/work/mh0033/m300883/North_Atlantic_SST_pattern/data/MPI_GE_CMIP6/pattern_corr/r${ens}i1p1f1/
    mkdir -p $corr_dir

    MPI_GE_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6/anomaly/r${ens}i1p1f1/

    mpi_ge_files=$(ls ${MPI_GE_dir}/*.nc)
    
    echo "Processing ensemble r${ens}i1p1f1 - Number of files: $(echo $mpi_ge_files | wc -w)"

    # Function to calculate pattern correlation for a single file
    pattern_corr() {
        infile=$1
        ens=$2
        corr_dir=$3
        era_file=$4
        
        outfile="${corr_dir}$(basename $infile)"
        cdo -r -f nc -fldcor -remapbil,$infile $era_file $infile $outfile
    }

    # Export the function for parallel
    export -f pattern_corr

    # Use parallel to process all files (7 parallel jobs per ensemble)
    parallel -j 7 pattern_corr {} $ens $corr_dir $era_file ::: $mpi_ge_files
}

# Export the function for parallel
export -f process_ensemble
export era_file

# Process all ensemble members in parallel (1 to 50)
echo "=========================================="
echo "Calculating pattern correlation for ensembles 1 to 50"
echo "=========================================="
parallel -j 1 process_ensemble ::: {1..50}