#!/bin/bash
#SBATCH --job-name=corr_MPI_ERA5
#SBATCH --output=corr_MPI_ERA5.%j.out
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=10
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --account=mh0033

# Load modules
module load python3/unstable
module load intel-oneapi-mpi/2021.5.0-intel-2021.5.0

# Activate conda
source ~/.bashrc
conda activate north_atlantic_sst

# Set PMI library
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

# Launch MPI job directly - let SLURM auto-detect MPI interface
# srun python MPI_ERA5_patterncorr.py
srun python MPI_ERA5_patterncorr_rmglm.py