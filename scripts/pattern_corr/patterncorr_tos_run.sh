#!/bin/bash
#SBATCH --job-name=corr_MPI_ERA5_tos
#SBATCH --output=corr_MPI_ERA5_tos.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=00:30:00
#SBATCH --partition=compute
#SBATCH --mem=0
#SBATCH --account=mh0033

# Corrected + simplified: the old job used 4 nodes / 40 MPI ranks for a
# workload that is I/O-bound, not CPU-bound, and looped year-by-year in
# Python. Now that the correlation itself is vectorized (one xr.corr call
# per ensemble covering its whole time series), a single node with a small
# process pool over the 50 ensembles is enough -- no MPI needed at all.

module load python3/unstable

source /sw/spack-levante/mambaforge-23.1.0-1-Linux-x86_64-3boc6i/etc/profile.d/conda.sh
conda activate north_atlantic_sst
export PYTHONNOUSERSITE=1

cd /work/mh0033/m300883/North_Atlantic_SST_pattern
python3 scripts/pattern_corr/MPI_ERA5_patterncorr_tos.py
