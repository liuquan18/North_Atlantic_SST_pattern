#!/bin/bash
#SBATCH --job-name=tos_anomaly
#SBATCH --output=tos_anomaly.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mem=0
#SBATCH --account=mh0033

# Corrected MPI-GE preprocessing: uses Omon/tos (actual SST, on the native
# curvilinear ocean grid) instead of Amon/ts + sftlf land mask. tos is
# already ocean-masked, so no land-mask step is needed; it does need
# remapping onto the same regular grid ERA5 was remapped to.
#
# Runs on a single compute node (not 4+ nodes like the correlation step
# used to): this workload is I/O-bound, not CPU-bound. The main speedup
# is reusing one set of remap weights for all 50 ensemble members x ~14
# files each, instead of regenerating bilinear weights per call.

module load cdo/2.5.0-gcc-11.2.0
module load parallel

# Same regular grid the observed ERA5 anomaly was already remapped to
# (permanent archive copy -- the original scratch copy of this file has
# since been purged, see CLAUDE.md).
target_grid="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Amon/ts/gn/v20190710/ts_Amon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-186912.nc"

base_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos
weights_file=${base_dir}/remap_weights_tos_to_gn.nc
mkdir -p "$base_dir"

# Every ensemble member shares the same native ocean grid, so generate the
# remap weights exactly once and reuse them for every subsequent remap call.
sample_tos="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/tos/gn/v20190710/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-186912.nc"
if [ ! -f "$weights_file" ]; then
    echo "=========================================="
    echo "Generating remap weights (once, reused for all 50 members)"
    echo "=========================================="
    cdo -s genbil,"$target_grid" "$sample_tos" "$weights_file"
fi

process_ensemble() {
    ens=$1
    target_grid=$2
    weights_file=$3

    hist_dir="/pool/data/CMIP6/data/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r${ens}i1p1f1/Omon/tos/gn/v*/"
    ssp245_dir="/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r${ens}i1p1f1/Omon/tos/gn/v*/"

    clim_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/climatology/r${ens}i1p1f1/
    ano_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/anomaly/r${ens}i1p1f1/
    rm_glbm_dir=/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/removed_glbm/r${ens}i1p1f1/

    mkdir -p "$clim_dir" "$ano_dir" "$rm_glbm_dir"

    hist_clim_files=$(ls ${hist_dir}/tos_*.nc)
    ssp_clim_files=$(ls ${ssp245_dir}/tos_Omon_MPI-ESM1-2-LR_ssp245_r${ens}i1p1f1_gn_2015*.nc)
    clim_files="${hist_clim_files} ${ssp_clim_files}"

    echo "=========================================="
    echo "tos JJA Climatology (1991-2020) for r${ens}i1p1f1"
    echo "=========================================="

    # Climatology stays on the native ocean grid -- it's only used to
    # subtract from native-grid monthly data below. Remapping happens once,
    # on the small JJA-only anomaly result, not on every raw monthly field.
    cdo -r -f nc -ymonmean -selmon,6/8 -selyear,1991/2020 -mergetime $clim_files \
        "${clim_dir}/clim_r${ens}i1p1f1_JJA_native.nc"

    ssp_all_files=$(ls ${ssp245_dir}/tos_Omon_MPI-ESM1-2-LR_ssp245_r${ens}i1p1f1_gn_*.nc)
    all_files="${hist_clim_files} ${ssp_all_files}"

    echo "Calculating tos JJA anomalies for r${ens}i1p1f1 (1850-2100)..."

    calculate_anomaly() {
        infile=$1
        ens=$2
        ano_dir=$3
        clim_dir=$4
        target_grid=$5
        weights_file=$6

        outfile="${ano_dir}$(basename "$infile")"
        cdo -r -f nc -remap,"$target_grid","$weights_file" \
            -ymonsub -selmon,6/8 "$infile" "${clim_dir}/clim_r${ens}i1p1f1_JJA_native.nc" \
            "$outfile"
    }
    export -f calculate_anomaly
    parallel -j 5 calculate_anomaly {} "$ens" "$ano_dir" "$clim_dir" "$target_grid" "$weights_file" ::: $all_files

    # Remove the North-Atlantic-box-relative global mean, same as before,
    # now operating on the already-remapped (regular-grid) anomaly.
    echo "Removing global mean for r${ens}i1p1f1..."
    anomaly_files=$(ls ${ano_dir}/tos_*.nc)
    remove_global_mean() {
        infile=$1
        rm_glbm_dir=$2

        outfile="${rm_glbm_dir}$(basename "$infile")"
        cdo -r -f nc \
            -sellonlatbox,280,360,0,70 -sub "$infile" -enlarge,"$infile" -fldmean "$infile" \
            "$outfile"
    }
    export -f remove_global_mean
    parallel -j 5 remove_global_mean {} "$rm_glbm_dir" ::: $anomaly_files
}
export -f process_ensemble

echo "=========================================="
echo "Processing all 50 ensemble members (tos, single node, weights reused)"
echo "=========================================="
parallel -j 10 process_ensemble {} "$target_grid" "$weights_file" ::: {1..50}

echo "Done."
