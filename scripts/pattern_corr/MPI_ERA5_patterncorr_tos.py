"""
Pattern correlation between ERA5 2023 JJA SST and MPI-GE tos anomalies.

Corrected version: uses the actual SST variable (Omon/tos, already remapped
onto ERA5's grid by scripts/MPI_GE_pattern/anomaly_tos.sh) instead of Amon/ts
+ land mask, and replaces the old 4-node/40-rank MPI job with a single
process on one node. The per-timestep loop that used to compute the
correlation year-by-year is gone too (src/pattern_correlation.py now
vectorizes it) -- the actual per-ensemble work is a few seconds of I/O plus
one xr.corr call, so a small process pool over the 50 members is enough.

Produces both methodologies from doc/pattern_corr_JJA.pptx:
  1. remove spatial mean before correlating   -> data/pattern_corr_JJA_tos/
  2. remove global mean before correlating     -> data/pattern_corr_JJA_tos_noglm/
"""

import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import xarray as xr

import src.pattern_correlation as pc

PROJECT_ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
ANOMALY_BASE = Path("/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos")

ERA5_JJA = PROJECT_ROOT / "data/ERA5/anomaly/anom_2023_jja.nc"
ERA5_JJA_RM_GLBM = PROJECT_ROOT / "data/ERA5/anomaly/anom_2023_jja_rm_glbm.nc"

OUT_SPATIAL_MEAN_REMOVED = PROJECT_ROOT / "data/pattern_corr_JJA_tos"
OUT_GLOBAL_MEAN_REMOVED = PROJECT_ROOT / "data/pattern_corr_JJA_tos_noglm"

N_WORKERS = 10


def read_ensemble_data(ens_num, variant_dir):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        model_data = xr.open_mfdataset(
            f"{variant_dir}/r{ens_num}i1p1f1/*.nc",
            combine="by_coords",
            data_vars="minimal",
        ).tos
    model_data = model_data.groupby("time.year").mean("time")
    model_data = model_data.rename({"year": "time"})
    return model_data.load()


def process_one(ens_num, era5_file, variant_dir, out_dir, rm_spatial_mean):
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"mpi_era5_patterncorr_ens{ens_num:02d}.nc"
    if output_file.exists():
        return ens_num, f"{output_file} (already done, skipped)"

    era5_data = xr.open_dataset(era5_file).squeeze()
    var_name = list(era5_data.data_vars)[-1]
    era5_data = era5_data[var_name]

    model_data = read_ensemble_data(ens_num, variant_dir)
    corr = pc.spatial_corr_model_era(model_data, era5_data, rm_spatial_mean=rm_spatial_mean)
    corr.to_netcdf(output_file)
    return ens_num, str(output_file)


def run_variant(label, era5_file, variant_dir, out_dir, rm_spatial_mean):
    print(f"=== {label}: 50 ensembles, {N_WORKERS} workers, 1 node ===", flush=True)
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {
            pool.submit(process_one, ens, era5_file, variant_dir, out_dir, rm_spatial_mean): ens
            for ens in range(1, 51)
        }
        for future in as_completed(futures):
            ens_num, output_file = future.result()
            print(f"  ensemble {ens_num:02d} -> {output_file}", flush=True)


if __name__ == "__main__":
    run_variant(
        "Method 1: remove spatial mean",
        ERA5_JJA,
        ANOMALY_BASE / "anomaly",
        OUT_SPATIAL_MEAN_REMOVED,
        rm_spatial_mean=True,
    )
    run_variant(
        "Method 2: remove global mean",
        ERA5_JJA_RM_GLBM,
        ANOMALY_BASE / "removed_glbm",
        OUT_GLOBAL_MEAN_REMOVED,
        rm_spatial_mean=False,
    )
    print("Done.", flush=True)
