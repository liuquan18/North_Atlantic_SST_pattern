#!/bin/bash
# Example preprocessing script for SST data using CDO
# This script demonstrates the workflow for preprocessing SST data
# Adjust paths and parameters according to your data location on Levante

set -e  # Exit on error

# Configuration
# TODO: Update these paths to match your data location on Levante
OBS_DATA_PATH="/path/to/observations/sst_monthly_2023.nc"
MODEL_DATA_PATH="/path/to/mpi_ge/member001/sst_monthly.nc"
MODEL_CLIM_PATH="/path/to/mpi_ge/member001/sst_monthly_all_years.nc"

# Output directories
DATA_DIR="$(dirname "$0")/../data"
INTERMEDIATE_DIR="${DATA_DIR}/intermediate"
PROCESSED_DIR="${DATA_DIR}/processed"

# Create output directories
mkdir -p "${INTERMEDIATE_DIR}"
mkdir -p "${PROCESSED_DIR}"

echo "=========================================="
echo "SST Data Preprocessing with CDO"
echo "=========================================="
echo ""

# Check if CDO is available
if ! command -v cdo &> /dev/null; then
    echo "Error: CDO is not available. Please load the CDO module or activate the conda environment."
    echo "Try: module load cdo"
    echo "Or: conda activate north_atlantic_sst"
    exit 1
fi

echo "CDO version: $(cdo --version | head -n1)"
echo ""

# ==========================================
# Process Observations
# ==========================================
echo "Step 1: Processing observation data..."
echo "----------------------------------------"

if [ ! -f "${OBS_DATA_PATH}" ]; then
    echo "Warning: Observation data not found at ${OBS_DATA_PATH}"
    echo "Please update the OBS_DATA_PATH variable in this script."
else
    # Select JJA months
    echo "  - Selecting JJA months (June, July, August)..."
    cdo -select,month=6,7,8 -selyear,2023 "${OBS_DATA_PATH}" "${INTERMEDIATE_DIR}/obs_jja2023.nc"
    
    # Calculate seasonal mean
    echo "  - Calculating JJA seasonal mean..."
    cdo -timmean "${INTERMEDIATE_DIR}/obs_jja2023.nc" "${INTERMEDIATE_DIR}/obs_jja2023_mean.nc"
    
    # Select North Atlantic region (0-70N, 280-360E)
    echo "  - Selecting North Atlantic region..."
    cdo -sellonlatbox,280,360,0,70 "${INTERMEDIATE_DIR}/obs_jja2023_mean.nc" "${PROCESSED_DIR}/obs_jja2023_natl.nc"
    
    echo "  ✓ Observation processing complete"
    echo "    Output: ${PROCESSED_DIR}/obs_jja2023_natl.nc"
fi

echo ""

# ==========================================
# Process Model Data
# ==========================================
echo "Step 2: Processing model data..."
echo "----------------------------------------"

if [ ! -f "${MODEL_DATA_PATH}" ]; then
    echo "Warning: Model data not found at ${MODEL_DATA_PATH}"
    echo "Please update the MODEL_DATA_PATH variable in this script."
elif [ ! -f "${MODEL_CLIM_PATH}" ]; then
    echo "Warning: Model climatology data not found at ${MODEL_CLIM_PATH}"
    echo "Please update the MODEL_CLIM_PATH variable in this script."
else
    # Calculate model climatology
    echo "  - Calculating model JJA climatology..."
    cdo -select,month=6,7,8 "${MODEL_CLIM_PATH}" "${INTERMEDIATE_DIR}/model_jja_all.nc"
    cdo -timmean "${INTERMEDIATE_DIR}/model_jja_all.nc" "${INTERMEDIATE_DIR}/model_climatology.nc"
    
    # Process 2023 model data
    echo "  - Processing 2023 model data..."
    cdo -select,month=6,7,8 -selyear,2023 "${MODEL_DATA_PATH}" "${INTERMEDIATE_DIR}/model_jja2023.nc"
    cdo -timmean "${INTERMEDIATE_DIR}/model_jja2023.nc" "${INTERMEDIATE_DIR}/model_jja2023_mean.nc"
    
    # Calculate anomaly
    echo "  - Calculating anomalies..."
    cdo -sub "${INTERMEDIATE_DIR}/model_jja2023_mean.nc" "${INTERMEDIATE_DIR}/model_climatology.nc" \
        "${INTERMEDIATE_DIR}/model_jja2023_anomaly.nc"
    
    # Select North Atlantic region
    echo "  - Selecting North Atlantic region..."
    cdo -sellonlatbox,280,360,0,70 "${INTERMEDIATE_DIR}/model_jja2023_anomaly.nc" \
        "${INTERMEDIATE_DIR}/model_jja2023_anomaly_natl.nc"
    
    # Regrid to observation grid (if obs data exists)
    if [ -f "${PROCESSED_DIR}/obs_jja2023_natl.nc" ]; then
        echo "  - Regridding to observation grid..."
        cdo remapbil,"${PROCESSED_DIR}/obs_jja2023_natl.nc" \
            "${INTERMEDIATE_DIR}/model_jja2023_anomaly_natl.nc" \
            "${PROCESSED_DIR}/model_jja2023_natl.nc"
    else
        cp "${INTERMEDIATE_DIR}/model_jja2023_anomaly_natl.nc" "${PROCESSED_DIR}/model_jja2023_natl.nc"
    fi
    
    echo "  ✓ Model processing complete"
    echo "    Output: ${PROCESSED_DIR}/model_jja2023_natl.nc"
fi

echo ""
echo "=========================================="
echo "Preprocessing Complete!"
echo "=========================================="
echo ""
echo "Processed files are in: ${PROCESSED_DIR}/"
echo ""
echo "Next steps:"
echo "  1. Verify the output files"
echo "  2. Run pattern correlation analysis"
echo ""
echo "Example verification commands:"
echo "  cdo sinfon ${PROCESSED_DIR}/obs_jja2023_natl.nc"
echo "  cdo fldmean ${PROCESSED_DIR}/obs_jja2023_natl.nc"
echo ""
