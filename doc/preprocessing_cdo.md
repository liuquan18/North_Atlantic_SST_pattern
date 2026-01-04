# Data Preprocessing with CDO

This document describes the data preprocessing workflow using CDO (Climate Data Operators).

## Prerequisites

CDO should be available in your conda environment. Verify installation:

```bash
cdo --version
```

## Data Sources

### Observations (JJA 2023)
- Source: [Specify observation dataset, e.g., HadISST, NOAA OI SST]
- Location on Levante: `/path/to/observations/`
- Variable: Sea Surface Temperature (SST)

### Model Data (MPI_GE)
- Source: MPI Grand Ensemble
- Location on Levante: `/path/to/mpi_ge/`
- Variable: Sea Surface Temperature
- Ensemble members: Multiple realizations

## Preprocessing Steps

### 1. Select JJA Season (June, July, August)

For observation data:
```bash
cdo -select,month=6,7,8 /path/to/obs/sst_monthly_2023.nc data/obs_jja2023.nc
```

For model data (example for one ensemble member):
```bash
cdo -select,month=6,7,8 /path/to/mpi_ge/member001/sst_monthly.nc data/mpi_ge_jja_member001.nc
```

### 2. Calculate Seasonal Mean

Average over the JJA season:
```bash
# Observations
cdo -timmean data/obs_jja2023.nc data/obs_jja2023_mean.nc

# Model
cdo -timmean data/mpi_ge_jja_member001.nc data/mpi_ge_jja_mean_member001.nc
```

### 3. Calculate Climatology (if needed)

For model climatology (using multiple years):
```bash
# Select JJA for all years
cdo -select,month=6,7,8 /path/to/mpi_ge/member001/sst_monthly_all_years.nc data/mpi_ge_jja_all.nc

# Calculate climatological mean
cdo -timmean data/mpi_ge_jja_all.nc data/mpi_ge_climatology.nc
```

### 4. Calculate Anomalies

Subtract climatology from the pattern:
```bash
# For observations (if you have obs climatology)
cdo -sub data/obs_jja2023_mean.nc data/obs_climatology.nc data/obs_jja2023_anomaly.nc

# For model
cdo -sub data/mpi_ge_jja_mean_member001.nc data/mpi_ge_climatology.nc data/mpi_ge_anomaly_member001.nc
```

### 5. Regional Selection (North Atlantic)

Select the North Atlantic region (example coordinates):
```bash
# Select latitude: 0-70N, longitude: 280-360E (80W-0E)
cdo -sellonlatbox,280,360,0,70 data/obs_jja2023_anomaly.nc data/obs_jja2023_anomaly_natl.nc
cdo -sellonlatbox,280,360,0,70 data/mpi_ge_anomaly_member001.nc data/mpi_ge_anomaly_natl_member001.nc
```

### 6. Regridding (if needed)

If observation and model grids differ:
```bash
# Regrid model to observation grid
cdo remapbil,data/obs_jja2023_anomaly_natl.nc data/mpi_ge_anomaly_natl_member001.nc data/mpi_ge_anomaly_natl_member001_regrid.nc
```

## Complete Preprocessing Pipeline

Here's a complete script example:

```bash
#!/bin/bash
# preprocess_sst.sh

# Define paths
OBS_RAW="/path/to/observations/sst_2023.nc"
MODEL_RAW="/path/to/mpi_ge/member001/sst.nc"
MODEL_CLIM="/path/to/mpi_ge/member001/sst_all_years.nc"

# Create output directories
mkdir -p data/processed
mkdir -p data/intermediate

# Process observations
echo "Processing observations..."
cdo -select,month=6,7,8 -selyear,2023 $OBS_RAW data/intermediate/obs_jja2023.nc
cdo -timmean data/intermediate/obs_jja2023.nc data/intermediate/obs_jja2023_mean.nc
cdo -sellonlatbox,280,360,0,70 data/intermediate/obs_jja2023_mean.nc data/processed/obs_jja2023_natl.nc

# Process model data
echo "Processing model data..."
# Calculate model climatology
cdo -select,month=6,7,8 $MODEL_CLIM data/intermediate/model_jja_all.nc
cdo -timmean data/intermediate/model_jja_all.nc data/intermediate/model_climatology.nc

# Process 2023 model data
cdo -select,month=6,7,8 -selyear,2023 $MODEL_RAW data/intermediate/model_jja2023.nc
cdo -timmean data/intermediate/model_jja2023.nc data/intermediate/model_jja2023_mean.nc

# Calculate anomaly
cdo -sub data/intermediate/model_jja2023_mean.nc data/intermediate/model_climatology.nc data/intermediate/model_jja2023_anomaly.nc

# Select North Atlantic region
cdo -sellonlatbox,280,360,0,70 data/intermediate/model_jja2023_anomaly.nc data/intermediate/model_jja2023_anomaly_natl.nc

# Regrid to observation grid
cdo remapbil,data/processed/obs_jja2023_natl.nc data/intermediate/model_jja2023_anomaly_natl.nc data/processed/model_jja2023_natl.nc

echo "Preprocessing complete!"
echo "Output files:"
echo "  - data/processed/obs_jja2023_natl.nc"
echo "  - data/processed/model_jja2023_natl.nc"
```

## Useful CDO Commands

### Information Commands
```bash
# Show file information
cdo sinfon file.nc

# Show variable names
cdo showname file.nc

# Show time information
cdo showtimestamp file.nc
```

### Quality Checks
```bash
# Check for missing values
cdo info file.nc

# Calculate statistics
cdo fldmean file.nc  # Field mean
cdo fldstd file.nc   # Field standard deviation
```

## Notes

- Always check the coordinate systems and units before processing
- Ensure time axes are properly formatted
- For large files, CDO operations can be chained with pipes for efficiency
- Consider using CDO's `-O` flag to overwrite existing files during development

## Resources

- [CDO Documentation](https://code.mpimet.mpg.de/projects/cdo)
- [CDO Reference Card](https://code.mpimet.mpg.de/projects/cdo/embedded/cdo_refcard.pdf)
