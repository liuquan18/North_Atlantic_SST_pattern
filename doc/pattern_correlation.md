# Pattern Correlation Analysis

## Overview

Pattern correlation is a statistical measure used to quantify the spatial similarity between two fields (e.g., SST patterns). It is essentially the Pearson correlation coefficient computed over spatial points, optionally weighted by grid cell area.

## Mathematical Definition

For two spatial patterns X and Y, the pattern correlation r is defined as:

```
r = Σᵢ wᵢ(xᵢ - x̄)(yᵢ - ȳ) / √[Σᵢ wᵢ(xᵢ - x̄)² × Σᵢ wᵢ(yᵢ - ȳ)²]
```

where:
- xᵢ, yᵢ are the values at grid point i
- wᵢ is the weight at grid point i (e.g., cos(latitude) for area weighting)
- x̄, ȳ are the weighted means
- The sum is over all valid (non-NaN) grid points

## Implementation Details

The `calculate_pattern_correlation` function in `src/pattern_correlation.py` implements this calculation with the following features:

### Features
1. **Flexible Input**: Accepts both numpy arrays and xarray DataArrays
2. **Automatic Masking**: NaN values are automatically excluded
3. **Weighted Correlation**: Optional weights (e.g., for area weighting)
4. **Centered/Uncentered**: Option to compute correlation with or without mean removal

### Area Weighting

For global or regional SST patterns, it's important to account for the fact that grid cells have different areas. Near the equator, grid cells are larger than near the poles. We use cosine of latitude as weights:

```python
weights = cos(latitude)
```

This ensures that each unit of Earth's surface area contributes equally to the correlation.

## Usage Examples

### Basic Usage

```python
import numpy as np
from src import calculate_pattern_correlation

# Simple patterns (no weighting)
pattern1 = np.random.randn(10, 20)
pattern2 = np.random.randn(10, 20)
corr = calculate_pattern_correlation(pattern1, pattern2)
print(f"Pattern correlation: {corr:.3f}")
```

### With Area Weighting

```python
import numpy as np
import xarray as xr
from src.pattern_correlation import calculate_spatial_correlation

# Load data
obs_data = xr.open_dataset('data/observations.nc')['sst']
model_data = xr.open_dataset('data/model.nc')['sst']

# Calculate with area weighting
corr = calculate_spatial_correlation(
    obs_data,
    model_data,
    area_weighted=True
)
print(f"Area-weighted pattern correlation: {corr:.3f}")
```

### Working with Anomalies

Pattern correlation is typically computed on anomalies (deviations from climatology):

```python
import xarray as xr
from src.pattern_correlation import calculate_spatial_correlation

# Load data
obs = xr.open_dataset('data/jja2023_sst.nc')['sst']
obs_clim = xr.open_dataset('data/sst_climatology.nc')['sst']

# Calculate anomalies
obs_anom = obs - obs_clim

# Similarly for model
model = xr.open_dataset('data/mpi_ge_member001.nc')['sst']
model_clim = xr.open_dataset('data/mpi_ge_climatology.nc')['sst']
model_anom = model - model_clim

# Calculate pattern correlation
corr = calculate_spatial_correlation(obs_anom, model_anom, area_weighted=True)
```

## Interpretation

- **r ≈ 1**: Strong positive similarity (patterns are very similar)
- **r ≈ 0**: No linear relationship (patterns are uncorrelated)
- **r ≈ -1**: Strong negative similarity (patterns are opposite)

For SST patterns:
- r > 0.7: Generally considered good agreement
- r > 0.8: Very good agreement
- r > 0.9: Excellent agreement

## References

1. Taylor, K. E. (2001). Summarizing multiple aspects of model performance in a single diagram. Journal of Geophysical Research, 106(D7), 7183-7192.
2. Wilks, D. S. (2011). Statistical methods in the atmospheric sciences (Vol. 100). Academic press.
