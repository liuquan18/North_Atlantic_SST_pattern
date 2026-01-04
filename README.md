# North Atlantic SST Pattern Analysis

This project compares the pattern correlation between the SST pattern in JJA (June-July-August) 2023 with the SST pattern in MPI_GE (Max Planck Institute Grand Ensemble). All data, in monthly time scales, are available on Levante.

## Project Structure

```
North_Atlantic_SST_pattern/
├── data/           # Data directory (processed and raw data)
├── scripts/        # Analysis and processing scripts
├── src/            # Source code (Python modules)
├── test/           # Unit tests
├── doc/            # Documentation
├── environment.yml # Conda environment specification
└── requirements.txt # Python package requirements
```

## Setup on Levante

### 1. Clone the Repository

```bash
git clone https://github.com/liuquan18/North_Atlantic_SST_pattern.git
cd North_Atlantic_SST_pattern
```

### 2. Create the Conda Environment

Run the setup script:

```bash
./scripts/setup_environment.sh
```

Or manually create the environment:

```bash
conda env create -f environment.yml
conda activate north_atlantic_sst
```

Alternatively, you can use pip:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 3. Verify Installation

```bash
conda activate north_atlantic_sst
python -c "import xarray; import numpy; import matplotlib; print('All packages installed successfully!')"
```

## Data

The project uses monthly SST data from:
- Observations: JJA 2023 SST patterns
- Model: MPI_GE ensemble members

Data should be stored in the `data/` directory with appropriate subdirectories for raw and processed data.

## Preprocessing with CDO

Climate Data Operators (CDO) is used for preprocessing. Example preprocessing steps:

```bash
# Select JJA months (June, July, August)
cdo -select,month=6,7,8 input.nc jja_data.nc

# Calculate seasonal mean
cdo -timmean jja_data.nc jja_mean.nc

# Calculate anomalies
cdo -sub jja_mean.nc climatology.nc jja_anomalies.nc
```

## Pattern Correlation

The core functionality is provided by the `src/pattern_correlation.py` module:

```python
from src import calculate_pattern_correlation
import xarray as xr

# Load SST patterns
pattern_obs = xr.open_dataset('data/jja2023_sst_anomaly.nc')['sst']
pattern_model = xr.open_dataset('data/mpi_ge_sst_anomaly.nc')['sst']

# Calculate pattern correlation
correlation = calculate_pattern_correlation(pattern_obs, pattern_model)
print(f"Pattern correlation: {correlation:.3f}")
```

For xarray DataArrays with lat/lon coordinates, you can use the convenience function:

```python
from src.pattern_correlation import calculate_spatial_correlation

correlation = calculate_spatial_correlation(
    pattern_obs, 
    pattern_model,
    area_weighted=True  # Weight by cosine of latitude
)
```

## Development

### Running Tests

```bash
pytest test/
```

### Code Style

The project uses `black` for code formatting and `flake8` for linting:

```bash
black src/ scripts/ test/
flake8 src/ scripts/ test/
```

## Contributing

1. Create a new branch for your feature
2. Make your changes
3. Run tests and linters
4. Submit a pull request

## License

See LICENSE file for details.

## Contact

Quan Liu - [liuquan18](https://github.com/liuquan18)