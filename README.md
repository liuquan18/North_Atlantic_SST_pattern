# North Atlantic SST Pattern Analysis

This project compares an observed summer (JJA) North Atlantic sea surface
temperature pattern against climate simulations, by spatial pattern
correlation. All data are on Levante.

Two studies live here:

- **JJA 2026, North Atlantic + Mediterranean** (current) — 30–60°N, 80°W–40°E,
  scored against MPI-GE, **km-scale ICON (EPOC, 10 km atm / 5 km ocean)** and
  **EERIE ICON-ESM-ER**. Sea surface temperature only (`tos` / `to`); no `ts`
  or `tas`. Pipeline, data traps and methodology:
  **[scripts/pattern_2026/README.md](scripts/pattern_2026/README.md)**.
  Outputs in `data/pattern_2026/` and `figures/pattern_2026/`.
- **JJA 2023, North Atlantic** (original) — 0–70°N, 80°W–0°, MPI-GE only.
  Scripts in `scripts/pattern_corr/`, helpers in `src/pattern_correlation.py`.

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

### 3. Install the Package (Optional but Recommended)

For easier imports and development, install the package in editable mode:

```bash
# After activating your environment
pip install -e .

# Or with development dependencies
pip install -e ".[dev,viz]"
```

This allows you to import the package from anywhere:
```python
from src import calculate_pattern_correlation
# or
from src.pattern_correlation import calculate_spatial_correlation
```

### 4. Verify Installation

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