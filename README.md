# North Atlantic SST Pattern Analysis

Compares the observed **JJA 2026** sea surface temperature pattern over the
North Atlantic + Mediterranean (20–60°N, 80°W–40°E) against climate
simulations by spatial pattern correlation: the MPI-ESM1-2-LR grand ensemble,
MPI-ESM1.2-ER, **km-scale ICON (EPOC, 10 km atm / 5 km ocean)** and EERIE
ICON-ESM-ER. Sea surface temperature only (`tos` / `to`). All data are on
Levante.

Methodology, data traps and results: **[doc/pattern_2026.md](doc/pattern_2026.md)**.

## Project structure

```
North_Atlantic_SST_pattern/
├── scripts/            pipeline (see scripts/README.md)
│   ├── config.sh       shared settings
│   ├── run_analysis.sh analysis + all figures
│   ├── pre_process/    cdo stages per dataset (sbatch)
│   ├── analysis/       pattern correlation, region means
│   └── plotting/       one script per figure
├── src/                pattern_2026.py (analysis helpers), viz2026.py (styling)
├── test/               unit tests for src/
├── doc/                write-ups
├── data/pattern_2026/  outputs (not tracked)
├── figures/pattern_2026/
└── logs/               SLURM logs (not tracked)
```

## Setup on Levante

```bash
source /sw/spack-levante/mambaforge-23.1.0-1-Linux-x86_64-3boc6i/etc/profile.d/conda.sh
conda env create -f environment.yml     # first time only
conda activate north_atlantic_sst
export PYTHONNOUSERSITE=1               # keeps ~/.local numpy from shadowing the env's
```

CDO is loaded as a module (`module load cdo/2.5.0-gcc-11.2.0`); `scripts/config.sh`
does this for the pipeline.

## Running

```bash
sbatch scripts/pre_process/01_era5.sh   # ... through 06, see scripts/README.md
bash scripts/run_analysis.sh
```

## Tests

```bash
pytest test/
```

## Contact

Quan Liu - [liuquan18](https://github.com/liuquan18)
