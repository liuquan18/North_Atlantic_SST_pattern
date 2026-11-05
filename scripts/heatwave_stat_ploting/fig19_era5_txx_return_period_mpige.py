"""
Figure 19 -- as figure 18, for TXx: the hottest day of each JJA in land-mean
daily maximum temperature over southern Europe (04_txx_table.py), as anomaly
from each dataset's own 1991-2020 mean, scored against the stationary GEV of
the 500 MPI-GE summers around each year (03_mpige_window_rp.py txx).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

fig18 = importlib.import_module("fig18_era5_return_period_mpige")

if __name__ == "__main__":
    fig18.main("txx")
