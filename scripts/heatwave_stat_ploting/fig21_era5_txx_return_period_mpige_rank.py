"""
Figure 21 -- as figure 19 (TXx), but without any fitted distribution or
return period: the percentile rank (0-100%) of each ERA5 summer among the 500
MPI-GE summers of its window (05_mpige_window_rank.py txx).
"""
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")

if __name__ == "__main__":
    fig18.main("txx", method="rank")
