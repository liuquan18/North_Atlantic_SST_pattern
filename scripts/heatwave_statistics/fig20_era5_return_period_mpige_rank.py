"""
Figure 20 -- as figure 18 (heatwave days), but without any fitted distribution:
the return period of each ERA5 summer from its rank among the 500 MPI-GE
summers of its window (05_mpige_window_rank.py hwd).
"""
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")

if __name__ == "__main__":
    fig18.main("hwd", method="rank")
