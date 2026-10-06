"""
Figure 1 (heatwave analogues) -- the same layout as fig1_patterns_simulations.py,
but each model's season is chosen to be high on BOTH axes of figure 8: a
2026-like SST pattern and many European heatwave days.

Candidates are the seasons (and members) of figure 8: JJA 1990-2025, with
both an SST pattern and daily Tmax. Within each dataset every season is
ranked on its pattern correlation with ERA5 JJA 2026 and on its land-mean
JJA heatwave days over the whole European map; its score is the LOWER of the
two percentile ranks, and the season with the highest score is shown. A
season therefore wins only by being near the top of both distributions --
unlike the highest r alone (fig1_patterns_simulations.py), or a sum of the
two, where the long hot tail of heatwave days dominates. This also holds
where the fitted slope of figure 8 is ~0 (EERIE).

MPI-GE is the exception: among its 1,800 seasons that rule lands on the long
hot tail (many heatwave days, moderate r), far above the relationship the
fitted line describes. Its season is instead the highest-r one among those
close to the line (fig8_hwd_vs_corr.SELECTION / along_fit).

  top row     heatwave days of that season over the whole of Europe, with
              their land mean
  bottom row  the SST pattern of the same season (box mean removed), with its
              pattern correlation against ERA5 JJA 2026

MPI-ESM1.2-ER has no daily Tmax, so it has no season to select and is left out.
"""
import sys
from pathlib import Path

import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.heatwave as hw
import src.pattern_2026 as p2
from fig1_patterns_simulations import VARIANT, column, land_mean_series, plot_figure
from scipy import stats
from fig8_hwd_vs_corr import PERIOD, SELECTION, season_table, select_season

KEYS = ["MPI-GE", "ICON-EPOC-hist", "EERIE"]


def top_percent(rank):
    """Percentile rank as 'top x %' (rank 1.0 is the single highest season)."""
    return f"top {max(100 * (1 - rank), 0.1):.0f} %" if rank < 0.995 else "top 1 %"


def main():
    corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{VARIANT}.nc").load()
    with xr.open_dataset(hw.metrics_file("ERA5")) as ds:
        obs = float(land_mean_series(ds.hwd.sel(year=p2.REF_YEAR).load()))

    cols = []
    for key in KEYS:
        t = season_table(key, corr)
        s = select_season(key, t, stats.linregress(t.r, t.hwd))
        year = int(s.year)
        mem = int(s.member) if s.member == s.member else None
        d = p2.DATASETS[key]
        when = f"JJA {year}" + (f"  ·  member r{mem}" if mem else "")
        print(f"{key:16s} {when}: r = {s.r:+.2f} ({top_percent(s.rank_r)}), "
              f"{s.hwd:.1f} heatwave days ({top_percent(s.rank_hwd)}) of {len(t)} seasons")
        c = column(key, year, mem, float(s.r), f"{d.label}\n{d.resolution}\n{when}")
        cols.append(c)

    plot_figure(
        cols,
        "European summer heatwaves and the North Atlantic + Mediterranean SST pattern: "
        "each model's season with both a 2026-like pattern and many heatwave days",
        p2.FIG_DIR / "fig1_heatwave_analogues_2026.png",
        hw_note=f"\nseason: highest lower percentile rank of (pattern correlation, land-mean "
                f"heatwave days) among the dataset's JJA {PERIOD[0]}–{PERIOD[1]} seasons and "
                f"members; {', '.join(k for k, v in SELECTION.items() if v == 'along_fit')}: "
                "highest pattern correlation close to the fitted line of figure 8\n"
                f"top-row value: land-mean heatwave days (ERA5 JJA {p2.REF_YEAR}: "
                f"{obs:.1f} days)")


if __name__ == "__main__":
    main()
