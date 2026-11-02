"""
Figure 22 circulation box difference -- as fig22_boxdiff_mpige.py, with a
circulation field in place of SST: the heatwave box minus the North Atlantic
box of the SST index, and the two box means, each observed summer against the
MPI-GE summers of its climate.

The field is chosen by the argument `z300` (default) or `slp` (CIRC), as fig
28's blocking index; the outputs are named fig22_<field>_boxdiff_mpige.

(a) ERA5 EUR − NAL: the JJA anomaly of the field averaged over the heatwave
    box (src.heatwave.MEAN_BOX, 35–60°N, 10°W–30°E) minus that over the North
    Atlantic box of the SST index (p2.DIPOLE_BOXES['nal'], 40–50°N, 60–20°W),
    on the 50th / 90th / 99th percentiles and mean of the 500 MPI-GE summers in
    the 10 years around each year. A ridge / high over Europe with a trough /
    low upstream, the pattern of a blocked summer, is positive.
(b) the same for the heatwave box mean alone and
(c) for the North Atlantic box mean alone (one y axis). Box means as fig 28's
    z300_index / slp_index: cos-lat weighted over all cells (land and sea),
    anomalies against 1991–2020 of the same dataset and member.

  z300  300 hPa geopotential height, m (data/sst_anomaly/<ds>_jja_z300_anom_1deg.nc,
        pre_process/09_era5_zg.sh, 10_model_zg.sh). Its rise with the warming
        column, common to both boxes, largely cancels in (a).
  slp   mean sea-level pressure, hPa (data/circulation/<ds>_jja_slp_2p5deg.nc,
        pre_process/11_slp.sh).

Marked and right-hand box plots as in fig 22; the runs with upper-air or SLP
output differ (CIRC 'runs').
"""
import importlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
fig18 = importlib.import_module("fig18_era5_return_period_mpige")
fig22 = importlib.import_module("fig22_era5_hwd_corr_mpige")
fig28 = importlib.import_module("fig28_hwd_sst_blocking_regression")
boxdiff = importlib.import_module("fig22_boxdiff_mpige")
import src.heatwave as hw
import src.pattern_2026 as p2

NAL = p2.DIPOLE_BOXES["nal"]
LABEL = fig22.Z300_LABEL
#: the two fields: loader (key, box) -> (member, year) Series, short name,
#: unit, long name, caption note, and the box plots' runs
CIRC = {
    "z300": dict(load=fig22.z300_m, name="Z300", unit="m",
                 long="300 hPa geopotential height", desc="a ridge over Europe with a trough upstream",
                 note=("the rise of Z300 with the warming column, which largely cancels in (a), "
                       "shows here"),
                 runs=("ICON-EPOC: 1 run, 1991–2025; ICON-EERIE: r1, 1991–2030; MPI-ESM1.2-ER: "
                       "r2, 1991–2030; the other members have no upper-air output")),
    "slp": dict(load=fig28.slp_index, name="SLP", unit="hPa", long="mean sea-level pressure",
                desc="high pressure over Europe with low pressure upstream",
                note="unlike Z300, SLP has little trend, so (a) mainly sharpens the contrast",
                runs=("ICON-EPOC: 1 run, 1991–2025; ICON-EERIE: 3, 1991–2030; MPI-ESM1.2-ER: "
                      "r2, 1991–2030; r1 and r3 have no SLP over 1991–2020")),
}


def rows(var):
    """(loader key -> (member, year) Series, y label) of each panel for field `var`."""
    c = CIRC[var]
    load, name, unit = c["load"], c["name"], c["unit"]
    return [(lambda k: load(k) - load(k, NAL), f"{LABEL} {name} anomaly ({unit})"),
            (lambda k: load(k), f"Heatwave-box {name}\nanomaly ({unit})"),
            (lambda k: load(k, NAL), f"North Atlantic {name}\nanomaly ({unit})")]


def caption(var):
    """Title and caption of the figure for field `var`."""
    c = CIRC[var]
    title = (f"The Europe–North Atlantic {c['name']} contrast and its two boxes against the "
             "simulated summers of its climate")
    text = (
        f"(a) ERA5 box-difference index {LABEL} (black): JJA {c['long']} anomaly averaged over "
        f"the heatwave box ({hw.MEAN_BOX_LABEL}) minus that over the North Atlantic box of the "
        f"SST index ({NAL['detail']}), in {c['unit']}, land and sea, anomalies against 1991–2020 "
        f"of the same dataset; positive: {c['desc']}, as in a blocked summer. Blue: the same "
        "index in MPI-ESM1.2-LR: the 50th, 90th and 99th percentiles (solid) and the mean "
        "(dashed) of its 50 members in the 10 years around each year (y−5 to y+4, 500 summers). "
        "(b) The same for the heatwave box mean alone and (c) for the North Atlantic box mean "
        f"alone (same y axis): {c['note']}. "
        "Marked: 2003, 2026 and every summer ranked at or above 90% of its window in heatwave "
        f"days. Right: the same quantities in ERA5, {fig22.BOX_YEARS[0]}–{p2.REF_YEAR}, and in "
        "the higher-resolution runs, all historical + scenario summers "
        f"{fig22.BOX_YEARS[0]}–{fig22.BOX_YEARS[1]} where available, pooled over members "
        f"({c['runs']}; sample size in brackets; box 25th–75th percentile, line median, "
        "whiskers 5th–95th percentile, dots beyond).")
    return title, text


def main(var="z300"):
    meth = fig18.METHODS["rank"]
    win_hwd = pd.read_csv(fig18.DIR / f"{meth['windows']}.csv")
    years, _, _ = fig18.load_series("hwd")
    marked = sorted(int(c.rsplit("_", 1)[1]) for c in win_hwd if c.startswith("pct_fixed_"))
    panels = rows(var)
    series = [load("ERA5").xs(1, level="member").to_xarray()
              .sel(year=slice(int(years[0]), p2.REF_YEAR)).dropna("year") for load, _ in panels]
    wins = [fig22.window_levels(load("MPI-GE").to_xarray(), win_hwd.year.values)
            for load, _ in panels]
    boxes = fig22.box_samples([load for load, _ in panels], hwd=False)
    title, text = caption(var)
    boxdiff.figure(win_hwd.year.values, series, wins, [y for _, y in panels], marked, boxes,
                   title, text, p2.FIG_DIR / f"fig22_{var}_boxdiff_mpige.png")


if __name__ == "__main__":
    main(*sys.argv[1:2])
