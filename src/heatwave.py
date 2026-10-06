"""
European summer heatwaves, following the definition of Xu et al. (2026,
Nature Climate Change, doi:10.1038/s41558-026-02762-2):

    "A heatwave was defined as a period of at least three consecutive days
    during which daily maximum temperature anomalies -- calculated relative to
    the multi-year mean for each Julian day -- exceeded the calendar-day-specific
    90th-percentile threshold calculated using a 15-day moving window."

Differences from the paper:

  * Reference period 1991-2020 (the paper uses 1979-2023, its full record).
    Neither the km-scale runs nor EPOC reach back to 1979, and 1991-2020 is the
    climatology the SST anomalies in src/pattern_2026.py are taken against.
  * Each dataset is referenced to its *own* 1991-2020 climatology and
    threshold, exactly as each dataset's SST anomaly is.

Everything else follows the paper: the climatology is the raw (unsmoothed)
calendar-day mean; the threshold is the 90th percentile of the *anomalies*
pooled over the reference years and a centred 15-day window (+-7 days, so
30 x 15 = 450 samples per calendar day); heatwaves are only detected inside
the Northern Hemisphere warm-season window, May-September, so an event is cut
at the window edges.

The JJA metrics are then taken from the days of those events that fall in
June-August.

All functions work on plain numpy arrays shaped (time, point) so that the
land-only grid can be packed into one dimension; the driver script handles the
xarray bookkeeping.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
DATA_DIR = PROJECT_ROOT / "data" / "heatwave_2026"

REF_PERIOD = (1991, 2020)
PERCENTILE = 90.0
HALF_WINDOW = 7          # 15-day centred window
MIN_DURATION = 3         # days

# day-of-year on a 365-day calendar (Feb 29 never enters: only Apr-Oct is used)
DETECT_DOY = (121, 273)  # 1 May .. 30 Sep
JJA_DOY = (152, 243)     # 1 Jun .. 31 Aug

#: Europe domain of the common 0.25 deg grid (nodes on multiples of 0.25 deg,
#: coincident with the EERIE `gr` grid). Must match GRID_EU025 in
#: scripts/european_heatwave/hw_config.sh.
REGION = {"lon": (-15.0, 45.0), "lat": (30.0, 72.0)}

#: dataset key (as in src/pattern_2026.py) -> file stem. MPI-ESM1.2-ER is absent:
#: its echam6 output was only written as monthly means, so there is no daily
#: Tmax to detect heatwaves in.
DATASETS = {
    "ERA5": "era5",
    "MPI-GE": "mpige",
    "ICON-EPOC-hist": "epoc2_020",
    "EERIE": "eerie_hist_ssp245",
}

#: Datasets with several realizations on the common grid, each read from its
#: own r<m>/ subdirectory of the daily Tmax work dir; their metrics carry a
#: `member` dimension, NaN in years a member was not run.
MEMBERS = {"EERIE": [1, 2, 3]}

#: Datasets kept on their native grid (cropped to Europe) instead of the common
#: 0.25 deg one, with (member, year, lat, lon) metrics and their own land mask.
NATIVE = {"MPI-GE"}


def grid_tag(key: str) -> str:
    return "native" if key in NATIVE else "eu025"


def metrics_file(key: str) -> Path:
    return DATA_DIR / f"{DATASETS[key]}_heatwave_{grid_tag(key)}.nc"


def landmask_file(key: str) -> Path:
    return DATA_DIR / (f"landmask_{DATASETS[key]}_native.nc" if key in NATIVE
                       else "landmask_eu025.nc")


# --------------------------------------------------------------------------
# calendar
# --------------------------------------------------------------------------
def noleap_doy(times) -> np.ndarray:
    """Day of year on a 365-day calendar: 1 Mar is always 60, 1 May always 121."""
    t = pd.DatetimeIndex(times)
    doy = t.dayofyear.values.copy()
    doy[t.is_leap_year & (t.month > 2)] -= 1
    return doy


# --------------------------------------------------------------------------
# climatology and threshold
# --------------------------------------------------------------------------
def daily_climatology(tmax: np.ndarray, doy: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """
    Mean Tmax for each calendar day over the reference years.

    Returns an array (366, point) indexed by doy; calendar days without
    reference data are NaN.
    """
    clim = np.full((366,) + tmax.shape[1:], np.nan, dtype=np.float64)
    for d in np.unique(doy[ref]):
        clim[d] = np.nanmean(tmax[ref & (doy == d)], axis=0)
    return clim


def percentile_threshold(anom: np.ndarray, doy: np.ndarray, ref: np.ndarray,
                         doys, q: float = PERCENTILE,
                         half_window: int = HALF_WINDOW) -> np.ndarray:
    """
    Calendar-day q-th percentile of the anomalies, pooled over the reference
    years and the days within +-half_window of that calendar day.

    Only the calendar days in `doys` are filled (the detection window); the
    window around them must be covered by the data, which is why the
    pre-processing keeps April-October rather than just May-September.
    Returns (366, point), NaN outside `doys`.
    """
    thr = np.full((366,) + anom.shape[1:], np.nan, dtype=np.float64)
    for d in doys:
        sel = ref & (np.abs(doy - d) <= half_window)
        thr[d] = np.nanpercentile(anom[sel], q, axis=0)
    return thr


# --------------------------------------------------------------------------
# event detection
# --------------------------------------------------------------------------
def run_length(mask: np.ndarray) -> np.ndarray:
    """
    For every True element, the length of the consecutive run (along axis 0)
    it belongs to; 0 where False.
    """
    mask = np.asarray(mask, dtype=bool)
    n = mask.shape[0]
    fwd = np.zeros(mask.shape, dtype=np.int32)
    bwd = np.zeros(mask.shape, dtype=np.int32)
    for t in range(n):
        fwd[t] = np.where(mask[t], (fwd[t - 1] if t else 0) + 1, 0)
    for t in range(n - 1, -1, -1):
        bwd[t] = np.where(mask[t], (bwd[t + 1] if t < n - 1 else 0) + 1, 0)
    return np.where(mask, fwd + bwd - 1, 0)


def season_metrics(anom: np.ndarray, thr: np.ndarray, doy: np.ndarray,
                   min_duration: int = MIN_DURATION) -> dict[str, np.ndarray]:
    """
    Heatwave metrics for one year.

    anom, thr : (day, point) daily Tmax anomaly and its threshold for every
                day of that year's data (any calendar days; only those inside
                the detection window are used).
    doy       : (day,) 365-day-calendar day of year, consecutive.

    Returns per point:
      hwd   number of JJA days that belong to a heatwave
      hwn   number of heatwave events with at least one day in JJA
      hwcum cumulative heat: sum of the Tmax anomaly over JJA heatwave days (degC day)
      hwmax longest heatwave (full length, days) with at least one day in JJA
      hwpeak hottest JJA heatwave day: its Tmax anomaly (degC; NaN if none)
      hwexcess excess heat: sum over JJA heatwave days of the anomaly ABOVE the
            threshold (degC day) -- the cumulative-heat definition of
            Perkins-Kirkpatrick & Lewis (2020), counting only heat beyond the
            90th percentile rather than beyond the climatological mean
      onset first heatwave day in the May-Sep detection window (day of year; NaN if none)
    """
    win = (doy >= DETECT_DOY[0]) & (doy <= DETECT_DOY[1])
    if np.any(np.diff(doy[win]) != 1):
        raise ValueError("detection window is not a consecutive daily series")
    a, th, d = anom[win], thr[win], doy[win]

    hot = a > th                          # NaN compares False
    rl = run_length(hot)
    hw = rl >= min_duration

    jja = (d >= JJA_DOY[0]) & (d <= JJA_DOY[1])
    hw_jja = hw[jja]

    # an event is counted for JJA if it starts in JJA, or is already running on 1 June
    start = hw & ~np.vstack([np.zeros((1,) + hw.shape[1:], bool), hw[:-1]])
    start_jja = start[jja]
    hwn = start_jja.sum(0) + (hw_jja[0] & ~start_jja[0])

    any_hw = hw.any(0)
    onset = np.where(any_hw, d[np.argmax(hw, axis=0)], np.nan)

    return {
        "hwd": hw_jja.sum(0).astype(np.int16),
        "hwn": hwn.astype(np.int16),
        "hwcum": np.where(hw_jja, a[jja], 0.0).sum(0),
        "hwmax": np.where(hw_jja, rl[jja], 0).max(0).astype(np.int16),
        "hwpeak": np.where(hw_jja.any(0), np.where(hw_jja, a[jja], -np.inf).max(0), np.nan),
        "hwexcess": np.where(hw_jja, a[jja] - th[jja], 0.0).sum(0),
        "onset": onset,
    }
