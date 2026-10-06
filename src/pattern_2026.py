"""
Shared machinery for the JJA 2026 North Atlantic + Mediterranean SST pattern study.

Reference pattern: observed (ERA5) JJA 2026 SST anomaly over 30-60N, 80W-40E,
relative to a 1991-2020 JJA climatology. Every model dataset is scored against
it with the same two metrics used in the earlier 2023 work:

  "spatial"  remove the area-weighted mean over the analysis box from both
             fields, then correlate -- a pure shape comparison, blind to how
             warm the basin is overall.
  "global"   remove each dataset's own global-ocean-mean JJA SST anomaly from
             both fields, then take the uncentred (cosine) correlation -- the
             basin-mean warmth that survives after the global warming signal is
             taken out still counts as part of the pattern.

The correlation takes the spatial dimensions explicitly, so a field that
also carries an ensemble `member` dimension is handled without member being
mistaken for a spatial axis.

All datasets arrive from the cdo stage on a common 0.25 deg regional grid. The
correlation itself is computed after coarsening to 1 deg: MPI-ESM1-2-LR's ocean
grid is ~1 deg, and scoring it on 0.25 deg structure it cannot represent would
flatter the km-scale runs for the wrong reason.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

PROJECT_ROOT = Path("/work/mh0033/m300883/North_Atlantic_SST_pattern")
DATA_DIR = PROJECT_ROOT / "data" / "pattern_2026"
RESULT_DIR = DATA_DIR / "results"
FIG_DIR = PROJECT_ROOT / "figures" / "pattern_2026"

REF_YEAR = 2026
CLIM_PERIOD = (1991, 2020)
REGION = {"lon": (-80.0, 40.0), "lat": (30.0, 60.0)}

#: Fraction of a 1 deg cell that must be valid ocean in *every* dataset for the
#: cell to enter the correlation.
MIN_OCEAN_FRACTION = 0.5


# --------------------------------------------------------------------------
# dataset registry
# --------------------------------------------------------------------------
class Dataset:
    """One SST source, already reduced to the two common products by the cdo stage."""

    def __init__(self, key, label, anom_file, gm_file, kind, resolution, period, note=""):
        self.key = key
        self.label = label
        self.anom_file = DATA_DIR / anom_file
        self.gm_file = DATA_DIR / gm_file
        self.kind = kind              # "obs" | "model"
        self.resolution = resolution  # human-readable native ocean resolution
        self.period = period
        self.note = note

    def __repr__(self):
        return f"<Dataset {self.key}>"


DATASETS = {
    d.key: d
    for d in [
        Dataset("ERA5", "ERA5 (observed)",
                "era5_jja_anom_na025.nc", "era5_jja_gmsst.nc",
                "obs", "0.25 deg", (1940, 2026),
                "JJA 2026 built from hourly ERA5T (Jul/Aug not yet in the monthly stream)"),
        Dataset("MPI-GE", "MPI-ESM1-2-LR grand ensemble",
                "mpige_jja_anom_na025.nc", "mpige_jja_gmsst.nc",
                "model", "~1 deg ocean", (1850, 2100),
                "50 members, historical + ssp245"),
        # epoc2_020: irad_co2=3 with greenhouse_ssp585.nc and bc_ozone_ssp585_<year>,
        # i.e. transient greenhouse-gas forcing -- but Kinne aerosols pinned at 1850.
        # MPI-ESM1.2-ER: same model family as MPI-GE above, but with the
        # eddy-resolving TP6M ocean instead of the ~1 deg one -- so the pair
        # isolates ocean resolution within a single model.
        Dataset("MPI-ER", "MPI-ESM1.2-ER",
                "mpier_jja_anom_na025.nc", "mpier_jja_gmsst.nc",
                "model", "~100 km atm / 10 km ocean", (1950, 2099),
                "3 realizations, historical + ssp585"),
        Dataset("ICON-EPOC-hist", "km-scale ICON (EPOC, transient GHG)",
                "epoc2_020_jja_anom_na025.nc", "epoc2_020_jja_gmsst.nc",
                "model", "10 km atm / 5 km ocean", (1990, 2025),
                "historical + SSP5-8.5 greenhouse gases, aerosols held at 1850"),
        # epoc2_010: irad_co2=2 with aerosols and ozone pinned at 1990 -- a
        # constant-forcing control, not a historical run.
        Dataset("ICON-EPOC-ctrl", "km-scale ICON (EPOC, control)",
                "epoc2_010_jja_anom_na025.nc", "epoc2_010_jja_gmsst.nc",
                "model", "10 km atm / 5 km ocean", (1990, 2024),
                "constant 1990 forcing"),
        Dataset("EERIE", "EERIE ICON-ESM-ER (hist+ssp245)",
                "eerie_hist_ssp245_jja_anom_na025.nc", "eerie_hist_ssp245_jja_gmsst.nc",
                "model", "10 km atm / 5 km ocean", (1950, 2050),
                "3 realizations: r1 1950-2050, r2 1975-2014, r3 1975-2020"),
        Dataset("EERIE-ctrl", "EERIE ICON-ESM-ER (control)",
                "eerie_control_jja_anom_na025.nc", "eerie_control_jja_gmsst.nc",
                "model", "10 km atm / 5 km ocean", (1950, 2050),
                "constant 1950 forcing"),
    ]
}

#: The columns of figure 1, in order. MPI-ER sits next to MPI-GE so the two
#: ocean resolutions of the same model read as a pair.
MAIN_KEYS = ["ERA5", "MPI-GE", "MPI-ER", "ICON-EPOC-hist", "EERIE"]


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def _single_var(ds: xr.Dataset) -> xr.DataArray:
    """The one SST field in the file (named var34 / tos / to depending on source)."""
    names = [v for v in ds.data_vars if ds[v].ndim >= 1 and not v.endswith("_bnds")]
    if len(names) != 1:
        names = [v for v in names if v not in ("time_bnds", "lat_bnds", "lon_bnds")]
    return ds[names[0]]


def _to_year_axis(da: xr.DataArray) -> xr.DataArray:
    """Replace the (mid-season) time axis with a plain integer `year` coordinate."""
    if "time" in da.dims:
        try:                       # datetime64 or cftime axis
            years = da["time"].dt.year.values
        except (TypeError, AttributeError):
            years = da["time"].values   # already a bare year index
        da = da.rename({"time": "year"}).assign_coords(
            year=("year", np.asarray(years, dtype=int)))
    da = da.drop_vars([c for c in da.coords if c not in da.dims], errors="ignore")
    return da


def load_anomaly(key: str, chunks="auto") -> xr.DataArray:
    """
    Regional JJA SST anomaly field, dims (year[, member], lat, lon), degC.

    Chunked by default: the 50-member MPI-GE file is ~2.9 GB on the 0.25 deg
    grid, and every downstream step (coarsening, pattern removal, correlation)
    reduces it, so there is no reason to hold it whole.
    """
    ds = DATASETS[key]
    da = _to_year_axis(
        _single_var(xr.open_dataset(ds.anom_file, chunks=chunks)).squeeze(drop=True))
    return da.rename("sst_anom").assign_attrs(dataset=ds.key, label=ds.label)


def load_gmsst(key: str) -> xr.DataArray:
    """Global-ocean-mean JJA SST anomaly, dims (year[, member]), degC."""
    ds = DATASETS[key]
    da = _to_year_axis(_single_var(xr.open_dataset(ds.gm_file)).squeeze(drop=True))
    return da.rename("gmsst").assign_attrs(dataset=ds.key).load()


# --------------------------------------------------------------------------
# grid handling
# --------------------------------------------------------------------------
def area_weights(da: xr.DataArray) -> xr.DataArray:
    """cos(latitude) area weights broadcast over the lat/lon grid."""
    return np.cos(np.deg2rad(da["lat"])).broadcast_like(da["lat"] * da["lon"]).rename("weights")


def coarsen_to_1deg(da: xr.DataArray, min_fraction: float = MIN_OCEAN_FRACTION) -> xr.DataArray:
    """
    Area-weighted 4x4 coarsening of the 0.25 deg regional grid onto 1 deg.

    Blocks that are less than `min_fraction` valid ocean are returned as NaN, so
    a 1 deg cell that is mostly land does not enter the correlation on the
    strength of one or two wet 0.25 deg cells.
    """
    w = np.cos(np.deg2rad(da["lat"])) * xr.ones_like(da["lon"])
    valid = da.notnull()
    wv = w.where(valid)
    num = (da.fillna(0) * wv.fillna(0)).coarsen(lat=4, lon=4, boundary="exact").sum()
    den = wv.fillna(0).coarsen(lat=4, lon=4, boundary="exact").sum()
    frac = (w * valid).coarsen(lat=4, lon=4, boundary="exact").sum() / \
           w.coarsen(lat=4, lon=4, boundary="exact").sum()
    # divide by a NaN-masked denominator rather than by zero: all-land blocks
    # would otherwise emit an invalid-value warning on every call
    out = (num / den.where(den > 0)).where(frac >= min_fraction)
    return out.assign_attrs(da.attrs)


def common_ocean_mask(fields: dict[str, xr.DataArray]) -> xr.DataArray:
    """
    Grid cells that are valid ocean in *every* dataset.

    Land masks differ between ERA5, the ~1 deg MPI ocean grid and the two
    km-scale ICON grids -- especially around the Mediterranean coast and the
    Black Sea. Without intersecting them the datasets would be correlated over
    slightly different domains and the numbers would not be comparable.
    """
    mask = None
    for da in fields.values():
        # a cell counts as ocean only if it is valid in every season and member
        # that exists: seasons a member was not run (EERIE r2/r3 cover fewer
        # years than r1) are all-NaN padding and must not mask every cell
        present = da.notnull().any(["lat", "lon"])
        m = (da.notnull() | ~present).all([d for d in da.dims if d not in ("lat", "lon")])
        mask = m if mask is None else (mask & m)
    return mask.rename("ocean_mask")


# --------------------------------------------------------------------------
# the two pattern definitions
# --------------------------------------------------------------------------
def remove_spatial_mean(da: xr.DataArray) -> xr.DataArray:
    """Subtract the area-weighted mean over the analysis box."""
    w = np.cos(np.deg2rad(da["lat"])) * xr.ones_like(da["lon"])
    w = w.where(da.notnull())
    return da - (da * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))


def remove_global_mean(da: xr.DataArray, gmsst: xr.DataArray) -> xr.DataArray:
    """Subtract the dataset's own global-ocean-mean JJA SST anomaly, year by year."""
    return da - gmsst


def make_pattern(da: xr.DataArray, gmsst: xr.DataArray, variant: str) -> xr.DataArray:
    if variant == "spatial":
        return remove_spatial_mean(da)
    if variant == "global":
        return remove_global_mean(da, gmsst)
    raise ValueError(f"unknown variant {variant!r}")


def pattern_corr(field: xr.DataArray, ref: xr.DataArray, *, centered: bool,
                 spatial_dims=("lat", "lon")) -> xr.DataArray:
    """
    Area-weighted spatial pattern correlation of `field` against the 2D `ref`.

    `centered=True`  removes each field's own weighted spatial mean first
                     (Pearson correlation over space) -- used with the
                     "spatial" variant.
    `centered=False` leaves the means in place (uncentred cosine similarity)
                     -- used with the "global" variant, where the basin-mean
                     anomaly left after removing global-mean warming is itself
                     part of the pattern being compared.

    Any dimensions of `field` beyond `spatial_dims` (year, member, ...) are
    broadcast over, so the whole ensemble-by-year matrix comes out of one call.
    """
    spatial_dims = list(spatial_dims)
    w = np.cos(np.deg2rad(field["lat"])) * xr.ones_like(field["lon"])

    valid = field.notnull() & ref.notnull()
    x = field.where(valid)
    y = ref.where(valid)
    wv = w.where(valid)

    if centered:
        x = x - (x * wv).sum(spatial_dims) / wv.sum(spatial_dims)
        y = y - (y * wv).sum(spatial_dims) / wv.sum(spatial_dims)

    num = (wv * x * y).sum(spatial_dims)
    den = np.sqrt((wv * x**2).sum(spatial_dims) * (wv * y**2).sum(spatial_dims))
    return (num / den).where(den > 0)


#: Variants drawn in the figures. The analysis computes and stores both, but
#: they give near-identical correlations here -- JJA 2026's basin-mean warmth is
#: almost all the global signal (box +0.66 degC against a global ocean mean of
#: +0.56), so removing the box mean and removing the global mean leave nearly
#: the same field. The "global mean removed" numbers stay in
#: results/corr_global.nc and best_analogues.csv.
PLOT_VARIANTS = ["spatial"]

VARIANTS = {
    "spatial": dict(centered=True,
                    title="spatial mean removed",
                    long="area-weighted box mean removed (centred pattern correlation)"),
    "global": dict(centered=False,
                   title="global mean removed",
                   long="global-ocean-mean JJA anomaly removed (uncentred pattern correlation)"),
}


# --------------------------------------------------------------------------
# regional averages
# --------------------------------------------------------------------------
#: Sub-regions of the analysis box, for area-mean time series.
#:
#: "natl" is the analysis box with the Mediterranean box cut out, so that
#: box = natl + med exactly and the three curves decompose rather than
#: duplicate each other. Note the Mediterranean box, taken literally, also
#: holds a strip of open Atlantic west of Gibraltar: that is 11% of its ocean
#: area and moves its JJA 2026 mean by 0.03 degC, so it is left in rather than
#: silently redefining the box.
MED_BOX = {"lat": (30.0, 50.0), "lon": (-10.0, 40.0)}

#: `detail` is kept short enough to sit in a figure's left margin.
REGIONS = {
    "global": dict(label="Global ocean",
                   detail="each dataset's own grid"),
    "box": dict(label="Analysis box",
                detail="30–60°N, 80°W–40°E"),
    "natl": dict(label="North Atlantic",
                 detail="box minus Mediterranean"),
    "med": dict(label="Mediterranean",
                detail="30–50°N, 10°W–40°E"),
}


def region_masks(ocean_mask: xr.DataArray) -> dict[str, xr.DataArray]:
    """Boolean masks for 'box', 'natl' and 'med' on the common 1 deg grid."""
    lat, lon = ocean_mask["lat"], ocean_mask["lon"]
    in_med = ((lat >= MED_BOX["lat"][0]) & (lat <= MED_BOX["lat"][1])
              & (lon >= MED_BOX["lon"][0]) & (lon <= MED_BOX["lon"][1]))
    return {
        "box": ocean_mask,
        "med": ocean_mask & in_med,
        "natl": ocean_mask & ~in_med,
    }


def member_dim(da: xr.DataArray) -> str | None:
    """
    Name of the ensemble dimension of `da`, or None if it has none.

    Ensembles of different size must not share a dimension *name* inside one
    Dataset: xarray would align MPI-GE's 50 members and MPI-ESM1.2-ER's 3 onto
    a common index, padding the smaller one to 50 with NaN and making it report
    50 members. Writers therefore store the dimension as `member_<key>`; this
    accepts either spelling so older files keep working.
    """
    for d in da.dims:
        if d == "member" or str(d).startswith("member_"):
            return str(d)
    return None


def area_mean(da: xr.DataArray, mask: xr.DataArray | None = None) -> xr.DataArray:
    """cos(lat)-weighted mean over lat/lon, skipping missing cells."""
    field = da.where(mask) if mask is not None else da
    w = (np.cos(np.deg2rad(field["lat"])) * xr.ones_like(field["lon"])).where(field.notnull())
    return (field * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))


# --------------------------------------------------------------------------
# event composites (figure 7)
# --------------------------------------------------------------------------
#: With at least this many members, the ensemble mean is a clean estimate of
#: the forced part of a correlation series; below it, it still carries a large
#: share of each member's own internal variability.
ENSEMBLE_MEAN_MIN_MEMBERS = 10

#: Records at least this long get a quadratic forced baseline, shorter ones a
#: linear one. A parabola fitted to a ~90-yr record can bend into a single
#: 60-70-yr swing, which is exactly the variability under study, and remove it.
QUADRATIC_MIN_YEARS = 100


def internal_component(r: xr.DataArray, *, control: bool = False) -> xr.DataArray:
    """
    A correlation series with its forced part removed.

    Composites of raw correlation slope upward through any forced record,
    because the seasons resembling a warm 2026 cluster late in the century; what
    is left once the forced part is gone is the internally generated part, which
    is what can carry decadal persistence.

      large ensemble    minus the ensemble mean, year by year
      control run       minus its own mean (there is no forcing to remove)
      anything else     minus a per-member polynomial in time: quadratic for
                        records of QUADRATIC_MIN_YEARS or more, linear otherwise
    """
    mdim = member_dim(r)
    if mdim and int(r.notnull().any("year").sum()) >= ENSEMBLE_MEAN_MIN_MEMBERS:
        return r - r.mean(mdim)

    def detrend(x: np.ndarray, years: np.ndarray) -> np.ndarray:
        ok = np.isfinite(x)
        out = np.full_like(x, np.nan)
        if ok.sum() < 3:
            return out
        if control:
            deg = 0
        else:
            deg = 2 if years[ok].max() - years[ok].min() + 1 >= QUADRATIC_MIN_YEARS else 1
        out[ok] = x[ok] - np.polyval(np.polyfit(years[ok], x[ok], deg), years[ok])
        return out

    return xr.apply_ufunc(detrend, r, r["year"], input_core_dims=[["year"], ["year"]],
                          output_core_dims=[["year"]], vectorize=True).transpose(*r.dims)


def event_onsets(x: np.ndarray, threshold: float) -> np.ndarray:
    """
    Indices where `x` reaches `threshold` from below.

    Only the first season of a run above the threshold counts, so a spell of
    several high seasons is one event, aligned on when it began. A run already
    above the threshold at the start of the record has no observed onset and is
    skipped.
    """
    above = np.asarray(x) >= threshold
    prev_below = np.zeros_like(above)
    prev_below[1:] = np.isfinite(x[:-1]) & ~above[:-1]
    return np.flatnonzero(above & prev_below)


def epoch_windows(x: np.ndarray, onsets: np.ndarray, lag: int) -> np.ndarray:
    """(n_events, 2*lag+1) windows of `x` centred on each onset, NaN past the record ends."""
    x = np.asarray(x, dtype=float)
    padded = np.concatenate([np.full(lag, np.nan), x, np.full(lag, np.nan)])
    idx = np.asarray(onsets)[:, None] + np.arange(2 * lag + 1)[None, :]
    return padded[idx] if len(onsets) else np.empty((0, 2 * lag + 1))


def running_mean(x: np.ndarray, width: int) -> np.ndarray:
    """Centred running mean; NaN wherever the window is not complete."""
    x = np.asarray(x, dtype=float)
    out = np.full(len(x), np.nan)
    if width <= 1:
        return x.copy()
    if len(x) >= width:
        out[width // 2: len(x) - width // 2] = np.convolve(x, np.ones(width) / width, "valid")
    return out


def random_onset_band(series: list[np.ndarray], n_events: int, lag: int, *,
                      n_draws: int = 2000, quantiles=(0.05, 0.95),
                      seed: int = 0) -> np.ndarray:
    """
    Quantiles, per lag, of the composite mean of `n_events` onsets drawn at random.

    The null for a composite: what averaging this many windows of the same
    series gives when nothing singles out the centre season. Onsets are drawn
    uniformly over every valid season of every series (member), so a long
    ensemble and a single short run are each judged against their own record.
    Returns an array of shape (len(quantiles), 2*lag+1).
    """
    rng = np.random.default_rng(seed)
    pool = [(i, j) for i, s in enumerate(series) for j in np.flatnonzero(np.isfinite(s))]
    pool = np.array(pool)
    means = np.empty((n_draws, 2 * lag + 1))
    for d in range(n_draws):
        pick = pool[rng.integers(len(pool), size=n_events)]
        w = np.concatenate([epoch_windows(series[i], pick[pick[:, 0] == i, 1], lag)
                            for i in np.unique(pick[:, 0])])
        with np.errstate(invalid="ignore"):
            means[d] = np.nanmean(w, axis=0)
    return np.nanquantile(means, quantiles, axis=0)



def decluster(idx: np.ndarray, x: np.ndarray, min_separation: int) -> np.ndarray:
    """
    Thin event indices `idx` of series `x` so that kept events are at least
    `min_separation` apart; where two are closer, the one with the higher
    value wins. Keeps one warm spell from entering a composite several times.
    """
    kept: list[int] = []
    for i in sorted(idx, key=lambda j: x[j], reverse=True):
        if all(abs(i - k) >= min_separation for k in kept):
            kept.append(int(i))
    return np.array(sorted(kept), dtype=int)
