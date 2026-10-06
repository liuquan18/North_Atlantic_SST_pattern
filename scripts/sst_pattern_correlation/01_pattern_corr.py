"""
Score every simulation against the observed JJA 2026 North Atlantic +
Mediterranean SST pattern.

Inputs  : data/pattern_2026/<ds>_jja_anom_na025.nc and <ds>_jja_gmsst.nc
Outputs : data/pattern_2026/results/
            corr_<variant>.nc        correlation time series per dataset
            reference_<variant>.nc   the ERA5 2026 target pattern (0.25 and 1 deg)
            best_analogues.csv       the top-scoring season of each simulation
            ocean_mask_1deg.nc       the common domain all scores are taken over
            summary.json             headline numbers quoted in the figures

Run after the pre_process/ stages have produced their outputs.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, "/work/mh0033/m300883/North_Atlantic_SST_pattern")
import src.pattern_2026 as p2


def main():
    p2.RESULT_DIR.mkdir(parents=True, exist_ok=True)

    available = [k for k, d in p2.DATASETS.items()
                 if d.anom_file.exists() and d.gm_file.exists()]
    missing = [k for k in p2.DATASETS if k not in available]
    if missing:
        print(f"!! not on disk yet, skipping: {missing}")
    print(f"datasets: {available}")

    anom = {k: p2.load_anomaly(k) for k in available}
    gm = {k: p2.load_gmsst(k) for k in available}
    for k in available:
        print(f"  {k:16s} {dict(anom[k].sizes)}  years "
              f"{int(anom[k].year.min())}-{int(anom[k].year.max())}")

    # ---- common 1 deg analysis grid and shared ocean domain ----------------
    coarse = {k: p2.coarsen_to_1deg(v) for k, v in anom.items()}
    mask = p2.common_ocean_mask(coarse)
    ncell = int(mask.sum())
    print(f"\ncommon ocean cells on the 1 deg grid: {ncell} "
          f"of {mask.size} ({100*ncell/mask.size:.0f}%)")
    mask.to_netcdf(p2.RESULT_DIR / "ocean_mask_1deg.nc")

    if p2.REF_YEAR not in anom["ERA5"].year.values:
        raise SystemExit(f"ERA5 has no JJA {p2.REF_YEAR}")

    summary = {"reference_year": p2.REF_YEAR,
               "region": p2.REGION,
               "climatology": list(p2.CLIM_PERIOD),
               "n_common_ocean_cells_1deg": ncell,
               "datasets": {}}
    rows = []

    for variant, vinfo in p2.VARIANTS.items():
        print(f"\n===== variant: {variant} ({vinfo['long']}) =====")
        centered = vinfo["centered"]

        # patterns on both grids: 1 deg for scoring, 0.25 deg for the maps
        pat_fine = {k: p2.make_pattern(anom[k], gm[k], variant) for k in available}
        pat_1deg = {k: p2.make_pattern(coarse[k], gm[k], variant).where(mask)
                    for k in available}

        ref_1deg = pat_1deg["ERA5"].sel(year=p2.REF_YEAR).compute()
        ref_fine = pat_fine["ERA5"].sel(year=p2.REF_YEAR)
        xr.Dataset({"ref_1deg": ref_1deg, "ref_025": ref_fine}).to_netcdf(
            p2.RESULT_DIR / f"reference_{variant}.nc")

        corrs = {}
        for k in available:
            # .compute() once: everything below (argmax, count, to_netcdf) would
            # otherwise re-walk the dask graph, and for MPI-GE that means
            # re-coarsening 2.9 GB from 0.25 deg on every single access.
            r = p2.pattern_corr(pat_1deg[k], ref_1deg, centered=centered).compute()
            corrs[k] = r.rename(k)

            # best analogue = the single simulated season closest to the
            # observed 2026 pattern (over members too, where there are any)
            if "member" in r.dims:
                stacked = r.stack(z=("member", "year"))
                zbest = stacked.isel(z=int(stacked.fillna(-9).argmax("z")))
                best_year, best_member, best_r = (
                    int(zbest["year"]), int(zbest["member"]), float(zbest))
            else:
                ybest = r.isel(year=int(r.fillna(-9).argmax("year")))
                best_year, best_member, best_r = int(ybest["year"]), None, float(ybest)

            # for the observations, the reference year itself is trivially r=1
            note = ""
            if k == "ERA5":
                excl = r.where(r["year"] != p2.REF_YEAR)
                ybest = excl.isel(year=int(excl.fillna(-9).argmax("year")))
                best_year, best_r = int(ybest["year"]), float(ybest)
                note = f"excluding {p2.REF_YEAR} itself"

            rows.append(dict(variant=variant, dataset=k, label=p2.DATASETS[k].label,
                             best_year=best_year, best_member=best_member,
                             best_r=round(best_r, 4),
                             n_seasons=int(r.count()), note=note))
            mem = f", member r{best_member}" if best_member else ""
            print(f"  {k:16s} best analogue: JJA {best_year}{mem}  r = {best_r:+.3f}"
                  f"   ({int(r.count())} seasons scored) {note}")

        # The models are free-running, so their own "2026" is an arbitrary draw
        # rather than a forecast -- but it is the literal comparison, so record it.
        for k in available:
            if k == "ERA5" or p2.REF_YEAR not in corrs[k].year.values:
                continue
            r26 = corrs[k].sel(year=p2.REF_YEAR)
            entry = {"variant": variant, "r_max": round(float(r26.max()), 4)}
            if "member" in r26.dims:
                entry.update(r_mean=round(float(r26.mean()), 4),
                             r_min=round(float(r26.min()), 4),
                             best_member=int(r26["member"][int(r26.argmax())]))
            summary.setdefault("nominal_2026", {}).setdefault(k, []).append(entry)

        # each ensemble gets its own member dimension: a shared "member" would
        # align MPI-GE's 50 and MPI-ESM1.2-ER's 3 onto one index and pad the
        # smaller with NaN (see src.pattern_2026.member_dim)
        ds_out = xr.Dataset({
            k: (v.rename({"member": f"member_{k}"}) if "member" in v.dims else v)
            for k, v in corrs.items()
        })
        ds_out.attrs.update(variant=variant, description=vinfo["long"],
                            reference=f"ERA5 JJA {p2.REF_YEAR}",
                            region=str(p2.REGION))
        ds_out.to_netcdf(p2.RESULT_DIR / f"corr_{variant}.nc")

    # ---- like-for-like best analogues --------------------------------------
    # The raw "best season" comparison is unfair: MPI-GE gets 50 members x 251
    # years = 12,550 draws at a good analogue, while the km-scale runs get 35.
    # The best of a large sample is higher than the best of a small one even if
    # the underlying distributions are identical. So also report, for every
    # dataset, the distribution of "best of n" where n is the shortest record.
    #
    # NB: resampling is with replacement, so for a dataset whose record is
    # already only n seasons long this is degenerate -- bounded by its own
    # maximum, and not an independent estimate. Read the output as a reference
    # distribution built from the *large*-sample datasets, against which the
    # short records' actual maxima are compared.
    rng = np.random.default_rng(0)
    n_draw = min(int(xr.open_dataset(p2.RESULT_DIR / "corr_spatial.nc")[k].count())
                 for k in available if k != "ERA5")
    summary["fair_sample_size"] = n_draw
    print(f"\n===== best-of-{n_draw} (equalised sample size, 2000 bootstrap draws) =====")
    fair_rows = []
    for variant in p2.VARIANTS:
        ds_corr = xr.open_dataset(p2.RESULT_DIR / f"corr_{variant}.nc")
        for k in available:
            vals = ds_corr[k].values.ravel()
            vals = vals[np.isfinite(vals)]
            if k == "ERA5" or vals.size < n_draw:
                continue
            draws = rng.choice(vals, size=(2000, n_draw), replace=True)
            maxima = draws.max(axis=1)
            lo, med, hi = np.percentile(maxima, [5, 50, 95])
            fair_rows.append(dict(variant=variant, dataset=k, n_draw=n_draw,
                                  best_of_n_median=round(float(med), 4),
                                  best_of_n_p5=round(float(lo), 4),
                                  best_of_n_p95=round(float(hi), 4),
                                  n_available=int(vals.size)))
            print(f"  {variant:8s} {k:16s} best-of-{n_draw}: "
                  f"r = {med:+.3f}  (5–95%: {lo:+.3f} to {hi:+.3f})")
    pd.DataFrame(fair_rows).to_csv(p2.RESULT_DIR / "best_of_n_bootstrap.csv", index=False)

    # ---- context: how warm and how unusual was the observed 2026 box --------
    w = np.cos(np.deg2rad(anom["ERA5"]["lat"])) * xr.ones_like(anom["ERA5"]["lon"])
    w = w.where(anom["ERA5"].isel(year=0).notnull())
    box_mean = (anom["ERA5"] * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))
    box_mean.rename("box_mean_anom").to_netcdf(p2.RESULT_DIR / "era5_box_mean_anom.nc")
    v2026 = float(box_mean.sel(year=p2.REF_YEAR))
    rank = int((box_mean >= box_mean.sel(year=p2.REF_YEAR)).sum())
    summary["era5_box_mean_anom_2026_degC"] = round(v2026, 3)
    summary["era5_box_mean_rank"] = rank
    summary["era5_n_years"] = int(box_mean.count())
    summary["era5_gmsst_2026_degC"] = round(float(gm["ERA5"].sel(year=p2.REF_YEAR)), 3)
    print(f"\nERA5 JJA {p2.REF_YEAR} box-mean anomaly: {v2026:+.2f} degC "
          f"(rank {rank} of {int(box_mean.count())} since {int(box_mean.year.min())})")
    print(f"ERA5 JJA {p2.REF_YEAR} global-ocean-mean anomaly: "
          f"{summary['era5_gmsst_2026_degC']:+.2f} degC")

    for k in available:
        d = p2.DATASETS[k]
        summary["datasets"][k] = dict(label=d.label, resolution=d.resolution,
                                      years=[int(anom[k].year.min()), int(anom[k].year.max())],
                                      n_members=int(anom[k].sizes.get("member", 1)),
                                      note=d.note)

    df = pd.DataFrame(rows)
    df.to_csv(p2.RESULT_DIR / "best_analogues.csv", index=False)
    (p2.RESULT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote results to {p2.RESULT_DIR}")


if __name__ == "__main__":
    main()
