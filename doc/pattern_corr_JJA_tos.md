# North Atlantic SST Pattern 2023 JJA — tos-corrected pattern correlation

Supersedes the analysis in `pattern_corr_JJA.pptx`, which used `Amon/ts` (surface
temperature, land+ocean) masked by land fraction as a stand-in for SST. This
version uses `Omon/tos` — the model's actual sea surface temperature — and
also cuts the compute footprint significantly. See [CLAUDE.md](../CLAUDE.md)
for the data paths and environment notes referenced below.

## What changed and why

**Variable.** `Amon/ts` sits on the atmosphere's regular Gaussian grid and
needs an `sftlf` land-fraction mask to approximate "ocean only." `Omon/tos` is
the actual SST diagnostic, already ocean-masked (missing over land), but it
lives on MPI-ESM1-2-LR's native curvilinear tripolar ocean grid (256×220
points) rather than a regular lat/lon grid. That grid has to be remapped
before it's comparable to ERA5's regular grid at all — this was silently
absent from the old pipeline because `ts` happened to already sit on the
target grid.

**Grid handling.** Every one of the 50 ensemble members shares the identical
native ocean grid, so the bilinear remap weights (`cdo genbil`) are generated
once and reused (`cdo remap`) for every subsequent file, instead of being
regenerated per file. Climatology is kept on the native grid for the
`ymonsub` step and only the (much smaller) resulting anomaly is remapped —
remapping happens on ~250 timesteps of JJA-only data per member, not on every
raw monthly field.

**Compute.** The old correlation job requested 4 nodes / 40 MPI ranks for a
workload that is I/O-bound, not CPU-bound, and looped year-by-year in Python.
`spatial_corr_model_era` now calls `xr.corr` once per ensemble over the whole
time series (`xr.corr` broadcasts over dims outside `dim=`, verified to give
identical results to the old loop) instead of looping. Combined with weight
reuse, both stages now run on a single compute node:

| Stage | Before | After | Measured |
|---|---|---|---|
| Preprocessing (climatology + anomaly + global-mean removal, 50 members, 1850–2100) | not on a compute node at all (interactive `parallel`, no `#SBATCH`) | 1 node, `#SBATCH --nodes=1` | 43s |
| Pattern correlation (50 members × both methods) | 4 nodes, 40 MPI ranks, up to 2h requested | 1 node, no MPI, `ProcessPoolExecutor` | 21s (+ a first attempt that surfaced a real alignment bug, see below) |

**A real bug this surfaced.** `calculate_spatial_correlation_keepareamean`
used `xr.align(..., join="exact")`, which requires bit-identical coordinate
labels. ERA5 and the model tos anomaly are both remapped onto the same
physical grid but by independent `cdo` pipelines, so their `lat` values agree
to ~1e-8 but not bit-for-bit — `join="exact"` rejected this outright
(`AlignmentError`). Fixed by snapping field1's spatial coordinates onto
field2's once they're confirmed to agree within tolerance
(`src/pattern_correlation.py`). This was a latent risk in the original
pipeline too, just never triggered.

**Scratch data doesn't persist.** All the old `ts`-based scratch intermediates
had already been purged by the time this work started — Levante's scratch
retention policy, not something this project did. The new `tos` intermediates
live under `/scratch/m/m300883/nalt/MPI_GE_CMIP6_tos/` and should be assumed
equally ephemeral; only the correlation results and figures are kept
long-term, under `data/` and `figures/`.

## Method 1: remove spatial mean before correlating

Centered (standard Pearson) pattern correlation between each ensemble
member's JJA `tos` anomaly and ERA5's observed 2023 JJA anomaly, both
area-weighted by cos(latitude), computed globally (not restricted to the
North Atlantic box).

![ensemble-mean correlation, both methods](../figures/pattern_corr_avg_ens_JJA_tos.png)

The ensemble-mean correlation (blue curve above) rises from around −0.1 in
the 19th century to around +0.2 by 2100. This tracks the global-warming
fingerprint: as forced warming comes to dominate the global SST field, every
ensemble member's pattern increasingly resembles 2023's (which is itself
warming-dominated), independent of any North-Atlantic-specific process.

![violin and threshold count, method 1](../figures/violin_count_pattern_corr_JJA_tos.png)

Only 8 (member, year) combinations across all 50 members and 251 years reach
r ≥ 0.4, and they cluster in the 2040s–2080s — consistent with the
warming-fingerprint explanation above rather than a distinctive North
Atlantic pattern.

![composite patterns, method 1](../figures/high_corr_patterns_JJA_tos.png)

The composite of those 8 events (left) is a near-global warm pattern, not
something specific to the North Atlantic — this method is dominated by the
forced trend, which is exactly why method 2 exists.

## Method 2: remove global mean before correlating (extended back to 1850)

Uncentered pattern correlation, computed only within the North Atlantic box
(280–360°E, 0–70°N) after subtracting each field's own global mean. This removes the shared forced-warming signal and
isolates whether the *regional* pattern — independent of how much the planet
has warmed overall — resembles 2023.

![violin and threshold count, method 2](../figures/violin_count_pattern_corr_JJA_tos_noglm.png)

No secular trend here (as expected — the global mean, and with it most of
the forced trend, has been removed). 494 of 12,550 (member, year) pairs
(~3.9%) reach r ≥ 0.4 across the full 1850–2100 record, so a pattern this
strong is not a rare event under internal variability alone.

![composite patterns, method 2](../figures/high_corr_patterns_JJA_tos_noglm.png)

Both the composite of all 494 events (left) and the single strongest match
(right — ensemble 18, **1917**, r = 0.72) reproduce the same qualitative
structure: a warm horseshoe around the subpolar North Atlantic with a cool
patch in the central subpolar gyre — the "North Atlantic warming hole"
pattern that's the actual object of interest for this project. That the
strongest match comes from 1917 — decades before any detectable Atlantic
forced-warming signal — supports reading 2023's pattern as consistent with
internal variability, not necessarily a forced response.

## Where things are

- Preprocessing: `scripts/MPI_GE_pattern/anomaly_tos.sh` (submit with `sbatch`)
- Correlation: `scripts/pattern_corr/MPI_ERA5_patterncorr_tos.py`, run via `scripts/pattern_corr/patterncorr_tos_run.sh`
- Plots: `scripts/pattern_corr/plots_pattern_corr_tos.py`
- Results: `data/pattern_corr_JJA_tos/`, `data/pattern_corr_JJA_tos_noglm/`
- Figures: `figures/pattern_corr_avg_ens_JJA_tos.png`, `figures/violin_count_pattern_corr_JJA_tos*.png`, `figures/high_corr_patterns_JJA_tos*.png`

The original `ts`+land-mask scripts (`scripts/MPI_GE_pattern/anomaly.sh`,
`scripts/pattern_corr/MPI_ERA5_patterncorr*.py`, the `_norm` plotting variant)
are left in place for provenance but are superseded by the above.
