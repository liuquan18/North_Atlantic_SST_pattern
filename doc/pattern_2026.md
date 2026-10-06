# JJA 2026 North Atlantic + Mediterranean SST pattern

Pattern-correlation study of the **observed JJA 2026** sea surface temperature
pattern over **30–60°N, 80°W–40°E** against four simulation sources, including
two km-scale ICON configurations that are new to this repository.

It succeeds an earlier JJA 2023 analysis (MPI-GE only), whose scripts were
removed from this branch; they remain in the git history on `main`.

## What is different from the 2023 analysis

| | 2023 pipeline | this pipeline |
|---|---|---|
| reference season | JJA 2023 | JJA 2026 |
| region | 0–70°N, 80°W–0° | 30–60°N, 80°W–40°E (adds the Mediterranean) |
| variable | `Amon/ts` then `Omon/tos` | `tos` / `to` only — never `ts`/`tas` |
| simulations | MPI-GE | MPI-GE, km-scale ICON (EPOC ×2), EERIE ICON-ESM-ER (×2, forced run 3 members) |
| scoring grid | ERA5's grid | common 1° grid, common ocean mask across all datasets |

## Running it

Script layout and run order: [scripts/README.md](../scripts/README.md).

```bash
cd /work/mh0033/m300883/North_Atlantic_SST_pattern
sbatch scripts/pre_process/01_era5.sh          # observations, 1940-2026
sbatch scripts/pre_process/02_mpige.sh         # MPI-ESM1-2-LR, 50 members
sbatch scripts/pre_process/03_epoc_icon.sh     # km-scale ICON (EPOC), 2 experiments
sbatch scripts/pre_process/04_eerie.sh         # EERIE ICON-ESM-ER, 2 experiments (forced: 3 members)
sbatch scripts/pre_process/05_mpi_er.sh        # MPI-ESM1.2-ER, 3 realizations
sbatch scripts/pre_process/06_kmscale_zoom.sh  # zoom data for figure 5 (needs 03)

# then the analysis and all six figures (handles the conda env itself):
bash scripts/run_analysis.sh
```

Stages 01–05 are independent and idempotent — each skips work whose output is
already on disk, so a partial run can be resubmitted as is.

## The common interface

Every cdo stage reduces its dataset to exactly two files in
`data/pattern_2026/`, which is all the Python knows about:

- `<ds>_jja_anom_na025.nc` — JJA-mean SST anomaly, one step per year (plus a
  `member` dimension for MPI-GE), on the common 0.25° regional grid
- `<ds>_jja_gmsst.nc` — global-ocean-mean JJA SST anomaly, one value per year

Anomalies are always against the dataset's **own** 1991–2020 JJA monthly
climatology, computed on its **native** grid; remapping happens once, on the
small anomaly field, never on raw monthly data.

## What the analysis found

Numbers below are from `data/pattern_2026/results/` (run of 2026-09-22); figures
in `figures/pattern_2026/`.

**The observed reference.** JJA 2026 over the box is +0.69 °C, the 3rd warmest
summer since 1940. The global ocean mean that summer is +0.56 °C, so only about
0.13 °C of that is regional rather than global warming — which is why the two
pattern variants give almost identical correlations throughout. The pattern
itself is a cold central/subpolar North Atlantic against a warm eastern
boundary and a very warm Mediterranean (locally beyond +3 °C), and it is present
in all three months separately, not carried by one (figure 4).

**Closest observed analogue: JJA 2003**, r = +0.63 — the European heatwave
summer. This is used as the benchmark everywhere rather than an arbitrary
threshold.

**Best analogue each simulation can produce** (spatial-mean-removed variant):

| simulation | ocean | seasons available | best season | r |
|---|---|---|---|---|
| MPI-ESM1-2-LR grand ensemble | ~1° | 12,550 | JJA 2016, member r24 | **+0.67** |
| MPI-ESM1.2-ER | ~10 km | 450 | JJA 2038, member r2 | +0.60 |
| km-scale ICON (EPOC, transient GHG) | 5 km | 36 | JJA 2002 | +0.59 |
| km-scale ICON (EPOC, control) | 5 km | 35 | JJA 2005 | +0.53 |
| EERIE ICON-ESM-ER (hist+ssp245) | 5 km, archived 0.25° | 187 (3 members) | JJA 2039, member r1 | +0.50 |
| EERIE ICON-ESM-ER (control) | 5 km, archived 0.25° | 101 | JJA 2041 | +0.61 |

**The raw ranking is an artefact of sample size.** MPI-GE gets 12,550 draws at a
good analogue; MPI-ESM1.2-ER gets 450 and the km-scale ICON runs get 36 (transient) and 35 (control).
Equalising to 35 draws (`best_of_n_bootstrap.csv`):

| | best-of-35 | 5–95% |
|---|---|---|
| MPI-ESM1-2-LR (~1° ocean) | +0.39 | +0.27 to +0.53 |
| MPI-ESM1.2-ER (~10 km ocean) | +0.46 | +0.31 to +0.60 |
| km-scale ICON (EPOC, 5 km) | +0.59 (its actual best of 36) | — |

**Ocean resolution helps, and the cleanest evidence is the MPI pair.** MPI-GE
and MPI-ESM1.2-ER are the same model family differing mainly in ocean
resolution (~1° vs ~10 km), and the finer one reaches +0.46 against +0.39 on
equal sampling. The km-scale ICON run's +0.59 out of its own 36 seasons sits
above MPI-GE's 95th percentile. The ordering is not monotonic in resolution
overall — EERIE is lowest at +0.38 — but EERIE is only archived at 0.25°, so
what is scored there is not its 5 km field.

*Caveat on `best_of_n_bootstrap.csv`:* resampling is with replacement, so for a
dataset whose record is already only 35–36 seasons the "best-of-35" estimate is
degenerate (bounded by its own maximum) and is not an independent estimate. Use
the bootstrap as a **reference distribution from the large-sample datasets**
(MPI-GE, and EERIE at n=187) and compare the short records' actual maxima
against it, as above.

**A 2026-like pattern is rare everywhere.** Measured against the observed 2003
analogue, only **1 of 12,550** MPI-GE seasons reaches it — and **0 of 450**
(MPI-ESM1.2-ER), **0 of 36** (EPOC) and **0 of 187** (EERIE) (figure 3).

**It looks like internal variability, not a forced signal.** The constant-forcing
control runs do as well as their forced counterparts — EERIE control +0.61 vs
forced +0.50, EPOC control +0.53 vs transient +0.59 — so nothing in the
greenhouse-gas forcing is generating this pattern. Consistently, each model's
own nominal JJA 2026 is unremarkable (MPI-GE best member +0.54, 50-member mean
+0.13; EERIE r1 +0.21 — r2/r3 end before 2026), confirming that in free-running simulations the calendar
year carries no forecast meaning.

**Figures show one pattern variant.** Removing the box mean and removing the
global mean give near-identical correlations (differences of <0.01 in r
throughout), because JJA 2026's basin-mean warmth is almost all the global
signal. Figures 1–3 therefore show only the box-mean-removed variant; both are
still computed and stored.

**Regional means (figure 6).** Area-mean JJA anomalies for four regions, where
`box = natl + med` exactly. Observed JJA 2026: global ocean **+0.56 °C**,
analysis box **+0.66**, North Atlantic **+0.44**, Mediterranean **+1.64**. The
box figure is pulled up almost entirely by the Mediterranean, which is only
~19% of the box area by weight. The Mediterranean anomaly is far outside its own
1991–2020 range (95th percentile +0.74 °C).

The Mediterranean also warms roughly twice as fast as the global ocean in the
observations, and the models disagree about by how much (K/decade over each
record):

| | global | box | North Atlantic | Mediterranean |
|---|---|---|---|---|
| ERA5 (1940–2026) | +0.085 | +0.123 | +0.111 | **+0.174** |
| MPI-GE (1850–2100) | +0.074 | +0.087 | +0.083 | +0.104 |
| MPI-ESM1.2-ER (1950–2099, ssp585) | +0.194 | +0.233 | +0.218 | +0.301 |
| EERIE (1950–2050) | +0.125 | +0.192 | +0.180 | +0.244 |
| km-scale ICON (1990–2025) | +0.133 | +0.346 | +0.323 | +0.446 |

Absolute trends are not comparable across these rows — the records cover
different periods and scenarios (ssp245 for MPI-GE and EERIE, ssp585 for
MPI-ESM1.2-ER and EPOC). The comparable quantity is the **Mediterranean
amplification**, the ratio of the Mediterranean trend to the global one:
observed **2.05**, MPI-GE 1.40, MPI-ESM1.2-ER 1.55, EERIE 1.95, EPOC 3.36 (over
36 years only). Every model but EPOC under-does the observed amplification, and
here too the finer MPI ocean is closer than the coarse one.

MPI-ESM1-2-LR warms the Mediterranean at only ~1.4× its global rate against
~2.0× observed, which is consistent with the cold Mediterranean bias visible in
figure 5. The km-scale and EERIE trends are over much shorter records that
include the strong recent warming, so they are not directly comparable to the
1940–2026 observed trend — treat the last two rows as record-specific, not as
climate sensitivities. Both control runs are flat (≤0.016 K/decade globally),
confirming these are forced trends.

**Resolution caveat (figure 5).** The correlations are computed at 1°, which is
blind to what a 5 km ocean resolves. At native resolution MPI-ESM1-2-LR's ~1°
ocean is visibly blocky in the Gulf Stream extension and barely resolves the
Adriatic and Aegean, and its Mediterranean is several °C too cold in JJA. Note
also that EERIE's 5 km ocean is only *archived* at 0.25°, so its km-scale
structure is not available here — EPOC is the only genuinely km-scale field in
this study.

## Data sources and the traps in each

**ERA5** `/pool/data/ERA5/E5/sf/an/1M/034` (var34, SSTK), 1940–2026.
The monthly stream currently ends at **2026-06**, so JJA 2026 cannot be built
from it. July and August 2026 come from the near-real-time **ERA5T** hourly
stream `/pool/data/ERA5/ET/sf/an/1H/034` (which runs to 2026-09-16). For June
2026, where both exist, the hourly-derived monthly mean matches the archived
monthly mean to 7×10⁻⁴ K, so the whole of 2026 is taken from ERA5T for
within-season consistency.

**MPI-ESM1-2-LR grand ensemble** `Omon/tos`, historical + ssp245, r1–r50,
1850–2100. All 50 members share one tripolar ocean grid, so remap weights are
generated once and reused — the whole 50-member stage runs in well under a
minute on one node.

**km-scale ICON (EPOC)** 10 km atmosphere / 5 km ocean.
`/work/bm1313/b383127/epoc-icon-2024.10{,_aerosols}/experiments/epoc2_0{10,20}/work/run_*`
Three things differ from a CMIP archive and each will silently corrupt results:

1. The field is `to` in the `oce_2d_1mth_mean` stream, on the **unstructured**
   ICON R02B09 ocean grid (14,914,033 cells). The grid is not in the data files
   and must be attached from
   `/pool/data/ICON/grids/public/mpim/0045/icon_grid_0045_R02B09_O.nc`
   (uuid `b7ad8f68-8dde-11ee-b8df-3133d6395582`).
2. **Land cells are written as exactly 0.0 with no `_FillValue`.** Left alone
   they bleed into every coastal target cell during remapping, which would
   wreck the Mediterranean — half of this study's region. Masking them moves
   the regional minimum from 0.0 °C to 3.3 °C in a test month. The zero set is
   bit-identical between 1990 and 2020, so a static mask built once is exact,
   and it keeps the remap linear (which is why anomalies can be computed on the
   small remapped fields rather than on 6 GB of native intermediates).
3. Monthly files are stamped with the **end** of their averaging period:
   `..._19900801T000000Z.nc` holds the **July 1990** mean.

Forcing, read off the run scripts rather than the pad (they differ):
`epoc2_010` has `irad_co2=2` with Kinne aerosols and ozone pinned at 1990 — a
**constant-1990-forcing control**, not a historical run. `epoc2_020` has
`irad_co2=3` with `greenhouse_ssp585.nc` and `bc_ozone_ssp585_<year>` — a
**transient historical + SSP5-8.5 greenhouse-gas run**, but with Kinne aerosols
held at 1850, which is what the `_aerosols` build name refers to.

Coverage: `epoc2_010` 1990–2024; `epoc2_020` is still running (at 2026-03 as
of 2026-10-05), so its last complete JJA is 2025. Neither has JJA 2026 yet.
From 2025-01 on, `epoc2_020` files are written uncompressed (~2.5 GB instead
of ~0.85 GB); variables, names and the three traps above are unchanged.

**MPI-ESM1.2-ER** T127 atmosphere (~100 km) + the eddy-resolving MPIOM TP6M
ocean (3602×2394 curvilinear, ~0.1° / ~10 km) — the same model family as MPI-GE
but with a ~10× finer ocean, so the pair isolates ocean resolution within one
model. Three realizations under
`/work/uo0122/u241089/MPIESM/{ER,ER3,ER5}-{hist,ssp585}/outdata/mpiom/`
(ER = member 1, ER3 = 2, ER5 = 3), historical 1950–2014 + ssp585 2015–2099, one
~2.2 GB year-file each. Traps are mild next to EPOC: the variable is lowercase
**`tos`** (the readme says `TOS`, and `-selname,TOS` aborts), and the real cost
is decompression — ~10 s to pull JJA out of a year-file, so the extracted field
is written once and reused for both the remap and the global mean. Monthly means
are stamped at the *end* of their month, which for `-selmon` is already correct.
`/work/mh1421/data/mpiesm/` holds only the hist legs, so use the `uo0122` path.

**EERIE ICON-ESM-ER** 10 km atmosphere / 5 km ocean, CMORized `Omon/tos` on a
0.25° regular `gr` grid — no unstructured handling needed.
`hist-1950` (1950–2014) + `highres-future-ssp245` (2015–2050) under
`/work/bm1344/DKRZ/CMOR/EERIE/HighResMIP/MPI-M/ICON-ESM-ER/`, and
`eerie-control-1950` (1950–2050) under `/pool/data/EERIE/EERIE/MPI-M/`.
Only `r1i1p1f1` is in the CMOR tree; members r2 (erc2023, SST 1975–2014) and r3 (erc2024, 1975–2020) are read from the raw 0.25° output under `/work/bm1344/k202193/ICON/` (see `04_eerie.sh`), r2 against a 1991–2014 climatology because its SST stops in 2014. Note the archive is stored at
0.25°, so the 5 km ocean structure is **not** available here — figure 5 makes
that explicit.

## Methodology notes

**Two pattern definitions** are computed and stored, matching the 2023 work.
They turn out to give near-identical correlations here, so the figures draw only
the first; `results/corr_global.nc` and `best_analogues.csv` keep both, and
`src.pattern_2026.PLOT_VARIANTS` controls what is drawn.


- *spatial mean removed* — subtract the area-weighted mean over the box from
  both fields, then correlate (centred Pearson correlation over space). A pure
  shape comparison, blind to how warm the basin is overall.
- *global mean removed* — subtract each dataset's own global-ocean-mean JJA SST
  anomaly, then take the uncentred (cosine) correlation. Basin-mean warmth that
  survives the removal of the global warming signal still counts as pattern.

**Scoring grid.** Correlations are computed after coarsening the common 0.25°
grid to 1°. MPI-ESM1-2-LR's ocean grid is ~1°, and scoring it on 0.25°
structure it cannot represent would flatter the km-scale runs for the wrong
reason. A 1° cell enters the correlation only if it is ≥50% valid ocean.

**Common ocean mask.** Land masks differ between ERA5, the ~1° MPI ocean grid
and the ICON ocean grids, especially around the Mediterranean coast. The mask
used is the **intersection** over all datasets, so every correlation is taken
over exactly the same set of cells.

**Best analogue.** Models are free-running, so their calendar year carries no
forecast meaning — "JJA 2026" in ssp245 is an arbitrary draw. Each model is
therefore shown at its best-matching season over every year and member it
simulates, chosen independently under each of the two variants.
