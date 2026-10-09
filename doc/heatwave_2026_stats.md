# European summer heatwaves 2003 and 2026: how extreme, how rare?

Statistical analysis of the two strongest observed southern-European summers
in heatwave terms, against the observed record (ERA5 1940–2026) and the
simulations: the MPI-ESM1.2-LR grand ensemble (MPI-GE, 50 members) and the
km-scale ICON runs, both transient (EPOC hist+ssp585, EERIE hist+ssp245) and
constant-forcing controls (EERIE, Sapphire sap0006, EPOC). Code:
`scripts/heatwave_statistics/`, `src/extremes.py`; tables:
`data/heatwave_stats/`; figures 13–16 in `figures/pattern_2026/`.

## Index

Land-mean JJA **heatwave days** (`hwd`) and **cumulative heat** (`hwcum`,
°C·day) over southern Europe, 35–60°N, 10°W–30°E (`src/heatwave.py`
`MEAN_BOX`), heatwaves after Xu et al. (2026): ≥3 days with Tmax above the
calendar-day 90th percentile of the dataset's **own** 1991–2020 climate.

| ERA5 | heatwave days | cumulative heat |
|---|---|---|
| JJA 2003 | 10.4 | 73 °C·day |
| JJA 2026 | 20.1 | 155 °C·day |
| 1991–2020 mean ± std | 4.9 ± 2.7 | 34 ± 20 |
| 1950–1989 mean, detrended std | 2.0, 1.1 | |

2003 was matched or beaten in 2015 and in every summer 2021–2026.

ERA5 heatwave metrics were extended from 1991 back to 1940 for this analysis
(1991–2026 values bit-identical to before). Fits use 1950–2026; 1940–49 are a
sensitivity test.

## Method

**Warming covariate.** T = JJA global-mean SST anomaly of the same dataset
(1991–2020 reference), 4-year trailing mean (MPI-GE: 50-member mean, i.e. the
forced signal). ERA5: T = −0.38 K in the 1950s, −0.06 K in 2003, +0.47 K in
2026 (the 2023–26 SST surge).

**Nonstationary GEV** (`src/extremes.py`), maximum likelihood, |ξ| ≤ 0.5:

- *shift*: μ = μ₀ + μ₁T, σ constant. This is the usual choice for temperature;
- *shift+scale*: μ = μ₀ + μ₁T, log σ = s₀ + s₁T.

Every season gets its return period and normal-equivalent z **in its own
climate**. 95% intervals come from 1000 bootstrap refits (300 for the
sensitivity variants).

**Why the scale term is needed.** With the threshold fixed, heatwave days
spread out as the mean rises. In ERA5 the detrended std is 1.1 days in
1950–89 and 3.0 in 2007–26; in MPI-GE it is 2.2 → 4.3 days. A shift-only fit
can absorb this only through a heavy tail (ξ ≈ +0.3), and it then leaves the
residuals heteroscedastic. The diagnostics, for 1950–2026 heatwave days:

| fit | AIC | ξ | z-spread, 1st / 2nd half | Levene p |
|---|---|---|---|---|
| shift | 331 | +0.32 | 0.64 / 1.22 | < 0.001 |
| shift+scale | 300 | +0.10 | 1.04 / 0.97 | 0.90 |
| shift+scale, 20-yr covariate | 292 | +0.03 | 1.02 / 1.00 | 0.92 |

Shift+scale is therefore the main model. The shift fit is shown as a
(misspecified) sensitivity case: it puts the 2026 value at only 1-in-180
years in the 1950s climate, when no 1950–89 summer exceeded 4.5 days.

**Simulations.**

- MPI-GE is scored empirically (no fit). The comparison windows are:
  - 1850–1900;
  - ±10 calendar years around the event;
  - the seasons at ERA5's warming level (forced T within ±0.05 K): 1999–2002
    for ERA5's 2003 level and 2039–2045 for ERA5's 2026 level. MPI-GE warms
    more slowly than the observed 2023–26 SST surge.
- The transient km-scale runs get the same NS-GEV.
- The controls get a stationary GEV (linear drift removed) plus empirical
  counts.

## Results

### 1. Observed record (fig 13)

How exceptional each summer was **for its time**, as rank among 77 seasons,
z, and return period (RP) in its own climate:

| method | 2003 | 2026 |
|---|---|---|
| linear trend in time, 1991–2026 (the original basis) | #2, z 1.8 | #1, z 2.6 |
| linear trend in time, 1950–2026 | #6, z 1.6 | #1, z 4.1 |
| GEV shift, 1950–2026 | #5, 1-in-21 yr | #1, 1-in-88 yr |
| **GEV shift+scale, 1950–2026** | **#3, 1-in-34 yr (18–145)** | **#14, 1-in-7 yr (4–19)** |
| shift+scale, 10 / 20-yr covariate | #2 / #1, 1-in-42 / 70 yr | #10 / #11, 1-in-9 / 8 yr |
| shift+scale, event left out | #3, 1-in-41 yr | #12, 1-in-7.5 yr |
| shift+scale, 1940–2026 | #4, 1-in-26 yr | #13, 1-in-7.5 yr |

Cumulative heat gives the same picture: 2003 is 1-in-32 years; 2026 is
1-in-8 years.

The two years behave differently:

- **2003 is robustly exceptional for its time.** It is in the top 1–6 under
  every method, typically 1-in-20 to 1-in-70 years. The only seasons that
  rival it in their own climate are 1957 and 1976 (and 1952), which were
  modest in absolute terms but stood out against the narrow 1950–80s
  distribution.
- **2026 is #1 only if the forced change is a straight line with constant
  spread.** Once the distribution may widen with warming, 2026 is about a
  1-in-7 to 1-in-9-year summer **in the 2026 climate**. This holds whether or
  not the 2023–26 SST surge counts as climate state (10- and 20-year
  covariates) and whether or not 2026 itself is in the fit. The 2022, 2024
  and 2025 summers (14–15 days) already sat in the same widened
  distribution.

**Attribution-style contrasts** (shift+scale, heatwave days, 95% CI):

| | 2003 value (10.4 d) | 2026 value (20.1 d) |
|---|---|---|
| RP in 2026 climate | 2.2 yr (1.8–3.1) | 6.8 yr (3.8–19) |
| RP in 2003 climate | 34 yr (18–145) | ~1,100 yr (>140) |
| RP in 1950s climate | ~4,000 yr (>280) | >2,600 yr (lower bound) |
| how much weaker an equally rare 1950s summer was | 6.2 days (5.2–7.1) | 17.4 days (15.9–18.3) |

A 2003-sized summer is now an every-other-year event. A summer as rare in the
2026 climate as 2003 was in its own would have about 34 heatwave days
(fig 15a).

### 2. Model variability against ERA5 (fig 14)

Heatwave-day std in 1991–2020: ERA5 2.7 days, EERIE 2.8, MPI-GE 3.7, EPOC
3.9, controls 3.3–4.4. Early on the gap is larger: MPI-GE has 2.2 days at the
cold 1950–69 warming level, against ERA5's 1.1. The models are 1.3–2× too
variable, except EERIE.

**Spread growth with warming** (GEV s₁ per K of covariate): ERA5 2.6
(1.8–3.5), MPI-GE 1.05, EPOC 1.1, EERIE 0.4. ERA5's distribution widens
faster than in any model, and the 95% interval excludes all three.
**Forced shift** (μ₁, heatwave days per K): ERA5 6.9, MPI-GE 6.1, EPOC 7.5,
EERIE 4.1.

**Anomalies from the forced mean** (2003: +4.8 days, 2026: +7.2 days). In
absolute terms 2026 is reached by:

- 2–5% of seasons in EERIE and the controls;
- 5% of MPI-GE 1991–2026;
- 13% of MPI-GE at ERA5's 2026 warming, whose spread there is 6 days.

Divided by ERA5's local spread (std within ±15 years: 2.4 days around 2003,
3.6 around 2026), **2003 and 2026 are both about +2.0σ**. Every simulation
produces a ≥2σ season 3–6% of the time. Neither summer lies outside the
simulated internal variability once the growing spread is accounted for.

### 3. Return periods in the simulated climates (figs 15, 16)

Heatwave days; every simulation is counted against its own 1991–2020
threshold.

| source / climate | 2003 (10.4 d) | 2026 (20.1 d) |
|---|---|---|
| MPI-GE 1850–1900 (2,550 seasons) | 54 yr (41–74) | 1 season: ~2,500 yr (460–10⁵) |
| MPI-GE at ERA5 2003 warming (200) | 12 yr (7.5–20) | 1 season: 200 yr (36–7,900) |
| MPI-GE 2016–2036 (1,050) | 2.8 yr | 20 yr (15–26) |
| MPI-GE at ERA5 2026 warming (350) | 1.8 yr | **9 yr (6.7–12)** |
| EPOC hist, own 2003 / 2025 | 8 yr | 15 yr (5–170) |
| EERIE hist, own 2003 / 2026 | 20 yr (12–56) | **190 yr (61–5,400)** |
| controls (GEV): EERIE / Sapphire / EPOC | 7 / 10 / 9 yr | 99 / 140 / 78 yr; 1 season ≥ 20 d in 225 control years (EPOC) |

- **MPI-GE agrees with the ERA5 shift+scale fit** for 2026: about 1-in-9
  years at the matched warming level, 1-in-20 around calendar 2026. It does
  so for the right reason, a large forced shift with growing spread, but its
  spread is too large overall.
- **EERIE**, the model with ERA5-like variability, has a weaker forced shift
  and almost no spread growth. In it, 2026 is a 1-in-200-year summer.
- **The km-scale controls** give an unforced 2026-sized summer about 1-in-100
  years. That is the upper end of what internal variability can do in a
  1991–2020-like climate.
- **2003** is a 1-in-10 to 1-in-20-year summer in its own climate in most
  simulations, more frequent than ERA5's 1-in-34, because the models are more
  variable. In 1850–1900, MPI-GE makes it 1-in-54 years.

### 4. Return period in the MPI-GE climate of each year (fig 18)

This section is a **distribution-only view, with no trend model of the
observations**. A GEV is a stationary distribution, so it is fitted where
stationarity holds:

- For every year y, the fit runs across the 50 MPI-GE members in the
  calendar years y−5 … y+4. That gives 500 seasons of one background climate,
  so "1-in-N years" means 1 in N seasons of that climate.
- Each ERA5 summer is then scored against the GEV of its own year
  (`03_mpige_window_rp.py`).
- MPI-GE is used as simulated (no variability correction) and matched by
  calendar year.
- The GEV and the plain count of the 500 seasons agree everywhere: within
  13% in every year, and within 11% for return periods of 5 years or more.

The rarest ERA5 summers in the climate of their year (95% bootstrap CI;
k = MPI-GE seasons out of 500 at or above the observed value):

| year | heatwave days | return period | k / 500 |
|---|---|---|---|
| **2026** | 20.1 | **22 yr (17–31)** | 25 |
| **2003** | 10.4 | **13 yr (10–16)** | 36 |
| 1947 | 6.9 | 12 yr (10–16) | 38 |
| 2022 | 14.6 | 10 yr (8–12) | 55 |
| 2025, 2024, 2015 | 14.7, 14.3, 10.4 | 7–8 yr | 69–76 |

These are the only three summers rarer than 1-in-10 years. In this framing,
**2026 and 2003 are the two most extreme summers for their time**, consistent
with the trend-removed ranking of fig 17, and with no assumption about the
shape of the trend in the observations.

The same values scored in the MPI-GE climates of other years (fig 18c,
dashed and dotted curves; 1950s = geometric mean of the 1950–1959 windows):

| | 1950s | 2003 | 2026 |
|---|---|---|---|
| 2003 value (10.4 d) | 67 yr | 13 yr | 3.0 yr |
| 2026 value (20.1 d) | ~1,100 yr | 147 yr | 22 yr |

The 2026-value curve is noisy before about 1990, between about 280 and
2,800 years from window to window: there it is the far tail of a narrow
distribution, so read it as "several hundred to a few thousand years".

Why this differs from the ERA5-only shift+scale fit, which puts 2026 at
1-in-7 years:

1. **Calendar matching.** MPI-GE's 2021–2030 climate is cooler than the
   observed 2026 climate; ERA5's 2026 global SST level is reached in MPI-GE
   only around 2039–2045. At that warming level 2026 is about 1-in-9 years
   (section 3).
2. **Variability.** MPI-GE varies about 1.4× as much as ERA5. Taken alone this
   would *shorten* the return periods, so it works against point 1, not with
   it.

Limits of the window approach:

- **Assumed stationarity.** The ensemble mean still rises inside a window, by
  up to 3.4 days (about 0.6σ) in the windows centred on 1997 and 2019–2023,
  versus about 1 day around 2003. That widens the recent windows slightly and
  lengthens recent return periods.
- **1947.** It ranks third because 1940s ERA5 sits against a low MPI-GE
  baseline (6.9 days against a window median of 2.7). Early ERA5 is the least
  certain part of the record.

## Bottom line

1. **2003** is the most robust "extreme for its time" summer in the record:
   1-in-30 to 1-in-70 years in its own climate, and about 1-in-4,000 years in
   a 1950s climate. Today it is an every-other-year summer.
2. **2026** is the largest heatwave summer in absolute terms by a wide
   margin: 20 days against the previous maximum of 14.7 in 2025. It is
   exceptional against any climate before about 2015: >1-in-1,000 years in
   the 2003 climate, practically impossible in the 1950s. **In the climate of
   2026, though, it is about a 1-in-7 to 1-in-9-year summer** (ERA5 and
   MPI-GE agree). "2026 is the strongest even with the trend removed" holds
   only for a linear trend with constant spread, which the data reject.
3. The deciding factor is how fast the heatwave-day distribution widens with
   warming. ERA5's widens faster than in any of the models, and EERIE (the
   model with realistic variability) widens slowest, which is why it rates
   2026 as a 1-in-200-year summer. That observed-vs-model gap is the open
   question worth pursuing.

## Caveats

- The covariate is the global-mean SST, not global-mean surface temperature.
  The 2023–26 SST surge sets ERA5's "2026 climate" at +0.47 K. Longer
  smoothing (+0.31 / +0.20 K) leaves the 2026 result unchanged (1-in-8 to
  1-in-9 years).
- Heatwave days are a bounded, discrete-ish count. The GEV fits pass KS
  checks (p = 0.3–0.95), but return periods beyond about 3× the record length
  are extrapolation. The upper CI bounds hit the 10⁶-year cap in the 1950s
  rows.
- Seasons are treated as independent: no serial correlation in the bootstrap
  or the Clopper–Pearson intervals.
- Simulations are scored against their own 1991–2020 threshold, so absolute
  comparisons mix in model variability bias. The standardised comparison
  (section 2) removes it.
- EPOC control heatwaves use daily-mean temperature (its Tmax is on tape
  only). The EPOC hist cumulative-heat fit sits at the ξ = 0.5 bound (36
  seasons).
- MPI-ESM1.2-ER has no daily output and is absent.
