# ACA and Sickle Cell Disease: Insurance Coverage and Healthcare Utilization in California vs. Georgia, 2010–2018

Analysis code and analytic dataset for:

> deBettencourt J, Hulbert M, Archer NM. *A Difference-in-Differences Study on the Impact of
> the Affordable Care Act on Insurance Coverage and Healthcare Utilization Among Patients with
> Sickle Cell Disease.*

<!-- Add after the first Zenodo release:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

---

## Overview

This repository contains everything needed to reproduce the analyses, tables, and figures in
the manuscript. The study uses an **interrupted time series (ITS)** design with a
**difference-in-differences (DiD)** estimator to compare California (adopted Medicaid
expansion January 1, 2014) with Georgia (did not adopt), using publicly reported aggregate
data from the CDC Sickle Cell Data Collection (SCDC) program.

| | |
|---|---|
| **Design** | Retrospective, population-level quasi-experimental ITS with DiD estimation |
| **Unit of analysis** | age group × state × year |
| **Treatment state** | California (Medicaid expansion) |
| **Control state** | Georgia (non-expansion) |
| **Pre-period** | 2010–2013 |
| **Post-period** | 2014–2018 |
| **Observations** | 252 rows total (126 ED, 126 inpatient) |
| **Data source** | CDC Sickle Cell Data Collection (SCDC) program |
| **Reporting guideline** | STROBE with the RECORD extension |

No individual-level data are used; all values are aggregate counts and proportions as
published by CDC/SCDC. No data use agreement is required and the study did not require IRB
review.

---

## Repository contents

```
.
├── data/
│   ├── master_scd_utilization_GA_CA_2010_2018.csv   analytic dataset (252 × 22)
│   └── DATA_DICTIONARY.md                           variable definitions and missingness
├── R/
│   ├── 00_setup_functions.R          import, preparation, shared ITS/DiD functions, plot theme
│   ├── 01_Secondary_ED_analysis.R      ED-stream payer mix and all utilization outcomes
│   ├── 02_payer_mix_analysis.R       inpatient-only and combined ED+inpatient payer mix
│   └── 03_age_stratified_analysis.R  age-band effect modification (GEE)
├── output/
│   ├── figures/   generated .png figures
│   └── tables/    generated .csv result tables
├── run_all.R      runs all three analysis scripts and writes sessionInfo
├── CITATION.cff   machine-readable citation metadata
└── LICENSE        MIT license (code)
```

`00_setup_functions.R` is sourced by scripts 01 and 02 and is not run on its own.
`03_age_stratified_analysis.R` is self-contained. The three analysis scripts are independent
of each other and can be run in any order.

### Which script produces which manuscript element

| Manuscript element | Script |
|---|---|
| Table 1 — study population characteristics | `01`, Section 2 (descriptives) |
| Table 2 — inpatient-only and combined payer mix (primary outcome) | `02`, Sections 4–6 |
| Table 3 / Table S1 — ED-stream payer mix and utilization outcomes | `01`, Sections 5–7, 11 |
| Table 4 — age-band contrasts | `03`, Section 8 |
| Table S2 — utilization level *and* slope parameters | `01`, Section 13 |
| Payer-mix trend figures | `01` Figure 1; `02` Figures 9–13 |
| Utilization trend figures | `01` Figures 2–8 |
| Age-band forest plot | `03`, Section 10 |

---

## Reproducing the analysis

**Requirements:** R ≥ 4.4 (analyses were run under R 4.5.3). Install dependencies:

```r
install.packages(c("tidyverse", "lmtest", "sandwich", "nlme", "emmeans",
                   "broom", "car", "ggplot2", "scales", "geepack", "MASS", "here"))
```

**Run everything:**

```bash
git clone https://github.com/<username>/scd-aca-did.git
cd scd-aca-did
Rscript run_all.R
```

Or, from R / RStudio with the repository root as the working directory:

```r
source("R/01_secondary_ED_analysis.R")
source("R/02_payer_mix_analysis.R")
source("R/03_age_stratified_analysis.R")
```

All paths are resolved relative to the repository root, so no paths need editing. Figures are
written to `output/figures/` and result tables to `output/tables/`.

---

## Statistical approach

Scripts 01 and 02 fit, for each outcome, a weighted segmented regression:

```
Y_it = b0 + b1·time + b2·post_aca + b3·time_post
     + b4·california + b5·(california × time)
     + b6·(california × post_aca)    <- DiD level change  [PRIMARY]
     + b7·(california × time_post)   <- DiD slope change  [SECONDARY]
     + b8·age_group + e_it
```

- **Standard errors.** Cluster-robust (Liang–Zeger) sandwich estimator, clustered on
  state × age group (14 clusters).
- **Joint testing.** `linearHypothesis()` F-test of H0: b6 = 0 *and* b7 = 0.
- **Multiple comparisons.** Benjamini–Hochberg FDR correction applied to the DiD level-change
  p-values, **independently within each analysis** rather than pooled across all outcomes:
  9 outcomes in `01`, 8 in `02`, 4 contrasts in `03`.

### Weights

Weights differ by outcome, matching each outcome's natural denominator. This is deliberate and
worth noting when comparing coefficients across scripts:

| Outcome family | Weight | Script |
|---|---|---|
| ED payer mix, ED utilization, per-person outcomes | `people_with_scd_total` | `01` |
| Inpatient payer mix, LOS, admissions-from-ED | `admissions_total` | `01`, `02` |
| **Combined ED + inpatient payer mix (primary outcome)** | `ed_visits_total + admissions_total` | `02` |
| Age-band interaction model | `people_with_scd_total` | `03` |

### Post-period slope coding

Scripts 00/01/02 code `time_post = year − 2013` (2014 = 1 … 2018 = 5). Script 03 codes
`time_after = max(year − 2014, 0)` (2014 = 0 … 2018 = 4). Both are standard segmented-regression
parameterizations, but they anchor the level-change coefficient at different points on the time
axis, so **DiD level estimates are not numerically interchangeable between script 03 and
scripts 01/02**; slope estimates are unaffected. Table 4 reports contrasts *between* age bands
within script 03's model, which are internally consistent. See the header of
`R/03_age_stratified_analysis.R`.

### Assumption checks and sensitivity analyses included

| Check | Script / section |
|---|---|
| Parallel pre-trends, overall (2010–2013) | `01`, Section 3 |
| Parallel pre-trends, within each age band | `03`, Section 4 |
| Placebo intervention at 2012 | `01`, Section 10B |
| Medicare as negative-control outcome | `01`, Section 10C; `02`, Section 4 |
| Exclusion of 2014 transition year | `01`, Section 10A |
| GLS with AR(1) errors | `01`, Section 8 (provided, commented out) |
| Residual normality for LOS, and log-LOS sensitivity | `01`, Sections 14 and 14B |
| Rank check on the GEE sandwich covariance | `03`, Section 5b |
| Separate per-band models as a robustness check | `03`, Section 9 |

**Length of stay.** Residuals from the LOS model depart from normality (Shapiro–Wilk
W = 0.97, p = 0.004); a log transformation does not resolve this (W = 0.97, p = 0.003) and
yields substantively unchanged DiD estimates. Raw-scale results are reported. Inference
throughout relies on cluster-robust standard errors, which do not require normally distributed
errors.

**Age-band model caveats.** The full interaction model estimates 24 parameters from 14
clusters, so its empirical sandwich covariance can be rank-deficient. Script 03 checks this
explicitly and flags any contrast falling in a near-null direction via a `reliable` column.
The Young Adult band maps to a single CDC stratum (20–29), giving only 2 of the 14 clusters;
standalone band-level DiD estimates for that band are therefore not reported. 

---

## Combining ED and inpatient payer mix

ED visits and inpatient admissions have different denominators, so their payer proportions
cannot simply be averaged. `02_payer_mix_analysis.R` instead (1) converts each stratum's payer
proportions back to absolute encounter counts, (2) sums counts across the two streams, and
(3) re-expresses the totals as proportions of all encounters in that stratum. A sum check
verifying that combined proportions total ~1.00 runs automatically. Script 03 constructs the
same combined Medicaid proportion independently, using the identical formula.

---

## Data provenance and reuse

The analytic dataset was compiled from the publicly posted CDC Sickle Cell Data Collection
demographic and healthcare utilization tables for California and Georgia
(<https://www.cdc.gov/sickle-cell-research/php/data/index.html>). Annual sheets were extracted
per state, combined into state-level files, and merged into the single analytic file here. No
additional abstraction or derivation was performed prior to analysis; payer categories and
outcomes are defined and coded in the original CDC dataset.

Cases were identified by SCDC using administrative criteria — newborn screening, laboratory-
confirmed diagnosis, or multiple healthcare encounters with SCD-related ICD-9-CM/ICD-10-CM
codes — as described in Snyder et al., *MMWR Surveill Summ* 2022;71(9):1–18.

If you use these data, please cite both this repository and the original CDC SCDC source.

---

## Known limitations of the dataset

- **California and Georgia only.** These were the first two SCDC states; data from additional
  states were not publicly available at the time of analysis.
- **Administrative case identification** may both miss patients with SCD and include patients
  without it. The period also spans the ICD-9 to ICD-10 transition.
- **Acute care only.** The dataset captures ED visits and inpatient admissions, not outpatient
  care, infusion-center visits, or the majority of people with SCD who have no acute care
  encounter in a given year.
- **No genotype**, a major determinant of disease severity, is recorded.
- **Coarse age strata.** The seven CDC bins do not isolate the 19–26 group targeted by the
  ACA dependent-coverage provision; 20–29 is used as the Young Adult proxy.
- **"Self-pay" is a proxy for uninsured status**, following the HCUP convention.
- **Missing cells.** The 7 Georgia 2015 inpatient strata have no `admissions_from_ed_prop` or
  `direct_admission_prop` value, so that model is fitted on 119 rows rather than 126. All
  other outcomes use the full 126 rows per stream. A sensitivity analysis excluding 2015 for
  both states leaves the admissions-from-ED slope result unchanged. See
  `data/DATA_DICTIONARY.md`.

---

## License

- **Code** (`R/`, `run_all.R`): MIT License — see [`LICENSE`](LICENSE).
- **Data** (`data/`): derived from U.S. federal government public-domain material; released
  here under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). Please cite the
  CDC SCDC program as the original source.

## Authors

Joseph deBettencourt, MD; Monica Hulbert, MD; Natasha M. Archer, MD —
Dana-Farber/Boston Children's Cancer and Blood Disorders Center.
Corresponding author: joseph_debettencourt@dfci.harvard.edu

## Funding and acknowledgements

Joseph deBettencourt is supported by an NIH T32 Training Grant (PI Bauer, 5T32HL007574-45).
Thanks to the hematology department of Dana-Farber/Boston Children's Center for Cancer and
Blood Disorders. As well as the Computational Genomics Summer Institute (CGSI) for advice and assistance with the development of this project. 

## Use of AI assistance

Claude (Anthropic; Sonnet 4.6) was used for R code debugging and code cleaning, and as a
data-quality check to confirm that the compiled .csv matched the source CDC sheets. The study
authors independently reviewed that comparison and all code, and take full responsibility for
its correctness.
