# Data dictionary — `master_scd_utilization_GA_CA_2010_2018.csv`

**Dimensions:** 252 rows × 22 columns
**Grain:** one row per `data_type` × `state` × `year` × `age_group`
(2 data types × 2 states × 9 years × 7 age groups = 252)

**Source:** CDC Sickle Cell Data Collection (SCDC) program, demographic and healthcare
utilization tables for California and Georgia, accessed March 15, 2026.
<https://www.cdc.gov/sickle-cell-research/php/data/index.html>

All values are aggregate, publicly reported counts and proportions. No individual-level data
are present.

---

## Key / stratifying variables

| Column | Type | Description |
|---|---|---|
| `record_id` | integer | Row identifier |
| `state` | string | `California` (Medicaid expansion) or `Georgia` (non-expansion) |
| `year` | integer | Calendar year, 2010–2018 |
| `age_group` | string | `<10`, `10–19`, `20–29`, `30–39`, `40–49`, `50–59`, `60+` (en-dashes, not hyphens) |
| `data_type` | string | `ed` = emergency department stream; `hospital` = inpatient stream |

**Important:** `data_type` determines which measure columns are populated. The two streams are
stacked in one file, so roughly half of every measure column is `NA` by design.

| | `data_type == "ed"` (126 rows) | `data_type == "hospital"` (126 rows) |
|---|---|---|
| ED measures | populated | all `NA` |
| Inpatient measures | all `NA` | populated |
| Payer proportions | populated (share of ED visits) | populated (share of admissions) |
| `people_with_scd_total` | populated | populated |

## Denominator

| Column | Type | Description |
|---|---|---|
| `people_with_scd_total` | integer | Number of people with SCD identified in that state × year × age group. Used as the population weight and as the denominator for all per-person rates. |

## Emergency department measures (`data_type == "ed"` only)

| Column | Type | Description |
|---|---|---|
| `ed_visits_total` | float | Total ED visits in the stratum |
| `ed_visits_per_person_mean` | float | Mean ED visits per person with SCD |
| `readmission_30day_prop` | float | Proportion of encounters followed by a readmission within 30 days |
| `no_readmission_prop` | float | Complement of `readmission_30day_prop` |

## Inpatient measures (`data_type == "hospital"` only)

| Column | Type | Description |
|---|---|---|
| `admissions_total` | float | Total inpatient admissions in the stratum |
| `admissions_per_person_mean` | float | Mean admissions per person with SCD |
| `hospital_days_total` | float | Total inpatient days |
| `hospital_days_per_person_mean` | float | Mean inpatient days per person with SCD |
| `los_per_admission_mean` | float | Mean length of stay per admission, days |
| `admissions_from_ed_prop` | float | Proportion of admissions originating in the ED |
| `direct_admission_prop` | float | Proportion of admissions not originating in the ED (complement of the above) |

## Payer mix (both streams; proportions of encounters in that stream)

| Column | Type | Description |
|---|---|---|
| `payer_medicaid_prop` | float | Share of encounters with Medicaid as primary payer |
| `payer_medicare_prop` | float | Share with Medicare — used as the negative-control outcome |
| `payer_private_prop` | float | Share with private insurance |
| `payer_selfpay_prop` | float | Share self-pay. **Used as the proxy for uninsured status**, following the HCUP convention. Mapped to `uninsured_prop` in `prepare_df()`. |
| `payer_other_prop` | float | Share with any other payer |

The five payer proportions sum to approximately 1.00 within each row;
`02_payer_mix_analysis.R` runs this as an automatic sanity check after combining streams.

---

## Missing data

Beyond the by-design stream missingness described above, there is **one real gap**:

| Stratum | Columns affected | Rows |
|---|---|---|
| Georgia, 2015, all 7 age groups, `data_type == "hospital"` | `admissions_from_ed_prop`, `direct_admission_prop` | 7 |

These values were not reported by SCDC. `lm()` drops them via its default `na.omit`, so the
**admissions-from-ED model is fitted on 119 rows rather than 126**; every other model uses the
full 126 rows for its stream. A sensitivity analysis excluding 2015 for *both* states (so the
panel stays balanced) leaves the joint DiD test for admissions-from-ED significant.

---

## Derived variables created in `00_setup_functions.R`

These are not in the CSV; `prepare_df()` adds them at load time.

| Variable | Definition |
|---|---|
| `post_aca` | `1` if `year >= 2014`, else `0` |
| `california` | `1` if `state == "California"`, else `0` |
| `time` | `year - 2010` (continuous time index, 0 at study start) |
| `time_post` | `year - 2013` if `year >= 2014`, else `0` (post-period slope counter) |
| `cluster_id` | `state_ageGroup` — the clustering unit for Liang–Zeger standard errors (14 clusters) |
| `uninsured_prop` | Alias for `payer_selfpay_prop` |
| `n_medicaid`, `n_medicare`, `n_private`, `n_selfpay` | Payer proportion × `people_with_scd_total`, rounded |

`02_payer_mix_analysis.R` additionally builds `total_encounters` and the `comb_*` combined
payer proportions by converting proportions to counts, summing across streams, and
re-dividing — see that script's Section 3.
