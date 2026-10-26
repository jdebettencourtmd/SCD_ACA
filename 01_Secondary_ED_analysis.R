# =============================================================================
# SECONDARY ANALYSIS: ED-BASED INSURANCE COVERAGE & UTILIZATION OUTCOMES
# SCD Patients, California (Medicaid expansion) vs. Georgia (non-expansion),
# 2010-2018
#
# This is the SECONDARY, pre-specified analysis restricted to ED encounters
# (df_ed) for insurance-coverage outcomes and to the appropriate ED/hospital
# stream for each utilization outcome. 
#
# A separate, independently BH-corrected set of models using inpatient-only
# and combined ED+inpatient encounters lives in
# 02_payer_mix_analysis.R and should be run separately.
#
# Requires 00_setup_functions.R to have been sourced first.
# =============================================================================

source(if (requireNamespace("here", quietly = TRUE)) {
  here::here("R", "00_setup_functions.R")
} else "R/00_setup_functions.R")

# -----------------------------------------------------------------------------
# 2. DESCRIPTIVE STATISTICS
# -----------------------------------------------------------------------------

# ── 2a. ED data: pre vs. post ACA by state ────────────────────────────────────
desc_ed <- df_ed %>%
  group_by(state, post_aca) %>%
  summarise(
    n_obs           = n(),
    total_persons   = sum(people_with_scd_total),
    total_ed_visits = sum(ed_visits_total),
    mean_ed_rate    = weighted.mean(ed_visits_per_person_mean, people_with_scd_total),
    mean_uninsured  = weighted.mean(uninsured_prop,            people_with_scd_total),
    mean_medicaid   = weighted.mean(payer_medicaid_prop,       people_with_scd_total),
    mean_private    = weighted.mean(payer_private_prop,        people_with_scd_total),
    mean_readmit30  = weighted.mean(readmission_30day_prop,    people_with_scd_total),
    .groups = "drop"
  ) %>%
  mutate(period = if_else(post_aca == 0,
                          "Pre-ACA (2010-2013)",
                          "Post-ACA (2014-2018)")) %>%
  select(-post_aca)

cat("\n=== DESCRIPTIVE: ED Data ===\n")
print(desc_ed)

# ── 2b. Hospital data: pre vs. post ACA by state ──────────────────────────────
desc_hosp <- df_hosp %>%
  group_by(state, post_aca) %>%
  summarise(
    n_obs              = n(),
    total_persons      = sum(people_with_scd_total),
    total_admissions   = sum(admissions_total),
    mean_adm_rate      = weighted.mean(admissions_per_person_mean, people_with_scd_total),
    mean_los           = weighted.mean(los_per_admission_mean,     admissions_total),
    mean_hosp_days_pp  = weighted.mean(hospital_days_per_person_mean, people_with_scd_total),
    mean_adm_from_ed   = weighted.mean(admissions_from_ed_prop,    admissions_total),
    mean_uninsured     = weighted.mean(uninsured_prop,             people_with_scd_total),
    mean_medicaid      = weighted.mean(payer_medicaid_prop,        people_with_scd_total),
    mean_private       = weighted.mean(payer_private_prop,         people_with_scd_total),
    .groups = "drop"
  ) %>%
  mutate(period = if_else(post_aca == 0,
                          "Pre-ACA (2010-2013)",
                          "Post-ACA (2014-2018)")) %>%
  select(-post_aca)

cat("\n=== DESCRIPTIVE: Hospital Data ===\n")
print(desc_hosp)

# ── 2c. Annual trend summaries (for plotting) ─────────────────────────────────
trend_ed <- df_ed %>%
  group_by(state, year) %>%
  summarise(
    wt_ed_rate   = weighted.mean(ed_visits_per_person_mean, people_with_scd_total),
    wt_uninsured = weighted.mean(uninsured_prop,            people_with_scd_total),
    wt_medicaid  = weighted.mean(payer_medicaid_prop,       people_with_scd_total),
    wt_private   = weighted.mean(payer_private_prop,        people_with_scd_total),
    wt_readmit   = weighted.mean(readmission_30day_prop,    people_with_scd_total),
    .groups = "drop"
  )

trend_hosp <- df_hosp %>%
  group_by(state, year) %>%
  summarise(
    wt_los          = weighted.mean(los_per_admission_mean,        admissions_total),
    wt_adm_rate     = weighted.mean(admissions_per_person_mean,    people_with_scd_total),
    wt_hosp_days_pp = weighted.mean(hospital_days_per_person_mean, people_with_scd_total),
    wt_adm_from_ed  = weighted.mean(admissions_from_ed_prop,       admissions_total),
    wt_uninsured    = weighted.mean(uninsured_prop,                 people_with_scd_total),
    wt_medicaid     = weighted.mean(payer_medicaid_prop,            people_with_scd_total),
    wt_private      = weighted.mean(payer_private_prop,             people_with_scd_total),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 3. PRE-TREND PARALLELISM CHECK (Key DiD Assumption)
# -----------------------------------------------------------------------------
# DiD assumes parallel counterfactual trends pre-intervention.
# H0: state x time interaction = 0 in pre-period (should NOT be significant).
# Reference: Angrist & Pischke (2009) "Mostly Harmless Econometrics"

pre_ed   <- df_ed   %>% filter(year < 2014)
pre_hosp <- df_hosp %>% filter(year < 2014)

run_parallel_test <- function(outcome, label, data) {
  cat("\n=== PARALLEL PRE-TREND TEST:", label, "===\n")
  mod <- do.call("lm", list(
    formula = as.formula(paste(outcome,
                               "~ time * california + age_group")),
    data    = data,
    weights = data$people_with_scd_total
  ))
  print(coeftest(mod, vcov = vcovCL(mod, cluster = ~cluster_id)))
}

run_parallel_test("ed_visits_per_person_mean", "ED visits/person",    pre_ed)
run_parallel_test("uninsured_prop",            "Uninsured proportion", pre_ed)
run_parallel_test("payer_medicaid_prop",        "Medicaid proportion",  pre_ed)
run_parallel_test("admissions_per_person_mean","Admissions/person",   pre_hosp)
run_parallel_test("los_per_admission_mean",    "LOS/admission",       pre_hosp)


# -----------------------------------------------------------------------------
# 5. PRIMARY ANALYSES: INSURANCE COVERAGE OUTCOMES (ED data)
# -----------------------------------------------------------------------------
# Hypothesis: Fewer self-pay post-ACA in both states; larger decrease in CA.
# CA: increase in Medicaid (expansion) + private (subsidies).
# GA: increase in private only; no Medicaid change.

cat("\n\n", strrep("#", 70), "\n")
cat("# SECTION 5: INSURANCE COVERAGE OUTCOMES\n")
cat(strrep("#", 70), "\n")

m_uninsured <- run_its("uninsured_prop",       "Uninsured (self-pay) proportion", df_ed)
m_medicaid  <- run_its("payer_medicaid_prop",  "Medicaid proportion",             df_ed)
m_private   <- run_its("payer_private_prop",   "Private insurance proportion",    df_ed)

ft_uninsured <- joint_did_test(m_uninsured, "Uninsured proportion")
ft_medicaid  <- joint_did_test(m_medicaid,  "Medicaid proportion")
ft_private   <- joint_did_test(m_private,   "Private insurance proportion")


# -----------------------------------------------------------------------------
# 6. PRIMARY ANALYSES: ED UTILIZATION OUTCOMES (ED data)
# -----------------------------------------------------------------------------
# Hypothesis: Decreased ED utilization post-ACA; larger decrease in CA.
# Mechanism: improved outpatient access with better insurance coverage.

cat("\n\n", strrep("#", 70), "\n")
cat("# SECTION 6: ED UTILIZATION OUTCOMES\n")
cat(strrep("#", 70), "\n")

m_ed_rate  <- run_its("ed_visits_per_person_mean", "ED visits per person",          df_ed)
m_readmit  <- run_its("readmission_30day_prop",    "30-day readmission proportion", df_ed)

ft_ed_rate <- joint_did_test(m_ed_rate, "ED visits per person")
ft_readmit <- joint_did_test(m_readmit, "30-day readmission proportion")


# -----------------------------------------------------------------------------
# 7. PRIMARY ANALYSES: INPATIENT UTILIZATION OUTCOMES (Hospital data)
# -----------------------------------------------------------------------------
# Hypothesis: Decreased hospitalizations post-ACA; larger decrease in CA.

cat("\n\n", strrep("#", 70), "\n")
cat("# SECTION 7: INPATIENT UTILIZATION OUTCOMES\n")
cat(strrep("#", 70), "\n")

m_adm_rate   <- run_its("admissions_per_person_mean",    "Admissions per person",
                         df_hosp)
m_adm_ed     <- run_its("admissions_from_ed_prop",       "Admissions from ED (proportion)",
                         df_hosp, weight_var = "admissions_total")
m_los        <- run_its("los_per_admission_mean",        "LOS per admission (days)",
                         df_hosp, weight_var = "admissions_total")
m_hosp_days  <- run_its("hospital_days_per_person_mean", "Hospital days per person",
                         df_hosp)

ft_adm_rate  <- joint_did_test(m_adm_rate,  "Admissions per person")
ft_adm_ed    <- joint_did_test(m_adm_ed,    "Admissions from ED")
ft_los       <- joint_did_test(m_los,       "LOS per admission")
ft_hosp_days <- joint_did_test(m_hosp_days, "Hospital days per person")


# -----------------------------------------------------------------------------
# 8. GLS WITH AR(1) CORRECTION (if autocorrelation present)
# -----------------------------------------------------------------------------
# Run ACF plots to check residual autocorrelation for each primary outcome.
# If systematic pattern exists, use GLS AR(1) as confirmatory analysis.
# Reference: Penfold & Zhang (2013) Prev Chronic Dis

check_acf <- function(model_result, label) {
  acf(residuals(model_result$model),
      main = paste("Residual ACF:", label))
}


run_gls_ar1 <- function(outcome, label, input_df,
                         weight_var = "people_with_scd_total") {
  cat("\n=== GLS AR(1):", label, "===\n")
  data_sorted <- input_df %>% arrange(cluster_id, year)
  wt_formula  <- as.formula(paste("~1/", weight_var))
  mod_gls <- gls(
    its_formula(outcome),
    data        = data_sorted,
    weights     = varFixed(wt_formula),
    correlation = corAR1(form = ~ time | cluster_id),
    method      = "ML"
  )
  print(summary(mod_gls))
  return(mod_gls)
}


# -----------------------------------------------------------------------------
# 9. SUBGROUP ANALYSES BY AGE GROUP
# -----------------------------------------------------------------------------
# Rationale: ACA Medicaid expansion eligibility differs substantially by age.
# Pediatric SCD (<19): largely covered by CHIP pre-ACA; dependent coverage
#   to age 26 provision relevant for 19-25 age band within 10-19/20-29 groups.
# Young adults (20-29): most likely to gain coverage from expansion/subsidies?


age_groups_of_interest <- c("<10", "10\u201319", "20\u201329", "30\u201339")

run_subgroup_its <- function(input_df) {
  map_dfr(age_groups_of_interest, function(ag) {
    sub_df <- input_df %>% filter(as.character(age_group) == ag)
    cat("\n--- Subgroup:", ag, "---\n")

    sub_formula <- function(outcome) {
      as.formula(paste(
        outcome,
        "~ time + post_aca + time_post +",
        "california + california:time +",
        "california:post_aca + california:time_post"
      ))
    }

    map_dfr(
      c("uninsured_prop", "payer_medicaid_prop",
        "payer_private_prop", "ed_visits_per_person_mean"),
      function(outcome) {
        mod <- do.call("lm", list(
          formula = sub_formula(outcome),
          data    = sub_df,
          weights = sub_df$people_with_scd_total
        ))
        tidy(mod) %>% mutate(outcome = outcome, age_group = ag)
      }
    )
  })
}

cat("\n\n=== SUBGROUP ANALYSES: ED DATA ===\n")
subgroup_ed_results <- run_subgroup_its(df_ed)
print(subgroup_ed_results, n = Inf)


# -----------------------------------------------------------------------------
# 10. SENSITIVITY ANALYSES
# -----------------------------------------------------------------------------

# ── 10A. Exclude 2014 (transition year) ──────────────────────────────────────
cat("\n=== SENSITIVITY: Excluding 2014 ===\n")
df_ed_no2014   <- df_ed   %>% filter(year != 2014)
df_hosp_no2014 <- df_hosp %>% filter(year != 2014)

m_sens_unins  <- run_its("uninsured_prop",            "Uninsured (excl. 2014)",   df_ed_no2014)
m_sens_ed     <- run_its("ed_visits_per_person_mean", "ED rate (excl. 2014)",     df_ed_no2014)
m_sens_adm    <- run_its("admissions_per_person_mean","Adm rate (excl. 2014)",    df_hosp_no2014)
m_sens_los    <- run_its("los_per_admission_mean",    "LOS (excl. 2014)",         df_hosp_no2014,
                          weight_var = "admissions_total")

# ── 10B. Placebo test: fake intervention at 2012 ─────────────────────────────

cat("\n=== SENSITIVITY: Placebo test (fake intervention 2012) ===\n")

make_placebo_df <- function(input_df) {
  input_df %>%
    filter(year <= 2013) %>%
    mutate(
      post_aca  = as.integer(year >= 2012),
      time_post = as.integer(if_else(year >= 2012, year - 2011, 0))
    )
}

df_ed_placebo   <- make_placebo_df(df_ed)
df_hosp_placebo <- make_placebo_df(df_hosp)

m_placebo_unins <- run_its("uninsured_prop",            "PLACEBO: Uninsured", df_ed_placebo)
m_placebo_ed    <- run_its("ed_visits_per_person_mean", "PLACEBO: ED rate",   df_ed_placebo)
m_placebo_adm   <- run_its("admissions_per_person_mean","PLACEBO: Adm rate",  df_hosp_placebo)

# ── 10C. Medicare as negative control outcome ─────────────────────────────────
# Medicare eligibility not affected by Medicaid expansion.

cat("\n=== SENSITIVITY: Negative control - Medicare proportion ===\n")
m_medicare_neg <- run_its("payer_medicare_prop", "Medicare (negative control)", df_ed)


# -----------------------------------------------------------------------------
# 11. MULTIPLE COMPARISONS CORRECTION & JOINT F-TEST SUMMARY
# -----------------------------------------------------------------------------
# Apply Benjamini-Hochberg FDR correction across all primary DiD level-change
# p-values

pval_table <- bind_rows(
  # Insurance outcomes
  get_did_pval(m_uninsured, "Uninsured proportion"),
  get_did_pval(m_medicaid,  "Medicaid proportion"),
  get_did_pval(m_private,   "Private insurance proportion"),
  # ED utilization
  get_did_pval(m_ed_rate,   "ED visits per person"),
  get_did_pval(m_readmit,   "30-day readmission proportion"),
  # Inpatient utilization
  get_did_pval(m_adm_rate,  "Admissions per person"),
  get_did_pval(m_adm_ed,    "Admissions from ED proportion"),
  get_did_pval(m_los,       "LOS per admission"),
  get_did_pval(m_hosp_days, "Hospital days per person")
) %>%
  mutate(p_adj_BH = round(p.adjust(p_raw, method = "BH"), 4))

cat("\n\n=== MULTIPLE COMPARISONS: BH-adjusted p-values ===\n")
print(pval_table, n = Inf)

# Joint F-test results (H0: level change = 0 AND slope change = 0), compiled
# from the tibbles returned by joint_did_test() at each call above.
ftest_table <- bind_rows(
  ft_uninsured, ft_medicaid, ft_private,
  ft_ed_rate, ft_readmit,
  ft_adm_rate, ft_adm_ed, ft_los, ft_hosp_days
)

cat("\n\n=== JOINT DiD F-TESTS: summary ===\n")
print(ftest_table, n = Inf)


# -----------------------------------------------------------------------------
# 12. VISUALIZATION
# -----------------------------------------------------------------------------

# ── Figure 1: Insurance payer mix trends ─────────────────────────────────────
payer_long <- trend_ed %>%
  select(state, year, wt_uninsured, wt_medicaid, wt_private) %>%
  pivot_longer(cols      = starts_with("wt_"),
               names_to  = "payer",
               names_prefix = "wt_",
               values_to = "proportion") %>%
  mutate(payer = case_match(payer,
    "uninsured" ~ "Self-pay (Uninsured)",
    "medicaid"  ~ "Medicaid",
    "private"   ~ "Private"
  ))

fig1 <- ggplot(payer_long,
               aes(x = year, y = proportion, color = payer, shape = payer)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state) +
  aca_line + aca_label +
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(breaks = 2010:2018) +
  labs(title    = "Insurance Payer Mix Among SCD Patients, 2010-2018",
       subtitle = "Weighted by SCD population size per stratum; ED encounters",
       x = "Year", y = "Proportion of ED Encounters",
       color = "Payer", shape = "Payer",
       caption = "Dashed line = ACA implementation (Jan 2014)") +
  theme_scd
print(fig1)
ggsave(repo_path("output", "figures", "fig1_payer_mix.png"), fig1, width = 10, height = 5, dpi = 300)

# ── Figure 2: ED visits per person ───────────────────────────────────────────
fig2 <- ggplot(trend_ed, aes(x = year, y = wt_ed_rate, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_ed %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_ed %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_color_manual(values = state_colors) +
  labs(title    = "ED Visits per Person with SCD, 2010-2018",
       subtitle = "Segmented linear fit; shaded area = 95% CI",
       x = "Year", y = "Mean ED visits per person", color = "State",
       caption = "Dashed line = ACA implementation (Jan 2014)") +
  theme_scd
print(fig2)
ggsave(repo_path("output", "figures", "fig2_ed_utilization.png"), fig2, width = 8, height = 5, dpi = 300)

# ── Figure 3: Admissions per person ──────────────────────────────────────────
fig3 <- ggplot(trend_hosp, aes(x = year, y = wt_adm_rate, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_hosp %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_hosp %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_color_manual(values = state_colors) +
  labs(title    = "Inpatient Admissions per Person with SCD, 2010-2018",
       subtitle = "Segmented linear fit; shaded area = 95% CI",
       x = "Year", y = "Mean admissions per person", color = "State") +
  theme_scd
print(fig3)
ggsave(repo_path("output", "figures", "fig3_admissions.png"), fig3, width = 8, height = 5, dpi = 300)

# ── Figure 4: LOS per admission ──────────────────────────────────────────────
fig4 <- ggplot(trend_hosp, aes(x = year, y = wt_los, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_hosp %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_hosp %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_color_manual(values = state_colors) +
  labs(title    = "Mean Length of Stay per Admission - SCD Patients, 2010-2018",
       subtitle = "Weighted by admissions total; segmented linear fit with 95% CI",
       x = "Year", y = "Mean LOS (days)", color = "State") +
  theme_scd
print(fig4)
ggsave(repo_path("output", "figures", "fig4_los.png"), fig4, width = 8, height = 5, dpi = 300)

# ── Figure 5: Admissions from ED proportion ──────────────────────────────────
fig5 <- ggplot(trend_hosp, aes(x = year, y = wt_adm_from_ed, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_hosp %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_hosp %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = state_colors) +
  labs(title    = "Proportion of Inpatient Admissions via ED - SCD, 2010-2018",
       x = "Year", y = "Proportion admitted through ED", color = "State") +
  theme_scd
print(fig5)
ggsave(repo_path("output", "figures", "fig5_admissions_from_ed.png"), fig5, width = 8, height = 5, dpi = 300)

# ── Figure 6: 30-day readmission ─────────────────────────────────────────────
fig6 <- ggplot(trend_ed, aes(x = year, y = wt_readmit, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_ed %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_ed %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = state_colors) +
  labs(title = "30-Day Readmission Proportion Among SCD Patients, 2010-2018",
       x = "Year", y = "Proportion with 30-day readmission", color = "State") +
  theme_scd
print(fig6)
ggsave(repo_path("output", "figures", "fig6_readmission.png"), fig6, width = 8, height = 5, dpi = 300)

# ── Figure 7: Age-stratified ED rates ────────────────────────────────────────
fig7_data <- df_ed %>%
  group_by(state, year, age_group) %>%
  summarise(ed_rate = mean(ed_visits_per_person_mean), .groups = "drop")

fig7 <- ggplot(fig7_data, aes(x = year, y = ed_rate, color = state)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~age_group, scales = "free_y") +
  aca_line +
  scale_color_manual(values = state_colors) +
  scale_x_continuous(breaks = c(2010, 2013, 2016)) +
  labs(title = "ED Visits per Person by Age Group, 2010-2018",
       x = "Year", y = "ED visits per person", color = "State") +
  theme_scd +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(fig7)
ggsave(repo_path("output", "figures", "fig7_age_ed_rates.png"), fig7, width = 12, height = 8, dpi = 300)

# ── Figure 8: Hospital days per person (overall burden) ──────────────────────
fig8 <- ggplot(trend_hosp, aes(x = year, y = wt_hosp_days_pp, color = state)) +
  geom_point(size = 3) +
  geom_smooth(data = trend_hosp %>% filter(year <= 2013),
              method = "lm", se = TRUE, linewidth = 1) +
  geom_smooth(data = trend_hosp %>% filter(year >= 2014),
              method = "lm", se = TRUE, linewidth = 1) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_color_manual(values = state_colors) +
  labs(title    = "Hospital Days per Person with SCD, 2010-2018",
       subtitle = "Composite of admission frequency and length of stay",
       x = "Year", y = "Mean hospital days per person", color = "State") +
  theme_scd
print(fig8)
ggsave(repo_path("output", "figures", "fig8_hospital_days.png"), fig8, width = 8, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# 13. FINAL SUMMARY TABLE
# -----------------------------------------------------------------------------
all_models <- list(
  list(m_uninsured, "Uninsured proportion"),
  list(m_medicaid,  "Medicaid proportion"),
  list(m_private,   "Private insurance proportion"),
  list(m_ed_rate,   "ED visits per person"),
  list(m_readmit,   "30-day readmission proportion"),
  list(m_adm_rate,  "Admissions per person"),
  list(m_adm_ed,    "Admissions from ED proportion"),
  list(m_los,       "LOS per admission (days)"),
  list(m_hosp_days, "Hospital days per person")
)

results_table <- map_dfr(all_models, function(x) {
  extract_did_coefs(x[[1]], x[[2]])
})

cat("\n\n=== FINAL RESULTS TABLE ===\n")
print(results_table, n = Inf)

write_csv(results_table, repo_path("output", "tables", "SCD_ACA_DiD_Results.csv"))
write_csv(pval_table,    repo_path("output", "tables", "SCD_ACA_BH_pvalues.csv"))
write_csv(ftest_table,   repo_path("output", "tables", "SCD_ACA_JointFtest_Results.csv"))
cat("\nResults saved to SCD_ACA_DiD_Results.csv, SCD_ACA_BH_pvalues.csv,",
    "and SCD_ACA_JointFtest_Results.csv\n")

# =============================================================================
# 14. RESIDUAL DISTRIBUTION CHECK
# LOS is prone to a rightward skewed distribution; the values provided are
# the mean LOS per stratum, which should improve this skew. Checking with a
# residual distribution / normality test below.
# =============================================================================

# Check residual distribution from current LOS model
par(mfrow = c(1, 2))
hist(residuals(m_los$model),
     main = "LOS Model Residuals",
     xlab = "Residual",
     breaks = 20)
qqnorm(residuals(m_los$model))
qqline(residuals(m_los$model), col = "firebrick")

shapiro.test(residuals(m_los$model))

# =============================================================================
# 14B. SENSITIVITY: LOG-TRANSFORMED LOS
# =============================================================================
# 1. Create the log-transformed LOS variable
# (LOS is always > 0 for an admission, so a plain log() is safe here —
#  no zeros to worry about, unlike e.g. self-pay proportions)
df_hosp <- df_hosp %>%
  mutate(log_los = log(los_per_admission_mean))

# 2. Refit using run_its() — same weighting as the original model
m_los_log <- run_its(
  "log_los",
  "Log(LOS per admission)",
  df_hosp,
  weight_var = "admissions_total"   # matches the original m_los call
)

# 3. Check whether the residuals actually look better on the log scale
shapiro.test(residuals(m_los_log$model))

hist(residuals(m_los_log$model),
     main = "Log(LOS) Model Residuals", xlab = "Residual", breaks = 20)
qqnorm(residuals(m_los_log$model))
qqline(residuals(m_los_log$model), col = "firebrick")

# 4. Joint F-test and coefficient extraction for the log-scale model
ft_los_log    <- joint_did_test(m_los_log, "LOS per admission (log)")
coefs_los_log <- extract_did_coefs(m_los_log, "LOS per admission (log)")
print(ft_los_log)
print(coefs_los_log)

# 5. Back-transform log-scale DiD coefficients to interpretable % change
#    (exp(estimate) - 1) * 100, with a Wald-style 95% CI
log_to_pct_change <- function(estimate, se, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  tibble(
    pct_change = round((exp(estimate) - 1) * 100, 2),
    ci_lower   = round((exp(estimate - z * se) - 1) * 100, 2),
    ci_upper   = round((exp(estimate + z * se) - 1) * 100, 2)
  )
}

los_log_pct_change <- bind_cols(
  coefs_los_log,
  pmap_dfr(list(coefs_los_log$estimate, coefs_los_log$se), log_to_pct_change)
)

cat("\n=== LOG(LOS) DiD COEFFICIENTS, BACK-TRANSFORMED TO % CHANGE ===\n")
print(los_log_pct_change, n = Inf)

write_csv(coefs_los_log,      repo_path("output", "tables", "SCD_ACA_LOS_Log_Sensitivity_Coefs.csv"))
write_csv(los_log_pct_change, repo_path("output", "tables", "SCD_ACA_LOS_Log_Sensitivity_PctChange.csv"))

# =============================================================================
# END OF SCRIPT
# =============================================================================
