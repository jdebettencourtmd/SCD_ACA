# =============================================================================
# SCD / ACA DiD STUDY — SHARED SETUP AND CORE FUNCTIONS
# Title: Impact of the Affordable Care Act on Insurance Coverage and
#        Healthcare Utilization in Patients with Sickle Cell Disease:
#        A Comparison of Medicaid Expansion (California) vs.
#        Non-Expansion (Georgia) States, 2010-2018
#
# This file is sourced by both downstream analysis scripts and is not run
# directly. It is the single source of truth for:
#   - Package loading and raw data import
#   - Data preparation (prepare_df)
#   - Core ITS/DiD model-fitting infrastructure (its_formula, run_its,
#     joint_did_test, extract_did_coefs, get_did_pval)
#   - Shared plotting constants (theme_scd, aca_line, aca_label, state_colors)
#
#   source("00_setup_functions.R")
#   -> 01_primary_ED_analysis.R      (Table 3 / Appendix A: primary ED-based
#                                      insurance coverage & utilization models)
#   -> 02_payer_mix_analysis.R       (Table 2: inpatient-only and combined
#                                      ED+inpatient payer-mix sensitivity models)
#
# Design: Interrupted Time Series (ITS) with Difference-in-Differences (DiD)
# Unit of analysis: Age-group x State x Year (n = 126 per data type)
# ACA Implementation: January 1, 2014 (intervention year)
#
# Dataset: data/master_scd_utilization_GA_CA_2010_2018.csv
# Dataset is included in this repository and derived from
# CDC Sickle Cell Disease Demographic and Healthcare Utilization Data
# https://www.cdc.gov/sickle-cell-research/php/data/demographics-and-healthcare-utilization-data.html
# accessed March 15, 2026
#   - data_type == "ed"       -> ED visit analyses (126 rows)
#   - data_type == "hospital" -> Inpatient admission/LOS analyses (126 rows)
#
# R Code was created by Joseph deBettencourt, MD with use of Claude Sonnet 4.6
# (Anthropic, 2026) for code organization, data cleaning, and debugging some 
# comments in R Code may appear which are edits/comments made by AI debugging.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGE INSTALLATION & LOADING
# -----------------------------------------------------------------------------
# install.packages(c("tidyverse","lmtest","sandwich","nlme",
#                    "emmeans","broom","car","ggplot2","scales","here"))

library(tidyverse)
library(lmtest)
library(sandwich)
library(nlme)
library(emmeans)
library(broom)
library(car)
library(ggplot2)
library(scales)

# -----------------------------------------------------------------------------
# 1. PATH HANDLING, DATA IMPORT & PREPARATION
# -----------------------------------------------------------------------------
repo_path <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) here::here(...) else file.path(...)
}

# Output directories (tracked in the repository, contents git-ignored)
dir.create(repo_path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(repo_path("output", "tables"),  recursive = TRUE, showWarnings = FALSE)

df_raw <- read_csv(repo_path("data", "master_scd_utilization_GA_CA_2010_2018.csv"))

# ── Shared preparation function ───────────────────────────────────────────────

prepare_df <- function(data) {
  data %>%
    mutate(
      # Intervention coding
      post_aca   = as.integer(year >= 2014),
      # State coding: California = expansion (treatment), Georgia = control
      california = as.integer(state == "California"),
      # Time index: t=0 at 2010
      time       = as.integer(year - 2010),
      # Post-ACA slope counter
      time_post  = as.integer(if_else(year >= 2014, year - 2013, 0)),
      # Age group as ordered factor
      age_group  = factor(age_group,
                          levels = c("<10","10–19","20–29","30–39",
                                     "40–49","50–59","60+")),
      # Cluster ID for clustered SEs
      cluster_id = paste(state, age_group, sep = "_"),
      # Uninsured proxy (HCUP standard: self-pay = uninsured)
      uninsured_prop = payer_selfpay_prop,
      # Derived payer counts
      n_medicaid = round(payer_medicaid_prop * people_with_scd_total),
      n_medicare = round(payer_medicare_prop * people_with_scd_total),
      n_private  = round(payer_private_prop  * people_with_scd_total),
      n_selfpay  = round(payer_selfpay_prop  * people_with_scd_total)
    )
}

# ── Split into ED and hospital data frames ────────────────────────────────────
df_ed   <- df_raw %>% filter(data_type == "ed")       %>% prepare_df()
df_hosp <- df_raw %>% filter(data_type == "hospital") %>% prepare_df()

# Sanity checks
cat("ED rows:       ", nrow(df_ed),   "\n")
cat("Hospital rows: ", nrow(df_hosp), "\n")
table(df_ed$state,   df_ed$year)
table(df_hosp$state, df_hosp$year)


# -----------------------------------------------------------------------------
# 2. CORE ITS/DiD FUNCTIONS
# -----------------------------------------------------------------------------
# Model structure (Segmented regression ITS with DiD):
#
#   Y_it = b0 + b1*time + b2*post_aca + b3*time_post
#          + b4*california + b5*california:time
#          + b6*california:post_aca + b7*california:time_post
#          + b8*age_group + e_it
#
# Key coefficients:
#   b1 = pre-ACA annual trend in Georgia (reference/control)
#   b2 = level change at ACA in Georgia
#   b3 = slope change post-ACA in Georgia
#   b4 = baseline CA vs. GA difference (fixed effect)
#   b5 = differential pre-ACA trend CA vs. GA
#   b6 = DiD level change: CA additional effect at ACA  <-- PRIMARY
#   b7 = DiD slope change: CA additional post-ACA trend <-- PRIMARY
#
# Weights: people_with_scd_total by default; individual calls may override
#          (e.g. admissions_total for inpatient/LOS models, total_encounters
#          for combined ED+inpatient payer-mix models)
# SEs: clustered by state x age_group (Liang-Zeger)
# References:
#   Bernal et al. (2017) Int J Epidemiol
#   Dimick & Ryan (2014) JAMA
#   Bertrand et al. (2004) QJE
# -----------------------------------------------------------------------------

# ── ITS formula builder ───────────────────────────────────────────────────────
its_formula <- function(outcome) {
  as.formula(paste(
    outcome,
    "~ time + post_aca + time_post +",
    "california + california:time +",
    "california:post_aca + california:time_post +",
    "age_group"
  ))
}

# ── Main ITS/DiD fitting function ─────────────────────────────────────────────

run_its <- function(outcome, label, input_df, weight_var = "people_with_scd_total") {
  cat("\n\n", strrep("=", 70), "\n")
  cat("OUTCOME:", label, "\n")
  cat(strrep("=", 70), "\n")

  mod <- do.call("lm", list(
    formula = its_formula(outcome),
    data    = input_df,
    weights = input_df[[weight_var]]
  ))

  vcov_cl <- vcovCL(mod, cluster = ~cluster_id)
  ct      <- coeftest(mod, vcov = vcov_cl)
  print(ct)

  return(list(model = mod, coeftest = ct))
}

# ── Joint DiD F-test ──────────────────────────────────────────────────────────

joint_did_test <- function(model_result, label) {
  cat("\n=== Joint DiD F-test:", label, "===\n")
  mod      <- model_result$model
  vcov_cl  <- vcovCL(mod, cluster = ~cluster_id)
  cn       <- names(coef(mod))
  did_lvl  <- cn[grep("california.*post_aca|post_aca.*california", cn)]
  did_slp  <- cn[grep("california.*time_post|time_post.*california", cn)]
  cat("DiD terms detected:", did_lvl, "&", did_slp, "\n")

  lh <- linearHypothesis(
    mod,
    c(paste(did_lvl, "= 0"), paste(did_slp, "= 0")),
    vcov = vcov_cl
  )
  print(lh)

  tibble(
    outcome = label,
    f_stat  = round(lh$F[2], 4),
    df1     = lh$Df[2],
    df2     = lh$Res.Df[2],
    p_value = round(lh$`Pr(>F)`[2], 4)
  )
}

# ── Result extraction helper ──────────────────────────────────────────────────
extract_did_coefs <- function(model_result, outcome_label) {
  ct <- model_result$coeftest
  cn <- rownames(ct)
  # Dynamically find DiD coefficient rows
  ga_lvl  <- cn[grep("^post_aca$", cn)]
  ga_slp  <- cn[grep("^time_post$", cn)]
  did_lvl <- cn[grep("california.*post_aca|post_aca.*california", cn)]
  did_slp <- cn[grep("california.*time_post|time_post.*california", cn)]
  key_terms <- c(ga_lvl, ga_slp, did_lvl, did_slp)
  labels    <- c("GA: Level change at ACA",
                 "GA: Slope change post-ACA",
                 "DiD: CA level change vs. GA",
                 "DiD: CA slope change vs. GA")
  map2_dfr(key_terms, labels, function(term, lbl) {
    idx <- which(cn == term)
    tibble(
      outcome   = outcome_label,
      parameter = lbl,
      estimate  = round(ct[idx, 1], 4),
      se        = round(ct[idx, 2], 4),
      p_value   = round(ct[idx, 4], 4)
    )
  })
}

# ── BH-correction p-value extraction helper ───────────────────────────────────

get_did_pval <- function(model_result, outcome_label) {
  ct  <- model_result$coeftest
  cn  <- rownames(ct)
  idx <- grep("california.*post_aca|post_aca.*california", cn)
  tibble(
    outcome  = outcome_label,
    estimate = round(ct[idx, 1], 4),
    se       = round(ct[idx, 2], 4),
    p_raw    = round(ct[idx, 4], 4)
  )
}


# -----------------------------------------------------------------------------
# 3. SHARED PLOTTING CONSTANTS
# -----------------------------------------------------------------------------
theme_scd <- theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

aca_line  <- geom_vline(xintercept = 2013.5, linetype = "dashed",
                        color = "firebrick", linewidth = 0.8)
aca_label <- annotate("text", x = 2013.7, y = Inf, label = "ACA (2014)",
                      hjust = 0, vjust = 1.5, color = "firebrick", size = 3.5)

state_colors <- c("California" = "#2166ac", "Georgia" = "#d73027")

# =============================================================================
# END OF SETUP — proceed to 01_primary_ED_analysis.R and/or
# 02_payer_mix_analysis.R
# =============================================================================
