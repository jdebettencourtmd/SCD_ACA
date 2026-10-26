## =============================================================================
## Age-Stratified Effect Modification Analysis
## SCD / ACA Medicaid Expansion DiD Study (CA vs GA, 2010-2018)
##
## Purpose: Test whether the DiD effect on Medicaid-covered encounter
## proportion differs by age band, with the a priori hypothesis that the
## "young adult" band (proxying the ACA dependent-coverage-eligible / newly
## Medicaid-eligible adult population) shows the largest relative increase.

## ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c("tidyverse", "geepack", "broom", "MASS")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(tidyverse)
library(geepack)  
library(broom)

## ---- 1. PATHS AND CONFIGURATION ---------------------------------------------

repo_path <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) here::here(...) else file.path(...)
}
dir.create(repo_path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(repo_path("output", "tables"),  recursive = TRUE, showWarnings = FALSE)

data_path   <- repo_path("data", "master_scd_utilization_GA_CA_2010_2018.csv")

col_state       <- "state"                      
col_year        <- "year"
col_age_group   <- "age_group"                  
col_outcome     <- "medicaid_prop_combined"      
col_weight      <- "people_with_scd_total"       

treatment_state <- "California"
aca_year        <- 2014   # first post-intervention year

## ---- 2. Load and inspect ---------------------------------------------------
dat_raw <- read_csv(data_path)
glimpse(dat_raw)   

dat_combined <- dat_raw %>%
  select(state, age_group, year, data_type, people_with_scd_total,
         admissions_total, ed_visits_total, payer_medicaid_prop) %>%
  pivot_wider(
    id_cols = c(state, age_group, year, people_with_scd_total),
    names_from = data_type,
    values_from = c(admissions_total, ed_visits_total, payer_medicaid_prop)
  ) %>%
  mutate(
    medicaid_prop_combined =
      (payer_medicaid_prop_hospital * admissions_total_hospital +
         payer_medicaid_prop_ed * ed_visits_total_ed) /
      (admissions_total_hospital + ed_visits_total_ed)
  )

stopifnot(!any(is.na(dat_combined$medicaid_prop_combined)))

dat <- dat_combined %>%
  rename(
    state     = all_of(col_state),
    year      = all_of(col_year),
    age_group = all_of(col_age_group),
    outcome   = all_of(col_outcome),
    weight    = all_of(col_weight)
  ) %>%
  mutate(
    state_CA   = as.integer(state == treatment_state),
    time       = year - min(year),                       
    post       = as.integer(year >= aca_year),
    time_after = pmax(year - aca_year, 0)                 # 0 pre-ACA, 0/1/2/3/4 post
  )

## ---- 3. Define age bands ----------------------------------------------------

dat <- dat %>%
  mutate(age_group = str_replace_all(as.character(age_group), "–", "-"))

age_band_map <- c(
  "<10"    = "Pediatric",
  "10-19"  = "Pediatric",
  "20-29"  = "Young Adult",
  "30-39"  = "Older Adult",
  "40-49"  = "Older Adult",
  "50-59"  = "Older Adult",
  "60+"    = "Older Adult"
)


dat <- dat %>%
  mutate(
    age_band = factor(
      recode(as.character(age_group), !!!age_band_map),
      levels = c("Pediatric", "Young Adult", "Older Adult")   
    )
  ) %>%
  arrange(state, age_group, year)   

stopifnot(!any(is.na(dat$age_band))) 

## ---- 4. Sanity check: parallel pre-trends WITHIN each age band? --------------

pretrend_check <- dat %>%
  filter(year < aca_year) %>%
  group_by(age_band) %>%
  group_modify(~ {
    m <- lm(outcome ~ time * state_CA, data = .x, weights = weight)
    tidy(m) %>% filter(term == "time:state_CA")
  }) %>%
  ungroup()

print(pretrend_check)


## ---- 5. Primary model: full age-band interaction ---------------------------


dat <- dat %>%
  mutate(cluster_id = interaction(state, age_group, drop = TRUE))

model_full <- geeglm(
  outcome ~ (time + post + time_after + state_CA +
               state_CA:time + state_CA:post + state_CA:time_after) * age_band,
  id        = cluster_id,
  waves     = year,
  weights   = weight,
  data      = dat,
  family    = gaussian(link = "identity"),
  corstr    = "exchangeable"     # try "ar1" as a sensitivity check, matching
)                                 # your existing GLS AR(1) sensitivity approach

summary(model_full)
names(coef(model_full))   # <- inspect exact coefficient names before building
#    the contrast vectors in Section 7

## ---- 5b. Sanity check: is vcov(model_full) actually full rank? -------------

V_full_sym <- (vcov(model_full) + t(vcov(model_full))) / 2
p_full <- ncol(V_full_sym)
n_clusters_full <- length(unique(dat$cluster_id))
rank_full <- sum(eigen(V_full_sym, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
cat(sprintf(
  "\n[Sanity check] model_full: %d parameters estimated from %d clusters; vcov(model_full) rank = %d of %d.\n",
  p_full, n_clusters_full, rank_full, p_full))
if (rank_full < p_full) {
  warning(sprintf(
    paste("model_full's sandwich covariance is rank-deficient (rank %d of %d parameters,",
          "only %d clusters). Any lincom() contrast that falls along a near-null direction of",
          "this matrix will report an artificially tiny or NaN SE -- such a row is NOT a precise",
          "estimate and must be flagged/excluded, not reported at face value. See the",
          "`reliable` column added to band_did/ya_contrasts below."),
    rank_full, p_full, n_clusters_full))
}


## ---- 6. Age-band-specific DiD estimates (level and slope) ------------------

lincom <- function(model, contrast_vec, conf_level = 0.95) {
  cf <- coef(model)
  V  <- vcov(model)
  full <- setNames(rep(0, length(cf)), names(cf))
  missing_terms <- setdiff(names(contrast_vec), names(cf))
  if (length(missing_terms) > 0) {
    stop("Contrast references terms not in model: ",
         paste(missing_terms, collapse = ", "),
         "\nCheck names(coef(model)) for exact spelling.")
  }
  full[names(contrast_vec)] <- contrast_vec
  est <- as.numeric(full %*% cf)
  var_est <- as.numeric(t(full) %*% V %*% full)
  reliable <- isTRUE(var_est > 1e-10)
  se  <- suppressWarnings(sqrt(var_est))
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  mult <- qnorm(1 - (1 - conf_level) / 2)
  tibble(estimate = est, se = se,
         conf.low = est - mult * se, conf.high = est + mult * se,
         statistic = z, p.value = p, reliable = reliable)
}


band_did <- bind_rows(
  Pediatric_level   = lincom(model_full, c("post:state_CA" = 1)),
  Pediatric_slope   = lincom(model_full, c("time_after:state_CA" = 1)),
  YoungAdult_level  = lincom(model_full, c("post:state_CA" = 1,
                                           "post:state_CA:age_bandYoung Adult" = 1)),
  YoungAdult_slope  = lincom(model_full, c("time_after:state_CA" = 1,
                                           "time_after:state_CA:age_bandYoung Adult" = 1)),
  OlderAdult_level  = lincom(model_full, c("post:state_CA" = 1,
                                           "post:state_CA:age_bandOlder Adult" = 1)),
  OlderAdult_slope  = lincom(model_full, c("time_after:state_CA" = 1,
                                           "time_after:state_CA:age_bandOlder Adult" = 1)),
  .id = "band_parameter"
)

print(band_did)

## ---- 7. A priori contrasts: is Young Adult the largest? --------------------


ya_contrasts <- bind_rows(
  YA_vs_Pediatric_level = lincom(model_full, c("post:state_CA:age_bandYoung Adult" = 1)),
  YA_vs_Pediatric_slope = lincom(model_full, c("time_after:state_CA:age_bandYoung Adult" = 1)),
  YA_vs_OlderAdult_level = lincom(model_full, c(
    "post:state_CA:age_bandYoung Adult" = 1,
    "post:state_CA:age_bandOlder Adult" = -1)),
  YA_vs_OlderAdult_slope = lincom(model_full, c(
    "time_after:state_CA:age_bandYoung Adult" = 1,
    "time_after:state_CA:age_bandOlder Adult" = -1)),
  .id = "contrast"
) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH"))  

print(ya_contrasts)

## ---- 8. Robustness check: fully separate per-band models -------------------


separate_band_models <- dat %>%
  group_by(age_band) %>%
  group_modify(~ {
    cluster_id_band <- droplevels(interaction(.x$state, .x$age_group, drop = TRUE))
    n_clusters_band <- length(unique(cluster_id_band))
    m <- geeglm(
      outcome ~ time + post + time_after + state_CA +
        state_CA:time + state_CA:post + state_CA:time_after,
      id      = cluster_id_band,
      waves   = .x$year,
      weights = .x$weight,
      data    = .x,
      family  = gaussian(link = "identity"),
      corstr  = "exchangeable"
    )
    # Same canonical-ordering note as Section 7: fitted names come out as
    # "post:state_CA" / "time_after:state_CA", not "state_CA:post"/"state_CA:time_after".
    tidy(m) %>%
      filter(term %in% c("post:state_CA", "time_after:state_CA")) %>%
      mutate(n_clusters = n_clusters_band,
             reliable = n_clusters_band >= 4)   # rule of thumb floor for GEE sandwich variance
  }) %>%
  ungroup()

print(separate_band_models)

## ---- 9. Forest plot: DiD estimates by age band -----------------------------
plot_dat <- band_did %>%
  mutate(
    band = str_extract(band_parameter, "^[A-Za-z]+"),
    param = str_extract(band_parameter, "level|slope"),
    band = factor(band, levels = c("Pediatric", "YoungAdult", "OlderAdult"))
  )

ggplot(plot_dat, aes(x = band, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = "#2c3e50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ param, scales = "free_y",
             labeller = as_labeller(c(level = "DiD Level Change",
                                      slope = "DiD Slope Change"))) +
  labs(
    title = "Medicaid Proportion: DiD Estimate by Age Band",
    subtitle = "California vs Georgia, relative to Pediatric reference band",
    x = NULL, y = "Estimate (percentage points, 95% CI)"
  ) +
  theme_minimal(base_size = 12)

ggsave(repo_path("output", "figures", "figS1_did_by_age_band_forest_plot.png"),
       width = 9, height = 6, dpi = 150)

## ---- 10. Export result tables -----------------------------------------------
write_csv(band_did,             repo_path("output", "tables", "SCD_AgeBand_DiD_by_band.csv"))
write_csv(ya_contrasts,         repo_path("output", "tables", "SCD_AgeBand_YoungAdult_contrasts.csv"))
write_csv(separate_band_models, repo_path("output", "tables", "SCD_AgeBand_separate_models.csv"))
write_csv(pretrend_check,       repo_path("output", "tables", "SCD_AgeBand_pretrend_check.csv"))


