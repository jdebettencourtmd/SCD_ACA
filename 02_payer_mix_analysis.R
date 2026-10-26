# =============================================================================
# PAYER MIX SENSITIVITY ANALYSIS: INPATIENT-ONLY AND COMBINED (ED + INPATIENT)
# SCD Patients, California vs. Georgia, 2010-2018
#
#
# Produces Weighted payer mix trends for:
#        a) Inpatient admissions only 
#        b) Combined ED + inpatient encounters
# AND ITS/DiD models for payer mix in each data stream, with joint F-tests
#
# Weighting rationale:
#   - Inpatient: weighted by admissions_total (encounter-level weight)
#   - Combined:  weighted by (ed_visits_total + admissions_total) per stratum
#     so that each encounter contributes proportionally regardless of data type.
# 
# Requires 00_setup_functions.R to have been sourced first (provides df_ed,
# df_hosp, its_formula(), run_its(), joint_did_test(), get_did_pval(),
# theme_scd, aca_line, aca_label, state_colors).
# =============================================================================

# Source the shared setup. Works whether R was started at the repository root
# or elsewhere (if the `here` package is installed).
source(if (requireNamespace("here", quietly = TRUE)) {
  here::here("R", "00_setup_functions.R")
} else "R/00_setup_functions.R")

# -----------------------------------------------------------------------------
# 2. INPATIENT PAYER MIX
# -----------------------------------------------------------------------------


trend_hosp_payer <- df_hosp %>%
  group_by(state, year) %>%
  summarise(
    wt_medicaid = weighted.mean(payer_medicaid_prop, admissions_total),
    wt_medicare = weighted.mean(payer_medicare_prop, admissions_total),
    wt_private  = weighted.mean(payer_private_prop,  admissions_total),
    wt_selfpay  = weighted.mean(payer_selfpay_prop,  admissions_total),
    wt_other    = weighted.mean(payer_other_prop,    admissions_total),
    total_admissions = sum(admissions_total),
    .groups = "drop"
  )

cat("\n=== INPATIENT PAYER MIX BY YEAR AND STATE ===\n")
print(trend_hosp_payer, n = Inf)

# Pre/post summary
desc_hosp_payer <- df_hosp %>%
  mutate(period = if_else(year < 2014,
                          "Pre-ACA (2010\u20132013)",
                          "Post-ACA (2014\u20132018)")) %>%
  group_by(state, period) %>%
  summarise(
    total_admissions = sum(admissions_total),
    medicaid_pct = weighted.mean(payer_medicaid_prop, admissions_total),
    medicare_pct = weighted.mean(payer_medicare_prop, admissions_total),
    private_pct  = weighted.mean(payer_private_prop,  admissions_total),
    selfpay_pct  = weighted.mean(payer_selfpay_prop,  admissions_total),
    other_pct    = weighted.mean(payer_other_prop,    admissions_total),
    .groups = "drop"
  )

cat("\n=== INPATIENT PAYER MIX: Pre vs. Post ACA Summary ===\n")
print(desc_hosp_payer)

# -----------------------------------------------------------------------------
# 3. COMBINED (ED + INPATIENT) PAYER MIX
# -----------------------------------------------------------------------------


df_combined <- df_ed %>%
  select(state, year, age_group, cluster_id, post_aca, california,
         time, time_post, people_with_scd_total,
         ed_visits_total,
         ed_medicaid = payer_medicaid_prop,
         ed_medicare = payer_medicare_prop,
         ed_private  = payer_private_prop,
         ed_selfpay  = payer_selfpay_prop,
         ed_other    = payer_other_prop) %>%
  left_join(
    df_hosp %>%
      select(state, year, age_group,
             admissions_total,
             hosp_medicaid = payer_medicaid_prop,
             hosp_medicare = payer_medicare_prop,
             hosp_private  = payer_private_prop,
             hosp_selfpay  = payer_selfpay_prop,
             hosp_other    = payer_other_prop),
    by = c("state", "year", "age_group")
  ) %>%
  mutate(
    # Step 1: absolute encounter counts per payer per data type
    n_ed_medicaid   = ed_medicaid * ed_visits_total,
    n_ed_medicare   = ed_medicare * ed_visits_total,
    n_ed_private    = ed_private  * ed_visits_total,
    n_ed_selfpay    = ed_selfpay  * ed_visits_total,
    n_ed_other      = ed_other    * ed_visits_total,

    n_hosp_medicaid = hosp_medicaid * admissions_total,
    n_hosp_medicare = hosp_medicare * admissions_total,
    n_hosp_private  = hosp_private  * admissions_total,
    n_hosp_selfpay  = hosp_selfpay  * admissions_total,
    n_hosp_other    = hosp_other    * admissions_total,

    # Step 2: combined counts
    total_encounters = ed_visits_total + admissions_total,
    n_tot_medicaid   = n_ed_medicaid + n_hosp_medicaid,
    n_tot_medicare   = n_ed_medicare + n_hosp_medicare,
    n_tot_private    = n_ed_private  + n_hosp_private,
    n_tot_selfpay    = n_ed_selfpay  + n_hosp_selfpay,
    n_tot_other      = n_ed_other    + n_hosp_other,

    # Step 3: combined proportions (encounter-weighted)
    comb_medicaid    = n_tot_medicaid / total_encounters,
    comb_medicare    = n_tot_medicare / total_encounters,
    comb_private     = n_tot_private  / total_encounters,
    comb_selfpay     = n_tot_selfpay  / total_encounters,
    comb_other       = n_tot_other    / total_encounters
  )

# Verify proportions sum to ~1 (sanity check)
payer_sums <- df_combined %>%
  mutate(check = comb_medicaid + comb_medicare + comb_private +
                 comb_selfpay  + comb_other) %>%
  summarise(min_sum = min(check), max_sum = max(check))
cat("\n=== PAYER PROPORTION SUM CHECK (should be ~1.00) ===\n")
print(payer_sums)

# Annual trend: weighted by total_encounters
trend_combined_payer <- df_combined %>%
  group_by(state, year) %>%
  summarise(
    wt_medicaid      = weighted.mean(comb_medicaid, total_encounters),
    wt_medicare      = weighted.mean(comb_medicare, total_encounters),
    wt_private       = weighted.mean(comb_private,  total_encounters),
    wt_selfpay       = weighted.mean(comb_selfpay,  total_encounters),
    wt_other         = weighted.mean(comb_other,    total_encounters),
    total_encounters = sum(total_encounters),
    .groups = "drop"
  )

cat("\n=== COMBINED (ED + INPATIENT) PAYER MIX BY YEAR AND STATE ===\n")
print(trend_combined_payer, n = Inf)

# Pre/post summary
desc_combined_payer <- df_combined %>%
  mutate(period = if_else(year < 2014,
                          "Pre-ACA (2010\u20132013)",
                          "Post-ACA (2014\u20132018)")) %>%
  group_by(state, period) %>%
  summarise(
    total_encounters = sum(total_encounters),
    medicaid_pct     = sum(comb_medicaid * total_encounters) / sum(total_encounters),
    medicare_pct     = sum(comb_medicare * total_encounters) / sum(total_encounters),
    private_pct      = sum(comb_private  * total_encounters) / sum(total_encounters),
    selfpay_pct      = sum(comb_selfpay  * total_encounters) / sum(total_encounters),
    other_pct        = sum(comb_other    * total_encounters) / sum(total_encounters),
    .groups = "drop"
  )

cat("\n=== COMBINED PAYER MIX: Pre vs. Post ACA Summary ===\n")
print(desc_combined_payer)

# -----------------------------------------------------------------------------
# 4. ITS/DiD MODELS: INPATIENT PAYER MIX
# -----------------------------------------------------------------------------
# Reference: Bernal et al. (2017) Int J Epidemiol; Dimick & Ryan (2014) JAMA
# NOTE: its_formula(), run_its(), joint_did_test() come from
# 00_setup_functions.R — not redefined here.

cat("\n\n", strrep("#", 70), "\n")
cat("# INPATIENT PAYER MIX ITS/DiD MODELS\n")
cat(strrep("#", 70), "\n")

m_hosp_medicaid <- run_its("payer_medicaid_prop", "Inpatient: Medicaid proportion",
                            df_hosp, "admissions_total")
m_hosp_private  <- run_its("payer_private_prop",  "Inpatient: Private insurance proportion",
                            df_hosp, "admissions_total")
m_hosp_selfpay  <- run_its("payer_selfpay_prop",  "Inpatient: Self-pay (uninsured) proportion",
                            df_hosp, "admissions_total")
m_hosp_medicare <- run_its("payer_medicare_prop", "Inpatient: Medicare proportion",
                            df_hosp, "admissions_total")

ft_hosp_medicaid <- joint_did_test(m_hosp_medicaid, "Inpatient Medicaid")
ft_hosp_private  <- joint_did_test(m_hosp_private,  "Inpatient Private insurance")
ft_hosp_selfpay  <- joint_did_test(m_hosp_selfpay,  "Inpatient Self-pay (uninsured)")
ft_hosp_medicare <- joint_did_test(m_hosp_medicare, "Inpatient Medicare (negative control)")

# -----------------------------------------------------------------------------
# 5. ITS/DiD MODELS: COMBINED PAYER MIX
# -----------------------------------------------------------------------------
cat("\n\n", strrep("#", 70), "\n")
cat("# COMBINED (ED + INPATIENT) PAYER MIX ITS/DiD MODELS\n")
cat(strrep("#", 70), "\n")

m_comb_medicaid <- run_its("comb_medicaid", "Combined: Medicaid proportion",
                            df_combined, "total_encounters")
m_comb_private  <- run_its("comb_private",  "Combined: Private insurance proportion",
                            df_combined, "total_encounters")
m_comb_selfpay  <- run_its("comb_selfpay",  "Combined: Self-pay (uninsured) proportion",
                            df_combined, "total_encounters")
m_comb_medicare <- run_its("comb_medicare", "Combined: Medicare proportion",
                            df_combined, "total_encounters")

ft_comb_medicaid <- joint_did_test(m_comb_medicaid, "Combined Medicaid")
ft_comb_private  <- joint_did_test(m_comb_private,  "Combined Private insurance")
ft_comb_selfpay  <- joint_did_test(m_comb_selfpay,  "Combined Self-pay (uninsured)")
ft_comb_medicare <- joint_did_test(m_comb_medicare, "Combined Medicare (negative control)")

# -----------------------------------------------------------------------------
# 6. MULTIPLE COMPARISONS: BH CORRECTION ACROSS PAYER-MIX OUTCOMES,
#    AND JOINT F-TEST SUMMARY EXPORT
# -----------------------------------------------------------------------------

pval_payer <- bind_rows(
  get_did_pval(m_hosp_medicaid, "Inpatient: Medicaid"),
  get_did_pval(m_hosp_private,  "Inpatient: Private insurance"),
  get_did_pval(m_hosp_selfpay,  "Inpatient: Self-pay (uninsured)"),
  get_did_pval(m_hosp_medicare, "Inpatient: Medicare"),
  get_did_pval(m_comb_medicaid, "Combined: Medicaid"),
  get_did_pval(m_comb_private,  "Combined: Private insurance"),
  get_did_pval(m_comb_selfpay,  "Combined: Self-pay (uninsured)"),
  get_did_pval(m_comb_medicare, "Combined: Medicare")
) %>%
  mutate(p_adj_BH = round(p.adjust(p_raw, method = "BH"), 4))

cat("\n=== PAYER MIX DiD: BH-adjusted p-values ===\n")
print(pval_payer, n = Inf)

ftest_payer <- bind_rows(
  ft_hosp_medicaid, ft_hosp_private, ft_hosp_selfpay, ft_hosp_medicare,
  ft_comb_medicaid, ft_comb_private, ft_comb_selfpay, ft_comb_medicare
)

cat("\n=== PAYER MIX DiD: Joint F-test summary ===\n")
print(ftest_payer, n = Inf)

write_csv(pval_payer,  repo_path("output", "tables", "SCD_PayerMix_Inpatient_Combined_DiD.csv"))
write_csv(ftest_payer, repo_path("output", "tables", "SCD_PayerMix_JointFtest_Results.csv"))

# -----------------------------------------------------------------------------
# 7. VISUALIZATION
# -----------------------------------------------------------------------------
# theme_scd, aca_line, aca_label come from 00_setup_functions.R. That base
# theme does not rotate x-axis labels by default, so figures below that need
# rotated year labels add it explicitly (mirrors Figure 7 in
# 01_primary_ED_analysis.R, which does the same).

payer_colors <- c(
  "Medicaid"            = "#2166ac",
  "Medicare"            = "#4dac26",
  "Private"             = "#d6604d",
  "Self-pay (Uninsured)"= "#f4a582",
  "Other"               = "#878787"
)

# ── Figure 9: Inpatient payer mix — stacked area ──────────────────────────────
hosp_payer_long <- trend_hosp_payer %>%
  select(state, year, wt_medicaid, wt_medicare,
         wt_private, wt_selfpay, wt_other) %>%
  pivot_longer(cols = starts_with("wt_"),
               names_to = "payer", names_prefix = "wt_",
               values_to = "proportion") %>%
  mutate(payer = case_match(payer,
    "medicaid" ~ "Medicaid",
    "medicare" ~ "Medicare",
    "private"  ~ "Private",
    "selfpay"  ~ "Self-pay (Uninsured)",
    "other"    ~ "Other"
  ),
  payer = factor(payer, levels = c("Medicaid","Private","Medicare",
                                    "Self-pay (Uninsured)","Other")))

fig9 <- ggplot(hosp_payer_long,
               aes(x = year, y = proportion, fill = payer)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.3) +
  facet_wrap(~state) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.001)) +
  scale_fill_manual(values = payer_colors) +
  labs(
    title    = "Inpatient Payer Mix Among SCD Patients, 2010\u20132018",
    subtitle = "Weighted by admissions total per stratum",
    x = "Year", y = "Proportion of Inpatient Admissions",
    fill    = "Payer",
    caption = "Dashed line = ACA implementation (Jan 2014)"
  ) +
  theme_scd +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1))
print(fig9)
ggsave(repo_path("output", "figures", "fig9_inpatient_payer_mix_stacked.png"),
       fig9, width = 10, height = 5, dpi = 300)

# ── Figure 10: Inpatient payer mix — line plot (mirrors ED figure style) ──────
fig10 <- ggplot(hosp_payer_long,
                aes(x = year, y = proportion,
                    color = payer, shape = payer)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = payer_colors) +
  labs(
    title    = "Inpatient Payer Mix Among SCD Patients, 2010\u20132018",
    subtitle = "Weighted by admissions total per stratum",
    x = "Year", y = "Proportion of Inpatient Admissions",
    color = "Payer", shape = "Payer",
    caption = "Dashed line = ACA implementation (Jan 2014)"
  ) +
  theme_scd +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))
print(fig10)
ggsave(repo_path("output", "figures", "fig10_inpatient_payer_mix_lines.png"),
       fig10, width = 10, height = 5, dpi = 300)

# ── Figure 11: Combined payer mix — stacked area ──────────────────────────────
comb_payer_long <- trend_combined_payer %>%
  select(state, year, wt_medicaid, wt_medicare,
         wt_private, wt_selfpay, wt_other) %>%
  pivot_longer(cols = starts_with("wt_"),
               names_to = "payer", names_prefix = "wt_",
               values_to = "proportion") %>%
  mutate(payer = case_match(payer,
    "medicaid" ~ "Medicaid",
    "medicare" ~ "Medicare",
    "private"  ~ "Private",
    "selfpay"  ~ "Self-pay (Uninsured)",
    "other"    ~ "Other"
  ),
  payer = factor(payer, levels = c("Medicaid","Private","Medicare",
                                    "Self-pay (Uninsured)","Other")))

fig11 <- ggplot(comb_payer_long,
                aes(x = year, y = proportion, fill = payer)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.3) +
  facet_wrap(~state) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.001)) +
  scale_fill_manual(values = payer_colors) +
  labs(
    title    = "Combined (ED + Inpatient) Payer Mix Among SCD Patients, 2010\u20132018",
    subtitle = "Weighted by total encounters (ED visits + admissions) per stratum",
    x = "Year", y = "Proportion of Total Encounters",
    fill    = "Payer",
    caption = "Dashed line = ACA implementation (Jan 2014)"
  ) +
  theme_scd +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1))
print(fig11)
ggsave(repo_path("output", "figures", "fig11_combined_payer_mix_stacked.png"),
       fig11, width = 10, height = 5, dpi = 300)

# ── Figure 12: Combined payer mix — line plot ─────────────────────────────────
fig12 <- ggplot(comb_payer_long,
                aes(x = year, y = proportion,
                    color = payer, shape = payer)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state) +
  aca_line + aca_label +
  scale_x_continuous(breaks = 2010:2018) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = payer_colors) +
  labs(
    title    = "Combined (ED + Inpatient) Payer Mix Among SCD Patients, 2010\u20132018",
    subtitle = "Weighted by total encounters (ED visits + admissions) per stratum",
    x = "Year", y = "Proportion of Total Encounters",
    color = "Payer", shape = "Payer",
    caption = "Dashed line = ACA implementation (Jan 2014)"
  ) +
  theme_scd +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))
print(fig12)
ggsave(repo_path("output", "figures", "fig12_combined_payer_mix_lines.png"),
       fig12, width = 10, height = 5, dpi = 300)

# ── Figure 13: Side-by-side comparison — ED vs. Inpatient vs. Combined ────────


three_stream <- bind_rows(
  # ED stream (rebuilt here from df_ed for self-containment of this figure)
  df_ed %>%
    group_by(state, year) %>%
    summarise(
      medicaid = weighted.mean(payer_medicaid_prop, people_with_scd_total),
      private  = weighted.mean(payer_private_prop,  people_with_scd_total),
      selfpay  = weighted.mean(payer_selfpay_prop,  people_with_scd_total),
      .groups  = "drop"
    ) %>%
    mutate(stream = "ED Encounters"),

  # Inpatient stream
  trend_hosp_payer %>%
    select(state, year,
           medicaid = wt_medicaid,
           private  = wt_private,
           selfpay  = wt_selfpay) %>%
    mutate(stream = "Inpatient Admissions"),

  # Combined stream
  trend_combined_payer %>%
    select(state, year,
           medicaid = wt_medicaid,
           private  = wt_private,
           selfpay  = wt_selfpay) %>%
    mutate(stream = "Combined (ED + Inpatient)")
) %>%
  pivot_longer(cols = c(medicaid, private, selfpay),
               names_to = "payer", values_to = "proportion") %>%
  mutate(
    payer  = case_match(payer,
      "medicaid" ~ "Medicaid",
      "private"  ~ "Private",
      "selfpay"  ~ "Self-pay (Uninsured)"
    ),
    stream = factor(stream,
                    levels = c("ED Encounters",
                               "Inpatient Admissions",
                               "Combined (ED + Inpatient)"))
  )

fig13 <- ggplot(three_stream,
                aes(x = year, y = proportion,
                    color = stream, linetype = stream)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.8) +
  facet_grid(payer ~ state, scales = "free_y") +
  aca_line +
  scale_x_continuous(breaks = c(2010, 2013, 2016)) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = c(
    "ED Encounters"            = "#2166ac",
    "Inpatient Admissions"     = "#d73027",
    "Combined (ED + Inpatient)"= "#4dac26"
  )) +
  scale_linetype_manual(values = c(
    "ED Encounters"            = "solid",
    "Inpatient Admissions"     = "dashed",
    "Combined (ED + Inpatient)"= "dotdash"
  )) +
  labs(
    title    = "Payer Mix Comparison Across Data Streams, SCD Patients 2010\u20132018",
    subtitle = "Rows = payer type; columns = state; lines = data stream",
    x = "Year", y = "Proportion",
    color = "Data Stream", linetype = "Data Stream",
    caption = "Dashed vertical line = ACA implementation (Jan 2014)"
  ) +
  theme_scd +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 10, face = "bold"),
    legend.position  = "bottom"
  ) +
  guides(color    = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))
print(fig13)
ggsave(repo_path("output", "figures", "fig13_payer_mix_stream_comparison.png"),
       fig13, width = 11, height = 9, dpi = 300)

cat("\n\nAll payer mix analyses and figures complete.\n")
cat("Figures saved: fig9 through fig13\n")
cat("Results saved: SCD_PayerMix_Inpatient_Combined_DiD.csv,",
    "SCD_PayerMix_JointFtest_Results.csv\n")

# =============================================================================
# END OF SCRIPT
# =============================================================================
