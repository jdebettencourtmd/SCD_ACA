# =============================================================================
# RUN ALL ANALYSES
#
# Reproduces every model, table, and figure in the manuscript, and records the
# exact software environment used.
#
#   Rscript run_all.R
#
# Console output is teed to output/run_log.txt so the full numeric output of
# each model is preserved alongside the summary tables.
# =============================================================================

repo_root <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()

dir.create(file.path(repo_root, "output"), showWarnings = FALSE)
log_con <- file(file.path(repo_root, "output", "run_log.txt"), open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("SCD / ACA DiD analysis — full run\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

source(file.path(repo_root, "R", "01_primary_ED_analysis.R"))
source(file.path(repo_root, "R", "02_payer_mix_analysis.R"))
source(file.path(repo_root, "R", "03_age_stratified_analysis.R"))

cat("\n\n", strrep("=", 70), "\n")
cat("SESSION INFO\n")
cat(strrep("=", 70), "\n")
print(sessionInfo())

sink(type = "message")
sink()
close(log_con)

# Also write sessionInfo on its own, for the manuscript's software statement.
writeLines(capture.output(sessionInfo()),
           file.path(repo_root, "output", "sessionInfo.txt"))

cat("\nDone. Figures in output/figures/, tables in output/tables/,",
    "full log in output/run_log.txt\n")
