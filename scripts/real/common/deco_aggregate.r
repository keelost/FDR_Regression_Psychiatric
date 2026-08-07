#!/usr/bin/env Rscript
# ==================================================================
# Cross-Target Aggregation Report
# ==================================================================
# Description:
#   Scans the output-base of deco_validation.R, collects each target's
#   key_metrics and statistical_tests tables, and produces a combined
#   cross-target report (CSV + plain-text summary).
#
# Usage:
#   Rscript deco_aggregate.R \
#     --output-base /path/to/results/01.extra.analysis/01.deco_diagnosis \
#     --report-dir  <optional, defaults to <output-base>/00.aggregate_report>
# ==================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---- CLI ----
option_list <- list(
  make_option(c("-o", "--output-base"), type = "character", dest = "output_base",
              default = NULL,
              help = "Base directory containing per-target sub-folders. [REQUIRED]"),
  make_option(c("-r", "--report-dir"), type = "character", dest = "report_dir",
              default = NULL,
              help = "Directory for aggregate outputs. [default: <output-base>/00.aggregate_report]"),
  make_option(c("-t", "--targets"), type = "character", default = NULL,
              help = "Optional comma-separated target subset. Default: auto-detect all.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$output_base)) stop("--output-base is required.")

output_base <- opt$output_base
report_dir  <- if (is.null(opt$report_dir)) {
  file.path(output_base, "00.aggregate_report")
} else opt$report_dir
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Detect targets ----
detect_targets <- function(base) {
  subs <- list.dirs(base, recursive = FALSE, full.names = FALSE)
  # A valid target folder must contain tables/*_key_metrics.csv
  keep <- vapply(subs, function(s) {
    file.exists(file.path(base, s, "tables",
                          "decorrelation_validation_key_metrics.csv"))
  }, logical(1))
  subs[keep]
}

targets <- if (!is.null(opt$targets)) {
  trimws(strsplit(opt$targets, ",", fixed = TRUE)[[1]])
} else detect_targets(output_base)

if (length(targets) == 0) {
  stop("No valid target folders with key_metrics found under: ", output_base)
}

cat("============================================================\n")
cat("Cross-Target Aggregation\n")
cat("============================================================\n")
cat("Output base :", output_base, "\n")
cat("Report dir  :", report_dir, "\n")
cat("Targets     :", paste(targets, collapse = ", "), "\n\n")

# ---- Helper: safe read a target table, tagging the target name ----
read_target_table <- function(target, fname) {
  fpath <- file.path(output_base, target, "tables", fname)
  if (!file.exists(fpath)) {
    warning("Missing file for ", target, ": ", fpath)
    return(NULL)
  }
  dt <- fread(fpath)
  dt[, target := target]
  setcolorder(dt, c("target", setdiff(names(dt), "target")))
  dt[]
}

# ---- 1. Combine key metrics across targets ----
key_all <- rbindlist(
  lapply(targets, read_target_table,
         fname = "decorrelation_validation_key_metrics.csv"),
  use.names = TRUE, fill = TRUE)
fwrite(key_all, file.path(report_dir, "combined_key_metrics.csv"))

# ---- 2. Combine statistical tests across targets ----
test_all <- rbindlist(
  lapply(targets, read_target_table,
         fname = "decorrelation_validation_statistical_tests.csv"),
  use.names = TRUE, fill = TRUE)
fwrite(test_all, file.path(report_dir, "combined_statistical_tests.csv"))

# ---- 3. Combine per-covariate summaries across targets ----
covsum_all <- rbindlist(
  lapply(targets, read_target_table,
         fname = "decorrelation_validation_covariate_summary.csv"),
  use.names = TRUE, fill = TRUE)
fwrite(covsum_all, file.path(report_dir, "combined_covariate_summary.csv"))

# ---- 4. Per-target headline metrics (use the p>0 / full subset row) ----
# Prefer subset "p>0" as the representative full-data row; fallback to first.
pick_row <- function(dt_target) {
  r <- dt_target[subset == "p>0"]
  if (nrow(r) == 0) r <- dt_target[1]
  r
}
headline <- rbindlist(lapply(targets, function(tg) {
  dt <- key_all[target == tg]
  if (nrow(dt) == 0) return(NULL)
  pick_row(dt)
}), use.names = TRUE, fill = TRUE)
fwrite(headline, file.path(report_dir, "headline_per_target.csv"))

# ---- 5. Decorrelation effectiveness ranking ----
# Lower mean_abs_decor_cor and higher mean_reduction_pct = better decorrelation.
if (nrow(headline) > 0) {
  rank_dt <- headline[, .(
    target,
    n_snps,
    mean_abs_orig_cor,
    mean_abs_decor_cor,
    mean_reduction_pct,
    orig_vs_decor_r,
    theory_r,
    abs_diff_from_theory
  )][order(-mean_reduction_pct)]
  fwrite(rank_dt, file.path(report_dir, "decorrelation_ranking.csv"))
}

# ---- 6. Significant residual correlations after decorrelation ----
# Flag covariate pairs whose decorrelated target-covariate correlation is
# still significant (potential incomplete decorrelation).
if (nrow(test_all) > 0) {
  residual_sig <- test_all[test_type == "decor_target_vs_cov" &
                             significant == "Yes"]
  fwrite(residual_sig,
         file.path(report_dir, "residual_significant_after_decor.csv"))
} else {
  residual_sig <- data.table()
}

# ---- 7. Plain-text report ----
report_path <- file.path(report_dir, "aggregate_report.txt")
con <- file(report_path, open = "wt")
writeLines(c(
  "============================================================",
  "CROSS-TARGET DECORRELATION VALIDATION REPORT",
  "============================================================",
  paste("Generated       :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Output base     :", output_base),
  paste("Targets analysed:", length(targets)),
  paste("Target list     :", paste(targets, collapse = ", ")),
  ""
), con)

if (exists("rank_dt") && nrow(rank_dt) > 0) {
  writeLines(c(
    "------------------------------------------------------------",
    "1. Decorrelation effectiveness ranking (by mean_reduction_pct)",
    "------------------------------------------------------------"
  ), con)
  writeLines(capture.output(print(rank_dt)), con)
  writeLines("", con)
}

writeLines(c(
  "------------------------------------------------------------",
  "2. Residual significant target-covariate correlations",
  "   (decor_target_vs_cov, p < 0.05 -- incomplete decorrelation)",
  "------------------------------------------------------------",
  paste("Total flagged rows:", nrow(residual_sig))
), con)
if (nrow(residual_sig) > 0) {
  brief <- residual_sig[, .(target, subset, covariate,
                            correlation, p_value)]
  writeLines(capture.output(print(brief)), con)
}
writeLines("", con)

writeLines(c(
  "------------------------------------------------------------",
  "3. Output files in this report directory",
  "------------------------------------------------------------",
  "- combined_key_metrics.csv",
  "- combined_statistical_tests.csv",
  "- combined_covariate_summary.csv",
  "- headline_per_target.csv",
  "- decorrelation_ranking.csv",
  "- residual_significant_after_decor.csv",
  "- aggregate_report.txt (this file)"
), con)
close(con)

# ---- Console echo ----
cat("Aggregation complete.\n")
cat("Combined tables and report saved to:", report_dir, "\n\n")
if (exists("rank_dt") && nrow(rank_dt) > 0) {
  cat("---------- Decorrelation ranking ----------\n")
  print(rank_dt)
}
cat(sprintf("\nResidual significant correlations flagged: %d rows\n",
            nrow(residual_sig)))
