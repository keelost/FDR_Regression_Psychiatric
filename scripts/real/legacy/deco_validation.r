#!/usr/bin/env Rscript
# ==================================================================
# Decorrelation Validation Script (multi-target, parallel)
# Adapted for the fdrreg.R v2.0 CLI pipeline (harmonised, nolasso naming)
# ==================================================================
# Usage:
#   Rscript deco_validation.R \
#     --targets "adhd2019,scz2018,bd2012" \
#     --output-base /exeh_4/jinghong_qiu/SO_Lab/18.fdrreg_rebuild/01.extra.analysis/01.deco_diagnosis \
#     --max 4
#
# Input  (per target, from fdrreg.R v2.0):
#   <pipeline_base>/<target>/01.magma_input/<target>.overlap.4magma.txt
#   <pipeline_base>/<target>/01.magma_input/<trait>.overlap.4magma.txt
#   <pipeline_base>/<target>/02.fdrreg_results/fdrreg_model_info_nolasso.rds
#
#   NOTE: fdrreg.R v2.0 harmonises every trait to the target effect allele
#   BEFORE de-correlation. The harmonised (pre-decorrelation) z-scores are
#   stored in the MAGMA input files under the column "z"; the de-correlated
#   z-scores are stored under "z.decor". This validation rebuilds the raw
#   z-matrix from the harmonised "z" column so it matches exactly what
#   fdrreg.R feeds into Matpow(cov, -0.5).
#
# Output (per target):
#   <output_base>/<target>/tables/*.csv
#   <output_base>/<target>/plots/*.pdf, *.png
# ==================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(powerplus)
  library(ggplot2)
  library(scales)
  library(parallel)
})

# ==================================================================
# 0. COMMAND-LINE ARGUMENT PARSING
# ==================================================================

option_list <- list(
  make_option(c("-t", "--targets"), type = "character", default = NULL,
              help = "Comma-separated target names (e.g. adhd2019,scz2018). [REQUIRED]"),
  make_option(c("-o", "--output-base"), type = "character", dest = "output_base",
              default = NULL,
              help = "Base output directory; each target gets its own sub-folder. [REQUIRED]"),
  make_option(c("-p", "--pipeline-base"), type = "character", dest = "pipeline_base",
              default = "/exeh_4/jinghong_qiu/SO_Lab/18.fdrreg_rebuild",
              help = "Root directory of the fdrreg.R pipeline output. [default: %default]"),
  make_option(c("-q", "--p-thresholds"), type = "character", dest = "p_thresholds",
              default = "0,0.05,0.1,0.2,0.3",
              help = "Comma-separated p-value thresholds for null-SNP subsets. [default: %default]"),
  make_option(c("-m", "--max"), type = "integer", default = 4,
              help = "Max cores for parallel target processing. [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Multi-target decorrelation validation for fdrreg.R v2.0 pipeline"))

if (is.null(opt$targets))     stop("--targets is required. Example: --targets adhd2019,scz2018")
if (is.null(opt$output_base)) stop("--output-base is required.")

# Helper: split comma-separated argument into a clean character vector
split_csv_arg <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nchar(vals) > 0]
}

target_list   <- split_csv_arg(opt$targets)
pipeline_base <- opt$pipeline_base
output_base   <- opt$output_base
p_thresholds  <- as.numeric(split_csv_arg(opt$p_thresholds))
max_cores     <- opt$max

cat("============================================================\n")
cat("Multi-Target Decorrelation Validation\n")
cat("============================================================\n")
cat("Targets      :", paste(target_list, collapse = ", "), "\n")
cat("Pipeline base:", pipeline_base, "\n")
cat("Output base  :", output_base, "\n")
cat("P-thresholds :", paste(p_thresholds, collapse = ", "), "\n")
cat("Max cores    :", max_cores, "\n\n")

# ==================================================================
# 2. VALIDATION FUNCTIONS
# ==================================================================

# 2a. Theoretical expectation: Cor(orig_target, decor_target) under null
get_theoretical_orig_vs_decor <- function(cor.matr, target_row = 1) {
  S_half <- Matpow(cor.matr, 0.5)
  S_half[target_row, target_row]
}

# 2b. Fisher transformation for correlation inference
fisher_test_corr <- function(r, n, rho0 = 0) {
  if (n <= 3) stop("n must be > 3")
  z_obs  <- atanh(r)
  z_exp  <- atanh(rho0)
  se     <- 1 / sqrt(n - 3)
  z_stat <- (z_obs - z_exp) / se
  p_val  <- 2 * pnorm(-abs(z_stat))
  ci_low  <- tanh(z_obs - 1.96 * se)
  ci_high <- tanh(z_obs + 1.96 * se)
  list(r = r, n = n, rho0 = rho0, z_stat = z_stat, p_val = p_val,
       ci_low = ci_low, ci_high = ci_high)
}

# 2c. Core validation on a SNP subset
compute_validation_stats <- function(z.mat, cor.matr,
                                     target_row = 1, sub_idx = NULL) {
  if (is.null(sub_idx)) sub_idx <- seq_len(ncol(z.mat))

  z_sub   <- z.mat[, sub_idx, drop = FALSE]
  n_snps  <- ncol(z_sub)

  invsqrtS <- Matpow(cor.matr, -0.5)
  z_decor  <- invsqrtS %*% z_sub

  all_rows <- seq_len(nrow(z.mat))
  cov_rows <- setdiff(all_rows, target_row)

  orig_target  <- z_sub[target_row, ]
  orig_covs    <- z_sub[cov_rows, , drop = FALSE]
  decor_target <- z_decor[target_row, ]
  decor_covs   <- z_decor[cov_rows, , drop = FALSE]

  orig_cor <- vapply(seq_len(nrow(orig_covs)), function(i) {
    cor(orig_target, orig_covs[i, ])
  }, numeric(1))
  names(orig_cor) <- rownames(z.mat)[cov_rows]

  decor_cor <- vapply(seq_len(nrow(decor_covs)), function(i) {
    cor(decor_target, decor_covs[i, ])
  }, numeric(1))
  names(decor_cor) <- rownames(z.mat)[cov_rows]

  orig_decor_cor   <- cor(orig_target, decor_target)
  theo_orig_decor  <- get_theoretical_orig_vs_decor(cor.matr, target_row)

  list(
    n_snps               = n_snps,
    orig_target_vs_cov   = orig_cor,
    decor_target_vs_cov  = decor_cor,
    orig_vs_decor_target = orig_decor_cor,
    theoretical_orig_decor = theo_orig_decor
  )
}

# ==================================================================
# 3. TABLE-GENERATION FUNCTIONS
# ==================================================================

create_decorrelation_summary_table <- function(results) {
  rbindlist(lapply(names(results), function(nm) {
    res <- results[[nm]]
    data.table(
      subset             = nm,
      n_snps             = res$n_snps,
      orig_vs_decor_target = round(res$orig_vs_decor_target, 4),
      theoretical_expected = round(res$theoretical_orig_decor, 4),
      diff_from_theory   = round(res$orig_vs_decor_target -
                                   res$theoretical_orig_decor, 4),
      abs_diff           = round(abs(res$orig_vs_decor_target -
                                       res$theoretical_orig_decor), 4)
    )
  }))
}

create_covariate_correlation_table <- function(results, trait_names) {
  cov_names <- trait_names[-1]
  rbindlist(lapply(names(results), function(nm) {
    res <- results[[nm]]
    rbindlist(lapply(seq_along(cov_names), function(i) {
      orig_abs  <- abs(res$orig_target_vs_cov[i])
      decor_abs <- abs(res$decor_target_vs_cov[i])
      reduction_pct <- if (orig_abs > 0) {
        round(100 * (1 - decor_abs / orig_abs), 2)
      } else NA_real_
      data.table(
        subset              = nm,
        covariate           = cov_names[i],
        n_snps              = res$n_snps,
        orig_target_vs_cov  = round(res$orig_target_vs_cov[i], 4),
        decor_target_vs_cov = round(res$decor_target_vs_cov[i], 4),
        abs_orig            = round(orig_abs, 4),
        abs_decor           = round(decor_abs, 4),
        reduction           = round(orig_abs - decor_abs, 4),
        reduction_pct       = reduction_pct
      )
    }))
  }))
}

create_statistical_test_table <- function(results, trait_names) {
  cov_names <- trait_names[-1]
  test_list <- list()

  for (nm in names(results)) {
    res <- results[[nm]]

    for (i in seq_along(cov_names)) {
      ft <- fisher_test_corr(res$decor_target_vs_cov[i], res$n_snps, rho0 = 0)
      test_list[[paste(nm, cov_names[i], "decor", sep = "_")]] <- data.table(
        subset      = nm,
        covariate   = cov_names[i],
        test_type   = "decor_target_vs_cov",
        correlation = round(ft$r, 4),
        expected    = 0,
        n           = ft$n,
        z_stat      = round(ft$z_stat, 3),
        p_value     = format(ft$p_val, scientific = TRUE),
        ci_low      = round(ft$ci_low, 4),
        ci_high     = round(ft$ci_high, 4),
        significant = ifelse(ft$p_val < 0.05, "Yes", "No")
      )
    }

    ft_od <- fisher_test_corr(res$orig_vs_decor_target, res$n_snps,
                              rho0 = res$theoretical_orig_decor)
    test_list[[paste(nm, "orig_decor", sep = "_")]] <- data.table(
      subset      = nm,
      covariate   = "original_vs_decor_target",
      test_type   = "orig_vs_decor_target",
      correlation = round(ft_od$r, 4),
      expected    = round(ft_od$rho0, 4),
      n           = ft_od$n,
      z_stat      = round(ft_od$z_stat, 3),
      p_value     = format(ft_od$p_val, scientific = TRUE),
      ci_low      = round(ft_od$ci_low, 4),
      ci_high     = round(ft_od$ci_high, 4),
      significant = ifelse(ft_od$p_val < 0.05, "Yes", "No")
    )
  }
  rbindlist(test_list, use.names = TRUE, fill = TRUE)
}

create_key_metrics_table <- function(results) {
  data.table(
    subset            = names(results),
    n_snps            = vapply(results, `[[`, integer(1), "n_snps"),
    mean_abs_orig_cor = vapply(results, function(x)
      round(mean(abs(x$orig_target_vs_cov)), 4), numeric(1)),
    mean_abs_decor_cor = vapply(results, function(x)
      round(mean(abs(x$decor_target_vs_cov)), 4), numeric(1)),
    median_abs_orig_cor = vapply(results, function(x)
      round(median(abs(x$orig_target_vs_cov)), 4), numeric(1)),
    median_abs_decor_cor = vapply(results, function(x)
      round(median(abs(x$decor_target_vs_cov)), 4), numeric(1)),
    max_abs_orig_cor  = vapply(results, function(x)
      round(max(abs(x$orig_target_vs_cov)), 4), numeric(1)),
    max_abs_decor_cor = vapply(results, function(x)
      round(max(abs(x$decor_target_vs_cov)), 4), numeric(1)),
    mean_reduction = vapply(results, function(x)
      round(mean(abs(x$orig_target_vs_cov) - abs(x$decor_target_vs_cov)), 4),
      numeric(1)),
    mean_reduction_pct = vapply(results, function(x) {
      o <- abs(x$orig_target_vs_cov); d <- abs(x$decor_target_vs_cov)
      v <- o > 0
      if (sum(v) > 0) round(100 * mean(1 - d[v] / o[v]), 2) else NA_real_
    }, numeric(1)),
    orig_vs_decor_r     = vapply(results, function(x)
      round(x$orig_vs_decor_target, 4), numeric(1)),
    theory_r            = vapply(results, function(x)
      round(x$theoretical_orig_decor, 4), numeric(1)),
    abs_diff_from_theory = vapply(results, function(x)
      round(abs(x$orig_vs_decor_target - x$theoretical_orig_decor), 4),
      numeric(1))
  )
}

create_covariate_summary_table <- function(results, trait_names) {
  cov_names <- trait_names[-1]
  rbindlist(lapply(seq_along(cov_names), function(i) {
    orig_vals  <- vapply(results, function(x) x$orig_target_vs_cov[i], numeric(1))
    decor_vals <- vapply(results, function(x) x$decor_target_vs_cov[i], numeric(1))
    v <- abs(orig_vals) > 0
    rpct <- if (sum(v) > 0)
      round(100 * mean(1 - abs(decor_vals[v]) / abs(orig_vals[v])), 2) else NA_real_
    data.table(
      covariate_index   = i,
      covariate_name    = cov_names[i],
      mean_orig_cor     = round(mean(orig_vals), 4),
      mean_decor_cor    = round(mean(decor_vals), 4),
      mean_abs_orig     = round(mean(abs(orig_vals)), 4),
      mean_abs_decor    = round(mean(abs(decor_vals)), 4),
      mean_reduction    = round(mean(abs(orig_vals) - abs(decor_vals)), 4),
      mean_reduction_pct = rpct
    )
  }))
}

# ==================================================================
# 4. COMPREHENSIVE REPORT
# ==================================================================

create_comprehensive_report <- function(results, trait_names, output_prefix) {

  cat("\n========================================\n")
  cat("Z-SCORE DECORRELATION VALIDATION REPORT\n")
  cat("========================================\n\n")

  cat("1. Summary Table (Original vs Decorrelated Target)\n")
  cat("----------------------------------------------------------\n")
  tbl_summary <- create_decorrelation_summary_table(results)
  print(tbl_summary)

  cat("\n2. Covariate Correlation Details\n")
  cat("----------------------------------------------------------\n")
  tbl_cov <- create_covariate_correlation_table(results, trait_names)
  print(tbl_cov)

  cat("\n3. Statistical Tests (Fisher)\n")
  cat("----------------------------------------------------------\n")
  tbl_test <- create_statistical_test_table(results, trait_names)
  print(tbl_test)

  cat("\n4. Key Metrics per Subset\n")
  cat("----------------------------------------------------------\n")
  tbl_key <- create_key_metrics_table(results)
  print(tbl_key)

  cat("\n5. Per-Covariate Summary (averaged over subsets)\n")
  cat("----------------------------------------------------------\n")
  tbl_covsum <- create_covariate_summary_table(results, trait_names)
  print(tbl_covsum)

  fwrite(tbl_summary, paste0(output_prefix, "_summary.csv"))
  fwrite(tbl_cov,     paste0(output_prefix, "_covariates.csv"))
  fwrite(tbl_test,    paste0(output_prefix, "_statistical_tests.csv"))
  fwrite(tbl_key,     paste0(output_prefix, "_key_metrics.csv"))
  fwrite(tbl_covsum,  paste0(output_prefix, "_covariate_summary.csv"))

  cat("\nTables saved with prefix:", output_prefix, "\n")

  invisible(list(
    summary_table         = tbl_summary,
    covariate_table       = tbl_cov,
    statistical_test_table = tbl_test,
    key_metrics           = tbl_key,
    covariate_summary     = tbl_covsum
  ))
}

# ==================================================================
# 5. PLOTTING FUNCTION
# ==================================================================

create_validation_plots <- function(results, full_stats, p_thresholds,
                                     output_prefix) {

  plot_data <- data.table(
    p_threshold = p_thresholds,
    subset_name = names(results),
    n_snps      = vapply(results, `[[`, integer(1), "n_snps"),
    obs_cor     = vapply(results, function(x)
      x$orig_vs_decor_target, numeric(1)),
    theory_cor  = vapply(results, function(x)
      x$theoretical_orig_decor, numeric(1)),
    abs_diff    = vapply(results, function(x)
      abs(x$orig_vs_decor_target - x$theoretical_orig_decor), numeric(1))
  )
  theoretical_value <- full_stats$theoretical_orig_decor

  p1 <- ggplot(plot_data, aes(x = p_threshold, y = obs_cor)) +
    geom_line(color = "#2E86AB", linewidth = 1.2) +
    geom_point(color = "#2E86AB", size = 3) +
    geom_hline(yintercept = theoretical_value,
               linetype = "dashed", color = "#A23B72", linewidth = 1) +
    geom_text(aes(label = paste0("n=", format(n_snps, big.mark = ","))),
              vjust = -1, size = 3, color = "gray50") +
    labs(title = "Validation: Original vs Decorrelated Target",
         subtitle = sprintf("Theoretical expectation under null: %.4f",
                            theoretical_value),
         x = "P-value threshold for null subset (p > cutoff)",
         y = "Correlation (original vs decorrelated target)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray50")) +
    scale_y_continuous(
      limits = c(min(plot_data$obs_cor) * 0.98,
                 max(c(plot_data$obs_cor, theoretical_value)) * 1.02)) +
    annotate("text", x = max(p_thresholds) * 0.95,
             y = theoretical_value * 1.005,
             label = "Theoretical expectation",
             color = "#A23B72", hjust = 1, vjust = -0.5, size = 3.5)

  p2 <- ggplot(plot_data, aes(x = p_threshold, y = abs_diff)) +
    geom_line(color = "#F18F01", linewidth = 1.2) +
    geom_point(color = "#F18F01", size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Absolute Deviation from Theoretical Expectation",
         x = "P-value threshold for null subset (p > cutoff)",
         y = "|observed - theoretical|") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  p3 <- ggplot(plot_data, aes(x = p_threshold, y = n_snps)) +
    geom_line(color = "#C73E1D", linewidth = 1.2) +
    geom_point(color = "#C73E1D", size = 3) +
    geom_text(aes(label = format(n_snps, big.mark = ",")),
              vjust = -1, size = 3, color = "gray50") +
    labs(title = "Sample Size by P-value Threshold",
         x = "P-value threshold for null subset (p > cutoff)",
         y = "Number of SNPs") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_y_continuous(labels = comma)

  cov_data <- rbindlist(lapply(names(results), function(nm) {
    res <- results[[nm]]
    n_covs <- length(res$orig_target_vs_cov)
    rbindlist(lapply(seq_len(n_covs), function(i) {
      data.table(
        subset          = nm,
        type            = c("Original", "Decorrelated"),
        abs_correlation = c(abs(res$orig_target_vs_cov[i]),
                            abs(res$decor_target_vs_cov[i]))
      )
    }))
  }))

  p4 <- ggplot(cov_data, aes(x = subset, y = abs_correlation, fill = type)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.7)) +
    labs(title = "Distribution of Absolute Correlations by Subset",
         x = "P-value threshold subset",
         y = "Absolute correlation",
         fill = "Correlation type") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("#2E86AB", "#A23B72"))

  pdf(paste0(output_prefix, ".pdf"), width = 10, height = 7)
  print(p1); print(p2); print(p3); print(p4)
  dev.off()

  png(paste0(output_prefix, "_correlation.png"),
      width = 1000, height = 700, res = 100)
  print(p1); dev.off()

  png(paste0(output_prefix, "_deviation.png"),
      width = 1000, height = 700, res = 100)
  print(p2); dev.off()

  png(paste0(output_prefix, "_sample_size.png"),
      width = 1000, height = 700, res = 100)
  print(p3); dev.off()

  png(paste0(output_prefix, "_boxplot.png"),
      width = 1000, height = 700, res = 100)
  print(p4); dev.off()

  cat("\nPlots saved with prefix:", output_prefix, "\n")

  invisible(list(correlation_plot = p1, deviation_plot = p2,
                 sample_size_plot = p3, boxplot = p4))
}

# ==================================================================
# 6. PER-TARGET DRIVER FUNCTION
# ==================================================================
# Encapsulates: load pipeline data -> run validation -> save tables/plots.
# Returns a short status list for the top-level summary.
#
# Data source change for fdrreg.R v2.0:
#   Harmonised z-scores now live in the MAGMA input files, not in
#   00.overlap_data/*.overlap.csv. We read the "z" column (harmonised,
#   pre-decorrelation) and the "pval" column from
#   01.magma_input/<name>.overlap.4magma.txt (space-separated).

run_target_validation <- function(target_name, pipeline_base,
                                   output_base, p_thresholds) {

  res_status <- list(target = target_name, status = "OK", message = "")

  tryCatch({
    # --- Derived paths ---
    target_pipeline_dir <- file.path(pipeline_base, target_name)
    magma_dir <- file.path(target_pipeline_dir, "01.magma_input")
    fdr_dir   <- file.path(target_pipeline_dir, "02.fdrreg_results")
    out_dir   <- file.path(output_base, target_name)
    table_dir <- file.path(out_dir, "tables")
    plot_dir  <- file.path(out_dir, "plots")
    dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(plot_dir,  showWarnings = FALSE, recursive = TRUE)

    # --- Load model info (nolasso naming) ---
    model_info_path <- file.path(fdr_dir, "fdrreg_model_info_nolasso.rds")
    if (!file.exists(model_info_path)) {
      stop("fdrreg_model_info_nolasso.rds not found at: ", model_info_path)
    }
    model_info <- readRDS(model_info_path)

    has_overlap         <- model_info$has_overlap
    traits_with_overlap <- model_info$traits_with_overlap
    cov_matrix          <- model_info$covariance_matrix

    if (!has_overlap || length(traits_with_overlap) == 0) {
      stop("No overlapping traits detected; decorrelation validation skipped.")
    }
    if (is.null(cov_matrix)) {
      stop("covariance_matrix is NULL in model_info. ",
           "Check the field name saved by fdrreg.R ",
           "(expected 'covariance_matrix').")
    }

    # --- Helper: read a MAGMA input file (space-separated) ---
    # The "z" column holds the harmonised (pre-decorrelation) z-score,
    # matching the raw_z_matrix that fdrreg.R feeds into Matpow(cov, -0.5).
    read_magma_input <- function(name) {
      fpath <- file.path(magma_dir, paste0(name, ".overlap.4magma.txt"))
      if (!file.exists(fpath)) stop("MAGMA input file not found: ", fpath)
      dt <- fread(fpath, sep = " ")
      required_cols <- c("z", "pval")
      missing_cols <- setdiff(required_cols, names(dt))
      if (length(missing_cols) > 0) {
        stop("Columns missing in ", fpath, ": ",
             paste(missing_cols, collapse = ", "))
      }
      dt
    }

    # --- Load target (harmonised z, original pval) ---
    target <- read_magma_input(target_name)

    # --- Load overlapping trait files (harmonised z) ---
    files_with <- list()
    for (trait in traits_with_overlap) {
      files_with[[trait]] <- read_magma_input(trait)
    }

    # --- Build z-matrix and correlation matrix ---
    # Uses the harmonised "z" column from every MAGMA input file. All files
    # share the same SNP ordering (fdrreg.R writes them aligned to the
    # overlap set), so cbind aligns rows correctly.
    z.stat.matr <- cbind(target$z, do.call(cbind, lapply(files_with, `[[`, "z")))
    colnames(z.stat.matr) <- c(target_name, names(files_with))
    cor.matr <- as.matrix(cov_matrix)
    trait_names <- colnames(z.stat.matr)
    target_rank <- 1

    # --- Run validation across p-value subsets ---
    full_stats <- compute_validation_stats(t(z.stat.matr), cor.matr,
                                            target_row = target_rank)
    results <- vector("list", length(p_thresholds))
    names(results) <- paste0("p>", p_thresholds)
    for (i in seq_along(p_thresholds)) {
      idx <- which(target$pval > p_thresholds[i])
      results[[i]] <- compute_validation_stats(t(z.stat.matr), cor.matr,
                                                target_row = target_rank,
                                                sub_idx = idx)
    }

    # --- Save tables and plots ---
    table_prefix <- file.path(table_dir, "decorrelation_validation")
    invisible(create_comprehensive_report(results, trait_names, table_prefix))

    plot_prefix <- file.path(plot_dir, "validation_plot")
    invisible(create_validation_plots(results, full_stats, p_thresholds, plot_prefix))

    res_status$message <- sprintf("SNPs=%d, overlap_traits=%d",
                                  nrow(target), length(traits_with_overlap))
    cat(sprintf("[DONE] %-15s : %s\n", target_name, res_status$message))

  }, error = function(e) {
    res_status$status  <<- "FAILED"
    res_status$message <<- conditionMessage(e)
    cat(sprintf("[FAIL] %-15s : %s\n", target_name, conditionMessage(e)))
  })

  res_status
}

# ==================================================================
# 7. PARALLEL EXECUTION OVER ALL TARGETS
# ==================================================================

cat("============================================================\n")
cat("Processing", length(target_list), "target(s)...\n")
cat("============================================================\n")

n_cores <- min(length(target_list), max_cores)

status_list <- mclapply(target_list, function(tg) {
  run_target_validation(tg, pipeline_base, output_base, p_thresholds)
}, mc.cores = n_cores)

# --- Top-level run summary ---
status_dt <- rbindlist(lapply(status_list, function(x) {
  data.table(target = x$target, status = x$status, message = x$message)
}))

cat("\n============================================================\n")
cat("All targets processed. Run status:\n")
cat("============================================================\n")
print(status_dt)

fwrite(status_dt, file.path(output_base, "run_status.csv"))
cat("\nRun status saved to:", file.path(output_base, "run_status.csv"), "\n")
