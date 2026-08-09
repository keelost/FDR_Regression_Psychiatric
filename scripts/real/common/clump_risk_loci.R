#!/usr/bin/env Rscript
# Project: FDRreg Analysis Pipeline
# File: clump_risk_loci.R
# Description: Identify independent genomic risk loci from FDRreg SNP-level
#              results via PLINK 1.9 clumping (r2=0.01, dist=1000kb), following
#              the FDRreg psychiatric GWAS paper. Supports multiple significance
#              metrics (BH q-value, FDRreg theoretical/empirical) and automatic
#              per-target population assignment (EUR/EAS) for batch processing.
# Date: 2026/06/02
# Version: 3.0

#-----------------------------#
#---- Argument Parsing -------#
#-----------------------------#

suppressPackageStartupMessages(library(optparse))

# Default reference panels (cleaned: '.'-ID variants removed, no duplicate rsIDs)
DEFAULT_BFILE_EUR <- Sys.getenv("FDRREG_BFILE_EUR", "")
DEFAULT_BFILE_EAS <- Sys.getenv("FDRREG_BFILE_EAS", "")

# Targets that should be clumped with the EAS panel; all others use EUR.
DEFAULT_EAS_TARGETS <- "scz.eas2019,mddco"

# Per-SNP results file name (relative to <target>/02.fdrreg_results/)
DEFAULT_INPUT_FILE <- "fdr_values_per_snp_nolasso.csv"

# Significance metric columns to clump on (each produces its own risk loci set)
DEFAULT_METRICS <- "qval_bh,fdr_theoretical,fdr_empirical"

option_list <- list(
  make_option(c("-t", "--target"), type = "character", default = "all",
              help = "Target name (e.g. scz2014), or 'all' to process every target dir. [default: %default]"),
  make_option(c("-m", "--metrics"), type = "character", default = DEFAULT_METRICS,
              help = "Comma-separated significance columns to clump on. [default: %default]"),
  make_option(c("-r", "--fdr-threshold"), type = "double", dest = "fdr_threshold",
              default = 0.05,
              help = "Primary significance threshold for defining risk loci. [default: %default]"),
  make_option(c("--extra-thresholds"), type = "character", dest = "extra_thresholds",
              default = "0.01",
              help = "Comma-separated extra thresholds to also report. [default: %default]"),
  make_option(c("--input-file"), type = "character", dest = "input_file",
              default = DEFAULT_INPUT_FILE,
              help = "Per-SNP results file name inside 02.fdrreg_results/. [default: %default]"),
    make_option(c("-o", "--output-dir"), type = "character", dest = "output_dir",
              default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
              help = "Root input directory holding per-target folders. [default: %default]"),
  make_option(c("--loci-out-dir"), type = "character", dest = "loci_out_dir",
              default = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "05.risk_loci"),
              help = "Root output directory for risk loci results. [default: %default]"),
  make_option(c("-p", "--plink"), type = "character",
              default = Sys.getenv("FDRREG_PLINK", "plink"),
              help = "Path to PLINK 1.9 executable. [default: %default]"),
  make_option(c("--bfile-eur"), type = "character", dest = "bfile_eur",
              default = DEFAULT_BFILE_EUR,
              help = "PLINK reference panel bfile prefix for EUR targets. [default: %default]"),
  make_option(c("--bfile-eas"), type = "character", dest = "bfile_eas",
              default = DEFAULT_BFILE_EAS,
              help = "PLINK reference panel bfile prefix for EAS targets. [default: %default]"),
  make_option(c("--eas-targets"), type = "character", dest = "eas_targets",
              default = DEFAULT_EAS_TARGETS,
              help = "Comma-separated targets that should use the EAS panel. [default: %default]"),
  make_option(c("-b", "--bfile"), type = "character", default = NULL,
              help = "Override: force this bfile for ALL targets (ignores auto EUR/EAS assignment). [default: none]"),
  make_option(c("--clump-r2"), type = "double", dest = "clump_r2", default = 0.01,
              help = "LD r2 threshold for clumping. [default: %default]"),
  make_option(c("--clump-kb"), type = "integer", dest = "clump_kb", default = 1000,
              help = "Physical distance threshold (kb) for clumping. [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages({
  library(data.table)
})

# Helper: split a comma-separated string into a clean character vector
split_csv_arg <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nchar(vals) > 0]
}

# Parse numeric thresholds
parse_thresholds <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) return(numeric(0))
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

metrics     <- split_csv_arg(opt$metrics)
if (length(metrics) == 0) stop("--metrics must list at least one column.")
extra_ths   <- parse_thresholds(opt$extra_thresholds)
all_ths     <- sort(unique(c(opt$fdr_threshold, extra_ths)), decreasing = TRUE)
eas_targets <- split_csv_arg(opt$eas_targets)

# Resolve which bfile a given target should use.
# Priority: explicit --bfile override > EAS list membership > EUR default.
resolve_bfile <- function(target, opt, eas_targets) {
  if (!is.null(opt$bfile)) return(list(bfile = opt$bfile, pop = "override"))
  if (target %in% eas_targets) return(list(bfile = opt$bfile_eas, pop = "EAS"))
  list(bfile = opt$bfile_eur, pop = "EUR")
}

cat("============================================================\n")
cat("FDRreg Risk Loci (PLINK clumping) Configuration\n")
cat("============================================================\n")
cat("Target(s)         :", opt$target, "\n")
cat("Loci output dir   :", opt$loci_out_dir, "\n")
cat("Input file        :", opt$input_file, "\n")
cat("Metrics           :", paste(metrics, collapse = ", "), "\n")
cat("Primary threshold :", opt$fdr_threshold, "\n")
cat("Reported thresholds:", paste(all_ths, collapse = ", "), "\n")
cat("PLINK             :", opt$plink, "\n")
if (!is.null(opt$bfile)) {
  cat("bfile (override)  :", opt$bfile, "(applied to ALL targets)\n")
} else {
  cat("bfile EUR         :", opt$bfile_eur, "\n")
  cat("bfile EAS         :", opt$bfile_eas, "\n")
  cat("EAS targets       :", paste(eas_targets, collapse = ", "), "\n")
}
cat("clump r2 / kb     :", opt$clump_r2, "/", opt$clump_kb, "\n\n")

#-----------------------------#
#---- Helper Functions -------#
#-----------------------------#

# Run PLINK clumping for one threshold and return number of loci + clumped table
run_clump <- function(assoc_file, plink, bfile, r2, kb, p_thresh, out_prefix) {
  # --clump-p1: index-variant threshold; --clump-p2: clumped-variant threshold.
  # Setting both to p_thresh ensures only significant SNPs form/enter loci.
  cmd <- sprintf(
    "%s --bfile %s --clump %s --clump-snp-field SNP --clump-field P --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out %s",
    plink, bfile, assoc_file, p_thresh, p_thresh, r2, kb, out_prefix
  )
  cat("    Running:", cmd, "\n")
  status <- system(cmd)

  clumped_file <- paste0(out_prefix, ".clumped")
  if (status != 0 || !file.exists(clumped_file)) {
    # PLINK produces no .clumped file when zero index SNPs pass the threshold.
    cat("    No .clumped output (likely 0 significant loci at this threshold).\n")
    return(list(n_loci = 0L, clumped = NULL))
  }

  clumped <- fread(clumped_file)
  # PLINK 1.9 .clumped columns: CHR, F, SNP, BP, P, TOTAL, NSIG, S05, S01, S001, S0001, SP2
  list(n_loci = nrow(clumped), clumped = clumped)
}

# Process one metric column for a single target
process_metric <- function(target, metric, fdr_dt, panel, loci_dir, opt, all_ths) {
  if (!metric %in% names(fdr_dt)) {
    warning("Column '", metric, "' not found for target ", target, ". Skipped.")
    return(NULL)
  }

  # Build PLINK assoc-style input: SNP + P (P = metric value used for ranking).
  # Floor exact zeros so PLINK does not choke on P = 0.
  clump_in <- data.table(
    SNP = fdr_dt$snpid,
    P   = pmax(fdr_dt[[metric]], .Machine$double.xmin)
  )
  clump_in <- clump_in[is.finite(P)]

  assoc_file <- file.path(loci_dir, paste0(target, ".", metric, ".clump_input.txt"))
  fwrite(clump_in, assoc_file, sep = " ")

  metric_rows <- list()
  for (th in all_ths) {
    n_sig <- sum(fdr_dt[[metric]] < th, na.rm = TRUE)
    cat(sprintf("\n  [%s < %g] significant SNPs: %d\n", metric, th, n_sig))

    th_tag     <- gsub("\\.", "p", format(th, scientific = FALSE, trim = TRUE))
    out_prefix <- file.path(loci_dir, sprintf("%s.%s.thr%s", target, metric, th_tag))

    res <- run_clump(assoc_file, opt$plink, panel$bfile,
                     opt$clump_r2, opt$clump_kb, th, out_prefix)

    if (!is.null(res$clumped)) {
      loci_out <- paste0(out_prefix, ".risk_loci.csv")
      fwrite(res$clumped, loci_out)
      cat(sprintf("    Risk loci: %d  ->  %s\n", res$n_loci, loci_out))
    } else {
      cat(sprintf("    Risk loci: %d\n", res$n_loci))
    }

    metric_rows[[length(metric_rows) + 1]] <- data.table(
      target      = target,
      population  = panel$pop,
      metric      = metric,
      threshold   = th,
      n_sig_snps  = n_sig,
      n_risk_loci = res$n_loci
    )
  }
  rbindlist(metric_rows)
}

# Process a single target across all metrics
process_target <- function(target, opt, metrics, all_ths, eas_targets) {
  cat("------------------------------------------------------------\n")
  cat("Processing target:", target, "\n")

  panel <- resolve_bfile(target, opt, eas_targets)
  cat("Population / panel:", panel$pop, "->", panel$bfile, "\n")
  cat("------------------------------------------------------------\n")

  fdr_dir  <- file.path(opt$output_dir, target, "02.fdrreg_results")
  fdr_file <- file.path(fdr_dir, opt$input_file)

  if (!file.exists(fdr_file)) {
    warning("Missing input file for target ", target, ": ", fdr_file, ". Skipped.")
    return(NULL)
  }

  loci_dir <- file.path(opt$loci_out_dir, target)
  dir.create(loci_dir, showWarnings = FALSE, recursive = TRUE)

  fdr_dt <- fread(fdr_file)

  target_rows <- list()
  for (metric in metrics) {
    res <- process_metric(target, metric, fdr_dt, panel, loci_dir, opt, all_ths)
    if (!is.null(res)) target_rows[[metric]] <- res
  }

  if (length(target_rows) == 0) return(NULL)

  summary_dt <- rbindlist(target_rows)
  fwrite(summary_dt, file.path(loci_dir, paste0(target, ".risk_loci_summary.csv")))
  cat("\nSummary for", target, ":\n")
  print(summary_dt)
  summary_dt
}

#-----------------------------#
#---- Main Execution ---------#
#-----------------------------#

# Resolve target list
if (identical(opt$target, "all")) {
  candidate_dirs <- list.dirs(opt$output_dir, full.names = FALSE, recursive = FALSE)
  targets <- candidate_dirs[
    file.exists(file.path(opt$output_dir, candidate_dirs,
                          "02.fdrreg_results", opt$input_file))
  ]
  if (length(targets) == 0) stop("No targets with ", opt$input_file, " found under ", opt$output_dir)
  cat("Found", length(targets), "targets:", paste(targets, collapse = ", "), "\n\n")
} else {
  targets <- opt$target
}

all_summaries <- list()
for (tg in targets) {
  s <- process_target(tg, opt, metrics, all_ths, eas_targets)
  if (!is.null(s)) all_summaries[[tg]] <- s
}

# Write a combined summary across all processed targets
if (length(all_summaries) > 0) {
  combined <- rbindlist(all_summaries)
  dir.create(opt$loci_out_dir, showWarnings = FALSE, recursive = TRUE)
  combined_file <- file.path(opt$loci_out_dir, "risk_loci_summary_ALL.csv")
  fwrite(combined, combined_file)
  cat("\n============================================================\n")
  cat("Combined summary written to:", combined_file, "\n")
  cat("============================================================\n")
  print(combined)

  # Also write a wide comparison table: rows = target x threshold, cols = metrics
  wide <- dcast(combined, target + population + threshold ~ metric,
                value.var = "n_risk_loci")
  wide_file <- file.path(opt$loci_out_dir, "risk_loci_comparison_wide.csv")
  fwrite(wide, wide_file)
  cat("\nWide comparison (loci per metric) written to:", wide_file, "\n")
  print(wide)
}

cat("\nDone.\n")
