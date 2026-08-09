# Project: FDRreg Sensitivity Analysis - Significant Count Summary
# File: summarize_z_abs_significance.R
# Description: Aggregate significant (fdr_the < 0.01) SNP and MetaXcan gene
#              counts across all targets and the two covariate modes (abs, split).
#              MetaXcan counts are reported per region (long) and summed (wide).
# Usage: Rscript summarize_z_abs_significance.R
# Date: 2026/07/06

suppressPackageStartupMessages(library(data.table))

base_dir <- file.path(Sys.getenv('FDRREG_RESULTS_DIR', ''), '01.extra.analysis', '08.z_abs')
cov_modes <- c('abs', 'split')
sig_threshold <- 0.01

# Discover target directories that contain at least one mode subfolder.
all_subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
targets <- all_subdirs[sapply(all_subdirs, function(d) {
  any(dir.exists(file.path(base_dir, d, cov_modes)))
})]
if (length(targets) == 0) stop("No target directories found under: ", base_dir)
cat(sprintf("Found %d target(s): %s\n", length(targets), paste(targets, collapse = ", ")))

main_records <- list()   # one row per target x mode
region_records <- list() # MetaXcan per-region detail

for (target in targets) {
  for (mode in cov_modes) {
    mode_dir <- file.path(base_dir, target, mode)

    # ---- SNP ---- #
    snp_file <- file.path(mode_dir, paste0(target, '_snp_fdr_the.csv'))
    snp_n <- NA_integer_; snp_tot <- NA_integer_
    if (file.exists(snp_file)) {
      snp_dt <- fread(snp_file)
      snp_n <- sum(snp_dt$fdr_the < sig_threshold, na.rm = TRUE)
      snp_tot <- nrow(snp_dt)
    }

    # ---- MetaXcan (one file per brain region) ---- #
    mx_dir <- file.path(mode_dir, 'metaxcan')
    mx_n <- NA_integer_; mx_tot <- NA_integer_
    if (dir.exists(mx_dir)) {
      mx_files <- list.files(mx_dir, pattern = '_fdr_the\\.csv$', full.names = TRUE)
      if (length(mx_files) > 0) {
        mx_n <- 0L; mx_tot <- 0L
        for (mx_file in mx_files) {
          mx_dt <- fread(mx_file)
          # Region is stored in the file; fall back to filename if absent.
          region_label <- if ('region' %in% names(mx_dt) && nrow(mx_dt) > 0) {
            mx_dt$region[1]
          } else {
            sub(paste0('^', target, '_'), '',
                sub('_fdr_the\\.csv$', '', basename(mx_file)))
          }
          r_sig <- sum(mx_dt$fdr_the < sig_threshold, na.rm = TRUE)
          r_tot <- nrow(mx_dt)
          mx_n <- mx_n + r_sig
          mx_tot <- mx_tot + r_tot
          region_records[[length(region_records) + 1]] <- data.table(
            target = target, cov_mode = mode, region = region_label,
            metaxcan_n_sig = r_sig, metaxcan_total = r_tot
          )
        }
      }
    }

    main_records[[length(main_records) + 1]] <- data.table(
      target = target, cov_mode = mode,
      snp_n_sig = snp_n, snp_total = snp_tot,
      metaxcan_n_sig = mx_n, metaxcan_total = mx_tot
    )
  }
}

main_dt <- rbindlist(main_records)
main_dt[, cov_mode := factor(cov_mode, levels = cov_modes)]
setorder(main_dt, target, cov_mode)
main_dt[, cov_mode := as.character(cov_mode)]

main_out <- file.path(base_dir, 'z_abs_significance_summary.csv')
fwrite(main_dt, main_out)

cat("\n---------- Significant counts (fdr_the < 0.01) ----------\n")
print(main_dt)
cat(sprintf("\nSaved main summary to: %s\n", main_out))

if (length(region_records) > 0) {
  region_dt <- rbindlist(region_records)
  region_dt[, cov_mode := factor(cov_mode, levels = cov_modes)]
  setorder(region_dt, target, cov_mode, region)
  region_dt[, cov_mode := as.character(cov_mode)]
  region_out <- file.path(base_dir, 'z_abs_metaxcan_significance_by_region.csv')
  fwrite(region_dt, region_out)
  cat(sprintf("Saved MetaXcan per-region detail to: %s\n", region_out))
}
