# Project: FDRreg Sensitivity Analysis - Covariate Sign Modes (SNP + MetaXcan)
# File: fdrreg_sensitivity_analysis_z_abs.R
# Description: Sensitivity analysis testing two covariate transformations for the
#              FDRreg feature matrix. Only the covariates change; the response
#              variable (target z) is unchanged. Modes:
#                1. "abs"   : abs(z)  (identical to the original pipeline)
#                2. "split" : directional split. For each trait z,
#                             z_pos = ifelse(z > 0, abs(z), 0)
#                             z_neg = ifelse(z <= 0, abs(z), 0)
#              Levels:
#                PART 1: SNP        - reuse the harmonised + decorrelated MAGMA
#                                     inputs written by fdrreg.R v2.0 (z.decor),
#                                     so the SNP set and sign convention exactly
#                                     match the original analysis. No re-harmonise,
#                                     no Matpow, no LDSC read.
#                PART 2: MetaXcan   - per brain region, no decorrelation, no
#                                     allele harmonisation (gene-level z-scores).
#              Only computes fdr.the (theoretical null).
# Usage:
#   Rscript fdrreg_sensitivity_analysis_z_abs.R \
#       --target wal.scz2018 \
#       --traits_with ad,adhd2019,... --traits_no asd2019,cannabis,... \
#       --lasso N --seed 100
# Date: 2026/07/06
# Version: 3.0 (reuse fdrreg.R MAGMA inputs; gene/SNP list matches original)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(FDRreg)
  library(stringr)
})

#-------------------------------#
#---- Command-Line Arguments ---#
#-------------------------------#

option_list <- list(
  make_option("--target", type = "character", default = NULL,
              help = "Target disease name (required)"),
  make_option("--traits_with", type = "character", default = "",
              help = "Comma-separated traits WITH sample overlap (decorrelated upstream)"),
  make_option("--traits_no", type = "character", default = "",
              help = "Comma-separated traits WITHOUT sample overlap"),
  make_option("--lasso", type = "character", default = "N",
              help = "LASSO feature selection (ignored). Default: N"),
  make_option("--seed", type = "integer", default = 100,
              help = "Random seed for reproducibility. Default: 100")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$target)) stop("Argument --target is required.")

parse_trait_list <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nchar(vals) > 0]
}

target_name <- trimws(opt$target)
traits_with <- parse_trait_list(opt$traits_with)
traits_no <- parse_trait_list(opt$traits_no)
all_traits <- c(traits_with, traits_no)
random_seed <- opt$seed
if (length(all_traits) == 0) stop("At least one trait must be specified.")

set.seed(random_seed)
start_time <- Sys.time()

# The two covariate transformation modes to evaluate.
COV_MODES <- c("abs", "split")

cat("============================================================\n")
cat("FDRreg Covariate Sign-Mode Sensitivity Analysis (SNP + MetaXcan)\n")
cat("  Target:                     ", target_name, "\n")
cat("  Traits with sample overlap: ",
    ifelse(length(traits_with) == 0, "(none)", paste(traits_with, collapse = ", ")), "\n")
cat("  Traits without overlap:     ",
    ifelse(length(traits_no) == 0, "(none)", paste(traits_no, collapse = ", ")), "\n")
cat("  Covariate modes:            ", paste(COV_MODES, collapse = ", "), "\n")
cat("  Random seed:                ", random_seed, "\n")
cat("============================================================\n")

#-------------------------------#
#---- Fixed Path Parameters ----#
#-------------------------------#

base_output_path <- Sys.getenv('FDRREG_RESULTS_DIR', '')
if (!nzchar(base_output_path)) stop('Set FDRREG_RESULTS_DIR before running v7 z-score sensitivity analysis.')
# v7 results are written to a dedicated directory to keep them separate from v8.
sensitivity_output_path <- file.path(base_output_path, '01.extra.analysis', '08.z_abs_v7')
target_output_root <- file.path(sensitivity_output_path, target_name)
dir.create(target_output_root, showWarnings = FALSE, recursive = TRUE)

# Directory where fdrreg.R (v2.0) wrote the harmonised + decorrelated MAGMA inputs.
magma_input_dir <- file.path(base_output_path, target_name, '01.magma_input')

sig_threshold <- 0.01

#-------------------------------#
#---- Covariate Transform Helper ----#
#-------------------------------#

# Build an FDRreg feature matrix from a named list of signed z-score vectors.
#   mode == "abs"   : abs(z) per trait (one column per trait)
#   mode == "split" : two columns per trait,
#                       <trait>_pos = ifelse(z > 0, abs(z), 0)
#                       <trait>_neg = ifelse(z <= 0, abs(z), 0)
build_features <- function(z_list, mode) {
  if (mode == "abs") {
    m <- do.call(cbind, lapply(z_list, abs))
    colnames(m) <- names(z_list)
  } else if (mode == "split") {
    cols <- list()
    for (tr in names(z_list)) {
      z <- z_list[[tr]]
      cols[[paste0(tr, "_pos")]] <- ifelse(z > 0, abs(z), 0)
      cols[[paste0(tr, "_neg")]] <- ifelse(z <= 0, abs(z), 0)
    }
    m <- do.call(cbind, cols)
    colnames(m) <- names(cols)
  } else {
    stop("Unknown covariate mode: ", mode)
  }
  as.matrix(m)
}

#=============================================================#
#---- PART 1: SNP SENSITIVITY (reuse fdrreg.R MAGMA inputs) ---
#=============================================================#
# fdrreg.R Phase 4 already wrote, per trait and target:
#   <name>.overlap.4magma.txt  with columns snpid, ..., z.decor, p.decor
# z.decor is exactly the signed z we need:
#   - target        : decorrelated target z (the FDRreg response)
#   - traits_with   : decorrelated trait z
#   - traits_no     : harmonised (non-decorrelated) trait z
# All files share the same overlapping-SNP set and row order, so no
# re-harmonisation or Matpow is required here. This guarantees the SNP set
# (and downstream gene list) matches the original analysis exactly.

cat("\n==================== PART 1: SNP ====================\n")

if (!dir.exists(magma_input_dir)) {
  stop("MAGMA input directory not found: ", magma_input_dir,
       "\nRun fdrreg.R for this target first.")
}

magma_file <- function(name) {
  file.path(magma_input_dir, paste0(name, '.overlap.4magma.txt'))
}

# ---- Load target: snpid defines canonical order; z.decor is the response ----
target_mf <- magma_file(target_name)
if (!file.exists(target_mf)) {
  stop("Target MAGMA input not found: ", target_mf,
       "\nRun fdrreg.R for this target first.")
}
target_magma <- fread(target_mf)
for (col in c("snpid", "z.decor")) {
  if (!col %in% names(target_magma)) {
    stop(sprintf("Column '%s' missing from %s", col, target_mf))
  }
}
overlapping_snps <- target_magma$snpid
target_z_snp <- as.numeric(target_magma[["z.decor"]])
cat(sprintf("Loaded target MAGMA input: %d SNPs\n", length(overlapping_snps)))

# ---- Validate all requested traits exist as MAGMA inputs ----
missing_traits <- all_traits[!file.exists(vapply(all_traits, magma_file, character(1)))]
if (length(missing_traits) > 0) {
  stop("MAGMA inputs missing for traits: ",
       paste(missing_traits, collapse = ", "),
       "\nThese were not produced by the fdrreg.R run for this target.\n",
       "Re-run fdrreg.R with the same trait set, or remove them here.")
}
valid_traits      <- all_traits
traits_with_valid <- intersect(traits_with, valid_traits)
traits_no_valid   <- intersect(traits_no,   valid_traits)
all_traits_valid  <- c(traits_with_valid, traits_no_valid)
has_overlap       <- length(traits_with_valid) > 0

if (has_overlap) {
  cat(sprintf("Using pre-decorrelated z.decor for %d overlapping traits.\n",
              length(traits_with_valid)))
} else {
  cat("No traits with sample overlap (decorrelation was not applied upstream).\n")
}

# ---- Build per-trait SIGNED z list, aligned to the target SNP order ----
# z.decor already encodes the correct sign convention and decorrelation state.
snp_z_list <- list()
for (trait in all_traits_valid) {
  tm <- fread(magma_file(trait))
  for (col in c("snpid", "z.decor")) {
    if (!col %in% names(tm)) {
      stop(sprintf("Column '%s' missing from %s", col, magma_file(trait)))
    }
  }
  # Align to the target's SNP order (defensive; files should already match).
  idx <- match(overlapping_snps, tm$snpid)
  if (anyNA(idx)) {
    stop(sprintf("Trait '%s' MAGMA input does not cover all target SNPs.", trait))
  }
  snp_z_list[[trait]] <- as.numeric(tm[["z.decor"]][idx])
}

for (mode in COV_MODES) {
  cat(sprintf("\n[SNP] mode = %s\n", mode))
  fm <- build_features(snp_z_list, mode)
  cat(sprintf("  feature matrix: %d x %d\n", nrow(fm), ncol(fm)))
  set.seed(random_seed)
  fit <- FDRreg(target_z_snp, fm, nulltype = 'theoretical', method = 'pr')
  fdr_snp <- fit$FDR
  cat(sprintf("  fdr.the < %.2f: %d / %d\n", sig_threshold,
              sum(fdr_snp < sig_threshold, na.rm = TRUE), length(fdr_snp)))
  mode_dir <- file.path(target_output_root, mode)
  dir.create(mode_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(data.table(snpid = overlapping_snps, fdr_the = fdr_snp),
         file.path(mode_dir, paste0(target_name, '_snp_fdr_the.csv')))
}

#=============================================================#
#---- PART 2: MetaXcan SENSITIVITY ANALYSIS (per brain region) ---
#=============================================================#
# MetaXcan (S-PrediXcan) per-region results carry a signed zscore column whose
# sign is already consistent across traits (relative to predicted expression),
# so there is no SNP-level allele harmonisation to apply here. No decorrelation
# (matching the original MetaXcan pipeline). The target response is the real
# signed zscore, unchanged across modes.

cat("\n================= PART 2: MetaXcan =================\n")

metaxcan_input_dir <- file.path(base_output_path, target_name, '06.metaxcan_v7')

brain_regions <- c(
  'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
  'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
  'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
  'Brain_Hippocampus', 'Brain_Hypothalamus',
  'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
  'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra'
)

extract_trait_from_filename <- function(filename, brain_region) {
  # v7 S-PrediXcan outputs are named: gtex_v7_<trait>_in_<region>.csv
  bn <- basename(filename)
  bn <- gsub("^gtex_v7_", "", bn)
  gsub(paste0("_in_", brain_region, "\\.csv$"), "", bn)
}

# Accumulate per-mode MetaXcan results across regions.
metaxcan_by_mode <- setNames(vector("list", length(COV_MODES)), COV_MODES)

if (!dir.exists(metaxcan_input_dir)) {
  warning("MetaXcan input directory not found: ", metaxcan_input_dir)
} else {
  for (brain_region in brain_regions) {
    pattern <- paste0("_in_", brain_region, "\\.csv$")
    filenames <- list.files(metaxcan_input_dir, pattern = pattern, full.names = TRUE)
    if (length(filenames) == 0) {
      warning(sprintf("  No files for region: %s", brain_region)); next
    }

    all_trait_names <- sapply(filenames, extract_trait_from_filename, brain_region = brain_region)
    target_idx <- which(all_trait_names == target_name)
    if (length(target_idx) != 1) {
      warning(sprintf("  Target file not uniquely found in region: %s", brain_region)); next
    }

    trait_idx <- which(all_trait_names %in% all_traits)
    trait_idx <- setdiff(trait_idx, target_idx)
    if (length(trait_idx) == 0) {
      warning(sprintf("  No requested trait files in region: %s", brain_region)); next
    }

    file_indices <- c(target_idx, trait_idx)
    trait_names_region <- all_trait_names[trait_idx]

    files <- lapply(filenames[file_indices], function(x) {
      dt <- fread(x, data.table = FALSE)
      dt <- dt[complete.cases(dt$pvalue), ]
      dt[order(dt$gene), ]
    })
    target_file <- files[[1]]
    trait_files <- files[-1]

    # Align to common genes across target + all traits.
    common_genes <- Reduce(intersect, c(list(target_file$gene),
                                         lapply(trait_files, `[[`, "gene")))
    if (length(common_genes) == 0) {
      warning(sprintf("  No common genes in region: %s", brain_region)); next
    }
    target_file <- target_file[target_file$gene %in% common_genes, ]
    target_file <- target_file[order(target_file$gene), ]
    trait_files <- lapply(trait_files, function(dt) {
      dt <- dt[dt$gene %in% common_genes, ]; dt[order(dt$gene), ]
    })

    target_z_mx <- target_file$zscore
    gene_ids_clean <- gsub("\\..*", "", target_file$gene)

    mx_z_list <- setNames(lapply(trait_files, `[[`, "zscore"), trait_names_region)

    for (mode in COV_MODES) {
      fm <- build_features(mx_z_list, mode)
      set.seed(random_seed)
      fit <- FDRreg(target_z_mx, fm, nulltype = 'theoretical', method = 'pr')
      metaxcan_by_mode[[mode]][[brain_region]] <- data.table(
        region = brain_region,
        ENSEMBL_GENE_ID = gene_ids_clean,
        fdr_the = fit$FDR
      )
    }
    cat(sprintf("  %s: %d genes, %d traits\n",
                brain_region, length(gene_ids_clean), length(trait_names_region)))
  }

  # Write MetaXcan results as ONE FILE PER BRAIN REGION (per mode).
  # Layout: <target>/<mode>/metaxcan/<target>_<region>_fdr_the.csv
  for (mode in COV_MODES) {
    parts <- metaxcan_by_mode[[mode]]
    if (length(parts) == 0) {
      warning(sprintf("[MetaXcan] no results for mode: %s", mode)); next
    }
    mx_dir <- file.path(target_output_root, mode, 'metaxcan')
    dir.create(mx_dir, showWarnings = FALSE, recursive = TRUE)

    for (brain_region in names(parts)) {
      region_dt <- parts[[brain_region]]
      n_sig <- sum(region_dt$fdr_the < sig_threshold, na.rm = TRUE)
      cat(sprintf("[MetaXcan] mode = %s | %s | fdr.the < %.2f: %d / %d\n",
                  mode, brain_region, sig_threshold, n_sig, nrow(region_dt)))
      fwrite(region_dt,
             file.path(mx_dir, paste0(target_name, '_', brain_region, '_fdr_the.csv')))
    }
  }
}

#-------------------------------#
#---- Config + Completion ------#
#-------------------------------#

saveRDS(list(
  target = target_name, traits_with = traits_with, traits_no = traits_no,
  all_traits = all_traits, seed = random_seed, threshold = sig_threshold,
  cov_modes = COV_MODES, has_snp_decorrelation = has_overlap,
  snp_source = "magma_input", snp_count = length(overlapping_snps),
  metaxcan_version = "v7",
  timestamp = Sys.time()
), file.path(target_output_root, paste0(target_name, '_sensitivity_config.rds')))

end_time <- Sys.time()
cat("\n============================================================\n")
cat("Sensitivity Analysis Complete!\n")
cat(sprintf("Target: %s\n", target_name))
cat(sprintf("Results root: %s\n", target_output_root))
cat(sprintf("Duration: %.2f seconds\n",
            as.numeric(difftime(end_time, start_time, units = "secs"))))
cat("============================================================\n")

q("no")
