# Project: FDRreg Analysis Pipeline
# File: FDRreg_pipeline.R
# Description: Complete pipeline for FDRreg analysis with theoretical and empirical null models.
#              Traits are harmonised to the target effect allele using a vectorized
#              data.table implementation (stack-safe, CPU-parallel, low-memory) BEFORE
#              overlap identification, so all z-scores share a consistent sign convention.
#              Handles targets with or without sample-overlapping traits.
#              Optional LASSO feature selection prior to FDRreg.
#              Saves de-correlated z-scores and p-values for MAGMA analysis.
#              Assesses covariate contributions via Hessian-based standard errors.
# Usage:
#   Rscript FDRreg_pipeline.R \
#       --target adhd2017 \
#       --traits_with alco2018,asd2019,bd2018 \
#       --traits_no antisocial,anx,cannabis \
#       --lasso N \
#       --seed 100
# Date: 2026/06/02
# Author: Jinghong QIU (Modified by Chief Scientist)
# Version: 2.0

#-----------------------------#
#---- Configuration Setup ----#
#-----------------------------#

# Clear environment and set options
rm(list = ls())
chooseCRANmirror(ind = 22)  # HK mirror setting
options(stringsAsFactors = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(FDRreg)
  library(powerplus)
  library(HelpersMG)
  library(glmnet)
  library(doParallel)
  library(parallel)
})

#-------------------------------------#
#---- Harmonisation configuration ----#
#-------------------------------------#

# Allele convention for the *.impute.map.txt files.
# Validated by smoke test: z is signed to a1 (effect allele).
TARGET_EA <- "a1"; TARGET_OA <- "a2"
TRAIT_EA  <- "a1"; TRAIT_OA  <- "a2"

# Drop strand-ambiguous palindromic SNPs (A/T, C/G).
# TRUE  -> matches TwoSampleMR action = 2 when no allele frequency is available.
# FALSE -> keep palindromes assuming forward strand (action = 1 style).
DROP_PALINDROMIC <- TRUE

# Cores for the parallel harmonisation passes.
# NOTE: higher values speed up harmonisation but increase peak memory in
# Pass B, where each core transiently holds one full trait table before it is
# subset to the overlap. Tune to your node's RAM.
HARMONISE_CORES <- max(1, min(4, detectCores() - 1))

#-------------------------------#
#---- Command-Line Arguments ---#
#-------------------------------#

option_list <- list(
  make_option("--target", type = "character", default = NULL,
              help = "Target disease name (required), e.g. adhd2017"),
  make_option("--traits_with", type = "character", default = "",
              help = "Comma-separated traits WITH sample overlap (de-correlation applied). Empty allowed."),
  make_option("--traits_no", type = "character", default = "",
              help = "Comma-separated traits WITHOUT sample overlap (optional)."),
  make_option("--lasso", type = "character", default = "N",
              help = "Apply LASSO feature selection before FDRreg (Y/N). Default: N"),
  make_option("--seed", type = "integer", default = 100,
              help = "Random seed for reproducibility. Default: 100"),
  make_option("--input_dir", type = "character",
              default = Sys.getenv("FDRREG_INPUT_DIR", ""),
              help = "Directory containing *.impute.map.txt files"),
  make_option("--output_dir", type = "character",
              default = Sys.getenv("FDRREG_OUTPUT_DIR", ""),
              help = "Pipeline output directory"),
  make_option("--ldsc_results", type = "character",
              default = Sys.getenv("FDRREG_LDSC_RESULTS", ""),
              help = "LDSC genetic-correlation result table")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required argument
if (is.null(opt$target)) {
  stop("Argument --target is required.")
}

# Helper: parse comma-separated string into a clean character vector
parse_trait_list <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nchar(vals) > 0]
}

# Assign parsed parameters
target_name                   <- trimws(opt$target)
traits_with_sample_overlap    <- parse_trait_list(opt$traits_with)
traits_without_sample_overlap <- parse_trait_list(opt$traits_no)
use_lasso                     <- toupper(trimws(opt$lasso)) == "Y"
random_seed                   <- opt$seed

# Set reproducible seed
set.seed(random_seed)

# Time tracking initialization
start_time <- Sys.time()

cat("============================================================\n")
cat("Run configuration:\n")
cat("  Target:                     ", target_name, "\n")
cat("  Traits with sample overlap: ", length(traits_with_sample_overlap), "\n")
cat("  Traits without overlap:     ", length(traits_without_sample_overlap), "\n")
cat("  LASSO feature selection:    ", ifelse(use_lasso, "Y", "N"), "\n")
cat("  Random seed:                ", random_seed, "\n")
cat("  Harmonise cores:            ", HARMONISE_CORES, "\n")
cat("  Drop palindromic SNPs:      ", DROP_PALINDROMIC, "\n")
cat("============================================================\n")

#-------------------------------#
#---- Fixed Path Parameters ----#
#-------------------------------#

base_path <- opt$input_dir
output_dir <- opt$output_dir
ldsc_results_path <- opt$ldsc_results

if (!nzchar(base_path) || !dir.exists(base_path)) {
  stop("--input_dir must point to an existing directory: ", base_path)
}
if (!nzchar(output_dir)) {
  stop("--output_dir is required (or set FDRREG_OUTPUT_DIR).")
}
if (!nzchar(ldsc_results_path) || !file.exists(ldsc_results_path)) {
  stop("--ldsc_results must point to an existing file: ", ldsc_results_path)
}

# Result-file mode label (keeps lasso / non-lasso outputs separated)
mode_label <- ifelse(use_lasso, "lasso", "nolasso")

# Create target-specific output directory structure
target_output_dir  <- file.path(output_dir, target_name)
magma_output_dir   <- file.path(target_output_dir, "01.magma_input")
fdr_output_dir     <- file.path(target_output_dir, "02.fdrreg_results")
overlap_data_dir   <- file.path(target_output_dir, "00.overlap_data")
contrib_output_dir <- file.path(target_output_dir, "03.contribution")

dir.create(overlap_data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(magma_output_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(fdr_output_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(contrib_output_dir, showWarnings = FALSE, recursive = TRUE)

#--------------------------------------#
#---- Harmonisation Helper Functions --#
#--------------------------------------#

# Reverse-complement of allele strings (vectorized, works on single-char SNPs
# and multi-char strings alike; chartr maps each base independently).
complement_allele <- function(a) chartr("ACGT", "TGCA", a)

# Strand-ambiguous palindromic SNPs: A/T or C/G pairs.
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
  (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# Classify allele orientation between target (e1/e2) and trait (o1/o2).
# Returns a numeric sign vector:
#    1  -> alleles already aligned (direct or strand-flipped)
#   -1  -> alleles swapped (swap or strand-flipped + swap), z must be negated
#   NA  -> incompatible, or palindromic when DROP_PALINDROMIC is TRUE (dropped)
classify_sign <- function(e1, e2, o1, o2) {
  co1 <- complement_allele(o1)
  co2 <- complement_allele(o2)

  aligned <- (o1 == e1 & o2 == e2) | (co1 == e1 & co2 == e2)
  swapped <- (o1 == e2 & o2 == e1) | (co1 == e2 & co2 == e1)

  s <- rep(NA_real_, length(e1))
  s[aligned] <- 1
  s[swapped] <- -1

  if (DROP_PALINDROMIC) {
    s[is_palindromic(o1, o2)] <- NA_real_
  }
  s
}

# Full harmonisation of a trait table against target alleles.
# target_alleles: data.table(snpid, e1, e2)
# trait_dt:       data.table with snpid, TRAIT_EA, TRAIT_OA, z, pval, ...
# Returns the trait table (all columns preserved) restricted to usable SNPs,
# with z realigned to the target effect allele and alleles set to e1/e2.
harmonise_trait_vec <- function(target_alleles, trait_dt) {
  m <- merge(target_alleles, trait_dt, by = "snpid")
  if (nrow(m) == 0L) return(NULL)

  s <- classify_sign(m$e1, m$e2, m[[TRAIT_EA]], m[[TRAIT_OA]])
  keep <- which(!is.na(s))
  if (length(keep) == 0L) return(NULL)

  m <- m[keep]
  s <- s[keep]

  # Realign to the target effect allele
  m[, z := z * s]
  if ("beta" %in% names(m)) {
    m[, beta := beta * s]
  }
  m[, (TRAIT_EA) := e1]
  m[, (TRAIT_OA) := e2]
  m[, c("e1", "e2") := NULL]
  m[]
}

#--------------------------------------#
#---- Phase 1: Harmonisation ----------#
#--------------------------------------#

cat("============================================================\n")
cat("Phase 1: Harmonising traits to the target effect allele...\n")
cat("============================================================\n")

# Resolve file paths
target_path <- file.path(base_path, paste0(target_name, '.impute.map.txt'))
all_trait_names <- c(traits_with_sample_overlap, traits_without_sample_overlap)

if (length(all_trait_names) == 0) {
  stop("At least one trait (with or without overlap) must be provided.")
}

trait_file_paths <- file.path(base_path, paste0(all_trait_names, '.impute.map.txt'))

# Check file existence
valid_idx <- file.exists(trait_file_paths)
if (any(!valid_idx)) {
  missing_traits <- all_trait_names[!valid_idx]
  warning("Missing files for traits: ", paste(missing_traits, collapse = ", "))
  trait_file_paths <- trait_file_paths[valid_idx]
  all_trait_names <- all_trait_names[valid_idx]
  traits_with_sample_overlap    <- intersect(traits_with_sample_overlap, all_trait_names)
  traits_without_sample_overlap <- intersect(traits_without_sample_overlap, all_trait_names)
}

# Flag: whether this target has sample-overlapping traits
has_overlap <- length(traits_with_sample_overlap) > 0
if (has_overlap) {
  cat(sprintf("Sample-overlapping traits detected (%d). De-correlation will be applied.\n",
              length(traits_with_sample_overlap)))
} else {
  cat("No sample-overlapping traits. De-correlation will be SKIPPED.\n")
}

# ---- Pass A (light): read alleles only, decide kept SNPs per trait ----
# Only snpid + effect/other alleles are needed to classify orientation, so this
# pass stays memory-cheap even at genome-wide scale.
cat("Pass A: reading alleles and classifying orientation (memory-light)...\n")

target_alleles_full <- fread(target_path, select = c("snpid", TARGET_EA, TARGET_OA))
setnames(target_alleles_full, c("snpid", TARGET_EA, TARGET_OA), c("snpid", "e1", "e2"))
setkey(target_alleles_full, snpid)

kept_snpid_list <- mclapply(seq_along(trait_file_paths), function(i) {
  ta <- fread(trait_file_paths[i], select = c("snpid", TRAIT_EA, TRAIT_OA))
  m  <- merge(target_alleles_full, ta, by = "snpid")
  if (nrow(m) == 0L) return(character(0))
  s <- classify_sign(m$e1, m$e2, m[[TRAIT_EA]], m[[TRAIT_OA]])
  m$snpid[!is.na(s)]
}, mc.cores = min(length(trait_file_paths), HARMONISE_CORES))
names(kept_snpid_list) <- all_trait_names

# Detect and report traits with no usable SNPs (hard failure, no silent drop)
empty_traits <- all_trait_names[vapply(kept_snpid_list, length, integer(1)) == 0]
if (length(empty_traits) > 0) {
  stop(sprintf(
    "Harmonisation produced zero usable SNPs for: %s\nCheck allele columns / DROP_PALINDROMIC setting.",
    paste(empty_traits, collapse = ", ")
  ))
}

for (nm in all_trait_names) {
  cat(sprintf("  %-22s usable SNPs after harmonisation: %d\n", nm, length(kept_snpid_list[[nm]])))
}

#--------------------------------------#
#---- Phase 2: Overlap & Filtering ----#
#--------------------------------------#

cat("\n============================================================\n")
cat("Phase 2: Identifying overlapping SNPs and loading data...\n")
cat("============================================================\n")

# Post-harmonisation overlap = intersection of every trait's usable SNPs.
# (kept SNPs already come from a merge with the target, so target is included.)
overlapping_snps <- Reduce(intersect, kept_snpid_list)
cat(sprintf("Final overlapping SNPs (post-harmonisation): %d\n\n", length(overlapping_snps)))

if (length(overlapping_snps) == 0) {
  stop("No overlapping SNPs found after harmonisation. Check data compatibility.")
}

# Save overlapping SNP IDs for reference (IDs only, not harmonised data)
fwrite(data.table(snpid = overlapping_snps),
       file.path(overlap_data_dir, "overlapping_snps_list.csv"))

# Free the kept-SNP lists; only the intersection is needed from here on
rm(kept_snpid_list); gc()

# Target alleles restricted to the overlap (reused for Pass B harmonisation)
target_alleles_overlap <- target_alleles_full[.(overlapping_snps)]
rm(target_alleles_full); gc()

# ---- Load target full table, subset to overlap (target z is the reference) ----
cat("Loading target data and subsetting to overlap...\n")
target_full <- fread(target_path)
setkey(target_full, snpid)
target_data <- target_full[.(overlapping_snps)]  # ordered to overlapping_snps
rm(target_full); gc()

# ---- Pass B (heavy): read each trait fully, subset to overlap, harmonise ----
# Each full table is materialized only transiently, then immediately reduced to
# the overlap, keeping peak memory to roughly HARMONISE_CORES full tables.
cat("Pass B: harmonising trait z-scores on the overlap set (parallel)...\n")

harmonised_traits <- mclapply(seq_along(trait_file_paths), function(i) {
  dt <- fread(trait_file_paths[i])
  setkey(dt, snpid)
  dt <- dt[.(overlapping_snps)]                 # subset to overlap, free the rest
  h  <- harmonise_trait_vec(target_alleles_overlap, dt)
  rm(dt)
  # Reorder to the canonical overlap order so all cbind operations align by row
  setkey(h, snpid)
  h[.(overlapping_snps)]
}, mc.cores = min(length(trait_file_paths), HARMONISE_CORES))
names(harmonised_traits) <- all_trait_names

# Sanity check: every trait must cover the full overlap in identical order
bad <- names(harmonised_traits)[vapply(harmonised_traits,
                                       function(x) nrow(x) != length(overlapping_snps),
                                       logical(1))]
if (length(bad) > 0) {
  stop("Row-count mismatch after harmonisation for: ", paste(bad, collapse = ", "))
}

rm(target_alleles_overlap); gc()

#---------------------------------------------#
#---- Phase 3: De-correlation & Feature Prep --#
#---------------------------------------------#

cat("============================================================\n")
cat("Phase 3: De-correlation and feature preparation...\n")
cat("============================================================\n")

# Traits without sample overlap (always needed) - harmonised tables
files_no_overlap <- harmonised_traits[traits_without_sample_overlap]

# Initialize containers conditionally
files_with_overlap <- list()
covariance_matrix  <- NULL
decorr_z_matrix    <- NULL

if (has_overlap) {
  # ---- Branch A: targets WITH sample-overlapping traits ----
  cat("De-correlation branch: building covariance matrix and applying Matpow...\n")

  # Load LDSC genetic correlation results
  ldsc_results <- fread(ldsc_results_path)

  # Harmonised traits with sample overlap
  files_with_overlap <- harmonised_traits[traits_with_sample_overlap]

  # Build covariance matrix from LDSC intercepts
  decorr_traits <- c(target_name, traits_with_sample_overlap)
  n_decorr <- length(decorr_traits)
  covariance_matrix <- matrix(NA, nrow = n_decorr, ncol = n_decorr)
  rownames(covariance_matrix) <- decorr_traits
  colnames(covariance_matrix) <- decorr_traits

  for (i in 1:n_decorr) {
    for (j in 1:n_decorr) {
      trait_i <- decorr_traits[i]
      trait_j <- decorr_traits[j]

      if (i == j) {
        covariance_matrix[i, j] <- 1
      } else {
        prefix_i <- ifelse(trait_i == target_name, "T:", "L:")
        prefix_j <- ifelse(trait_j == target_name, "T:", "L:")

        keys_to_try <- list(
          list(p1 = paste0(prefix_i, trait_i), p2 = paste0(prefix_j, trait_j)),
          list(p1 = paste0(prefix_j, trait_j), p2 = paste0(prefix_i, trait_i))
        )

        found <- FALSE
        value <- NA

        for (k in 1:length(keys_to_try)) {
          current_key <- keys_to_try[[k]]
          matched_row <- ldsc_results[p1 == current_key$p1 & p2 == current_key$p2]

          if (nrow(matched_row) == 1) {
            value <- matched_row$gcov_intercept
            found <- TRUE
            break
          }
        }

        if (found) {
          covariance_matrix[i, j] <- value
          covariance_matrix[j, i] <- value
        } else {
          warning(paste("No matching data found for trait pair:", trait_i, "and", trait_j))
          covariance_matrix[i, j] <- 0
          covariance_matrix[j, i] <- 0
        }
      }
    }
  }

  # Calculate de-correlated z-statistics using Matpow.
  # Trait z-scores are already harmonised to the target effect allele.
  raw_z_matrix <- cbind(target_data$z,
                        do.call(cbind, lapply(files_with_overlap, `[[`, "z")))
  decorr_z_matrix <- Matpow(covariance_matrix, -0.5) %*% t(raw_z_matrix)

  # Target signed z (de-correlated)
  target_z <- decorr_z_matrix[1, ]
  cat("De-correlation complete.\n")

} else {
  # ---- Branch B: targets WITHOUT sample-overlapping traits ----
  cat("No-overlap branch: using original target z-scores (no de-correlation).\n")
  target_z <- target_data$z
}

#----------------------------------------------#
#---- Phase 4: Save MAGMA Input Files ----------#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 4: Saving de-correlated z-scores and p-values for MAGMA...\n")
cat("============================================================\n")

# (a) Target MAGMA input
cat("Saving target MAGMA input...\n")
if (has_overlap) {
  target_magma <- cbind(decorr_z_matrix[1, ], target_data)
  colnames(target_magma)[1] <- 'z.decor'
  target_magma$p.decor <- 2 * pnorm(abs(target_magma$z.decor), lower.tail = FALSE)
} else {
  target_magma <- copy(target_data)
  target_magma$z.decor <- target_magma$z
  target_magma$p.decor <- target_magma$pval
}
fwrite(as.data.frame(target_magma),
       file.path(magma_output_dir, paste0(target_name, '.overlap.4magma.txt')),
       sep = ' ')
cat(sprintf("  Written: %s.overlap.4magma.txt (%d SNPs)\n", target_name, nrow(target_magma)))

# (b) Traits with sample overlap: de-correlated z and p (only if has_overlap)
if (has_overlap) {
  cat("Saving traits with sample overlap MAGMA input...\n")
  for (i in 1:length(files_with_overlap)) {
    with_magma <- cbind(decorr_z_matrix[i + 1, ], files_with_overlap[[i]])
    colnames(with_magma)[1] <- 'z.decor'
    with_magma$p.decor <- 2 * pnorm(abs(with_magma$z.decor), lower.tail = FALSE)
    fwrite(as.data.frame(with_magma),
           file.path(magma_output_dir, paste0(traits_with_sample_overlap[i], '.overlap.4magma.txt')),
           sep = ' ')
    cat(sprintf("  Written: %s.overlap.4magma.txt (%d SNPs)\n",
                traits_with_sample_overlap[i], nrow(with_magma)))
  }
}

# (c) Traits without sample overlap: harmonised z and original p
cat("Saving traits without sample overlap MAGMA input...\n")
for (i in seq_along(files_no_overlap)) {
  no_magma <- files_no_overlap[[i]]
  no_magma$p.decor <- no_magma$pval
  no_magma$z.decor <- no_magma$z
  fwrite(as.data.frame(no_magma),
         file.path(magma_output_dir, paste0(traits_without_sample_overlap[i], '.overlap.4magma.txt')),
         sep = ' ')
  cat(sprintf("  Written: %s.overlap.4magma.txt (%d SNPs)\n",
              traits_without_sample_overlap[i], nrow(no_magma)))
}

cat("All MAGMA input files saved.\n")

#----------------------------------------------#
#---- Phase 5: Feature Matrix Construction -----#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 5: Building feature matrix (absolute z-scores)...\n")
cat("============================================================\n")

# Absolute z-scores. drop = FALSE preserves matrix structure when only ONE
# overlapping trait exists (otherwise the row collapses to a vector).
if (has_overlap) {
  feature_matrix <- cbind(
    abs(t(decorr_z_matrix[-1, , drop = FALSE])),
    do.call(cbind, lapply(files_no_overlap, function(x) abs(x$z)))
  )
  feature_names <- c(traits_with_sample_overlap, traits_without_sample_overlap)
} else {
  feature_matrix <- do.call(cbind, lapply(files_no_overlap, function(x) abs(x$z)))
  feature_names <- traits_without_sample_overlap
}
colnames(feature_matrix) <- feature_names

#----------------------------------------------#
#---- Phase 5b: Optional LASSO Selection -------#
#----------------------------------------------#

if (use_lasso) {
  cat("\n------------------------------------------------------------\n")
  cat("LASSO feature selection (cv.glmnet, alpha = 1)...\n")
  cat("------------------------------------------------------------\n")

  # Register parallel backend for cv.glmnet
  n_cores <- max(1, min(detectCores() - 1, 4))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Cross-validated LASSO: response is |target z|
  lasso_cv <- cv.glmnet(
    feature_matrix,
    abs(target_z),
    family     = 'gaussian',
    nlambda    = 50,
    alpha      = 1,
    standardize = TRUE,
    parallel   = TRUE
  )

  stopCluster(cl)

  # Coefficients at lambda.min (drop intercept)
  lasso_coef <- as.matrix(coef(lasso_cv, s = "lambda.min"))
  lasso_coef <- unlist(lasso_coef)[-1]

  selected_mask   <- lasso_coef != 0
  fdr_feature_mat <- feature_matrix[, selected_mask, drop = FALSE]

  dropped_features <- setdiff(colnames(feature_matrix), colnames(fdr_feature_mat))
  cat('Variables not selected from gwas:',
      ifelse(length(dropped_features) == 0, "(none)", paste(dropped_features, collapse = ", ")),
      '\n')
  cat(sprintf('Variables retained: %d / %d\n',
              ncol(fdr_feature_mat), ncol(feature_matrix)))

  if (ncol(fdr_feature_mat) == 0) {
    stop("LASSO removed all features. Cannot proceed with FDRreg.")
  }
} else {
  fdr_feature_mat <- feature_matrix
}

#----------------------------------------------#
#---- Phase 6: FDRreg Analysis ----------------#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 6: Performing FDRreg analysis (target signed + features abs)...\n")
cat("============================================================\n")

# Run FDRreg with theoretical null
cat("Running FDRreg (nulltype = 'theoretical', method = 'pr')...\n")
fdr_theoretical <- FDRreg(
  target_z,
  fdr_feature_mat,
  nulltype = "theoretical",
  method   = "pr"
)

# Run FDRreg with empirical null
cat("Running FDRreg (nulltype = 'empirical', method = 'pr')...\n")
fdr_empirical <- FDRreg(
  target_z,
  fdr_feature_mat,
  nulltype = "empirical",
  method   = "pr"
)

#----------------------------------------------#
#---- Phase 7: Covariate Contribution ---------#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 7: Assessing covariate contributions...\n")
cat("============================================================\n")

# Helper: extract contribution table from a fitted FDRreg model
extract_contribution <- function(fdr_model, feature_mat) {
  model_se   <- SEfromHessian(fdr_model$model$hessian)
  model_coef <- fdr_model$model$coef
  covariate_idx <- 2:length(model_coef)  # skip intercept (index 1)

  covariate_z    <- model_coef[covariate_idx] / model_se[covariate_idx]
  covariate_pval <- 2 * pnorm(abs(covariate_z), lower.tail = FALSE)

  assessment <- data.frame(
    covariate   = colnames(feature_mat),
    pvalue      = covariate_pval,
    coefficient = model_coef[covariate_idx],
    std_error   = model_se[covariate_idx],
    stringsAsFactors = FALSE
  )

  assessment[order(assessment$pvalue), ]
}

# (a) Contribution from theoretical null model
cat("Extracting contribution from theoretical null model...\n")
contrib_theoretical <- extract_contribution(fdr_theoretical, fdr_feature_mat)
fwrite(contrib_theoretical,
       file.path(contrib_output_dir,
                 paste0("contribution_", target_name, "_", mode_label, "_theoretical.csv")),
       sep = ',')
cat(sprintf("  Written: contribution_%s_%s_theoretical.csv (%d covariates)\n",
            target_name, mode_label, nrow(contrib_theoretical)))
cat("\n---------- Top Covariates (Theoretical Null) ----------\n")
print(head(contrib_theoretical, 10))

# (b) Contribution from empirical null model
cat("\nExtracting contribution from empirical null model...\n")
contrib_empirical <- extract_contribution(fdr_empirical, fdr_feature_mat)
fwrite(contrib_empirical,
       file.path(contrib_output_dir,
                 paste0("contribution_", target_name, "_", mode_label, "_empirical.csv")),
       sep = ',')
cat(sprintf("  Written: contribution_%s_%s_empirical.csv (%d covariates)\n",
            target_name, mode_label, nrow(contrib_empirical)))
cat("\n---------- Top Covariates (Empirical Null) ----------\n")
print(head(contrib_empirical, 10))

#----------------------------------------------#
#---- Phase 8: Baseline qvalue (BH) -----------#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 8: Computing baseline q-values (BH on raw pval)...\n")
cat("============================================================\n")

# Baseline q-value from the ORIGINAL (non-FDRreg) target pval via BH adjustment.
# This provides a standard multiple-testing comparison against FDRreg output.
target_qval <- p.adjust(target_data$pval, method = "BH")

#----------------------------------------------#
#---- Phase 9: Results Summary ----------------#
#----------------------------------------------#

cat("============================================================\n")
cat("Phase 9: Summarizing results...\n")
cat("============================================================\n")

# Significance thresholds
thresholds <- c(
  0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01,
  0.001, 5e-04, 5e-06, 5e-08
)

# Count values below each threshold
count_below <- function(value_vec, th_vec) {
  vapply(th_vec, function(th) sum(value_vec < th, na.rm = TRUE), integer(1))
}

combined_summary <- data.frame(
  threshold        = paste0("<", thresholds),
  fdr_theoretical  = count_below(fdr_theoretical$FDR, thresholds),
  fdr_empirical    = count_below(fdr_empirical$FDR, thresholds),
  qval_bh          = count_below(target_qval, thresholds),
  stringsAsFactors = FALSE
)

cat("\n========== FDRreg Analysis Results (target signed + features abs) ==========\n")
print(combined_summary)

#----------------------------------------------#
#---- Phase 10: Save All Results --------------#
#----------------------------------------------#

cat("\n============================================================\n")
cat("Phase 10: Saving all results...\n")
cat("============================================================\n")

# (a) Individual SNP results: FDRreg values + baseline q-value
fdr_results_dt <- data.table(
  snpid           = target_data$snpid,
  z_score_signed  = target_z,
  pval_raw        = target_data$pval,
  qval_bh         = target_qval,
  fdr_theoretical = fdr_theoretical$FDR,
  fdr_empirical   = fdr_empirical$FDR
)
fwrite(fdr_results_dt,
       file.path(fdr_output_dir, paste0("fdr_values_per_snp_", mode_label, ".csv")))
cat(sprintf("  Written: fdr_values_per_snp_%s.csv (%d SNPs)\n",
            mode_label, nrow(fdr_results_dt)))

# (b) Summary table
fwrite(combined_summary,
       file.path(fdr_output_dir, paste0("fdr_summary_table_", mode_label, ".csv")))
cat(sprintf("  Written: fdr_summary_table_%s.csv\n", mode_label))

# (c) Covariance matrix (only if de-correlation was applied)
if (has_overlap) {
  fwrite(as.data.frame(covariance_matrix),
         file.path(fdr_output_dir, "covariance_matrix.csv"),
         row.names = TRUE)
  cat("  Written: covariance_matrix.csv\n")
} else {
  cat("  Skipped: covariance_matrix.csv (no de-correlation performed)\n")
}

# (d) Complete model information
model_info <- list(
  target                 = target_name,
  has_overlap            = has_overlap,
  use_lasso              = use_lasso,
  drop_palindromic       = DROP_PALINDROMIC,
  traits_with_overlap    = traits_with_sample_overlap,
  traits_without_overlap = traits_without_sample_overlap,
  selected_features      = colnames(fdr_feature_mat),
  set_seed               = random_seed,
  covariance_matrix      = covariance_matrix,  # NULL if no overlap
  fdr_theoretical        = fdr_theoretical,
  fdr_empirical          = fdr_empirical,
  contrib_theoretical    = contrib_theoretical,
  contrib_empirical      = contrib_empirical
)
saveRDS(model_info,
        file.path(fdr_output_dir, paste0("fdrreg_model_info_", mode_label, ".rds")))
cat(sprintf("  Written: fdrreg_model_info_%s.rds\n", mode_label))

#-------------------------------#
#---- Completion Summary ------#
#-------------------------------#

end_time <- Sys.time()
cat('\n============================================================\n')
cat('Analysis complete!\n')
cat('============================================================\n')
cat('Target:', target_name, '\n')
cat('De-correlation applied:', has_overlap, '\n')
cat('LASSO applied:', use_lasso, '\n')
cat('Features used in FDRreg:', ncol(fdr_feature_mat), '\n')
cat('Overlapping SNPs:', length(overlapping_snps), '\n')
cat('\nTime Start:', format(start_time, "%a %b %d %X %Y"), '\n')
cat('Time End:',   format(end_time,   "%a %b %d %X %Y"), '\n')
cat('Time consuming:', end_time - start_time, '\n')

q("no")
