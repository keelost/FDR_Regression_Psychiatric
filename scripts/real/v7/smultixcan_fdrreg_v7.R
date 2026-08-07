#!/usr/bin/env Rscript
#*************************************************************#
# SmultiXcan FDRreg Analysis Pipeline
#*************************************************************#
# Usage:
#   Rscript smultixcan_fdrreg.R \
#     --target bd2018,adhd2016 \
#     --bio Y --lasso Y \
#     --traits college2013,asd2019
#
# R-4.0.2 or later
#*************************************************************#

suppressPackageStartupMessages({
  library(optparse)
  library(FDRreg)
  library(data.table)
  library(powerplus)
  library(stringr)
  library(dplyr)
  library(glmnet)
  library(HelpersMG)
})

# ======================== #
# ---- Helper Functions -- #
# ======================== #

# Safe p-value to z-score conversion
p_to_z_safe <- function(p, cap = 8.2) {
  p_clamped <- pmax(pmin(p, 1 - .Machine$double.eps), .Machine$double.eps)
  z <- qnorm(1 - p_clamped / 2)
  z <- pmin(pmax(z, -cap), cap)
  return(z)
}

# Count significant genes at various thresholds
count_sig <- function(fdr_vec, thresholds) {
  vapply(thresholds, function(th) sum(fdr_vec < th, na.rm = TRUE), integer(1))
}

# Safe FDRreg wrapper with error handling
safe_fdrreg <- function(target_z, predictors, nulltype, method) {
  tryCatch({
    FDRreg(target_z, predictors, nulltype = nulltype, method = method)
  }, error = function(e) {
    cat(sprintf("    FDRreg failed (%s null): %s\n", nulltype, conditionMessage(e)))
    NULL
  })
}

# Extract model assessment (feature p-values, betas, SEs)
extract_assessment <- function(fdr_model, feature_names) {
  if (is.null(fdr_model)) {
    return(data.frame(
      feature = feature_names, p = NA, beta = NA, se = NA,
      stringsAsFactors = FALSE
    ))
  }
  model_se <- SEfromHessian(fdr_model$model$hessian)
  model_coef <- fdr_model$model$coef
  model_z <- model_coef[-1] / model_se[-1]
  model_p <- 2 * pnorm(abs(model_z), lower.tail = FALSE)
  data.frame(
    feature = feature_names, p = model_p,
    beta = model_coef[-1], se = model_se[-1],
    stringsAsFactors = FALSE
  )
}

# Merge assessment data frames, filling missing columns with NA
merge_assessments <- function(base_df, new_df, suffix) {
  if (is.null(new_df)) {
    base_df[[paste0("p", suffix)]] <- NA
    base_df[[paste0("beta", suffix)]] <- NA
    base_df[[paste0("se", suffix)]] <- NA
    return(base_df)
  }
  setnames(new_df, c("p", "beta", "se"),
           c(paste0("p", suffix), paste0("beta", suffix), paste0("se", suffix)))
  merge(base_df, new_df, by = "feature", all = TRUE)
}

# ================================== #
# ---- Core Processing Function ---- #
# ================================== #

process_target <- function(target, user_traits = NULL, run_bio = TRUE,
                           run_lasso = TRUE, base_dir, bio_library = NULL,
                           seed = 100) {

  cat(sprintf("\n========== Target: %s ==========\n", target))

  # ---- Setup paths (v7) ---- #
  input_dir <- file.path(base_dir, target, "08.smultixcan_output_v7")
  output_dir <- file.path(base_dir, target, "09.smultixcan_fdrreg_v7")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # ---- Discover input files ---- #
  all_files <- list.files(input_dir, pattern = "\\.allbrain\\.txt$", full.names = TRUE)
  if (length(all_files) == 0) {
    stop("No .allbrain.txt files found in: ", input_dir)
  }

  # Identify target file by exact basename match
  target_basename <- paste0(target, ".allbrain.txt")
  target_file_idx <- which(basename(all_files) == target_basename)
  if (length(target_file_idx) == 0) {
    stop("Target file not found: ", target_basename)
  }

  # Identify trait files (exclude target)
  trait_file_idx <- setdiff(seq_along(all_files), target_file_idx)

  # Filter traits if user specified
  if (!is.null(user_traits)) {
    trait_basenames <- paste0(user_traits, ".allbrain.txt")
    matched <- which(basename(all_files[trait_file_idx]) %in% trait_basenames)
    if (length(matched) == 0) {
      stop("No matching trait files for: ", paste(user_traits, collapse = ", "))
    }
    trait_file_idx <- trait_file_idx[matched]
    cat(sprintf("Traits specified: %d matched\n", length(trait_file_idx)))
  } else {
    cat(sprintf("Traits: all %d in folder\n", length(trait_file_idx)))
  }

  # ---- Load data ---- #
  cat("Loading data...\n")
  load_indices <- c(target_file_idx, trait_file_idx)
  smultixcan_data <- lapply(load_indices, function(i) {
    dt <- fread(all_files[i], stringsAsFactors = FALSE, showProgress = FALSE)
    essential_cols <- intersect(c("gene", "gene_name", "pvalue"), names(dt))
    dt[, ..essential_cols]
  })

  trait_names <- str_remove(basename(all_files[trait_file_idx]), "\\.allbrain\\.txt$")

  # ---- Preprocess: remove NAs, sort, compute z-scores ---- #
  cat("Preprocessing...\n")
  for (i in seq_along(smultixcan_data)) {
    smultixcan_data[[i]] <- smultixcan_data[[i]][complete.cases(pvalue)]
    smultixcan_data[[i]] <- smultixcan_data[[i]][order(gene)]
    smultixcan_data[[i]][, zscore := p_to_z_safe(pvalue)]
  }

  # Find common genes across ALL files
  common_genes <- Reduce(intersect, lapply(smultixcan_data, `[[`, "gene"))
  cat(sprintf("Common genes: %d\n", length(common_genes)))
  if (length(common_genes) == 0) stop("No common genes found across files")

  for (i in seq_along(smultixcan_data)) {
    smultixcan_data[[i]] <- smultixcan_data[[i]][gene %in% common_genes][order(gene)]
  }

  # Verify consistency
  n_rows <- sapply(smultixcan_data, nrow)
  if (length(unique(n_rows)) != 1) {
    stop("Row count mismatch after filtering: ", paste(n_rows, collapse = ", "))
  }

  # ---- Build matrices ---- #
  target_data <- smultixcan_data[[1]]
  predictor_zscores <- sapply(2:length(smultixcan_data), function(i) {
    smultixcan_data[[i]]$zscore
  })
  colnames(predictor_zscores) <- trait_names

  # Target z-scores with random sign for null distribution
  set.seed(seed)
  target_zscores <- target_data$zscore * sign(rnorm(nrow(target_data)))

  # Extract Ensembl gene IDs (strip version suffix)
  gene_ids <- str_split_fixed(target_data$gene, "\\.", 2)[, 1]

  # Save gene list
  fwrite(data.frame(gene_id = gene_ids),
         file = file.path(output_dir, paste0(target, ".gene.list.txt")))

  # ---- Initialize result containers ---- #
  target_results <- copy(target_data)
  target_results[, gene := gene_ids]
  setnames(target_results, c("gene", "gene_name"),
           c("ENSEMBL_GENE_ID", "ENSEMBL_GENE_ID_name"))
  target_results[, qval := p.adjust(pvalue, method = "fdr")]

  # Initialize all FDR columns with NA
  target_results[, c("FDR.the", "FDR.emp",
                     "bio.FDR.the", "bio.FDR.emp",
                     "bio.FDR.the.lasso", "bio.FDR.emp.lasso") := NA]

  # Initialize assessment summary
  assessment_summary <- data.frame(feature = trait_names, stringsAsFactors = FALSE)
  assessment_summary[c("p1", "beta1", "se1",
                       "p2", "beta2", "se2",
                       "p3", "beta3", "se3")] <- NA

  # Track all models for saving
  model_list <- list(
    fdr_theoretical = NULL, fdr_empirical = NULL,
    fdr_bio_theoretical = NULL, fdr_bio_empirical = NULL,
    fdr_lasso_theoretical = NULL, fdr_lasso_empirical = NULL
  )

  # ========================================== #
  # ---- Step 1: Basic FDRreg (traits only) ----
  # ========================================== #
  cat("Step 1: Basic FDRreg (traits only)...\n")
  model_list$fdr_theoretical <- safe_fdrreg(target_zscores, predictor_zscores, 'theoretical', 'pr')
  model_list$fdr_empirical <- safe_fdrreg(target_zscores, predictor_zscores, 'empirical', 'pr')

  target_results[, FDR.the := if (!is.null(model_list$fdr_theoretical)) model_list$fdr_theoretical$FDR else NA]
  target_results[, FDR.emp := if (!is.null(model_list$fdr_empirical)) model_list$fdr_empirical$FDR else NA]

  assessment_basic <- extract_assessment(model_list$fdr_theoretical, trait_names)
  setnames(assessment_basic, c("p", "beta", "se"), c("p1", "beta1", "se1"))
  assessment_summary <- merge(assessment_summary, assessment_basic, by = "feature")

  # ========================================== #
  # ---- Step 2: FDRreg with bio annotations ----
  # ========================================== #
  bio_features_matrix <- NULL
  if (run_bio && !is.null(bio_library)) {
    cat("Step 2: FDRreg with bio annotations...\n")

    target_bio <- merge(
      target_results[, .(ENSEMBL_GENE_ID, ENSEMBL_GENE_ID_name)],
      bio_library, by = "ENSEMBL_GENE_ID", all.x = TRUE
    )
    target_bio[is.na(target_bio)] <- 0

    bio_feature_cols <- setdiff(names(target_bio),
                                c("ENSEMBL_GENE_ID", "ENSEMBL_GENE_ID_name"))
    bio_features_matrix <- as.matrix(target_bio[, ..bio_feature_cols])

    if (ncol(bio_features_matrix) > 0) {
      predictor_bio_zscores <- cbind(predictor_zscores, abs(bio_features_matrix))

      model_list$fdr_bio_theoretical <- safe_fdrreg(
        target_zscores, predictor_bio_zscores, 'theoretical', 'pr')
      model_list$fdr_bio_empirical <- safe_fdrreg(
        target_zscores, predictor_bio_zscores, 'empirical', 'pr')

      target_results[, bio.FDR.the := if (!is.null(model_list$fdr_bio_theoretical))
        model_list$fdr_bio_theoretical$FDR else NA]
      target_results[, bio.FDR.emp := if (!is.null(model_list$fdr_bio_empirical))
        model_list$fdr_bio_empirical$FDR else NA]

      assessment_bio <- extract_assessment(model_list$fdr_bio_theoretical,
                                           colnames(predictor_bio_zscores))
      assessment_summary <- merge_assessments(assessment_summary, assessment_bio, "2")
    } else {
      cat("  No bio features available\n")
    }
  } else {
    cat("Step 2: Skipped (bio annotations disabled)\n")
  }

  # ========================================== #
  # ---- Step 3: FDRreg with Lasso features ----
  # ========================================== #
  if (run_lasso && !is.null(bio_features_matrix) && ncol(bio_features_matrix) > 0) {
    cat("Step 3: FDRreg with Lasso-selected features...\n")
    tryCatch({
      lasso_cv <- cv.glmnet(
        bio_features_matrix, abs(target_zscores),
        family = "gaussian", nlambda = 50, alpha = 1,
        standardize = TRUE, parallel = TRUE
      )
      lasso_coef <- as.numeric(coef(lasso_cv, s = "lambda.min"))[-1]
      nonzero_idx <- which(lasso_coef != 0)

      if (length(nonzero_idx) > 0) {
        lasso_features <- bio_features_matrix[, nonzero_idx, drop = FALSE]
        predictor_lasso_zscores <- cbind(predictor_zscores, abs(lasso_features))
        cat(sprintf("  Lasso selected %d features\n", length(nonzero_idx)))

        model_list$fdr_lasso_theoretical <- safe_fdrreg(
          target_zscores, predictor_lasso_zscores, 'theoretical', 'pr')
        model_list$fdr_lasso_empirical <- safe_fdrreg(
          target_zscores, predictor_lasso_zscores, 'empirical', 'pr')

        target_results[, bio.FDR.the.lasso := if (!is.null(model_list$fdr_lasso_theoretical))
          model_list$fdr_lasso_theoretical$FDR else NA]
        target_results[, bio.FDR.emp.lasso := if (!is.null(model_list$fdr_lasso_empirical))
          model_list$fdr_lasso_empirical$FDR else NA]

        assessment_lasso <- extract_assessment(model_list$fdr_lasso_theoretical,
                                               colnames(predictor_lasso_zscores))
        assessment_summary <- merge_assessments(assessment_summary, assessment_lasso, "3")
      } else {
        cat("  Lasso selected no features\n")
      }
    }, error = function(e) {
      cat("  Lasso failed:", conditionMessage(e), "\n")
    })
  } else {
    cat("Step 3: Skipped (lasso disabled or no bio features)\n")
  }

  # ---- Fill NA for assessment columns that may be missing ---- #
  for (col in c("p1","beta1","se1","p2","beta2","se2","p3","beta3","se3")) {
    if (!col %in% names(assessment_summary)) assessment_summary[[col]] <- NA
  }

  # ====================================== #
  # ---- Summary Statistics ---- #
  # ====================================== #
  thresholds <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01,
                  0.001, 5e-04, 5e-06, 5e-08)

  sum_data <- data.frame(
    threshold = thresholds,
    qval = count_sig(target_results$qval, thresholds),
    fdr_the = count_sig(target_results$FDR.the, thresholds),
    fdr_emp = count_sig(target_results$FDR.emp, thresholds),
    bio_fdr_the = count_sig(target_results$bio.FDR.the, thresholds),
    bio_fdr_emp = count_sig(target_results$bio.FDR.emp, thresholds),
    bio_fdr_the_lasso = count_sig(target_results$bio.FDR.the.lasso, thresholds),
    bio_fdr_emp_lasso = count_sig(target_results$bio.FDR.emp.lasso, thresholds)
  )

  # ====================================== #
  # ---- Save Results ---- #
  # ====================================== #
  fwrite(target_results,
         file = file.path(output_dir, paste0(target, ".gene.bio.fdrreg.txt")),
         sep = " ")

  fwrite(assessment_summary,
         file = file.path(output_dir, paste0(target, ".contribution.csv")),
         sep = ",")

  fwrite(sum_data,
         file = file.path(output_dir, paste0(target, ".summary.csv")),
         sep = ",")

  save(list = names(model_list), envir = list2env(model_list),
       file = file.path(output_dir, paste0(target, ".fdrreg_models.RData")))

  cat(sprintf("\n--- %s Summary ---\n", target))
  cat(sprintf("  Genes: %d | Traits: %d | Bio: %s | Lasso: %s\n",
              length(gene_ids), length(trait_names),
              !is.null(model_list$fdr_bio_theoretical),
              !is.null(model_list$fdr_lasso_theoretical)))
  print(sum_data)
  cat(sprintf("  Saved to: %s\n", output_dir))

  return(invisible(NULL))
}

# ======================== #
# ---- Main Entry Point -- #
# ======================== #

main <- function() {
  # ---- Argument Specification ---- #
  option_list <- list(
    make_option(c("--target"), type = "character", default = NULL,
                help = "Target(s) to analyze, comma-separated [required]",
                metavar = "TARGET1,TARGET2"),
    make_option(c("--base_dir"), type = "character",
                default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
                help = "Base directory [default: %default]"),
    make_option(c("--bio_library"), type = "character",
                default = Sys.getenv("FDRREG_BIO_ENSEMBL", ""),
                help = "Biological annotation file [default: %default]"),
    make_option(c("--bio"), type = "character", default = "Y",
                help = "Run bio annotation analysis? Y/N [default: %default]"),
    make_option(c("--lasso"), type = "character", default = "Y",
                help = "Run Lasso feature selection? Y/N [default: %default]"),
    make_option(c("--traits"), type = "character", default = NULL,
                help = "Trait(s) as predictors, comma-separated [default: all]",
                metavar = "TRAIT1,TRAIT2"),
    make_option(c("--seed"), type = "integer", default = 100,
                help = "Random seed [default: %default]"),
    make_option(c("--jobs"), type = "integer", default = 1,
                help = "Number of targets to run in parallel [default: %default]")
  )

  opt_parser <- OptionParser(
    option_list = option_list,
    usage = "Usage: %prog --target TARGET1[,TARGET2] [options]",
    description = paste(
      "SmultiXcan FDRreg Analysis Pipeline.",
      "Example: Rscript %prog --target bd2018,adhd2016 --bio Y --lasso Y --traits college2013"
    )
  )

  opt <- parse_args(opt_parser)

  # ---- Validate ---- #
  if (is.null(opt$target)) {
    print_help(opt_parser)
    stop("--target is required", call. = FALSE)
  }

  targets <- trimws(strsplit(opt$target, ",")[[1]])

  # Expand the special keyword "all": every target that has a v7 S-MultiXcan
  # output folder is picked up automatically.
  if (length(targets) == 1 && tolower(targets) == "all") {
    candidate_dirs <- list.dirs(opt$base_dir, recursive = FALSE)
    has_v7 <- dir.exists(file.path(candidate_dirs, "08.smultixcan_output_v7"))
    targets <- basename(candidate_dirs)[has_v7]
    if (length(targets) == 0) {
      stop("No targets with 08.smultixcan_output_v7 found under: ", opt$base_dir)
    }
    cat(sprintf("--target all expanded to %d targets: %s\n",
                length(targets), paste(targets, collapse = ", ")))
  }

  # Clamp parallelism to a sane range.
  n_jobs <- max(1L, min(opt$jobs, length(targets)))

  run_bio <- toupper(opt$bio) == "Y"
  run_lasso <- toupper(opt$lasso) == "Y"
  user_traits <- if (!is.null(opt$traits)) {
    trimws(strsplit(opt$traits, ",")[[1]])
  } else {
    NULL
  }

  # Validate base directory
  if (!dir.exists(opt$base_dir)) {
    stop("Base directory does not exist: ", opt$base_dir)
  }

  # Validate bio library
  bio_library <- NULL
  if (run_bio) {
    if (!file.exists(opt$bio_library)) {
      warning("Bio library not found, disabling bio analysis: ", opt$bio_library)
      run_bio <- FALSE
    } else {
      cat("Loading biological annotation...\n")
      bio_library <- fread(opt$bio_library, showProgress = FALSE)
    }
  }

  # ---- Print Configuration ---- #
  cat("\n")
  cat("=========================================\n")
  cat(" SmultiXcan FDRreg Pipeline\n")
  cat("=========================================\n")
  cat(sprintf("  Targets      : %s\n", paste(targets, collapse = ", ")))
  cat(sprintf("  Base dir     : %s\n", opt$base_dir))
  cat(sprintf("  Bio library  : %s\n", opt$bio_library))
  cat(sprintf("  Bio analysis : %s\n", run_bio))
  cat(sprintf("  Lasso        : %s\n", run_lasso))
  cat(sprintf("  Traits       : %s\n",
              ifelse(is.null(user_traits), "all (auto-detect)", paste(user_traits, collapse = ", "))))
  cat(sprintf("  Seed         : %d\n", opt$seed))
  cat("=========================================\n")

  # ---- Run each target (optionally in parallel) ---- #
  # One worker per target. Errors are caught per target so one failure does
  # not abort the rest of the pool.
  run_one <- function(target) {
    tryCatch({
      process_target(
        target = target,
        user_traits = user_traits,
        run_bio = run_bio,
        run_lasso = run_lasso,
        base_dir = opt$base_dir,
        bio_library = bio_library,
        seed = opt$seed
      )
    }, error = function(e) {
      cat(sprintf("\n[ERROR] Target %s failed: %s\n", target, conditionMessage(e)))
    })
    invisible(NULL)
  }

  cat(sprintf("  Parallel jobs: %d\n", n_jobs))
  cat("=========================================\n")

  if (n_jobs > 1 && length(targets) > 1) {
    suppressPackageStartupMessages(library(parallel))
    # mc.preschedule = FALSE gives more even load when target runtimes vary.
    parallel::mclapply(targets, run_one,
                       mc.cores = n_jobs, mc.preschedule = FALSE)
  } else {
    for (target in targets) run_one(target)
  }

  cat("\n========== All targets completed ==========\n")
}

# Run
main()
