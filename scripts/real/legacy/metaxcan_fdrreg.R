#!/usr/bin/env Rscript

#**********************************#
#-----FDR_REG: metaxcan-fdrreg-----#
#**********************************#
# R-4.0.2
# Optimized: 2026/06/02
# Modified for GTEx v7 + multi-target parallel execution.

#---- Argument Parsing ----#
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-t", "--target"),
    type = "character",
    default = NULL,
    help = paste(
      "Target trait name(s). Comma-separated for several,",
      "or 'all' to auto-detect every target under input_dir. [REQUIRED]"
    )
  ),
  make_option(
    c("-j", "--jobs"),
    type = "integer",
    default = 1,
    help = "Number of targets to process in parallel. [default: %default]"
  ),
  make_option(
    c("--lasso"),
    type = "character",
    default = "Y",
    help = "Perform LASSO feature selection (Y/N). [default: %default]"
  ),
  make_option(
    c("--traits"),
    type = "character",
    default = NULL,
    help = "Comma-separated list of variable traits. If not specified, auto-detect. [default: NULL]"
  ),
  make_option(
    c("--bio_annotation"),
    type = "character",
    default = "Y",
    help = "Use biological annotations (Y/N). [default: %default]"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 100,
    help = "Random seed for reproducibility. [default: %default]"
  ),
  make_option(
    c("--input_dir"),
    type = "character",
    default = "/path/to/SO_Lab/18.fdrreg_rebuild",
    help = "Base input directory. [default: %default]"
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = NULL,
    help = "Base output directory. [default: NULL]"
  ),
  make_option(
    c("--bio_file"),
    type = "character",
    default = "/path/to/global.files/magma-library-uniq-ensembl.csv",
    help = "Path to biological annotation file. [default: %default]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "FDRreg analysis for MetaXcan v7 results"
)

opt <- parse_args(opt_parser)

if (is.null(opt$target)) {
  stop("--target is required. See: Rscript metaxcan_fdrreg.R --help")
}

perform_lasso <- toupper(opt$lasso) == "Y"
use_bio_annotation <- toupper(opt$bio_annotation) == "Y"

# Subdirectory that holds the S-PrediXcan v7 outputs for each target.
METAXCAN_SUBDIR <- "06.metaxcan"

# ---- Resolve the list of targets (supports comma-separated and 'all') ----
raw_targets <- trimws(strsplit(opt$target, ",")[[1]])
raw_targets <- raw_targets[nzchar(raw_targets)]

if (length(raw_targets) == 1 && tolower(raw_targets) == "all") {
  candidate_dirs <- list.dirs(
    opt$input_dir,
    full.names = TRUE,
    recursive = FALSE
  )

  TARGETS <- basename(
    candidate_dirs[
      dir.exists(file.path(candidate_dirs, METAXCAN_SUBDIR))
    ]
  )

  if (length(TARGETS) == 0) {
    stop(
      "No targets with a '", METAXCAN_SUBDIR,
      "' folder were found under: ", opt$input_dir
    )
  }
} else {
  TARGETS <- unique(raw_targets)
}

cat("Targets to process:", paste(TARGETS, collapse = ", "), "\n")
cat("Parallel jobs:", opt$jobs, "\n\n")

#---- Load packages ----#
suppressPackageStartupMessages({
  library(FDRreg)
  library(data.table)
  library(powerplus)
  library(stringr)
  library(dplyr)
  library(glmnet)
  library(HelpersMG)
  library(parallel)
})

#---- Helper Functions ----#

count_below_threshold <- function(values, threshold) {
  sum(values < threshold, na.rm = TRUE)
}

extract_fdr_contribution_safe <- function(
  fdr_model,
  feature_names,
  model_name
) {
  tryCatch({
    if (is.null(fdr_model)) {
      warning("FDR model is NULL for: ", model_name)
      return(NULL)
    }

    if (is.null(fdr_model$model$hessian)) {
      warning("Model hessian is NULL for: ", model_name)
      return(NULL)
    }

    features_se <- SEfromHessian(fdr_model$model$hessian)
    features_coef <- fdr_model$model$coef

    if (length(features_coef) <= 1) {
      warning("No covariates found in model: ", model_name)
      return(NULL)
    }

    if (length(feature_names) == 0) {
      warning("No feature names available for: ", model_name)
      return(NULL)
    }

    coef_idx <- seq.int(2, length(features_coef))

    n_features <- min(
      length(coef_idx),
      length(features_se) - 1L,
      length(feature_names)
    )

    if (n_features <= 0) {
      warning("No valid feature coefficients found for: ", model_name)
      return(NULL)
    }

    selected_coef_idx <- coef_idx[seq_len(n_features)]

    features_z_score <- (
      features_coef[selected_coef_idx] /
        features_se[selected_coef_idx]
    )

    features_pvalue <- 2 * pnorm(
      abs(features_z_score),
      lower.tail = FALSE
    )

    assessment <- data.frame(
      feature = feature_names[seq_len(n_features)],
      pvalue = features_pvalue,
      beta = features_coef[selected_coef_idx],
      se = features_se[selected_coef_idx],
      stringsAsFactors = FALSE
    )

    assessment$model <- model_name

    assessment <- assessment[
      order(assessment$pvalue),
      ,
      drop = FALSE
    ]

    return(assessment)

  }, error = function(e) {
    warning(
      "Failed to extract contribution for ",
      model_name,
      ": ",
      conditionMessage(e)
    )
    return(NULL)
  })
}

extract_trait_from_filename <- function(filename, brain_region) {
  trait_name <- basename(filename)

  # v8 S-PrediXcan outputs are named:
  #   gtex_v8_<trait>.overlap.4magma_in_<region>.csv
  trait_name <- sub("^gtex_v8_", "", trait_name)

  suffix_pattern <- paste0(
    "\\.overlap\\.4magma_in_",
    brain_region,
    "\\.csv$"
  )

  trait_name <- sub(
    suffix_pattern,
    "",
    trait_name
  )

  return(trait_name)
}

read_metaxcan_file <- function(filename, trait_name) {
  dt <- fread(
    filename,
    data.table = FALSE
  )

  required_columns <- c(
    "gene",
    "gene_name",
    "zscore",
    "pvalue"
  )

  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in file: ",
      filename,
      ". Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dt$gene <- as.character(dt$gene)
  dt$zscore <- suppressWarnings(
    as.numeric(dt$zscore)
  )
  dt$pvalue <- suppressWarnings(
    as.numeric(dt$pvalue)
  )

  valid_rows <- (
    !is.na(dt$gene) &
      nzchar(dt$gene) &
      is.finite(dt$zscore) &
      is.finite(dt$pvalue)
  )

  n_removed <- sum(!valid_rows)

  if (n_removed > 0) {
    cat(
      "  Removed",
      n_removed,
      "invalid rows from trait:",
      trait_name,
      "\n"
    )
  }

  dt <- dt[
    valid_rows,
    ,
    drop = FALSE
  ]

  if (nrow(dt) == 0) {
    stop(
      "No valid rows remain after filtering file: ",
      filename
    )
  }

  duplicated_gene_count <- sum(
    duplicated(dt$gene)
  )

  if (duplicated_gene_count > 0) {
    duplicated_genes <- unique(
      dt$gene[duplicated(dt$gene)]
    )

    stop(
      "Duplicated gene IDs detected in file: ",
      filename,
      ". Number of duplicated gene IDs: ",
      length(duplicated_genes),
      ". Examples: ",
      paste(
        head(duplicated_genes, 10),
        collapse = ", "
      )
    )
  }

  dt <- dt[
    order(dt$gene),
    ,
    drop = FALSE
  ]

  return(dt)
}

align_to_common_genes <- function(
  dt,
  common_gene_order,
  file_label
) {
  matched_indices <- match(
    common_gene_order,
    dt$gene
  )

  if (anyNA(matched_indices)) {
    stop(
      "Internal alignment failure for ",
      file_label,
      ": ",
      sum(is.na(matched_indices)),
      " common genes were not found."
    )
  }

  aligned_dt <- dt[
    matched_indices,
    ,
    drop = FALSE
  ]

  if (!identical(
    as.character(aligned_dt$gene),
    as.character(common_gene_order)
  )) {
    stop(
      "Gene order validation failed for: ",
      file_label
    )
  }

  return(aligned_dt)
}

generate_summary_row <- function(
  TARGET,
  brain_region,
  n_genes,
  n_traits,
  n_bio_features,
  n_lasso_features,
  has_bio_model,
  has_lasso_model,
  gene_results,
  thresholds
) {
  summary_row <- data.frame(
    target = TARGET,
    region = brain_region,
    total_genes = n_genes,
    n_traits = n_traits,
    n_bio_features = n_bio_features,
    n_lasso_features = n_lasso_features,
    has_bio_model = has_bio_model,
    has_lasso_model = has_lasso_model,
    stringsAsFactors = FALSE
  )

  all_fdr_methods <- c(
    "qval",
    "FDR_theoretical",
    "FDR_empirical",
    "bio_FDR_theoretical",
    "bio_FDR_empirical",
    "lasso_FDR_theoretical",
    "lasso_FDR_empirical"
  )

  for (method in all_fdr_methods) {
    for (thr in thresholds) {
      col_name <- paste(
        method,
        format(thr, scientific = FALSE),
        sep = "_"
      )

      if (
        method %in% names(gene_results) &&
          !all(is.na(gene_results[[method]]))
      ) {
        summary_row[[col_name]] <- count_below_threshold(
          gene_results[[method]],
          thr
        )
      } else {
        summary_row[[col_name]] <- NA_integer_
      }
    }
  }

  return(summary_row)
}

#---- Per-target analysis function ----#

run_target_analysis <- function(
  TARGET,
  opt,
  perform_lasso,
  use_bio_annotation,
  scz.gene.david,
  brain_regions,
  thresholds,
  metaxcan_subdir
) {
  # Per-fork reproducibility.
  set.seed(opt$seed)

  #---- Configuration (per target) ----#
  input_base_path <- file.path(
    opt$input_dir,
    TARGET,
    metaxcan_subdir
  )

  if (is.null(opt$output_dir)) {
    output_base <- file.path(opt$input_dir, TARGET, "07.metaxcan_fdrreg")
  } else {
    output_base <- file.path(opt$output_dir, TARGET, "07.metaxcan_fdrreg")
  }

  dir.create(file.path(output_base, "01.fdrreg_results"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_base, "02.contribution"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_base, "03.models"),
             showWarnings = FALSE, recursive = TRUE)

  cat("============================================================\n")
  cat("FDRreg Analysis Configuration\n")
  cat("============================================================\n")
  cat("Target:", TARGET, "\n")
  cat("Perform LASSO:", perform_lasso, "\n")
  cat("Use Biological Annotations:", use_bio_annotation, "\n")
  cat("Random Seed:", opt$seed, "\n")
  cat("Input Directory:", input_base_path, "\n")
  cat("Output Directory:", output_base, "\n\n")

  final_sum <- NULL
  final_contrib <- NULL

  cat(
    "Starting FDRreg analysis for",
    TARGET,
    "across",
    length(brain_regions),
    "brain regions...\n\n"
  )

  #---- Process each brain region ----#

  for (brain_region in brain_regions) {
    cat("\n========================================\n")
    cat("Processing brain region:", brain_region, "\n")
    cat("========================================\n")

    region_dir <- file.path(
      output_base,
      "01.fdrreg_results",
      brain_region
    )

    dir.create(
      region_dir,
      showWarnings = FALSE,
      recursive = TRUE
    )

    pattern <- paste0(
      "_in_",
      brain_region,
      "\\.csv$"
    )

    filenames <- list.files(
      path = input_base_path,
      pattern = pattern,
      full.names = TRUE
    )

    if (length(filenames) == 0) {
      warning(
        "No files found for brain region: ",
        brain_region
      )
      next
    }

    all_trait_names <- vapply(
      filenames,
      extract_trait_from_filename,
      character(1),
      brain_region = brain_region
    )

    target_idx <- which(
      all_trait_names == TARGET
    )

    if (length(target_idx) != 1) {
      warning(
        "Could not find unique target file for ",
        TARGET,
        " in brain region: ",
        brain_region
      )

      cat(
        "Available traits:",
        paste(all_trait_names, collapse = ", "),
        "\n"
      )

      next
    }

    if (!is.null(opt$traits)) {
      specified_traits <- trimws(
        strsplit(opt$traits, ",")[[1]]
      )

      specified_traits <- specified_traits[
        specified_traits != ""
      ]

      trait_indices <- setdiff(
        which(all_trait_names %in% specified_traits),
        target_idx
      )

      if (length(trait_indices) == 0) {
        warning(
          "None of the specified non-target traits found in brain region: ",
          brain_region
        )
        next
      }

      file_indices <- c(
        target_idx,
        trait_indices
      )

      trait_names_without_target <- all_trait_names[
        trait_indices
      ]
    } else {
      trait_indices <- setdiff(
        seq_along(all_trait_names),
        target_idx
      )

      if (length(trait_indices) == 0) {
        warning(
          "No variable traits found for brain region: ",
          brain_region
        )
        next
      }

      file_indices <- c(
        target_idx,
        trait_indices
      )

      trait_names_without_target <- all_trait_names[
        trait_indices
      ]
    }

    if (anyDuplicated(trait_names_without_target)) {
      duplicated_traits <- unique(
        trait_names_without_target[
          duplicated(trait_names_without_target)
        ]
      )

      stop(
        "Duplicated trait names detected in ",
        brain_region,
        ": ",
        paste(duplicated_traits, collapse = ", ")
      )
    }

    selected_filenames <- filenames[
      file_indices
    ]

    selected_trait_names <- all_trait_names[
      file_indices
    ]

    cat(
      "Loading",
      length(selected_filenames),
      "files (1 target +",
      length(trait_indices),
      "traits)...\n"
    )

    files <- lapply(
      seq_along(selected_filenames),
      function(i) {
        read_metaxcan_file(
          selected_filenames[i],
          selected_trait_names[i]
        )
      }
    )

    target_file <- files[[1]]
    trait_files <- files[-1]

    if (length(trait_files) == 0) {
      warning(
        "No trait files found for brain region: ",
        brain_region
      )
      next
    }

    # Align every file using the intersection across all files.
    cat("Aligning gene lists across all files...\n")

    all_gene_lists <- c(
      list(target_file$gene),
      lapply(
        trait_files,
        function(x) x$gene
      )
    )

    common_genes <- Reduce(
      intersect,
      all_gene_lists
    )

    # Preserve the target file gene order.
    common_genes <- target_file$gene[
      target_file$gene %in% common_genes
    ]

    cat(
      "  Target genes before alignment:",
      nrow(target_file),
      "\n"
    )

    cat(
      "  Common genes across all files:",
      length(common_genes),
      "\n"
    )

    if (length(common_genes) == 0) {
      warning(
        "No common genes found across all files for brain region: ",
        brain_region
      )
      next
    }

    if (length(common_genes) < 2) {
      warning(
        "Fewer than two common genes remain for brain region: ",
        brain_region
      )
      next
    }

    original_gene_counts <- c(
      target = nrow(target_file),
      vapply(
        trait_files,
        nrow,
        integer(1)
      )
    )

    if (any(original_gene_counts != length(common_genes))) {
      warning(
        "Gene lists were not identical before alignment in: ",
        brain_region,
        ". All files will be restricted to ",
        length(common_genes),
        " common genes."
      )
    }

    target_file <- align_to_common_genes(
      target_file,
      common_genes,
      paste0(TARGET, " target")
    )

    trait_files <- lapply(
      seq_along(trait_files),
      function(i) {
        align_to_common_genes(
          trait_files[[i]],
          common_genes,
          trait_names_without_target[i]
        )
      }
    )

    aligned_ok <- vapply(
      trait_files,
      function(x) {
        identical(
          as.character(x$gene),
          as.character(target_file$gene)
        )
      },
      logical(1)
    )

    if (!all(aligned_ok)) {
      failed_traits <- trait_names_without_target[
        !aligned_ok
      ]

      stop(
        "Gene alignment failed for traits: ",
        paste(failed_traits, collapse = ", ")
      )
    }

    # Extract target z-scores.
    target_z <- as.numeric(
      target_file$zscore
    )

    if (length(target_z) != length(common_genes)) {
      stop(
        "Target z-score length does not match the common gene count in: ",
        brain_region
      )
    }

    # Construct a guaranteed two-dimensional feature matrix.
    feature_vectors <- lapply(
      trait_files,
      function(x) {
        abs(as.numeric(x$zscore))
      }
    )

    feature_lengths <- vapply(
      feature_vectors,
      length,
      integer(1)
    )

    if (any(feature_lengths != length(common_genes))) {
      stop(
        "Feature vector lengths are inconsistent after alignment in: ",
        brain_region
      )
    }

    features_abs <- do.call(
      cbind,
      feature_vectors
    )

    if (is.null(dim(features_abs))) {
      features_abs <- matrix(
        features_abs,
        nrow = length(common_genes),
        ncol = length(feature_vectors)
      )
    }

    storage.mode(features_abs) <- "double"

    colnames(features_abs) <- trait_names_without_target
    rownames(features_abs) <- common_genes

    if (!all(
      dim(features_abs) ==
        c(
          length(common_genes),
          length(trait_names_without_target)
        )
    )) {
      stop(
        "Unexpected feature matrix dimension in ",
        brain_region,
        ". Observed: ",
        paste(dim(features_abs), collapse = " x ")
      )
    }

    if (any(!is.finite(features_abs))) {
      stop(
        "Non-finite values detected in feature matrix for: ",
        brain_region
      )
    }

    # Remove constant trait features before FDRreg.
    trait_sd <- apply(
      features_abs,
      2,
      sd,
      na.rm = TRUE
    )

    keep_trait_features <- (
      is.finite(trait_sd) &
        trait_sd > 0
    )

    if (!all(keep_trait_features)) {
      removed_traits <- trait_names_without_target[
        !keep_trait_features
      ]

      warning(
        "Removing constant trait features in ",
        brain_region,
        ": ",
        paste(removed_traits, collapse = ", ")
      )

      features_abs <- features_abs[
        ,
        keep_trait_features,
        drop = FALSE
      ]

      trait_names_without_target <- trait_names_without_target[
        keep_trait_features
      ]
    }

    if (ncol(features_abs) == 0) {
      warning(
        "No valid non-constant trait features remain for: ",
        brain_region
      )
      next
    }

    cat(
      "  Final feature matrix:",
      nrow(features_abs),
      "genes x",
      ncol(features_abs),
      "traits\n"
    )

    # Clean up gene IDs.
    gene_ids <- target_file$gene

    gene_ids_clean <- gsub(
      "\\..*",
      "",
      gene_ids
    )

    # Save gene list.
    fwrite(
      data.frame(
        gene_list = gene_ids_clean,
        stringsAsFactors = FALSE
      ),
      file.path(
        region_dir,
        paste0(
          brain_region,
          ".gene.list.txt"
        )
      )
    )

    #----------------------------------------------#
    #---- Step 1: Basic FDRreg Analysis ----------#
    #----------------------------------------------#

    cat("Running basic FDRreg analysis...\n")

    fdr_theoretical <- FDRreg(
      target_z,
      features_abs,
      nulltype = "theoretical",
      method = "pr"
    )

    fdr_empirical <- FDRreg(
      target_z,
      features_abs,
      nulltype = "empirical",
      method = "pr"
    )

    #----------------------------------------------#
    #---- Step 2: Biological Information Analysis #
    #----------------------------------------------#

    bio_features <- NULL
    bio_features_selected <- NULL
    bio_feature_names <- character(0)

    if (
      use_bio_annotation &&
        !is.null(scz.gene.david)
    ) {
      cat("Integrating biological information...\n")

      annotation_ids <- as.character(
        scz.gene.david$ENSEMBL_GENE_ID
      )

      annotation_match <- match(
        gene_ids_clean,
        annotation_ids
      )

      bio_info <- scz.gene.david[
        annotation_match,
        ,
        drop = FALSE
      ]

      bio_cols <- setdiff(
        names(scz.gene.david),
        c(
          "ENSEMBL_GENE_ID",
          "ID",
          "gene_name"
        )
      )

      numeric_bio_cols <- vapply(
        bio_cols,
        function(col_name) {
          raw_values <- bio_info[[col_name]]

          converted_values <- suppressWarnings(
            as.numeric(as.character(raw_values))
          )

          nonmissing_values <- (
            !is.na(raw_values) &
              nzchar(trimws(as.character(raw_values)))
          )

          all(
            !nonmissing_values |
              !is.na(converted_values)
          )
        },
        logical(1)
      )

      bio_cols <- bio_cols[
        numeric_bio_cols
      ]

      if (length(bio_cols) > 0) {
        bio_feature_vectors <- lapply(
          bio_cols,
          function(col_name) {
            values <- suppressWarnings(
              as.numeric(
                as.character(
                  bio_info[[col_name]]
                )
              )
            )

            values[is.na(values)] <- 0
            return(values)
          }
        )

        bio_features <- do.call(
          cbind,
          bio_feature_vectors
        )

        if (is.null(dim(bio_features))) {
          bio_features <- matrix(
            bio_features,
            nrow = length(gene_ids_clean),
            ncol = length(bio_cols)
          )
        }

        storage.mode(bio_features) <- "double"

        colnames(bio_features) <- bio_cols
        rownames(bio_features) <- gene_ids_clean

        bio_feature_names <- bio_cols

        # Remove constant biological features.
        bio_sd <- apply(
          bio_features,
          2,
          sd,
          na.rm = TRUE
        )

        keep_bio_features <- (
          is.finite(bio_sd) &
            bio_sd > 0
        )

        if (!all(keep_bio_features)) {
          removed_bio_features <- bio_feature_names[
            !keep_bio_features
          ]

          warning(
            "Removing constant biological features in ",
            brain_region,
            ": ",
            paste(
              removed_bio_features,
              collapse = ", "
            )
          )

          bio_features <- bio_features[
            ,
            keep_bio_features,
            drop = FALSE
          ]

          bio_feature_names <- bio_feature_names[
            keep_bio_features
          ]
        }

        if (ncol(bio_features) == 0) {
          bio_features <- NULL
          bio_feature_names <- character(0)
          cat(
            "  No non-constant biological features remain\n"
          )
        } else {
          cat(
            "  Found",
            ncol(bio_features),
            "biological features\n"
          )
        }
      } else {
        cat(
          "  No numeric biological features found\n"
        )
      }
    }

    #----------------------------------------------#
    #---- Step 3: LASSO Feature Selection --------#
    #----------------------------------------------#

    selected_lasso_features <- character(0)

    if (
      perform_lasso &&
        !is.null(bio_features) &&
        ncol(bio_features) > 0
    ) {
      cat("Performing LASSO feature selection...\n")

      tryCatch({
        lasso_cv <- cv.glmnet(
          x = bio_features,
          y = abs(target_z),
          family = "gaussian",
          nlambda = 50,
          alpha = 1,
          standardize = TRUE,
          parallel = FALSE
        )

        lasso_coef <- as.matrix(
          coef(
            lasso_cv,
            s = "lambda.min"
          )
        )

        coef_values <- lasso_coef[
          -1,
          1
        ]

        selected_idx <- which(
          coef_values != 0
        )

        if (length(selected_idx) > 0) {
          selected_lasso_features <- bio_feature_names[
            selected_idx
          ]

          bio_features_selected <- bio_features[
            ,
            selected_idx,
            drop = FALSE
          ]

          fwrite(
            data.frame(
              feature = selected_lasso_features,
              coefficient = coef_values[selected_idx],
              stringsAsFactors = FALSE
            ),
            file.path(
              region_dir,
              paste0(
                brain_region,
                ".lasso_selected_features.csv"
              )
            )
          )

          cat(
            "  LASSO selected",
            length(selected_lasso_features),
            "biological features\n"
          )
        } else {
          cat(
            "  No biological features selected by LASSO\n"
          )
        }
      }, error = function(e) {
        cat(
          "  LASSO failed:",
          conditionMessage(e),
          "\n"
        )
      })
    }

    #----------------------------------------------#
    #---- Step 4: Run Additional FDRreg Models ---#
    #----------------------------------------------#

    has_bio_model <- FALSE
    has_lasso_model <- FALSE

    fdr_bio_theoretical <- NULL
    fdr_bio_empirical <- NULL
    fdr_lasso_theoretical <- NULL
    fdr_lasso_empirical <- NULL

    bio_model_failed <- FALSE
    lasso_model_failed <- FALSE

    if (
      !is.null(bio_features) &&
        ncol(bio_features) > 0
    ) {
      cat(
        "Running FDRreg with all biological features...\n"
      )

      combined_features <- cbind(
        abs(features_abs),
        abs(bio_features)
      )

      combined_feature_names <- c(
        trait_names_without_target,
        bio_feature_names
      )

      colnames(combined_features) <- combined_feature_names

      tryCatch({
        capture.output({
          fdr_bio_theoretical <- FDRreg(
            target_z,
            combined_features,
            nulltype = "theoretical",
            method = "pr"
          )

          fdr_bio_empirical <- FDRreg(
            target_z,
            combined_features,
            nulltype = "empirical",
            method = "pr"
          )
        }, file = nullfile())

        if (
          length(fdr_bio_theoretical$model$coef) > 1
        ) {
          has_bio_model <- TRUE
          cat(
            "  Bio model fitted successfully\n"
          )
        } else {
          cat(
            "  Bio model reverted to no-covariates model\n"
          )
          bio_model_failed <- TRUE
        }
      }, error = function(e) {
        cat(
          "  Bio model failed:",
          conditionMessage(e),
          "\n"
        )
        bio_model_failed <- TRUE
      })
    }

    if (
      !is.null(bio_features_selected) &&
        ncol(bio_features_selected) > 0
    ) {
      cat(
        "Running FDRreg with LASSO-selected biological features...\n"
      )

      lasso_features <- cbind(
        abs(features_abs),
        abs(bio_features_selected)
      )

      lasso_feature_names <- c(
        trait_names_without_target,
        selected_lasso_features
      )

      colnames(lasso_features) <- lasso_feature_names

      tryCatch({
        capture.output({
          fdr_lasso_theoretical <- FDRreg(
            target_z,
            lasso_features,
            nulltype = "theoretical",
            method = "pr"
          )

          fdr_lasso_empirical <- FDRreg(
            target_z,
            lasso_features,
            nulltype = "empirical",
            method = "pr"
          )
        }, file = nullfile())

        if (
          length(fdr_lasso_theoretical$model$coef) > 1
        ) {
          has_lasso_model <- TRUE
          cat(
            "  LASSO model fitted successfully\n"
          )
        } else {
          cat(
            "  LASSO model reverted to no-covariates model\n"
          )
          lasso_model_failed <- TRUE
        }
      }, error = function(e) {
        cat(
          "  LASSO model failed:",
          conditionMessage(e),
          "\n"
        )
        lasso_model_failed <- TRUE
      })
    }

    #----------------------------------------------#
    #---- Step 5: Extract Model Contributions ----#
    #----------------------------------------------#

    cat("Extracting model contributions...\n")

    all_contributions <- extract_fdr_contribution_safe(
      fdr_theoretical,
      trait_names_without_target,
      "basic_theoretical"
    )

    if (
      has_bio_model &&
        !bio_model_failed
    ) {
      contrib_bio <- extract_fdr_contribution_safe(
        fdr_bio_theoretical,
        c(
          trait_names_without_target,
          bio_feature_names
        ),
        "bio_theoretical"
      )

      if (!is.null(contrib_bio)) {
        all_contributions <- rbind(
          all_contributions,
          contrib_bio
        )
      }
    }

    if (
      has_lasso_model &&
        !lasso_model_failed
    ) {
      contrib_lasso <- extract_fdr_contribution_safe(
        fdr_lasso_theoretical,
        c(
          trait_names_without_target,
          selected_lasso_features
        ),
        "lasso_theoretical"
      )

      if (!is.null(contrib_lasso)) {
        all_contributions <- rbind(
          all_contributions,
          contrib_lasso
        )
      }
    }

    #----------------------------------------------#
    #---- Step 6: Create Results DataFrame -------#
    #----------------------------------------------#

    cat("Creating results summary...\n")

    gene_results <- data.frame(
      ENSEMBL_GENE_ID = gene_ids_clean,
      gene_name = target_file$gene_name,
      zscore = target_z,
      pvalue = target_file$pvalue,
      stringsAsFactors = FALSE
    )

    gene_results$qval <- p.adjust(
      gene_results$pvalue,
      method = "fdr"
    )

    gene_results$FDR_theoretical <- fdr_theoretical$FDR
    gene_results$FDR_empirical <- fdr_empirical$FDR

    if (
      has_bio_model &&
        !bio_model_failed
    ) {
      gene_results$bio_FDR_theoretical <-
        fdr_bio_theoretical$FDR

      gene_results$bio_FDR_empirical <-
        fdr_bio_empirical$FDR
    } else {
      gene_results$bio_FDR_theoretical <- NA_real_
      gene_results$bio_FDR_empirical <- NA_real_
    }

    if (
      has_lasso_model &&
        !lasso_model_failed
    ) {
      gene_results$lasso_FDR_theoretical <-
        fdr_lasso_theoretical$FDR

      gene_results$lasso_FDR_empirical <-
        fdr_lasso_empirical$FDR
    } else {
      gene_results$lasso_FDR_theoretical <- NA_real_
      gene_results$lasso_FDR_empirical <- NA_real_
    }

    #----------------------------------------------#
    #---- Step 7: Generate Summary Statistics ----#
    #----------------------------------------------#

    summary_row <- generate_summary_row(
      TARGET = TARGET,
      brain_region = brain_region,
      n_genes = nrow(gene_results),
      n_traits = length(trait_names_without_target),
      n_bio_features = ifelse(
        !is.null(bio_features),
        ncol(bio_features),
        0
      ),
      n_lasso_features = length(
        selected_lasso_features
      ),
      has_bio_model = (
        has_bio_model &&
          !bio_model_failed
      ),
      has_lasso_model = (
        has_lasso_model &&
          !lasso_model_failed
      ),
      gene_results = gene_results,
      thresholds = thresholds
    )

    if (is.null(final_sum)) {
      final_sum <- summary_row
    } else {
      final_sum <- rbind(
        final_sum,
        summary_row
      )
    }

    #----------------------------------------------#
    #---- Step 8: Save Results -------------------#
    #----------------------------------------------#

    cat(
      "Saving results for",
      brain_region,
      "...\n"
    )

    fwrite(
      gene_results,
      file.path(
        region_dir,
        paste0(
          brain_region,
          ".gene.bio.fdrreg.txt"
        )
      ),
      sep = " "
    )

    if (
      !is.null(all_contributions) &&
        nrow(all_contributions) > 0
    ) {
      fwrite(
        all_contributions,
        file.path(
          region_dir,
          paste0(
            brain_region,
            ".contribution.csv"
          )
        ),
        sep = ","
      )
    }

    fwrite(
      summary_row,
      file.path(
        region_dir,
        paste0(
          brain_region,
          ".summary.csv"
        )
      ),
      sep = ","
    )

    save(
      fdr_theoretical,
      fdr_empirical,
      fdr_bio_theoretical,
      fdr_bio_empirical,
      fdr_lasso_theoretical,
      fdr_lasso_empirical,
      has_bio_model,
      has_lasso_model,
      bio_model_failed,
      lasso_model_failed,
      file = file.path(
        region_dir,
        paste0(
          brain_region,
          ".fdrreg_models.RData"
        )
      )
    )

    if (
      !is.null(all_contributions) &&
        nrow(all_contributions) > 0
    ) {
      all_contributions$target <- TARGET
      all_contributions$region <- brain_region

      if (is.null(final_contrib)) {
        final_contrib <- all_contributions
      } else {
        final_contrib <- rbind(
          final_contrib,
          all_contributions
        )
      }
    }

    cat(
      "  Completed processing:",
      brain_region,
      "\n"
    )
  }

  #----------------------------------------------#
  #---- Step 9: Save Global Summary ------------#
  #----------------------------------------------#

  cat("\n========================================\n")
  cat(
    "Saving global summary results for",
    TARGET,
    "...\n"
  )
  cat("========================================\n")

  summary_file <- file.path(
    output_base,
    paste0(
      TARGET,
      ".all.summary.csv"
    )
  )

  contribution_file <- file.path(
    output_base,
    paste0(
      TARGET,
      ".all.contribution.csv"
    )
  )

  if (!is.null(final_sum)) {
    fwrite(
      final_sum,
      summary_file
    )

    cat(
      "  Saved summary:",
      nrow(final_sum),
      "brain regions\n"
    )
  } else {
    cat(
      "  No brain regions were successfully processed.\n"
    )
  }

  if (!is.null(final_contrib)) {
    fwrite(
      final_contrib,
      contribution_file
    )

    cat(
      "  Saved contributions:",
      nrow(final_contrib),
      "rows\n"
    )
  }

  processed_regions <- if (
    is.null(final_sum)
  ) {
    0L
  } else {
    nrow(final_sum)
  }

  cat("\n============================================================\n")
  cat(
    "Analysis complete for",
    TARGET,
    "!\n"
  )
  cat("============================================================\n")
  cat(
    "Processed brain regions:",
    processed_regions,
    "\n"
  )

  if (file.exists(summary_file)) {
    cat("  Summary:", summary_file, "\n")
  }

  if (file.exists(contribution_file)) {
    cat("  Contributions:", contribution_file, "\n")
  }

  cat(
    "Individual region results in:",
    file.path(output_base, "01.fdrreg_results"),
    "\n\n"
  )

  return(invisible(TRUE))
}

#---- Static definitions (built once) ----#

brain_regions <- c(
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra"
)

thresholds <- c(
  0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02,
  0.01, 0.001, 5e-04, 5e-06, 5e-08
)

#---- Load the biological annotation table once (shared across forks) ----#

scz.gene.david <- NULL

if (use_bio_annotation) {
  cat("Loading biological information database...\n")

  if (!file.exists(opt$bio_file)) {
    stop(
      "Biological annotation file does not exist: ",
      opt$bio_file
    )
  }

  scz.gene.david <- fread(
    opt$bio_file,
    data.table = FALSE
  )

  if (!"ENSEMBL_GENE_ID" %in% names(scz.gene.david)) {
    stop(
      "Biological annotation file must contain column: ENSEMBL_GENE_ID"
    )
  }

  scz.gene.david$ENSEMBL_GENE_ID <- gsub(
    "\\..*",
    "",
    as.character(scz.gene.david$ENSEMBL_GENE_ID)
  )

  duplicated_annotation <- duplicated(
    scz.gene.david$ENSEMBL_GENE_ID
  )

  if (any(duplicated_annotation)) {
    warning(
      "Duplicated ENSEMBL_GENE_ID values detected. ",
      "The first record for each gene will be used."
    )

    scz.gene.david <- scz.gene.david[
      !duplicated_annotation,
      ,
      drop = FALSE
    ]
  }
}

#---- Run targets, optionally in parallel ----#

run_one <- function(tgt) {
  tryCatch({
    run_target_analysis(
      TARGET = tgt,
      opt = opt,
      perform_lasso = perform_lasso,
      use_bio_annotation = use_bio_annotation,
      scz.gene.david = scz.gene.david,
      brain_regions = brain_regions,
      thresholds = thresholds,
      metaxcan_subdir = METAXCAN_SUBDIR
    )

    data.frame(
      target = tgt,
      status = "OK",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      target = tgt,
      status = paste0("FAILED: ", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })
}

n_cores <- max(1L, min(opt$jobs, length(TARGETS)))

if (n_cores > 1L) {
  results <- mclapply(
    TARGETS,
    run_one,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
} else {
  results <- lapply(TARGETS, run_one)
}

status_df <- do.call(rbind, results)

cat("\n==================== Run summary ====================\n")
print(status_df, row.names = FALSE)

quit(
  save = "no",
  status = 0,
  runLast = FALSE
)
