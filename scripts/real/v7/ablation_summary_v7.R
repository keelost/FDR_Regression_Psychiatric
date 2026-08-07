#!/usr/bin/env Rscript
#*************************************************************#
# Ablation sensitivity analysis summary
#
# This script extracts the number of significant genes at:
#   FDR < 0.01
#
# Original analyses:
#   1. q-value
#   2. Theoretical FDRreg without bio annotations
#   3. Theoretical FDRreg with bio annotations
#
# Sensitivity analyses:
#   1. SA1: all biological annotations only
#   2. SA2: selected disorders plus partial biological annotations
#
# Supported pipelines:
#   1. MAGMA
#   2. MetaXcan / S-PrediXcan, per brain region
#   3. S-MultiXcan
#
# Missing models or all-NA results are reported as NA.
#*************************************************************#

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ============================================================ #
# ---- Argument parsing -------------------------------------- #
# ============================================================ #

option_list <- list(
  make_option(
    c("--base_dir"),
    type = "character",
    default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
    help = "Base directory containing original analyses. [default: %default]"
  ),
  make_option(
    c("--ablation_dir"),
    type = "character",
    default = paste0(
      file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "10.ablation")
    ),
    help = "Ablation analysis directory. [default: %default]"
  ),
  make_option(
    c("--threshold"),
    type = "double",
    default = 0.01,
    help = "Significance threshold. The default and recommended value is 0.01. [default: %default]"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = paste0(
      file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "10.ablation", "ablation.summary_v7.csv")
    ),
    help = "Output summary CSV file. [default: %default]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = paste(
    "Generate a combined summary table for the original FDRreg analyses",
    "and the two ablation sensitivity analyses."
  )
)

opt <- parse_args(opt_parser)

base_dir <- normalizePath(
  opt$base_dir,
  winslash = "/",
  mustWork = FALSE
)

ablation_dir <- normalizePath(
  opt$ablation_dir,
  winslash = "/",
  mustWork = FALSE
)

threshold <- opt$threshold
output_file <- opt$output

if (!dir.exists(base_dir)) {
  stop("Base directory does not exist: ", base_dir)
}

if (!dir.exists(ablation_dir)) {
  stop("Ablation directory does not exist: ", ablation_dir)
}

if (!is.finite(threshold) || threshold <= 0 || threshold >= 1) {
  stop("--threshold must be a finite number between 0 and 1.")
}

# ============================================================ #
# ---- Helper functions -------------------------------------- #
# ============================================================ #

# Read a comma-separated summary file.
#
# fill = TRUE is explicitly specified to avoid fread warnings if a model
# generated an incomplete row. Missing cells will be filled with NA.
read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    suppressWarnings(
      fread(
        file = path,
        sep = ",",
        header = TRUE,
        fill = TRUE,
        na.strings = c("NA", "NaN", "NULL", ""),
        showProgress = FALSE,
        data.table = TRUE
      )
    ),
    error = function(e) {
      message(
        "[READ ERROR] ",
        path,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# Read an ablation gene-level result file.
#
# The diagnostic results show that these files are valid:
#   MAGMA:     5 columns
#   MetaXcan:  6 columns
#   S-MultiXcan: expected 6 columns
read_ablation_safe <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    suppressWarnings(
      fread(
        file = path,
        header = TRUE,
        fill = TRUE,
        na.strings = c("NA", "NaN", "NULL", ""),
        showProgress = FALSE,
        data.table = TRUE
      )
    ),
    error = function(e) {
      message(
        "[READ ERROR] ",
        path,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# Convert values to numeric without generating warnings.
to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

# Extract the first available column from candidate column names.
get_first_available_column <- function(dt, candidate_names) {
  if (is.null(dt) || nrow(dt) == 0L) {
    return(NULL)
  }

  available <- intersect(candidate_names, names(dt))

  if (length(available) == 0L) {
    return(NULL)
  }

  dt[[available[1L]]]
}

# Extract a single integer count from a summary table.
#
# Return NA when:
#   1. The table is missing.
#   2. The requested column is missing.
#   3. The value is NA or non-numeric.
extract_integer_value <- function(dt, candidate_names) {
  values <- get_first_available_column(dt, candidate_names)

  if (is.null(values) || length(values) == 0L) {
    return(NA_integer_)
  }

  value <- to_numeric_safe(values[1L])

  if (length(value) == 0L || is.na(value) || !is.finite(value)) {
    return(NA_integer_)
  }

  as.integer(round(value))
}

# Count values strictly below the requested threshold.
#
# This follows the original code:
#   sum(FDR < threshold, na.rm = TRUE)
#
# If the entire model result is missing or all values are NA, return NA
# rather than zero.
count_significant <- function(values, threshold) {
  if (is.null(values) || length(values) == 0L) {
    return(NA_integer_)
  }

  values <- to_numeric_safe(values)

  if (all(is.na(values))) {
    return(NA_integer_)
  }

  as.integer(sum(values < threshold, na.rm = TRUE))
}

# Check whether an existing model column contains any non-missing value.
has_valid_values <- function(dt, column_name) {
  if (is.null(dt) || !column_name %in% names(dt)) {
    return(FALSE)
  }

  values <- to_numeric_safe(dt[[column_name]])
  any(!is.na(values))
}

# Find the requested threshold row in a long-format summary.
#
# This is used for:
#   1. MAGMA
#   2. S-MultiXcan
find_threshold_row <- function(dt, threshold, tolerance = 1e-12) {
  if (is.null(dt) || nrow(dt) == 0L || !"threshold" %in% names(dt)) {
    return(NULL)
  }

  threshold_values <- to_numeric_safe(dt$threshold)

  matched <- which(
    !is.na(threshold_values) &
      is.finite(threshold_values) &
      abs(threshold_values - threshold) <= tolerance
  )

  if (length(matched) == 0L) {
    return(NULL)
  }

  dt[matched[1L]]
}

# Create possible representations of a threshold for wide-format columns.
#
# The diagnostic result confirms that MetaXcan uses names such as:
#   qval_0.01
#   FDR_theoretical_0.01
#   bio_FDR_theoretical_0.01
make_threshold_labels <- function(threshold) {
  labels <- c(
    as.character(threshold),
    format(
      threshold,
      scientific = FALSE,
      trim = TRUE,
      digits = 15
    ),
    sprintf("%.2f", threshold),
    sprintf("%.3f", threshold),
    sprintf("%.4f", threshold),
    sprintf("%.5f", threshold),
    sprintf("%.6f", threshold),
    sprintf("%.8f", threshold)
  )

  # Remove unnecessary trailing zeros, while retaining valid decimal forms.
  trimmed_labels <- sub(
    "\\.?0+$",
    "",
    labels
  )

  unique(c(labels, trimmed_labels))
}

# Build candidate column names for a wide-format summary.
make_wide_column_names <- function(prefixes, threshold) {
  threshold_labels <- make_threshold_labels(threshold)

  unique(
    unlist(
      lapply(
        prefixes,
        function(prefix) {
          paste(prefix, threshold_labels, sep = "_")
        }
      ),
      use.names = FALSE
    )
  )
}

# ============================================================ #
# ---- Extract ablation counts ------------------------------- #
# ============================================================ #

extract_ablation_results <- function(path, threshold) {
  dt <- read_ablation_safe(path)

  file_exists <- file.exists(path)
  file_readable <- !is.null(dt)

  if (is.null(dt)) {
    return(list(
      file_exists = file_exists,
      file_readable = FALSE,
      total_genes = NA_integer_,
      qval_count = NA_integer_,
      sa1_count = NA_integer_,
      sa2_count = NA_integer_,
      sa1_column_exists = FALSE,
      sa2_column_exists = FALSE,
      sa1_all_na = NA,
      sa2_all_na = NA
    ))
  }

  sa1_exists <- "SA1_FDR_the" %in% names(dt)
  sa2_exists <- "SA2_FDR_the" %in% names(dt)

  sa1_values <- if (sa1_exists) dt[["SA1_FDR_the"]] else NULL
  sa2_values <- if (sa2_exists) dt[["SA2_FDR_the"]] else NULL
  qval_values <- if ("qval" %in% names(dt)) dt[["qval"]] else NULL

  list(
    file_exists = file_exists,
    file_readable = file_readable,
    total_genes = as.integer(nrow(dt)),
    qval_count = count_significant(qval_values, threshold),
    sa1_count = count_significant(sa1_values, threshold),
    sa2_count = count_significant(sa2_values, threshold),
    sa1_column_exists = sa1_exists,
    sa2_column_exists = sa2_exists,
    sa1_all_na = if (!sa1_exists) NA else
      all(is.na(to_numeric_safe(sa1_values))),
    sa2_all_na = if (!sa2_exists) NA else
      all(is.na(to_numeric_safe(sa2_values)))
  )
}

# ============================================================ #
# ---- Extract original MAGMA/S-MultiXcan results ------------ #
# ============================================================ #

extract_original_long_summary <- function(path, threshold) {
  dt <- read_csv_safe(path)

  file_exists <- file.exists(path)
  file_readable <- !is.null(dt)

  if (is.null(dt)) {
    return(list(
      file_exists = file_exists,
      file_readable = FALSE,
      threshold_found = FALSE,
      total_genes = NA_integer_,
      qval_count = NA_integer_,
      fdr_the_count = NA_integer_,
      bio_fdr_the_count = NA_integer_,
      qval_column_exists = FALSE,
      fdr_the_column_exists = FALSE,
      bio_fdr_the_column_exists = FALSE
    ))
  }

  threshold_row <- find_threshold_row(dt, threshold)

  qval_exists <- "qval" %in% names(dt)
  fdr_the_exists <- any(
    c("fdr_the", "FDR_theoretical") %in% names(dt)
  )
  bio_fdr_the_exists <- any(
    c("bio_fdr_the", "bio_FDR_theoretical") %in% names(dt)
  )

  list(
    file_exists = file_exists,
    file_readable = file_readable,
    threshold_found = !is.null(threshold_row),

    # Long summary files do not store total_genes.
    total_genes = NA_integer_,

    qval_count = extract_integer_value(
      threshold_row,
      c("qval")
    ),

    fdr_the_count = extract_integer_value(
      threshold_row,
      c("fdr_the", "FDR_theoretical")
    ),

    bio_fdr_the_count = extract_integer_value(
      threshold_row,
      c("bio_fdr_the", "bio_FDR_theoretical")
    ),

    qval_column_exists = qval_exists,
    fdr_the_column_exists = fdr_the_exists,
    bio_fdr_the_column_exists = bio_fdr_the_exists
  )
}

# ============================================================ #
# ---- Extract original MetaXcan results --------------------- #
# ============================================================ #

extract_original_metaxcan_summary <- function(path, threshold) {
  dt <- read_csv_safe(path)

  file_exists <- file.exists(path)
  file_readable <- !is.null(dt)

  if (is.null(dt)) {
    return(list(
      file_exists = file_exists,
      file_readable = FALSE,
      threshold_found = FALSE,
      total_genes = NA_integer_,
      qval_count = NA_integer_,
      fdr_the_count = NA_integer_,
      bio_fdr_the_count = NA_integer_,
      qval_column_exists = FALSE,
      fdr_the_column_exists = FALSE,
      bio_fdr_the_column_exists = FALSE,
      has_bio_model = NA
    ))
  }

  qval_candidates <- make_wide_column_names(
    prefixes = c("qval"),
    threshold = threshold
  )

  fdr_the_candidates <- make_wide_column_names(
    prefixes = c(
      "FDR_theoretical",
      "fdr_the"
    ),
    threshold = threshold
  )

  bio_fdr_the_candidates <- make_wide_column_names(
    prefixes = c(
      "bio_FDR_theoretical",
      "bio_fdr_the"
    ),
    threshold = threshold
  )

  qval_exists <- any(qval_candidates %in% names(dt))
  fdr_the_exists <- any(fdr_the_candidates %in% names(dt))
  bio_fdr_the_exists <- any(
    bio_fdr_the_candidates %in% names(dt)
  )

  # The requested threshold is considered available if all three requested
  # result columns exist. A model-specific value can still be NA.
  threshold_found <- (
    qval_exists &&
      fdr_the_exists &&
      bio_fdr_the_exists
  )

  total_genes <- extract_integer_value(
    dt,
    c("total_genes")
  )

  has_bio_model <- NA
  if ("has_bio_model" %in% names(dt)) {
    bio_flag <- dt[["has_bio_model"]][1L]

    if (is.logical(bio_flag)) {
      has_bio_model <- bio_flag
    } else {
      bio_flag_text <- toupper(trimws(as.character(bio_flag)))
      has_bio_model <- if (bio_flag_text %in% c("TRUE", "T", "1", "Y")) {
        TRUE
      } else if (bio_flag_text %in% c("FALSE", "F", "0", "N")) {
        FALSE
      } else {
        NA
      }
    }
  }

  list(
    file_exists = file_exists,
    file_readable = file_readable,
    threshold_found = threshold_found,
    total_genes = total_genes,

    qval_count = extract_integer_value(
      dt,
      qval_candidates
    ),

    fdr_the_count = extract_integer_value(
      dt,
      fdr_the_candidates
    ),

    bio_fdr_the_count = extract_integer_value(
      dt,
      bio_fdr_the_candidates
    ),

    qval_column_exists = qval_exists,
    fdr_the_column_exists = fdr_the_exists,
    bio_fdr_the_column_exists = bio_fdr_the_exists,
    has_bio_model = has_bio_model
  )
}

# ============================================================ #
# ---- Create one output row --------------------------------- #
# ============================================================ #

create_summary_row <- function(
    pipeline,
    target,
    region,
    original_results,
    ablation_results,
    original_path,
    ablation_path,
    threshold) {

  # Prefer total_genes from original summary when available.
  # Otherwise use the number of genes in the ablation result.
  original_total_genes <- original_results$total_genes
  ablation_total_genes <- ablation_results$total_genes

  data.table(
    pipeline = pipeline,
    target = target,
    region = region,
    threshold = threshold,

    original_total_genes = original_total_genes,
    ablation_total_genes = ablation_total_genes,

    # Original main-analysis results
    original_qval = original_results$qval_count,
    original_fdr_the = original_results$fdr_the_count,
    original_bio_fdr_the = original_results$bio_fdr_the_count,

    # Sensitivity-analysis results
    SA1_all_bio_only_fdr_the = ablation_results$sa1_count,
    SA2_disorders_partial_bio_fdr_the = ablation_results$sa2_count,

    # File availability
    original_file_exists = original_results$file_exists,
    original_file_readable = original_results$file_readable,
    original_threshold_found = original_results$threshold_found,

    ablation_file_exists = ablation_results$file_exists,
    ablation_file_readable = ablation_results$file_readable,

    # Original requested columns
    original_qval_column_exists =
      original_results$qval_column_exists,
    original_fdr_the_column_exists =
      original_results$fdr_the_column_exists,
    original_bio_fdr_the_column_exists =
      original_results$bio_fdr_the_column_exists,

    # Ablation model status
    SA1_column_exists = ablation_results$sa1_column_exists,
    SA2_column_exists = ablation_results$sa2_column_exists,
    SA1_all_na = ablation_results$sa1_all_na,
    SA2_all_na = ablation_results$sa2_all_na,

    # Source files
    original_summary_path = original_path,
    ablation_result_path = ablation_path
  )
}

# ============================================================ #
# ---- Locate files and generate summary --------------------- #
# ============================================================ #

summary_rows <- list()

# ------------------------------------------------------------ #
# ---- MAGMA ------------------------------------------------- #
# ------------------------------------------------------------ #

magma_ablation_root <- file.path(
  ablation_dir,
  "magma"
)

if (dir.exists(magma_ablation_root)) {
  magma_targets <- list.dirs(
    path = magma_ablation_root,
    full.names = FALSE,
    recursive = FALSE
  )

  magma_targets <- sort(
    magma_targets[magma_targets != ""]
  )

  cat(
    sprintf(
      "MAGMA: detected %d ablation target(s)\n",
      length(magma_targets)
    )
  )

  for (target in magma_targets) {
    ablation_path <- file.path(
      magma_ablation_root,
      target,
      paste0(target, ".ablation.fdrreg.txt")
    )

    original_path <- file.path(
      base_dir,
      target,
      "05.magma_fdrreg",
      paste0(target, ".summary.csv")
    )

    ablation_results <- extract_ablation_results(
      path = ablation_path,
      threshold = threshold
    )

    original_results <- extract_original_long_summary(
      path = original_path,
      threshold = threshold
    )

    summary_rows[[length(summary_rows) + 1L]] <-
      create_summary_row(
        pipeline = "magma",
        target = target,
        region = NA_character_,
        original_results = original_results,
        ablation_results = ablation_results,
        original_path = original_path,
        ablation_path = ablation_path,
        threshold = threshold
      )
  }
} else {
  message(
    "[WARNING] MAGMA ablation directory not found: ",
    magma_ablation_root
  )
}

# ------------------------------------------------------------ #
# ---- MetaXcan / S-PrediXcan -------------------------------- #
# ------------------------------------------------------------ #

metaxcan_ablation_root <- file.path(
  ablation_dir,
  "metaxcan_v7"
)

if (dir.exists(metaxcan_ablation_root)) {
  metaxcan_targets <- list.dirs(
    path = metaxcan_ablation_root,
    full.names = FALSE,
    recursive = FALSE
  )

  metaxcan_targets <- sort(
    metaxcan_targets[metaxcan_targets != ""]
  )

  cat(
    sprintf(
      "MetaXcan: detected %d ablation target(s)\n",
      length(metaxcan_targets)
    )
  )

  for (target in metaxcan_targets) {
    target_ablation_root <- file.path(
      metaxcan_ablation_root,
      target
    )

    brain_regions <- list.dirs(
      path = target_ablation_root,
      full.names = FALSE,
      recursive = FALSE
    )

    brain_regions <- sort(
      brain_regions[brain_regions != ""]
    )

    for (region in brain_regions) {
      ablation_path <- file.path(
        target_ablation_root,
        region,
        paste0(region, ".ablation.fdrreg.txt")
      )

      original_path <- file.path(
        base_dir,
        target,
        "07.metaxcan_fdrreg",
        "01.fdrreg_results",
        region,
        paste0(region, ".summary.csv")
      )

      ablation_results <- extract_ablation_results(
        path = ablation_path,
        threshold = threshold
      )

      original_results <- extract_original_metaxcan_summary(
        path = original_path,
        threshold = threshold
      )

      summary_rows[[length(summary_rows) + 1L]] <-
        create_summary_row(
          pipeline = "metaxcan",
          target = target,
          region = region,
          original_results = original_results,
          ablation_results = ablation_results,
          original_path = original_path,
          ablation_path = ablation_path,
          threshold = threshold
        )
    }
  }
} else {
  message(
    "[WARNING] MetaXcan ablation directory not found: ",
    metaxcan_ablation_root
  )
}

# ------------------------------------------------------------ #
# ---- S-MultiXcan ------------------------------------------- #
# ------------------------------------------------------------ #

smultixcan_ablation_root <- file.path(
  ablation_dir,
  "smultixcan_v7"
)

if (dir.exists(smultixcan_ablation_root)) {
  smultixcan_targets <- list.dirs(
    path = smultixcan_ablation_root,
    full.names = FALSE,
    recursive = FALSE
  )

  smultixcan_targets <- sort(
    smultixcan_targets[smultixcan_targets != ""]
  )

  cat(
    sprintf(
      "S-MultiXcan: detected %d ablation target(s)\n",
      length(smultixcan_targets)
    )
  )

  for (target in smultixcan_targets) {
    ablation_path <- file.path(
      smultixcan_ablation_root,
      target,
      paste0(target, ".ablation.fdrreg.txt")
    )

    original_path <- file.path(
      base_dir,
      target,
      "09.smultixcan_fdrreg",
      paste0(target, ".summary.csv")
    )

    ablation_results <- extract_ablation_results(
      path = ablation_path,
      threshold = threshold
    )

    original_results <- extract_original_long_summary(
      path = original_path,
      threshold = threshold
    )

    summary_rows[[length(summary_rows) + 1L]] <-
      create_summary_row(
        pipeline = "smultixcan",
        target = target,
        region = NA_character_,
        original_results = original_results,
        ablation_results = ablation_results,
        original_path = original_path,
        ablation_path = ablation_path,
        threshold = threshold
      )
  }
} else {
  message(
    "[WARNING] S-MultiXcan ablation directory not found: ",
    smultixcan_ablation_root
  )
}

# ============================================================ #
# ---- Combine and save results ------------------------------ #
# ============================================================ #

if (length(summary_rows) == 0L) {
  stop(
    "No valid ablation result directories were found under: ",
    ablation_dir
  )
}

summary_table <- rbindlist(
  summary_rows,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  summary_table,
  pipeline,
  target,
  region
)

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  summary_table,
  file = output_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

# ============================================================ #
# ---- Save compact requested table -------------------------- #
# ============================================================ #

if (grepl("\\.csv$", output_file, ignore.case = TRUE)) {
  compact_output_file <- sub(
    "\\.csv$",
    ".compact.csv",
    output_file,
    ignore.case = TRUE
  )
} else {
  compact_output_file <- paste0(
    output_file,
    ".compact.csv"
  )
}

compact_table <- summary_table[, .(
  pipeline,
  target,
  region,
  threshold,
  original_total_genes,
  ablation_total_genes,
  original_qval,
  original_fdr_the,
  original_bio_fdr_the,
  SA1_all_bio_only_fdr_the,
  SA2_disorders_partial_bio_fdr_the
)]

fwrite(
  compact_table,
  file = compact_output_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

# ============================================================ #
# ---- Save issue report ------------------------------------- #
# ============================================================ #

if (grepl("\\.csv$", output_file, ignore.case = TRUE)) {
  issue_output_file <- sub(
    "\\.csv$",
    ".issues.csv",
    output_file,
    ignore.case = TRUE
  )
} else {
  issue_output_file <- paste0(
    output_file,
    ".issues.csv"
  )
}

issue_table <- summary_table[
  original_file_exists == FALSE |
    original_file_readable == FALSE |
    original_threshold_found == FALSE |
    ablation_file_exists == FALSE |
    ablation_file_readable == FALSE |
    SA1_column_exists == FALSE |
    SA2_column_exists == FALSE |
    SA1_all_na == TRUE |
    SA2_all_na == TRUE |
    original_qval_column_exists == FALSE |
    original_fdr_the_column_exists == FALSE |
    original_bio_fdr_the_column_exists == FALSE
]

fwrite(
  issue_table,
  file = issue_output_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

# ============================================================ #
# ---- Completion report ------------------------------------- #
# ============================================================ #

cat("\n")
cat("============================================================\n")
cat("Ablation summary completed\n")
cat("============================================================\n")
cat("Significance criterion: value < ", threshold, "\n", sep = "")
cat("Total summary rows: ", nrow(summary_table), "\n", sep = "")
cat("\n")

cat("Rows by pipeline:\n")
print(
  summary_table[, .N, by = pipeline][order(pipeline)]
)

cat("\n")
cat("Full summary:\n")
cat("  ", output_file, "\n", sep = "")

cat("Compact summary:\n")
cat("  ", compact_output_file, "\n", sep = "")

cat("Issue report:\n")
cat("  ", issue_output_file, "\n", sep = "")

cat("\n")
cat("Number of rows with detected issues: ")
cat(nrow(issue_table), "\n")

if (nrow(issue_table) > 0L) {
  cat("\nDetected issues by pipeline:\n")
  print(
    issue_table[, .N, by = pipeline][order(pipeline)]
  )

  cat("\nReview the issue report for details.\n")
} else {
  cat("\nNo missing files, missing columns, or all-NA ablation models were detected.\n")
}

cat("============================================================\n")

quit(
  save = "no",
  status = 0,
  runLast = FALSE
)
