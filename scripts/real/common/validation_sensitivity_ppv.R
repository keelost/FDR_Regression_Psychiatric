#!/usr/bin/env Rscript

#============================================================#
# Temporal Validation of FDRreg Results
#============================================================#
# This script calculates:
#
# 1. Sensitivity
#    Sensitivity = TP / (TP + FN)
#
# 2. Positive Predictive Value
#    PPV = TP / (TP + FP)
#
# 3. McNemar's paired test of detection sensitivity
#
# Gold standard:
# Identifiers significant in the standard, unassisted BH-FDR
# analysis of the later, larger GWAS.
#
# McNemar table, restricted to gold-standard-positive IDs:
#
#                         Earlier standard detected
#                         Yes                 No
# Earlier FDRreg Yes      a                   b
# Earlier FDRreg No       c                   d
#
# a: Detected by both FDRreg and standard analysis
# b: Detected uniquely by FDRreg
# c: Detected uniquely by standard analysis
# d: Detected by neither method
#
# McNemar's test evaluates asymmetry between b and c.
# A significant result with b > c supports greater recovery
# of gold-standard signals by FDRreg.
#
# Usage:
#   Rscript validation_sensitivity_ppv.R --threshold 0.01
#   Rscript validation_sensitivity_ppv.R --threshold 0.05
#============================================================#

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

#============================================================#
# Argument parsing
#============================================================#

option_list <- list(
  make_option(
    c("-t", "--threshold"),
    type = "double",
    default = 0.01,
    dest = "threshold",
    help = paste0(
      "FDR significance threshold used for the later gold standard ",
      "and earlier-study discoveries [default: %default]"
    ),
    metavar = "THRESHOLD"
  ),
  make_option(
    c("-b", "--base-path"),
    type = "character",
    default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
    dest = "base_path",
    help = "Base directory containing all study folders [default: %default]",
    metavar = "DIR"
  ),
  make_option(
    c("-o", "--base-output"),
    type = "character",
    default = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "02.sen_ppv"),
    dest = "base_output",
    help = "Base output directory [default: %default]",
    metavar = "DIR"
  )
)

option_parser <- OptionParser(
  option_list = option_list,
  description = paste(
    "Temporal validation using sensitivity, PPV,",
    "and gold-standard-positive McNemar tests."
  )
)

opt <- parse_args(option_parser)

SIGNIFICANCE_THRESHOLD <- opt$threshold
base_path <- opt$base_path
base_output_dir <- opt$base_output

if (!is.finite(SIGNIFICANCE_THRESHOLD) ||
    SIGNIFICANCE_THRESHOLD <= 0 ||
    SIGNIFICANCE_THRESHOLD >= 1) {
  stop("--threshold must be a finite number between 0 and 1.")
}

if (!dir.exists(base_path)) {
  stop("Base path does not exist: ", base_path)
}

threshold_label <- format(
  SIGNIFICANCE_THRESHOLD,
  scientific = FALSE,
  trim = TRUE
)

output_dir <- file.path(
  base_output_dir,
  paste0("threshold_", threshold_label)
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop("Failed to create output directory: ", output_dir)
}

cat("============================================================\n")
cat("Temporal Validation Configuration\n")
cat("============================================================\n")
cat("Significance threshold :", threshold_label, "\n")
cat("Base input directory  :", base_path, "\n")
cat("Output directory      :", output_dir, "\n")
cat("============================================================\n\n")

#============================================================#
# Study comparisons
#============================================================#

trait_comparisons <- list(
  list(
    smaller = "scz2012",
    larger = "wal.scz2018",
    label = "SCZ"
  ),
  list(
    smaller = "scz2012",
    larger = "scz2014",
    label = "SCZ2014"
  ),
  list(
    smaller = "bd2012",
    larger = "bd2018",
    label = "BD"
  ),
  list(
    smaller = "mdd2013",
    larger = "mdd2019",
    label = "MDD"
  ),
  list(
    smaller = "scz2014",
    larger = "wal.scz2018",
    label = "SCZ_BIG"
  )
)

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

#============================================================#
# Analysis specifications
#============================================================#

analysis_specifications <- list(
  SNP = list(
    identifier_column = "snpid",
    standard_fdr_column = "qval_bh",
    method_columns = c(
      "qval_bh",
      "fdr_theoretical",
      "fdr_empirical"
    ),
    expected_columns = c(
      "snpid",
      "z_score_signed",
      "pval_raw",
      "qval_bh",
      "fdr_theoretical",
      "fdr_empirical"
    ),
    character_columns = c("snpid")
  ),
  MAGMA = list(
    identifier_column = "ID",
    standard_fdr_column = "qval",
    method_columns = c(
      "qval",
      "FDR.the",
      "FDR.emp",
      "bio.FDR.the",
      "bio.FDR.emp",
      "bio.FDR.the.lasso",
      "bio.FDR.emp.lasso"
    ),
    expected_columns = c(
      "ID",
      "CHR",
      "START",
      "STOP",
      "NSNPS",
      "NPARAM",
      "N",
      "ZSTAT",
      "P",
      "FDR.the",
      "FDR.emp",
      "qval",
      "bio.FDR.the",
      "bio.FDR.emp",
      "bio.FDR.the.lasso",
      "bio.FDR.emp.lasso"
    ),
    character_columns = c("ID")
  ),
  MetaXcan = list(
    identifier_column = "ENSEMBL_GENE_ID",
    standard_fdr_column = "qval",
    method_columns = c(
      "qval",
      "FDR_theoretical",
      "FDR_empirical",
      "bio_FDR_theoretical",
      "bio_FDR_empirical",
      "lasso_FDR_theoretical",
      "lasso_FDR_empirical"
    ),
    expected_columns = c(
      "ENSEMBL_GENE_ID",
      "gene_name",
      "zscore",
      "pvalue",
      "qval",
      "FDR_theoretical",
      "FDR_empirical",
      "bio_FDR_theoretical",
      "bio_FDR_empirical",
      "lasso_FDR_theoretical",
      "lasso_FDR_empirical"
    ),
    character_columns = c(
      "ENSEMBL_GENE_ID",
      "gene_name"
    )
  ),
  SmultiXcan = list(
    identifier_column = "ENSEMBL_GENE_ID",
    standard_fdr_column = "qval",
    method_columns = c(
      "qval",
      "FDR.the",
      "FDR.emp",
      "bio.FDR.the",
      "bio.FDR.emp",
      "bio.FDR.the.lasso",
      "bio.FDR.emp.lasso"
    ),
    expected_columns = c(
      "ENSEMBL_GENE_ID",
      "ENSEMBL_GENE_ID_name",
      "pvalue",
      "zscore",
      "qval",
      "FDR.the",
      "FDR.emp",
      "bio.FDR.the",
      "bio.FDR.emp",
      "bio.FDR.the.lasso",
      "bio.FDR.emp.lasso"
    ),
    character_columns = c(
      "ENSEMBL_GENE_ID",
      "ENSEMBL_GENE_ID_name"
    )
  )
)

#============================================================#
# Helper functions
#============================================================#

safe_ratio <- function(numerator, denominator) {
  if (is.na(numerator) || is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

is_model_available <- function(data, column_name) {
  if (!column_name %in% names(data)) {
    return(FALSE)
  }

  values <- suppressWarnings(as.numeric(data[[column_name]]))

  any(is.finite(values))
}

read_fdrreg_result <- function(file_path, specification) {
  if (!file.exists(file_path)) {
    warning("Result file not found: ", file_path)
    return(NULL)
  }

  result_data <- tryCatch(
    fread(
      file_path,
      fill = TRUE,
      na.strings = c("", "NA", "NaN", "NULL"),
      strip.white = TRUE,
      showProgress = FALSE
    ),
    error = function(e) {
      warning(
        "Failed to read file: ",
        file_path,
        "; reason: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(result_data)) {
    return(NULL)
  }

  identifier_column <- specification$identifier_column

  if (!identifier_column %in% names(result_data)) {
    warning(
      "Identifier column '",
      identifier_column,
      "' is missing from: ",
      file_path
    )
    return(NULL)
  }

  missing_columns <- setdiff(
    specification$expected_columns,
    names(result_data)
  )

  for (column_name in missing_columns) {
    if (column_name %in% specification$character_columns) {
      result_data[, (column_name) := NA_character_]
    } else {
      result_data[, (column_name) := NA_real_]
    }
  }

  numeric_columns <- setdiff(
    specification$expected_columns,
    specification$character_columns
  )

  numeric_columns <- intersect(
    numeric_columns,
    names(result_data)
  )

  for (column_name in numeric_columns) {
    result_data[
      ,
      (column_name) := suppressWarnings(
        as.numeric(get(column_name))
      )
    ]
  }

  result_data[
    ,
    (identifier_column) := trimws(
      as.character(get(identifier_column))
    )
  ]

  result_data <- result_data[
    !is.na(get(identifier_column)) &
      get(identifier_column) != ""
  ]

  duplicated_identifier_count <- sum(
    duplicated(result_data[[identifier_column]])
  )

  if (duplicated_identifier_count > 0) {
    warning(
      sprintf(
        paste0(
          "%d duplicated identifiers found in %s. ",
          "Only the first occurrence will be retained."
        ),
        duplicated_identifier_count,
        file_path
      )
    )

    result_data <- unique(
      result_data,
      by = identifier_column
    )
  }

  result_data
}

standardize_identifier <- function(data, identifier_column) {
  standardized_data <- copy(data)

  standardized_data[
    ,
    identifier := trimws(
      as.character(get(identifier_column))
    )
  ]

  standardized_data <- standardized_data[
    !is.na(identifier) & identifier != ""
  ]

  unique(
    standardized_data,
    by = "identifier"
  )
}

extract_discovery_ids <- function(
    data,
    method_column,
    comparable_ids,
    threshold) {

  if (!is_model_available(data, method_column)) {
    return(NULL)
  }

  discovery_ids <- data[
    identifier %in% comparable_ids &
      is.finite(get(method_column)) &
      get(method_column) < threshold,
    unique(identifier)
  ]

  discovery_ids
}

calculate_classification_metrics <- function(
    discovery_ids,
    gold_standard_ids,
    comparable_ids) {

  if (is.null(discovery_ids)) {
    return(list(
      true_positives = NA_integer_,
      false_positives = NA_integer_,
      false_negatives = NA_integer_,
      discovery_count = NA_integer_,
      gold_standard_count = length(gold_standard_ids),
      sensitivity = NA_real_,
      ppv = NA_real_
    ))
  }

  discovery_ids <- intersect(
    unique(discovery_ids),
    comparable_ids
  )

  gold_standard_ids <- intersect(
    unique(gold_standard_ids),
    comparable_ids
  )

  true_positives <- length(
    intersect(discovery_ids, gold_standard_ids)
  )

  false_positives <- length(
    setdiff(discovery_ids, gold_standard_ids)
  )

  false_negatives <- length(
    setdiff(gold_standard_ids, discovery_ids)
  )

  discovery_count <- length(discovery_ids)
  gold_standard_count <- length(gold_standard_ids)

  sensitivity <- safe_ratio(
    true_positives,
    true_positives + false_negatives
  )

  ppv <- safe_ratio(
    true_positives,
    true_positives + false_positives
  )

  list(
    true_positives = true_positives,
    false_positives = false_positives,
    false_negatives = false_negatives,
    discovery_count = discovery_count,
    gold_standard_count = gold_standard_count,
    sensitivity = sensitivity,
    ppv = ppv
  )
}

perform_gold_standard_mcnemar_test <- function(
    gold_standard_ids,
    standard_discovery_ids,
    fdrreg_discovery_ids) {

  if (is.null(standard_discovery_ids) ||
      is.null(fdrreg_discovery_ids)) {
    return(list(
      a_both_detected = NA_integer_,
      b_fdrreg_only = NA_integer_,
      c_standard_only = NA_integer_,
      d_neither_detected = NA_integer_,
      discordant_pairs = NA_integer_,
      mcnemar_chi_square = NA_real_,
      mcnemar_two_sided_p = NA_real_,
      exact_two_sided_p = NA_real_,
      exact_one_sided_p_fdrreg_greater = NA_real_,
      direction = NA_character_
    ))
  }

  gold_standard_ids <- unique(gold_standard_ids)

  if (length(gold_standard_ids) == 0) {
    return(list(
      a_both_detected = 0L,
      b_fdrreg_only = 0L,
      c_standard_only = 0L,
      d_neither_detected = 0L,
      discordant_pairs = 0L,
      mcnemar_chi_square = NA_real_,
      mcnemar_two_sided_p = NA_real_,
      exact_two_sided_p = NA_real_,
      exact_one_sided_p_fdrreg_greater = NA_real_,
      direction = "No gold-standard-positive identifiers"
    ))
  }

  standard_detected <- (
    gold_standard_ids %in% standard_discovery_ids
  )

  fdrreg_detected <- (
    gold_standard_ids %in% fdrreg_discovery_ids
  )

  a_both_detected <- sum(
    fdrreg_detected & standard_detected
  )

  b_fdrreg_only <- sum(
    fdrreg_detected & !standard_detected
  )

  c_standard_only <- sum(
    !fdrreg_detected & standard_detected
  )

  d_neither_detected <- sum(
    !fdrreg_detected & !standard_detected
  )

  mcnemar_table <- matrix(
    c(
      a_both_detected,
      b_fdrreg_only,
      c_standard_only,
      d_neither_detected
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      earlier_fdrreg = c("Detected", "Not_detected"),
      earlier_standard = c("Detected", "Not_detected")
    )
  )

  discordant_pairs <- b_fdrreg_only + c_standard_only

  if (discordant_pairs == 0) {
    return(list(
      a_both_detected = a_both_detected,
      b_fdrreg_only = b_fdrreg_only,
      c_standard_only = c_standard_only,
      d_neither_detected = d_neither_detected,
      discordant_pairs = 0L,
      mcnemar_chi_square = 0,
      mcnemar_two_sided_p = 1,
      exact_two_sided_p = 1,
      exact_one_sided_p_fdrreg_greater = 1,
      direction = "No discordant gold-standard identifiers"
    ))
  }

  mcnemar_result <- suppressWarnings(
    mcnemar.test(
      mcnemar_table,
      correct = TRUE
    )
  )

  exact_two_sided_result <- binom.test(
    x = b_fdrreg_only,
    n = discordant_pairs,
    p = 0.5,
    alternative = "two.sided"
  )

  exact_one_sided_result <- binom.test(
    x = b_fdrreg_only,
    n = discordant_pairs,
    p = 0.5,
    alternative = "greater"
  )

  direction <- if (b_fdrreg_only > c_standard_only) {
    "FDRreg recovered more gold-standard identifiers"
  } else if (b_fdrreg_only < c_standard_only) {
    "Standard analysis recovered more gold-standard identifiers"
  } else {
    "Equal discordant recovery"
  }

  list(
    a_both_detected = a_both_detected,
    b_fdrreg_only = b_fdrreg_only,
    c_standard_only = c_standard_only,
    d_neither_detected = d_neither_detected,
    discordant_pairs = discordant_pairs,
    mcnemar_chi_square = unname(
      mcnemar_result$statistic
    ),
    mcnemar_two_sided_p = mcnemar_result$p.value,
    exact_two_sided_p = exact_two_sided_result$p.value,
    exact_one_sided_p_fdrreg_greater =
      exact_one_sided_result$p.value,
    direction = direction
  )
}

empty_metric_row <- function(
    analysis_type,
    comparison_label,
    smaller_study,
    larger_study,
    brain_region,
    method_column,
    standard_fdr_column,
    comparable_count,
    gold_standard_count,
    is_standard_method) {

  data.table(
    analysis_type = analysis_type,
    comparison = comparison_label,
    smaller_study = smaller_study,
    larger_study = larger_study,
    brain_region = brain_region,
    method = method_column,
    method_role = if (is_standard_method) {
      "Standard_BH_FDR"
    } else {
      "FDRreg"
    },
    standard_fdr_column = standard_fdr_column,
    threshold = SIGNIFICANCE_THRESHOLD,
    model_available = FALSE,
    comparable_identifier_count = comparable_count,
    gold_standard_count = gold_standard_count,
    discovery_count = NA_integer_,
    true_positives = NA_integer_,
    false_positives = NA_integer_,
    false_negatives = NA_integer_,
    sensitivity = NA_real_,
    ppv = NA_real_
  )
}

process_temporal_comparison <- function(
    earlier_data,
    later_data,
    specification,
    analysis_type,
    comparison_label,
    smaller_study,
    larger_study,
    brain_region = NA_character_) {

  identifier_column <- specification$identifier_column
  standard_fdr_column <- specification$standard_fdr_column
  method_columns <- specification$method_columns

  earlier_data <- standardize_identifier(
    earlier_data,
    identifier_column
  )

  later_data <- standardize_identifier(
    later_data,
    identifier_column
  )

  comparable_ids <- intersect(
    earlier_data$identifier,
    later_data$identifier
  )

  if (length(comparable_ids) == 0) {
    warning(
      sprintf(
        paste0(
          "No common identifiers for %s, %s, ",
          "%s versus %s, region=%s"
        ),
        analysis_type,
        comparison_label,
        smaller_study,
        larger_study,
        brain_region
      )
    )

    return(list(
      metrics = data.table(),
      mcnemar = data.table()
    ))
  }

  earlier_data <- earlier_data[
    identifier %in% comparable_ids
  ]

  later_data <- later_data[
    identifier %in% comparable_ids
  ]

  if (!is_model_available(
    later_data,
    standard_fdr_column
  )) {
    warning(
      sprintf(
        paste0(
          "Later-study standard BH-FDR column is unavailable: ",
          "%s, %s, %s"
        ),
        analysis_type,
        larger_study,
        standard_fdr_column
      )
    )

    return(list(
      metrics = data.table(),
      mcnemar = data.table()
    ))
  }

  gold_standard_ids <- later_data[
    is.finite(get(standard_fdr_column)) &
      get(standard_fdr_column) < SIGNIFICANCE_THRESHOLD,
    unique(identifier)
  ]

  standard_discovery_ids <- extract_discovery_ids(
    data = earlier_data,
    method_column = standard_fdr_column,
    comparable_ids = comparable_ids,
    threshold = SIGNIFICANCE_THRESHOLD
  )

  metric_results <- list()
  mcnemar_results <- list()

  metric_index <- 1L
  mcnemar_index <- 1L

  for (method_column in method_columns) {
    is_standard_method <- (
      method_column == standard_fdr_column
    )

    method_available <- is_model_available(
      earlier_data,
      method_column
    )

    if (!method_available) {
      metric_results[[metric_index]] <- empty_metric_row(
        analysis_type = analysis_type,
        comparison_label = comparison_label,
        smaller_study = smaller_study,
        larger_study = larger_study,
        brain_region = brain_region,
        method_column = method_column,
        standard_fdr_column = standard_fdr_column,
        comparable_count = length(comparable_ids),
        gold_standard_count = length(gold_standard_ids),
        is_standard_method = is_standard_method
      )

      metric_index <- metric_index + 1L

      if (!is_standard_method) {
        mcnemar_results[[mcnemar_index]] <- data.table(
          analysis_type = analysis_type,
          comparison = comparison_label,
          smaller_study = smaller_study,
          larger_study = larger_study,
          brain_region = brain_region,
          standard_method = standard_fdr_column,
          fdrreg_method = method_column,
          threshold = SIGNIFICANCE_THRESHOLD,
          model_available = FALSE,
          comparable_identifier_count = length(comparable_ids),
          gold_standard_count = length(gold_standard_ids),
          a_both_detected = NA_integer_,
          b_fdrreg_only = NA_integer_,
          c_standard_only = NA_integer_,
          d_neither_detected = NA_integer_,
          discordant_pairs = NA_integer_,
          mcnemar_chi_square = NA_real_,
          mcnemar_two_sided_p = NA_real_,
          exact_two_sided_p = NA_real_,
          exact_one_sided_p_fdrreg_greater = NA_real_,
          direction = "FDRreg model unavailable"
        )

        mcnemar_index <- mcnemar_index + 1L
      }

      next
    }

    discovery_ids <- extract_discovery_ids(
      data = earlier_data,
      method_column = method_column,
      comparable_ids = comparable_ids,
      threshold = SIGNIFICANCE_THRESHOLD
    )

    classification_metrics <- calculate_classification_metrics(
      discovery_ids = discovery_ids,
      gold_standard_ids = gold_standard_ids,
      comparable_ids = comparable_ids
    )

    metric_results[[metric_index]] <- data.table(
      analysis_type = analysis_type,
      comparison = comparison_label,
      smaller_study = smaller_study,
      larger_study = larger_study,
      brain_region = brain_region,
      method = method_column,
      method_role = if (is_standard_method) {
        "Standard_BH_FDR"
      } else {
        "FDRreg"
      },
      standard_fdr_column = standard_fdr_column,
      threshold = SIGNIFICANCE_THRESHOLD,
      model_available = TRUE,
      comparable_identifier_count = length(comparable_ids),
      gold_standard_count =
        classification_metrics$gold_standard_count,
      discovery_count =
        classification_metrics$discovery_count,
      true_positives =
        classification_metrics$true_positives,
      false_positives =
        classification_metrics$false_positives,
      false_negatives =
        classification_metrics$false_negatives,
      sensitivity =
        classification_metrics$sensitivity,
      ppv =
        classification_metrics$ppv
    )

    metric_index <- metric_index + 1L

    if (!is_standard_method) {
      mcnemar_result <- perform_gold_standard_mcnemar_test(
        gold_standard_ids = gold_standard_ids,
        standard_discovery_ids = standard_discovery_ids,
        fdrreg_discovery_ids = discovery_ids
      )

      mcnemar_results[[mcnemar_index]] <- data.table(
        analysis_type = analysis_type,
        comparison = comparison_label,
        smaller_study = smaller_study,
        larger_study = larger_study,
        brain_region = brain_region,
        standard_method = standard_fdr_column,
        fdrreg_method = method_column,
        threshold = SIGNIFICANCE_THRESHOLD,
        model_available = TRUE,
        comparable_identifier_count = length(comparable_ids),
        gold_standard_count = length(gold_standard_ids),
        a_both_detected =
          mcnemar_result$a_both_detected,
        b_fdrreg_only =
          mcnemar_result$b_fdrreg_only,
        c_standard_only =
          mcnemar_result$c_standard_only,
        d_neither_detected =
          mcnemar_result$d_neither_detected,
        discordant_pairs =
          mcnemar_result$discordant_pairs,
        mcnemar_chi_square =
          mcnemar_result$mcnemar_chi_square,
        mcnemar_two_sided_p =
          mcnemar_result$mcnemar_two_sided_p,
        exact_two_sided_p =
          mcnemar_result$exact_two_sided_p,
        exact_one_sided_p_fdrreg_greater =
          mcnemar_result$exact_one_sided_p_fdrreg_greater,
        direction =
          mcnemar_result$direction
      )

      mcnemar_index <- mcnemar_index + 1L
    }
  }

  list(
    metrics = rbindlist(
      metric_results,
      fill = TRUE
    ),
    mcnemar = if (length(mcnemar_results) > 0) {
      rbindlist(
        mcnemar_results,
        fill = TRUE
      )
    } else {
      data.table()
    }
  )
}

#============================================================#
# Path construction functions
#============================================================#

get_snp_file <- function(study_name) {
  file.path(
    base_path,
    study_name,
    "02.fdrreg_results",
    "fdr_values_per_snp_nolasso.csv"
  )
}

get_magma_file <- function(study_name) {
  file.path(
    base_path,
    study_name,
    "05.magma_fdrreg",
    paste0(study_name, ".gene.bio.fdrreg.txt")
  )
}

get_metaxcan_file <- function(study_name, brain_region) {
  file.path(
    base_path,
    study_name,
    "07.metaxcan_fdrreg",
    "01.fdrreg_results",
    brain_region,
    paste0(brain_region, ".gene.bio.fdrreg.txt")
  )
}

get_smultixcan_file <- function(study_name) {
  file.path(
    base_path,
    study_name,
    "09.smultixcan_fdrreg",
    paste0(study_name, ".gene.bio.fdrreg.txt")
  )
}

#============================================================#
# Result containers
#============================================================#

all_metric_results <- list()
all_mcnemar_results <- list()

metric_result_index <- 1L
mcnemar_result_index <- 1L

append_comparison_results <- function(comparison_results) {
  if (nrow(comparison_results$metrics) > 0) {
    all_metric_results[[metric_result_index]] <<-
      comparison_results$metrics

    metric_result_index <<- metric_result_index + 1L
  }

  if (nrow(comparison_results$mcnemar) > 0) {
    all_mcnemar_results[[mcnemar_result_index]] <<-
      comparison_results$mcnemar

    mcnemar_result_index <<- mcnemar_result_index + 1L
  }

  invisible(NULL)
}

#============================================================#
# SNP analysis
#============================================================#

cat("========== SNP Analysis ==========\n")

snp_specification <- analysis_specifications$SNP

for (comparison in trait_comparisons) {
  cat(
    sprintf(
      "Processing SNP: %s vs %s\n",
      comparison$smaller,
      comparison$larger
    )
  )

  earlier_file <- get_snp_file(comparison$smaller)
  later_file <- get_snp_file(comparison$larger)

  earlier_data <- read_fdrreg_result(
    earlier_file,
    snp_specification
  )

  later_data <- read_fdrreg_result(
    later_file,
    snp_specification
  )

  if (is.null(earlier_data) || is.null(later_data)) {
    warning(
      "Skipping SNP comparison: ",
      comparison$smaller,
      " vs ",
      comparison$larger
    )
    next
  }

  comparison_results <- process_temporal_comparison(
    earlier_data = earlier_data,
    later_data = later_data,
    specification = snp_specification,
    analysis_type = "SNP",
    comparison_label = comparison$label,
    smaller_study = comparison$smaller,
    larger_study = comparison$larger
  )

  append_comparison_results(comparison_results)

  rm(earlier_data, later_data, comparison_results)
  gc(verbose = FALSE)
}

#============================================================#
# MAGMA analysis
#============================================================#

cat("\n========== MAGMA Analysis ==========\n")

magma_specification <- analysis_specifications$MAGMA

for (comparison in trait_comparisons) {
  cat(
    sprintf(
      "Processing MAGMA: %s vs %s\n",
      comparison$smaller,
      comparison$larger
    )
  )

  earlier_file <- get_magma_file(comparison$smaller)
  later_file <- get_magma_file(comparison$larger)

  earlier_data <- read_fdrreg_result(
    earlier_file,
    magma_specification
  )

  later_data <- read_fdrreg_result(
    later_file,
    magma_specification
  )

  if (is.null(earlier_data) || is.null(later_data)) {
    warning(
      "Skipping MAGMA comparison: ",
      comparison$smaller,
      " vs ",
      comparison$larger
    )
    next
  }

  comparison_results <- process_temporal_comparison(
    earlier_data = earlier_data,
    later_data = later_data,
    specification = magma_specification,
    analysis_type = "MAGMA",
    comparison_label = comparison$label,
    smaller_study = comparison$smaller,
    larger_study = comparison$larger
  )

  append_comparison_results(comparison_results)

  rm(earlier_data, later_data, comparison_results)
  gc(verbose = FALSE)
}

#============================================================#
# MetaXcan analysis
#============================================================#

cat("\n========== MetaXcan Analysis ==========\n")

metaxcan_specification <- analysis_specifications$MetaXcan

for (comparison in trait_comparisons) {
  cat(
    sprintf(
      "Processing MetaXcan: %s vs %s\n",
      comparison$smaller,
      comparison$larger
    )
  )

  for (brain_region in brain_regions) {
    cat("  Region: ", brain_region, "\n", sep = "")

    earlier_file <- get_metaxcan_file(
      comparison$smaller,
      brain_region
    )

    later_file <- get_metaxcan_file(
      comparison$larger,
      brain_region
    )

    earlier_data <- read_fdrreg_result(
      earlier_file,
      metaxcan_specification
    )

    later_data <- read_fdrreg_result(
      later_file,
      metaxcan_specification
    )

    if (is.null(earlier_data) || is.null(later_data)) {
      warning(
        "Skipping MetaXcan comparison: ",
        comparison$smaller,
        " vs ",
        comparison$larger,
        ", region=",
        brain_region
      )
      next
    }

    comparison_results <- process_temporal_comparison(
      earlier_data = earlier_data,
      later_data = later_data,
      specification = metaxcan_specification,
      analysis_type = "MetaXcan",
      comparison_label = comparison$label,
      smaller_study = comparison$smaller,
      larger_study = comparison$larger,
      brain_region = brain_region
    )

    append_comparison_results(comparison_results)

    rm(earlier_data, later_data, comparison_results)
    gc(verbose = FALSE)
  }
}

#============================================================#
# S-MultiXcan analysis
#============================================================#

cat("\n========== S-MultiXcan Analysis ==========\n")

smultixcan_specification <- analysis_specifications$SmultiXcan

for (comparison in trait_comparisons) {
  cat(
    sprintf(
      "Processing S-MultiXcan: %s vs %s\n",
      comparison$smaller,
      comparison$larger
    )
  )

  earlier_file <- get_smultixcan_file(
    comparison$smaller
  )

  later_file <- get_smultixcan_file(
    comparison$larger
  )

  earlier_data <- read_fdrreg_result(
    earlier_file,
    smultixcan_specification
  )

  later_data <- read_fdrreg_result(
    later_file,
    smultixcan_specification
  )

  if (is.null(earlier_data) || is.null(later_data)) {
    warning(
      "Skipping S-MultiXcan comparison: ",
      comparison$smaller,
      " vs ",
      comparison$larger
    )
    next
  }

  comparison_results <- process_temporal_comparison(
    earlier_data = earlier_data,
    later_data = later_data,
    specification = smultixcan_specification,
    analysis_type = "SmultiXcan",
    comparison_label = comparison$label,
    smaller_study = comparison$smaller,
    larger_study = comparison$larger
  )

  append_comparison_results(comparison_results)

  rm(earlier_data, later_data, comparison_results)
  gc(verbose = FALSE)
}

#============================================================#
# Combine results
#============================================================#

metric_results <- if (length(all_metric_results) > 0) {
  rbindlist(
    all_metric_results,
    fill = TRUE
  )
} else {
  data.table()
}

mcnemar_results <- if (length(all_mcnemar_results) > 0) {
  rbindlist(
    all_mcnemar_results,
    fill = TRUE
  )
} else {
  data.table()
}

#============================================================#
# Multiple-testing adjustment for McNemar results
#============================================================#

if (nrow(mcnemar_results) > 0) {
  mcnemar_results[
    ,
    mcnemar_two_sided_p_bh := p.adjust(
      mcnemar_two_sided_p,
      method = "BH"
    )
  ]

  mcnemar_results[
    ,
    exact_two_sided_p_bh := p.adjust(
      exact_two_sided_p,
      method = "BH"
    )
  ]

  mcnemar_results[
    ,
    exact_one_sided_p_fdrreg_greater_bh := p.adjust(
      exact_one_sided_p_fdrreg_greater,
      method = "BH"
    )
  ]

  mcnemar_results[
    ,
    supports_fdrreg_improvement :=
      model_available &
      !is.na(exact_one_sided_p_fdrreg_greater) &
      b_fdrreg_only > c_standard_only &
      exact_one_sided_p_fdrreg_greater < 0.05
  ]

  mcnemar_results[
    ,
    supports_fdrreg_improvement_bh :=
      model_available &
      !is.na(exact_one_sided_p_fdrreg_greater_bh) &
      b_fdrreg_only > c_standard_only &
      exact_one_sided_p_fdrreg_greater_bh < 0.05
  ]
}

#============================================================#
# Save results
#============================================================#

cat("\n========== Saving Results ==========\n")

metric_output_file <- file.path(
  output_dir,
  "temporal_validation_sensitivity_ppv.csv"
)

mcnemar_output_file <- file.path(
  output_dir,
  "temporal_validation_mcnemar_gold_positive.csv"
)

if (nrow(metric_results) > 0) {
  fwrite(
    metric_results,
    metric_output_file,
    na = "NA"
  )

  cat(
    "Saved sensitivity and PPV results: ",
    metric_output_file,
    "\n",
    sep = ""
  )
}

if (nrow(mcnemar_results) > 0) {
  fwrite(
    mcnemar_results,
    mcnemar_output_file,
    na = "NA"
  )

  cat(
    "Saved McNemar results: ",
    mcnemar_output_file,
    "\n",
    sep = ""
  )
}

# Save analysis-specific files
if (nrow(metric_results) > 0) {
  for (analysis_name in unique(metric_results$analysis_type)) {
    analysis_metric_file <- file.path(
      output_dir,
      paste0(
        tolower(analysis_name),
        "_sensitivity_ppv.csv"
      )
    )

    fwrite(
      metric_results[
        analysis_type == analysis_name
      ],
      analysis_metric_file,
      na = "NA"
    )
  }
}

if (nrow(mcnemar_results) > 0) {
  for (analysis_name in unique(mcnemar_results$analysis_type)) {
    analysis_mcnemar_file <- file.path(
      output_dir,
      paste0(
        tolower(analysis_name),
        "_mcnemar_gold_positive.csv"
      )
    )

    fwrite(
      mcnemar_results[
        analysis_type == analysis_name
      ],
      analysis_mcnemar_file,
      na = "NA"
    )
  }
}

#============================================================#
# Save compact combined report
#============================================================#

if (nrow(metric_results) > 0 &&
    nrow(mcnemar_results) > 0) {

  standard_metrics <- metric_results[
    method_role == "Standard_BH_FDR",
    .(
      analysis_type,
      comparison,
      smaller_study,
      larger_study,
      brain_region,
      threshold,
      standard_discovery_count = discovery_count,
      standard_true_positives = true_positives,
      standard_false_positives = false_positives,
      standard_false_negatives = false_negatives,
      standard_sensitivity = sensitivity,
      standard_ppv = ppv
    )
  ]

  fdrreg_metrics <- metric_results[
    method_role == "FDRreg",
    .(
      analysis_type,
      comparison,
      smaller_study,
      larger_study,
      brain_region,
      fdrreg_method = method,
      threshold,
      fdrreg_model_available = model_available,
      comparable_identifier_count,
      gold_standard_count,
      fdrreg_discovery_count = discovery_count,
      fdrreg_true_positives = true_positives,
      fdrreg_false_positives = false_positives,
      fdrreg_false_negatives = false_negatives,
      fdrreg_sensitivity = sensitivity,
      fdrreg_ppv = ppv
    )
  ]

  combined_report <- merge(
    fdrreg_metrics,
    standard_metrics,
    by = c(
      "analysis_type",
      "comparison",
      "smaller_study",
      "larger_study",
      "brain_region",
      "threshold"
    ),
    all.x = TRUE
  )

  combined_report <- merge(
    combined_report,
    mcnemar_results,
    by = c(
      "analysis_type",
      "comparison",
      "smaller_study",
      "larger_study",
      "brain_region",
      "fdrreg_method",
      "threshold"
    ),
    all.x = TRUE,
    suffixes = c("", "_mcnemar")
  )

  combined_report[
    ,
    sensitivity_difference :=
      fdrreg_sensitivity - standard_sensitivity
  ]

  combined_report[
    ,
    ppv_difference :=
      fdrreg_ppv - standard_ppv
  ]

  combined_report_file <- file.path(
    output_dir,
    "temporal_validation_combined_report.csv"
  )

  fwrite(
    combined_report,
    combined_report_file,
    na = "NA"
  )

  cat(
    "Saved combined report: ",
    combined_report_file,
    "\n",
    sep = ""
  )
}

#============================================================#
# Completion summary
#============================================================#

cat("\n============================================================\n")
cat("Temporal Validation Completed\n")
cat("============================================================\n")
cat("Threshold:", threshold_label, "\n")
cat("Metric result rows:", nrow(metric_results), "\n")
cat("McNemar result rows:", nrow(mcnemar_results), "\n")

if (nrow(metric_results) > 0) {
  cat(
    "Unavailable model rows:",
    sum(!metric_results$model_available),
    "\n"
  )
}

if (nrow(mcnemar_results) > 0) {
  cat(
    "Nominal directional improvements:",
    sum(
      mcnemar_results$supports_fdrreg_improvement,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "BH-adjusted directional improvements:",
    sum(
      mcnemar_results$supports_fdrreg_improvement_bh,
      na.rm = TRUE
    ),
    "\n"
  )
}

cat("Results directory:", output_dir, "\n")
cat("============================================================\n")
