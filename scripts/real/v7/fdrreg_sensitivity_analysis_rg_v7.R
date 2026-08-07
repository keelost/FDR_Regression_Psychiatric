#!/usr/bin/env Rscript

# Project: FDRreg Sensitivity Analysis
# File: fdrreg_sensitivity_analysis_rg.R
#
# Description:
#   Sensitivity analysis for:
#     1. SNP-level FDRreg
#     2. MAGMA gene-level FDRreg
#     3. SmultiXcan gene-level FDRreg
#
#   The SmultiXcan input files use the v7 directory structure:
#     08.smultixcan_output_v7
#     09.smultixcan_fdrreg_v7
#
#   All sensitivity-analysis output files are written to:
#     01.extra.analysis/06.rg_v7
#
# Usage:
#   Rscript fdrreg_sensitivity_analysis_rg.R \
#       --target scz2012 \
#       --traits_with adhd2019,alco2018,bd2018,ed2019,mdd2019,ocd2018 \
#       --traits_no asd2019,cannabis,dx.dep,ever.csh,ever.sh,ever.st,inte,life.threat.accident,mddco,neo.o,neuroti2018,no.dep,physical.crime,ptsd2019,sleep.dura,st.dep \
#       --lasso N \
#       --seed 100
#
# All code comments are written in English.

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(FDRreg)
  library(stringr)
})

#=============================================================#
#---- Command-line arguments ---------------------------------#
#=============================================================#

option_list <- list(
  make_option(
    "--target",
    type = "character",
    default = NULL,
    help = "Target disease or trait name. Required."
  ),
  make_option(
    "--traits_with",
    type = "character",
    default = "",
    help = "Comma-separated traits with sample overlap."
  ),
  make_option(
    "--traits_no",
    type = "character",
    default = "",
    help = "Comma-separated traits without sample overlap."
  ),
  make_option(
    "--lasso",
    type = "character",
    default = "N",
    help = "Retained for compatibility. LASSO is not used here."
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 100,
    help = "Random seed. Default: 100."
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

if (is.null(opt$target) || nchar(trimws(opt$target)) == 0L) {
  stop("Argument --target is required.")
}

#=============================================================#
#---- Helper functions ---------------------------------------#
#=============================================================#

parse_trait_list <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0L) {
    return(character(0))
  }

  values <- trimws(
    strsplit(x, ",", fixed = TRUE)[[1]]
  )

  values <- values[nchar(values) > 0L]

  unique(values)
}

strip_gene_version <- function(gene_ids) {
  sub(
    "\\..*$",
    "",
    trimws(as.character(gene_ids))
  )
}

p_to_z_safe <- function(p, cap = 8.2) {
  p <- suppressWarnings(as.numeric(p))

  if (any(!is.finite(p))) {
    stop("Non-finite p-values were found.")
  }

  if (any(p < 0 | p > 1)) {
    stop("P-values outside the interval [0, 1] were found.")
  }

  p_clamped <- pmax(
    pmin(p, 1 - .Machine$double.eps),
    .Machine$double.eps
  )

  z <- qnorm(1 - p_clamped / 2)

  pmin(
    pmax(z, 0),
    cap
  )
}

validate_numeric_vector <- function(x, label) {
  x <- suppressWarnings(as.numeric(x))

  if (length(x) == 0L) {
    stop(label, " is empty.")
  }

  if (any(!is.finite(x))) {
    stop(
      sprintf(
        "%s contains %d non-finite values.",
        label,
        sum(!is.finite(x))
      )
    )
  }

  x
}

align_to_reference <- function(
    dt,
    key_col,
    reference_keys,
    data_label) {

  if (!key_col %in% names(dt)) {
    stop(
      sprintf(
        "Required key column '%s' is missing from %s.",
        key_col,
        data_label
      )
    )
  }

  data_keys <- as.character(dt[[key_col]])
  reference_keys <- as.character(reference_keys)

  if (anyNA(data_keys) || any(data_keys == "")) {
    stop(
      "Missing or empty identifiers were found in ",
      data_label,
      "."
    )
  }

  duplicated_data_keys <- unique(
    data_keys[duplicated(data_keys)]
  )

  if (length(duplicated_data_keys) > 0L) {
    stop(
      data_label,
      " contains duplicated identifiers in column '",
      key_col,
      "'. Examples: ",
      paste(head(duplicated_data_keys, 10L), collapse = ", ")
    )
  }

  duplicated_reference_keys <- unique(
    reference_keys[duplicated(reference_keys)]
  )

  if (length(duplicated_reference_keys) > 0L) {
    stop(
      "The reference list contains duplicated identifiers. Examples: ",
      paste(head(duplicated_reference_keys, 10L), collapse = ", ")
    )
  }

  index <- match(
    reference_keys,
    data_keys
  )

  if (anyNA(index)) {
    missing_keys <- reference_keys[is.na(index)]

    stop(
      data_label,
      " does not contain all identifiers from the reference list. ",
      "Missing: ",
      length(missing_keys),
      ". Examples: ",
      paste(head(missing_keys, 10L), collapse = ", ")
    )
  }

  aligned_dt <- dt[index]

  if (!identical(
    as.character(aligned_dt[[key_col]]),
    reference_keys
  )) {
    stop(
      "Internal alignment failure for ",
      data_label
    )
  }

  aligned_dt
}

safe_ratio <- function(numerator, denominator) {
  if (
    is.na(numerator) ||
      is.na(denominator) ||
      denominator <= 0
  ) {
    return(NA_real_)
  }

  numerator / denominator
}

#=============================================================#
#---- Configuration ------------------------------------------#
#=============================================================#

base_output_path <- Sys.getenv("FDRREG_RESULTS_DIR", "")
if (!nzchar(base_output_path)) stop("Set FDRREG_RESULTS_DIR before running v7 rg sensitivity analysis.")

target_name <- trimws(opt$target)
traits_with <- parse_trait_list(opt$traits_with)
traits_no <- parse_trait_list(opt$traits_no)
all_traits <- unique(c(traits_with, traits_no))
random_seed <- opt$seed
sig_threshold <- 0.01

if (length(all_traits) == 0L) {
  stop(
    "At least one trait must be specified through ",
    "--traits_with or --traits_no."
  )
}

if (target_name %in% all_traits) {
  stop(
    "The target trait must not be included among the covariate traits."
  )
}

duplicated_across_groups <- intersect(
  traits_with,
  traits_no
)

if (length(duplicated_across_groups) > 0L) {
  stop(
    "The following traits occur in both --traits_with and --traits_no: ",
    paste(duplicated_across_groups, collapse = ", ")
  )
}

set.seed(random_seed)

start_time <- Sys.time()

#=============================================================#
#---- Version-specific paths ---------------------------------#
#=============================================================#

# Sensitivity-analysis output directory.
sensitivity_output_path <- file.path(
  base_output_path,
  "01.extra.analysis",
  "06.rg_v7"
)

target_output_dir <- file.path(
  sensitivity_output_path,
  target_name
)

dir.create(
  target_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# SNP and MAGMA paths.
magma_input_dir <- file.path(
  base_output_path,
  target_name,
  "01.magma_input"
)

magma_output_dir <- file.path(
  base_output_path,
  target_name,
  "04.magma_output"
)

magma_fdrreg_dir <- file.path(
  base_output_path,
  target_name,
  "05.magma_fdrreg"
)

# SmultiXcan v7 paths.
smultixcan_output_dir <- file.path(
  base_output_path,
  target_name,
  "08.smultixcan_output_v7"
)

smultixcan_fdrreg_dir <- file.path(
  base_output_path,
  target_name,
  "09.smultixcan_fdrreg_v7"
)

cat("============================================================\n")
cat("FDRreg Sensitivity Analysis Configuration\n")
cat("============================================================\n")
cat("Target:                         ", target_name, "\n", sep = "")
cat(
  "Traits with sample overlap:     ",
  ifelse(
    length(traits_with) == 0L,
    "(none)",
    paste(traits_with, collapse = ", ")
  ),
  "\n",
  sep = ""
)
cat(
  "Traits without sample overlap:  ",
  ifelse(
    length(traits_no) == 0L,
    "(none)",
    paste(traits_no, collapse = ", ")
  ),
  "\n",
  sep = ""
)
cat(
  "Total covariate traits:         ",
  length(all_traits),
  "\n",
  sep = ""
)
cat("Random seed:                    ", random_seed, "\n", sep = "")
cat("Significance threshold:         ", sig_threshold, "\n", sep = "")
cat("\nVersion-specific paths:\n")
cat(
  "Sensitivity output:             ",
  target_output_dir,
  "\n",
  sep = ""
)
cat(
  "SmultiXcan v7 raw output:       ",
  smultixcan_output_dir,
  "\n",
  sep = ""
)
cat(
  "SmultiXcan v7 FDRreg results:   ",
  smultixcan_fdrreg_dir,
  "\n",
  sep = ""
)
cat("============================================================\n")

# Summary variables.
snp_fdr_n_sig <- NA_integer_
snp_total <- NA_integer_
snp_qval_n_sig <- NA_integer_
snp_qval_total <- NA_integer_

magma_fdr_n_sig <- NA_integer_
magma_total <- NA_integer_
magma_qval_n_sig <- NA_integer_
magma_qval_total <- NA_integer_

smultixcan_fdr_n_sig <- NA_integer_
smultixcan_total <- NA_integer_
smultixcan_qval_n_sig <- NA_integer_
smultixcan_qval_total <- NA_integer_

overlapping_snps <- character(0)

has_overlap <- length(traits_with) > 0L

#=============================================================#
#---- Part 1: SNP-level sensitivity analysis -----------------#
#=============================================================#

cat("\n============================================================\n")
cat("PART 1: SNP Sensitivity Analysis\n")
cat("============================================================\n")

if (!dir.exists(magma_input_dir)) {
  stop(
    "MAGMA input directory not found: ",
    magma_input_dir
  )
}

get_magma_input_file <- function(trait_name) {
  file.path(
    magma_input_dir,
    paste0(trait_name, ".overlap.4magma.txt")
  )
}

read_magma_z <- function(trait_name) {
  input_file <- get_magma_input_file(trait_name)

  if (!file.exists(input_file)) {
    stop(
      "MAGMA input file not found for trait '",
      trait_name,
      "': ",
      input_file
    )
  }

  header <- names(
    fread(input_file, nrows = 0L)
  )

  required_columns <- c(
    "snpid",
    "z.decor"
  )

  missing_columns <- setdiff(
    required_columns,
    header
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Required columns are missing from ",
      input_file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dt <- fread(
    input_file,
    select = required_columns
  )

  if (nrow(dt) == 0L) {
    stop(
      "MAGMA input file is empty: ",
      input_file
    )
  }

  dt[, snpid := trimws(as.character(snpid))]

  if (
    anyNA(dt$snpid) ||
      any(dt$snpid == "")
  ) {
    stop(
      "Missing or empty snpid values found in: ",
      input_file
    )
  }

  if (anyDuplicated(dt$snpid)) {
    duplicated_ids <- unique(
      dt$snpid[duplicated(dt$snpid)]
    )

    stop(
      "Duplicated snpid values found in ",
      input_file,
      ". Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }

  dt[
    ,
    z.decor := validate_numeric_vector(
      z.decor,
      paste0("z.decor for ", trait_name)
    )
  ]

  dt
}

cat("Loading target MAGMA input...\n")

target_snp_data <- read_magma_z(target_name)

overlapping_snps <- target_snp_data$snpid
target_z_snp <- target_snp_data$z.decor

cat(
  "Loaded target SNPs: ",
  length(overlapping_snps),
  "\n",
  sep = ""
)

trait_snp_data <- setNames(
  vector("list", length(all_traits)),
  all_traits
)

for (trait in all_traits) {
  trait_dt <- read_magma_z(trait)

  trait_dt <- align_to_reference(
    dt = trait_dt,
    key_col = "snpid",
    reference_keys = overlapping_snps,
    data_label = paste0(
      "SNP MAGMA input for ",
      trait
    )
  )

  trait_snp_data[[trait]] <- trait_dt

  cat(
    "Loaded and aligned SNP input: ",
    trait,
    " (",
    nrow(trait_dt),
    " SNPs)\n",
    sep = ""
  )
}

cat("Building SNP feature matrix...\n")

feature_matrix_snp <- do.call(
  cbind,
  lapply(
    all_traits,
    function(trait) {
      abs(trait_snp_data[[trait]]$z.decor)
    }
  )
)

if (is.null(dim(feature_matrix_snp))) {
  feature_matrix_snp <- matrix(
    feature_matrix_snp,
    ncol = 1L
  )
}

colnames(feature_matrix_snp) <- all_traits
storage.mode(feature_matrix_snp) <- "double"

if (
  nrow(feature_matrix_snp) !=
    length(target_z_snp)
) {
  stop(
    "SNP feature matrix row count does not match ",
    "the target SNP z-score length."
  )
}

if (any(!is.finite(feature_matrix_snp))) {
  stop(
    "The SNP feature matrix contains non-finite values."
  )
}

cat(
  "SNP feature matrix dimensions: ",
  nrow(feature_matrix_snp),
  " x ",
  ncol(feature_matrix_snp),
  "\n",
  sep = ""
)

cat("Running FDRreg for SNPs...\n")

set.seed(random_seed)

fdr_snp_theoretical <- FDRreg(
  target_z_snp,
  feature_matrix_snp,
  nulltype = "theoretical",
  method = "pr"
)

fdr_snp <- fdr_snp_theoretical$FDR

snp_fdr_n_sig <- sum(
  fdr_snp < sig_threshold,
  na.rm = TRUE
)

snp_total <- length(fdr_snp)

cat(
  "SNP fdr.the < ",
  sig_threshold,
  ": ",
  snp_fdr_n_sig,
  " / ",
  snp_total,
  "\n",
  sep = ""
)

original_snp_fdr_file <- file.path(
  base_output_path,
  target_name,
  "02.fdrreg_results",
  "fdr_values_per_snp_nolasso.csv"
)

if (file.exists(original_snp_fdr_file)) {
  original_snp_results <- fread(
    original_snp_fdr_file
  )

  if ("qval_bh" %in% names(original_snp_results)) {
    snp_qval_n_sig <- sum(
      original_snp_results$qval_bh < sig_threshold,
      na.rm = TRUE
    )

    snp_qval_total <- nrow(
      original_snp_results
    )

    cat(
      "Original SNP qval < ",
      sig_threshold,
      ": ",
      snp_qval_n_sig,
      " / ",
      snp_qval_total,
      "\n",
      sep = ""
    )
  } else {
    warning(
      "qval_bh was not found in ",
      original_snp_fdr_file
    )
  }
} else {
  warning(
    "Original SNP result file was not found: ",
    original_snp_fdr_file
  )
}

snp_results_dt <- data.table(
  snpid = overlapping_snps,
  z_score_signed = target_z_snp,
  fdr_the = fdr_snp
)

snp_output_file <- file.path(
  target_output_dir,
  paste0(
    target_name,
    "_snp_fdr_the.csv"
  )
)

fwrite(
  snp_results_dt,
  snp_output_file
)

cat(
  "Written SNP sensitivity results: ",
  snp_output_file,
  "\n",
  sep = ""
)

rm(
  target_snp_data,
  trait_snp_data,
  feature_matrix_snp,
  fdr_snp_theoretical
)

gc()

#=============================================================#
#---- Part 2: MAGMA gene-level sensitivity analysis -----------#
#=============================================================#

cat("\n============================================================\n")
cat("PART 2: MAGMA Sensitivity Analysis\n")
cat("============================================================\n")

detect_magma_gene_column <- function(dt, data_label) {
  candidates <- c(
    "GENE",
    "ID"
  )

  matched <- candidates[
    candidates %in% names(dt)
  ]

  if (length(matched) == 0L) {
    stop(
      "No MAGMA gene identifier column was found in ",
      data_label,
      ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }

  matched[1]
}

detect_magma_p_column <- function(dt, data_label) {
  candidates <- c(
    "P",
    "p",
    "pvalue",
    "PVAL",
    "pval"
  )

  matched <- candidates[
    candidates %in% names(dt)
  ]

  if (length(matched) == 0L) {
    stop(
      "No MAGMA p-value column was found in ",
      data_label,
      ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }

  matched[1]
}

standardize_magma_gene_id <- function(dt, data_label) {
  gene_column <- detect_magma_gene_column(
    dt,
    data_label
  )

  dt[
    ,
    GENE_KEY := trimws(as.character(get(gene_column)))
  ]

  if (
    anyNA(dt$GENE_KEY) ||
      any(dt$GENE_KEY == "")
  ) {
    stop(
      "Missing or empty MAGMA gene IDs were found in ",
      data_label
    )
  }

  duplicated_ids <- unique(
    dt$GENE_KEY[
      duplicated(dt$GENE_KEY)
    ]
  )

  if (length(duplicated_ids) > 0L) {
    stop(
      "Duplicated MAGMA gene IDs were found in ",
      data_label,
      ". Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }

  attr(
    dt,
    "original_gene_column"
  ) <- gene_column

  dt
}

original_magma_gene_list <- NULL

original_magma_file <- file.path(
  magma_fdrreg_dir,
  paste0(
    target_name,
    ".gene.bio.fdrreg.txt"
  )
)

if (file.exists(original_magma_file)) {
  original_magma_results <- fread(
    original_magma_file
  )

  original_magma_results <- standardize_magma_gene_id(
    original_magma_results,
    original_magma_file
  )

  original_magma_gene_list <- original_magma_results$GENE_KEY

  magma_qval_total <- nrow(
    original_magma_results
  )

  if ("qval" %in% names(original_magma_results)) {
    magma_qval_n_sig <- sum(
      original_magma_results$qval < sig_threshold,
      na.rm = TRUE
    )
  }

  cat(
    "Original MAGMA reference genes: ",
    magma_qval_total,
    "\n",
    sep = ""
  )
} else {
  warning(
    "Original MAGMA FDRreg file was not found: ",
    original_magma_file
  )
}

target_magma_file <- file.path(
  magma_output_dir,
  paste0(
    target_name,
    ".genes.out"
  )
)

if (
  !file.exists(target_magma_file)
) {
  warning(
    "Target MAGMA file was not found: ",
    target_magma_file
  )
} else {
  target_magma <- fread(
    target_magma_file
  )

  target_magma <- standardize_magma_gene_id(
    target_magma,
    target_magma_file
  )

  target_p_column <- detect_magma_p_column(
    target_magma,
    target_magma_file
  )

  target_magma[
    ,
    P_VALUE := suppressWarnings(
      as.numeric(get(target_p_column))
    )
  ]

  if (any(!is.finite(target_magma$P_VALUE))) {
    stop(
      "Non-finite MAGMA p-values were found in ",
      target_magma_file
    )
  }

  if (
    any(
      target_magma$P_VALUE < 0 |
        target_magma$P_VALUE > 1
    )
  ) {
    stop(
      "MAGMA p-values outside [0, 1] were found in ",
      target_magma_file
    )
  }

  trait_magma_data <- list()
  valid_magma_traits <- character(0)

  for (trait in all_traits) {
    trait_file <- file.path(
      magma_output_dir,
      paste0(
        trait,
        ".genes.out"
      )
    )

    if (!file.exists(trait_file)) {
      warning(
        "MAGMA file not found for trait ",
        trait,
        ": ",
        trait_file
      )
      next
    }

    dt <- fread(
      trait_file
    )

    dt <- tryCatch(
      {
        standardize_magma_gene_id(
          dt,
          trait_file
        )
      },
      error = function(e) {
        warning(conditionMessage(e))
        NULL
      }
    )

    if (is.null(dt)) {
      next
    }

    p_column <- tryCatch(
      {
        detect_magma_p_column(
          dt,
          trait_file
        )
      },
      error = function(e) {
        warning(conditionMessage(e))
        NULL
      }
    )

    if (is.null(p_column)) {
      next
    }

    dt[
      ,
      P_VALUE := suppressWarnings(
        as.numeric(get(p_column))
      )
    ]

    if (any(!is.finite(dt$P_VALUE))) {
      warning(
        "Non-finite p-values found for MAGMA trait ",
        trait
      )
      next
    }

    if (
      any(
        dt$P_VALUE < 0 |
          dt$P_VALUE > 1
      )
    ) {
      warning(
        "P-values outside [0, 1] found for MAGMA trait ",
        trait
      )
      next
    }

    trait_magma_data[[trait]] <- dt
    valid_magma_traits <- c(
      valid_magma_traits,
      trait
    )
  }

  if (length(valid_magma_traits) > 0L) {
    if (!is.null(original_magma_gene_list)) {
      reference_magma_genes <- original_magma_gene_list
    } else {
      reference_magma_genes <- Reduce(
        intersect,
        c(
          list(target_magma$GENE_KEY),
          lapply(
            trait_magma_data[valid_magma_traits],
            function(x) x$GENE_KEY
          )
        )
      )

      reference_magma_genes <- sort(
        unique(reference_magma_genes)
      )
    }

    target_magma <- align_to_reference(
      dt = target_magma,
      key_col = "GENE_KEY",
      reference_keys = reference_magma_genes,
      data_label = paste0(
        "Target MAGMA output for ",
        target_name
      )
    )

    aligned_trait_magma_data <- setNames(
      vector(
        "list",
        length(valid_magma_traits)
      ),
      valid_magma_traits
    )

    for (trait in valid_magma_traits) {
      aligned_trait_magma_data[[trait]] <- align_to_reference(
        dt = trait_magma_data[[trait]],
        key_col = "GENE_KEY",
        reference_keys = reference_magma_genes,
        data_label = paste0(
          "MAGMA output for ",
          trait
        )
      )
    }

    trait_magma_data <- aligned_trait_magma_data

    feature_matrix_magma <- do.call(
      cbind,
      lapply(
        trait_magma_data,
        function(dt) {
          p_to_z_safe(dt$P_VALUE)
        }
      )
    )

    if (is.null(dim(feature_matrix_magma))) {
      feature_matrix_magma <- matrix(
        feature_matrix_magma,
        ncol = 1L
      )
    }

    colnames(feature_matrix_magma) <- valid_magma_traits
    storage.mode(feature_matrix_magma) <- "double"

    target_z_magma_unsigned <- p_to_z_safe(
      target_magma$P_VALUE
    )

    set.seed(random_seed)

    target_z_magma <- target_z_magma_unsigned *
      sign(
        rnorm(
          length(target_z_magma_unsigned)
        )
      )

    fdr_magma_theoretical <- FDRreg(
      target_z_magma,
      feature_matrix_magma,
      nulltype = "theoretical",
      method = "pr"
    )

    fdr_magma <- fdr_magma_theoretical$FDR

    magma_fdr_n_sig <- sum(
      fdr_magma < sig_threshold,
      na.rm = TRUE
    )

    magma_total <- length(fdr_magma)

    magma_results_dt <- data.table(
      GENE = target_magma$GENE_KEY,
      target_z = target_z_magma,
      fdr_the = fdr_magma
    )

    magma_output_file <- file.path(
      target_output_dir,
      paste0(
        target_name,
        "_magma_fdr_the.csv"
      )
    )

    fwrite(
      magma_results_dt,
      magma_output_file
    )

    cat(
      "MAGMA fdr.the < ",
      sig_threshold,
      ": ",
      magma_fdr_n_sig,
      " / ",
      magma_total,
      "\n",
      sep = ""
    )
  } else {
    warning(
      "No valid MAGMA trait files were found."
    )
  }
}

#=============================================================#
#---- Part 3: SmultiXcan v7 sensitivity analysis --------------#
#=============================================================#

cat("\n============================================================\n")
cat("PART 3: SmultiXcan v7 Sensitivity Analysis\n")
cat("============================================================\n")

cat(
  "SmultiXcan raw v7 directory: ",
  smultixcan_output_dir,
  "\n",
  sep = ""
)

cat(
  "SmultiXcan FDRreg v7 directory: ",
  smultixcan_fdrreg_dir,
  "\n",
  sep = ""
)

if (!dir.exists(smultixcan_output_dir)) {
  stop(
    "SmultiXcan v7 output directory does not exist: ",
    smultixcan_output_dir
  )
}

if (!dir.exists(smultixcan_fdrreg_dir)) {
  stop(
    "SmultiXcan v7 FDRreg directory does not exist: ",
    smultixcan_fdrreg_dir
  )
}

#-------------------------------------------------------------#
#---- SmultiXcan helper functions ----------------------------#
#-------------------------------------------------------------#

is_invalid_pvalue <- function(x) {
  x <- suppressWarnings(as.numeric(x))

  is.na(x) |
    !is.finite(x) |
    x < 0 |
    x > 1
}

read_and_align_smultixcan_v7 <- function(
    input_file,
    data_label,
    reference_genes,
    require_gene_name = FALSE,
    diagnostic_prefix) {

  if (!file.exists(input_file)) {
    stop(
      data_label,
      " file does not exist: ",
      input_file
    )
  }

  input_header <- names(
    fread(
      input_file,
      nrows = 0L
    )
  )

  required_columns <- c(
    "gene",
    "pvalue"
  )

  if (require_gene_name) {
    required_columns <- c(
      required_columns,
      "gene_name"
    )
  }

  missing_columns <- setdiff(
    required_columns,
    input_header
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Required columns are missing from ",
      data_label,
      ": ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  optional_columns <- intersect(
    c(
      "status",
      "n",
      "n_indep",
      "tmi"
    ),
    input_header
  )

  selected_columns <- unique(
    c(
      required_columns,
      optional_columns
    )
  )

  # Read pvalue as character to preserve NA and other raw values.
  dt <- fread(
    input_file,
    select = selected_columns,
    colClasses = c(
      pvalue = "character"
    )
  )

  if (nrow(dt) == 0L) {
    stop(
      data_label,
      " is empty: ",
      input_file
    )
  }

  dt[
    ,
    gene := trimws(as.character(gene))
  ]

  dt[
    ,
    gene_stripped := strip_gene_version(gene)
  ]

  if (
    anyNA(dt$gene_stripped) ||
      any(dt$gene_stripped == "")
  ) {
    stop(
      "Missing or empty gene identifiers were found in ",
      data_label,
      "."
    )
  }

  duplicated_genes <- unique(
    dt$gene_stripped[
      duplicated(dt$gene_stripped)
    ]
  )

  if (length(duplicated_genes) > 0L) {
    stop(
      "Duplicated gene identifiers were found in ",
      data_label,
      " after Ensembl version stripping. Examples: ",
      paste(
        head(duplicated_genes, 10L),
        collapse = ", "
      )
    )
  }

  dt[
    ,
    pvalue_raw := as.character(pvalue)
  ]

  dt[
    ,
    pvalue_numeric := suppressWarnings(
      as.numeric(pvalue_raw)
    )
  ]

  invalid_all <- is_invalid_pvalue(
    dt$pvalue_numeric
  )

  input_gene_set <- unique(
    dt$gene_stripped
  )

  missing_reference_genes <- setdiff(
    reference_genes,
    input_gene_set
  )

  extra_input_genes <- setdiff(
    input_gene_set,
    reference_genes
  )

  invalid_outside_reference <- invalid_all &
    !dt$gene_stripped %chin% reference_genes

  invalid_inside_reference_before_alignment <- invalid_all &
    dt$gene_stripped %chin% reference_genes

  cat("\n")
  cat(
    data_label,
    "\n",
    sep = ""
  )
  cat(
    "  Input file: ",
    input_file,
    "\n",
    sep = ""
  )
  cat(
    "  Raw genes: ",
    nrow(dt),
    "\n",
    sep = ""
  )
  cat(
    "  Reference genes: ",
    length(reference_genes),
    "\n",
    sep = ""
  )
  cat(
    "  Reference genes missing from raw output: ",
    length(missing_reference_genes),
    "\n",
    sep = ""
  )
  cat(
    "  Raw genes outside reference list: ",
    length(extra_input_genes),
    "\n",
    sep = ""
  )
  cat(
    "  Invalid p-values in complete raw output: ",
    sum(invalid_all),
    "\n",
    sep = ""
  )
  cat(
    "  Invalid p-values outside reference list: ",
    sum(invalid_outside_reference),
    "\n",
    sep = ""
  )
  cat(
    "  Invalid p-values inside reference list: ",
    sum(invalid_inside_reference_before_alignment),
    "\n",
    sep = ""
  )

  # Save all invalid raw records for documentation.
  if (any(invalid_all)) {
    invalid_all_file <- file.path(
      target_output_dir,
      paste0(
        diagnostic_prefix,
        "_invalid_pvalues_all_raw.csv"
      )
    )

    fwrite(
      dt[
        which(invalid_all)
      ],
      invalid_all_file
    )

    cat(
      "  Invalid raw p-value records written to: ",
      invalid_all_file,
      "\n",
      sep = ""
    )
  }

  if (length(missing_reference_genes) > 0L) {
    missing_reference_file <- file.path(
      target_output_dir,
      paste0(
        diagnostic_prefix,
        "_missing_reference_genes.csv"
      )
    )

    fwrite(
      data.table(
        ENSEMBL_GENE_ID = missing_reference_genes
      ),
      missing_reference_file
    )

    stop(
      data_label,
      " does not contain ",
      length(missing_reference_genes),
      " genes from the SmultiXcan v7 FDRreg reference list. ",
      "Examples: ",
      paste(
        head(missing_reference_genes, 10L),
        collapse = ", "
      ),
      ". Details were written to: ",
      missing_reference_file
    )
  }

  # Align first. Invalid p-values outside the FDRreg reference list
  # are deliberately excluded from the sensitivity analysis.
  aligned_dt <- align_to_reference(
    dt = dt,
    key_col = "gene_stripped",
    reference_keys = reference_genes,
    data_label = data_label
  )

  aligned_dt[
    ,
    pvalue := pvalue_numeric
  ]

  invalid_reference <- is_invalid_pvalue(
    aligned_dt$pvalue
  )

  cat(
    "  Invalid p-values after reference alignment: ",
    sum(invalid_reference),
    "\n",
    sep = ""
  )

  if (any(invalid_reference)) {
    invalid_reference_file <- file.path(
      target_output_dir,
      paste0(
        diagnostic_prefix,
        "_invalid_pvalues_in_reference.csv"
      )
    )

    fwrite(
      aligned_dt[
        which(invalid_reference)
      ],
      invalid_reference_file
    )

    stop(
      data_label,
      " contains ",
      sum(invalid_reference),
      " invalid p-values among the aligned FDRreg reference genes. ",
      "Details were written to: ",
      invalid_reference_file
    )
  }

  if (!identical(
    aligned_dt$gene_stripped,
    as.character(reference_genes)
  )) {
    stop(
      "Final SmultiXcan gene alignment failed for ",
      data_label,
      "."
    )
  }

  cat(
    "  Successfully aligned valid reference genes: ",
    nrow(aligned_dt),
    "\n",
    sep = ""
  )

  aligned_dt
}

#-------------------------------------------------------------#
#---- Load the SmultiXcan v7 FDRreg reference results --------#
#-------------------------------------------------------------#

original_smultixcan_file <- file.path(
  smultixcan_fdrreg_dir,
  paste0(
    target_name,
    ".gene.bio.fdrreg.txt"
  )
)

if (!file.exists(original_smultixcan_file)) {
  stop(
    "SmultiXcan v7 FDRreg file does not exist: ",
    original_smultixcan_file
  )
}

original_smultixcan_header <- names(
  fread(
    original_smultixcan_file,
    nrows = 0L
  )
)

if (
  !"ENSEMBL_GENE_ID" %in%
    original_smultixcan_header
) {
  stop(
    "ENSEMBL_GENE_ID is missing from the SmultiXcan v7 ",
    "FDRreg file: ",
    original_smultixcan_file
  )
}

original_smultixcan_results <- fread(
  original_smultixcan_file
)

original_smultixcan_results[
  ,
  gene_stripped := strip_gene_version(
    ENSEMBL_GENE_ID
  )
]

if (
  anyNA(original_smultixcan_results$gene_stripped) ||
    any(original_smultixcan_results$gene_stripped == "")
) {
  stop(
    "Missing or empty ENSEMBL_GENE_ID values were found in ",
    original_smultixcan_file
  )
}

duplicated_reference_genes <- unique(
  original_smultixcan_results$gene_stripped[
    duplicated(
      original_smultixcan_results$gene_stripped
    )
  ]
)

if (length(duplicated_reference_genes) > 0L) {
  stop(
    "Duplicated ENSEMBL_GENE_ID values were found in the ",
    "SmultiXcan v7 FDRreg reference file after version stripping. ",
    "Examples: ",
    paste(
      head(duplicated_reference_genes, 10L),
      collapse = ", "
    )
  )
}

reference_smultixcan_genes <-
  original_smultixcan_results$gene_stripped

smultixcan_qval_total <- nrow(
  original_smultixcan_results
)

if ("qval" %in% names(original_smultixcan_results)) {
  original_smultixcan_results[
    ,
    qval_numeric := suppressWarnings(
      as.numeric(qval)
    )
  ]

  invalid_reference_qval <- (
    is.na(original_smultixcan_results$qval_numeric) |
      !is.finite(
        original_smultixcan_results$qval_numeric
      ) |
      original_smultixcan_results$qval_numeric < 0 |
      original_smultixcan_results$qval_numeric > 1
  )

  if (any(invalid_reference_qval)) {
    warning(
      "The original SmultiXcan v7 FDRreg result contains ",
      sum(invalid_reference_qval),
      " invalid q-values. These values are excluded from the ",
      "qval significance count."
    )
  }

  smultixcan_qval_n_sig <- sum(
    original_smultixcan_results$qval_numeric <
      sig_threshold,
    na.rm = TRUE
  )
} else {
  smultixcan_qval_n_sig <- NA_integer_

  warning(
    "qval column was not found in ",
    original_smultixcan_file
  )
}

cat(
  "Original SmultiXcan v7 reference genes: ",
  length(reference_smultixcan_genes),
  "\n",
  sep = ""
)

cat(
  "Original SmultiXcan v7 qval < ",
  sig_threshold,
  ": ",
  ifelse(
    is.na(smultixcan_qval_n_sig),
    "NA",
    as.character(smultixcan_qval_n_sig)
  ),
  " / ",
  smultixcan_qval_total,
  "\n",
  sep = ""
)

#-------------------------------------------------------------#
#---- Load and align the target SmultiXcan v7 output ----------#
#-------------------------------------------------------------#

target_smultixcan_file <- file.path(
  smultixcan_output_dir,
  paste0(
    target_name,
    ".allbrain.txt"
  )
)

target_smultixcan <- read_and_align_smultixcan_v7(
  input_file = target_smultixcan_file,
  data_label = paste0(
    "Target SmultiXcan v7 output for ",
    target_name
  ),
  reference_genes = reference_smultixcan_genes,
  require_gene_name = TRUE,
  diagnostic_prefix = paste0(
    target_name,
    "_smultixcan_v7_target"
  )
)

#-------------------------------------------------------------#
#---- Load and align covariate SmultiXcan v7 outputs ----------#
#-------------------------------------------------------------#

trait_smultixcan_data <- setNames(
  vector(
    "list",
    length(all_traits)
  ),
  all_traits
)

for (trait in all_traits) {
  trait_file <- file.path(
    smultixcan_output_dir,
    paste0(
      trait,
      ".allbrain.txt"
    )
  )

  trait_smultixcan_data[[trait]] <-
    read_and_align_smultixcan_v7(
      input_file = trait_file,
      data_label = paste0(
        "SmultiXcan v7 output for ",
        trait
      ),
      reference_genes = reference_smultixcan_genes,
      require_gene_name = FALSE,
      diagnostic_prefix = paste0(
        target_name,
        "_",
        trait,
        "_smultixcan_v7"
      )
    )
}

valid_smultixcan_traits <- names(
  trait_smultixcan_data
)

if (length(valid_smultixcan_traits) == 0L) {
  stop(
    "No SmultiXcan v7 covariate traits are available."
  )
}

#-------------------------------------------------------------#
#---- Final alignment validation -----------------------------#
#-------------------------------------------------------------#

if (!identical(
  target_smultixcan$gene_stripped,
  reference_smultixcan_genes
)) {
  stop(
    "Target SmultiXcan v7 genes are not in the exact ",
    "FDRreg reference order."
  )
}

trait_alignment_ok <- vapply(
  trait_smultixcan_data,
  function(dt) {
    identical(
      dt$gene_stripped,
      reference_smultixcan_genes
    )
  },
  logical(1)
)

if (!all(trait_alignment_ok)) {
  stop(
    "Final SmultiXcan v7 gene alignment failed for: ",
    paste(
      names(trait_alignment_ok)[
        !trait_alignment_ok
      ],
      collapse = ", "
    )
  )
}

cat("\nAll SmultiXcan v7 datasets passed reference alignment.\n")
cat(
  "Aligned genes per dataset: ",
  length(reference_smultixcan_genes),
  "\n",
  sep = ""
)
cat(
  "Aligned covariate traits: ",
  length(valid_smultixcan_traits),
  "\n",
  sep = ""
)

#-------------------------------------------------------------#
#---- Build the SmultiXcan feature matrix --------------------#
#-------------------------------------------------------------#

feature_matrix_smultixcan <- do.call(
  cbind,
  lapply(
    trait_smultixcan_data,
    function(dt) {
      p_to_z_safe(
        dt$pvalue
      )
    }
  )
)

if (is.null(dim(feature_matrix_smultixcan))) {
  feature_matrix_smultixcan <- matrix(
    feature_matrix_smultixcan,
    ncol = 1L
  )
}

colnames(feature_matrix_smultixcan) <-
  valid_smultixcan_traits

storage.mode(feature_matrix_smultixcan) <-
  "double"

if (
  nrow(feature_matrix_smultixcan) !=
    length(reference_smultixcan_genes)
) {
  stop(
    "SmultiXcan feature matrix row count does not match ",
    "the FDRreg reference gene count."
  )
}

if (
  ncol(feature_matrix_smultixcan) !=
    length(valid_smultixcan_traits)
) {
  stop(
    "SmultiXcan feature matrix column count does not match ",
    "the number of covariate traits."
  )
}

if (any(!is.finite(feature_matrix_smultixcan))) {
  stop(
    "The SmultiXcan feature matrix contains non-finite values ",
    "after reference alignment."
  )
}

cat(
  "SmultiXcan feature matrix dimensions: ",
  nrow(feature_matrix_smultixcan),
  " x ",
  ncol(feature_matrix_smultixcan),
  "\n",
  sep = ""
)

#-------------------------------------------------------------#
#---- Construct the target z-score vector --------------------#
#-------------------------------------------------------------#

target_z_smultixcan_unsigned <- p_to_z_safe(
  target_smultixcan$pvalue
)

set.seed(random_seed)

random_signs_smultixcan <- ifelse(
  rnorm(
    length(target_z_smultixcan_unsigned)
  ) >= 0,
  1,
  -1
)

target_z_smultixcan <-
  target_z_smultixcan_unsigned *
  random_signs_smultixcan

if (
  length(target_z_smultixcan) !=
    nrow(feature_matrix_smultixcan)
) {
  stop(
    "The target SmultiXcan z-score length does not match ",
    "the feature matrix row count."
  )
}

if (any(!is.finite(target_z_smultixcan))) {
  stop(
    "The target SmultiXcan z-score vector contains ",
    "non-finite values."
  )
}

#-------------------------------------------------------------#
#---- Run SmultiXcan FDRreg ----------------------------------#
#-------------------------------------------------------------#

cat(
  "Running FDRreg for SmultiXcan v7 using ",
  length(reference_smultixcan_genes),
  " reference genes and ",
  length(valid_smultixcan_traits),
  " covariate traits...\n",
  sep = ""
)

set.seed(random_seed)

fdr_smultixcan_theoretical <- FDRreg(
  target_z_smultixcan,
  feature_matrix_smultixcan,
  nulltype = "theoretical",
  method = "pr"
)

fdr_smultixcan <-
  fdr_smultixcan_theoretical$FDR

if (
  length(fdr_smultixcan) !=
    length(reference_smultixcan_genes)
) {
  stop(
    "The SmultiXcan FDR vector length does not match ",
    "the reference gene count."
  )
}

smultixcan_fdr_n_sig <- sum(
  fdr_smultixcan < sig_threshold,
  na.rm = TRUE
)

smultixcan_total <- length(
  fdr_smultixcan
)

cat(
  "SmultiXcan v7 fdr.the < ",
  sig_threshold,
  ": ",
  smultixcan_fdr_n_sig,
  " / ",
  smultixcan_total,
  "\n",
  sep = ""
)

#-------------------------------------------------------------#
#---- Save SmultiXcan sensitivity results --------------------#
#-------------------------------------------------------------#

smultixcan_results_dt <- data.table(
  ENSEMBL_GENE_ID =
    target_smultixcan$gene_stripped,
  gene_name =
    target_smultixcan$gene_name,
  pvalue =
    target_smultixcan$pvalue,
  target_z =
    target_z_smultixcan,
  fdr_the =
    fdr_smultixcan
)

smultixcan_output_file <- file.path(
  target_output_dir,
  paste0(
    target_name,
    "_smultixcan_fdr_the.csv"
  )
)

fwrite(
  smultixcan_results_dt,
  smultixcan_output_file
)

cat(
  "Written SmultiXcan sensitivity results: ",
  smultixcan_output_file,
  "\n",
  sep = ""
)

# Save a dataset-level alignment summary.
smultixcan_alignment_summary <- data.table(
  dataset = c(
    target_name,
    all_traits
  ),
  role = c(
    "target",
    rep(
      "covariate",
      length(all_traits)
    )
  ),
  reference_gene_count = length(
    reference_smultixcan_genes
  ),
  aligned_gene_count = c(
    nrow(target_smultixcan),
    vapply(
      trait_smultixcan_data,
      nrow,
      integer(1)
    )
  ),
  all_reference_pvalues_valid = TRUE
)

smultixcan_alignment_summary_file <- file.path(
  target_output_dir,
  paste0(
    target_name,
    "_smultixcan_v7_alignment_summary.csv"
  )
)

fwrite(
  smultixcan_alignment_summary,
  smultixcan_alignment_summary_file
)

cat(
  "Written SmultiXcan alignment summary: ",
  smultixcan_alignment_summary_file,
  "\n",
  sep = ""
)

#=============================================================#
#---- Combined summary table ---------------------------------#
#=============================================================#

cat("\n============================================================\n")
cat("Creating Combined Summary Table\n")
cat("============================================================\n")

combined_summary <- data.frame(
  analysis = c(
    "snp",
    "magma",
    "smultixcan"
  ),
  fdr_the_n_sig = c(
    snp_fdr_n_sig,
    magma_fdr_n_sig,
    smultixcan_fdr_n_sig
  ),
  fdr_the_total = c(
    snp_total,
    magma_total,
    smultixcan_total
  ),
  qval_n_sig = c(
    snp_qval_n_sig,
    magma_qval_n_sig,
    smultixcan_qval_n_sig
  ),
  qval_total = c(
    snp_qval_total,
    magma_qval_total,
    smultixcan_qval_total
  ),
  sensitivity_ratio = c(
    safe_ratio(
      snp_fdr_n_sig,
      snp_qval_n_sig
    ),
    safe_ratio(
      magma_fdr_n_sig,
      magma_qval_n_sig
    ),
    safe_ratio(
      smultixcan_fdr_n_sig,
      smultixcan_qval_n_sig
    )
  ),
  stringsAsFactors = FALSE
)

print(combined_summary)

combined_summary_file <- file.path(
  target_output_dir,
  paste0(
    target_name,
    "_combined_sensitivity_summary.csv"
  )
)

fwrite(
  combined_summary,
  combined_summary_file
)

#=============================================================#
#---- Save analysis configuration ----------------------------#
#=============================================================#

config_info <- list(
  target = target_name,
  traits_with = traits_with,
  traits_no = traits_no,
  all_traits = all_traits,
  seed = random_seed,
  threshold = sig_threshold,
  has_snp_decorrelation = has_overlap,
  snp_source = "01.magma_input/*.overlap.4magma.txt",
  snp_z_column = "z.decor",
  magma_source = "04.magma_output/*.genes.out",
  magma_fdrreg_source = paste0(
    "05.magma_fdrreg/",
    target_name,
    ".gene.bio.fdrreg.txt"
  ),
  smultixcan_output_source = paste0(
    "08.smultixcan_output_v7/",
    target_name,
    ".allbrain.txt"
  ),
  smultixcan_fdrreg_source = paste0(
    "09.smultixcan_fdrreg_v7/",
    target_name,
    ".gene.bio.fdrreg.txt"
  ),
  sensitivity_output_directory = target_output_dir,
  snp_count = length(overlapping_snps),
  magma_gene_count = magma_total,
  smultixcan_gene_count = smultixcan_total,
  timestamp = Sys.time()
)

config_file <- file.path(
  target_output_dir,
  paste0(
    target_name,
    "_sensitivity_config_v7.rds"
  )
)

saveRDS(
  config_info,
  config_file
)

#=============================================================#
#---- Completion summary -------------------------------------#
#=============================================================#

end_time <- Sys.time()

cat("\n============================================================\n")
cat("Sensitivity Analysis Complete\n")
cat("============================================================\n")
cat("Target:          ", target_name, "\n", sep = "")
cat(
  "Results saved to: ",
  target_output_dir,
  "\n",
  sep = ""
)
cat(
  "Summary file:     ",
  combined_summary_file,
  "\n",
  sep = ""
)
cat(
  "Configuration:    ",
  config_file,
  "\n",
  sep = ""
)
cat(
  "Time start:       ",
  format(
    start_time,
    "%a %b %d %X %Y"
  ),
  "\n",
  sep = ""
)
cat(
  "Time end:         ",
  format(
    end_time,
    "%a %b %d %X %Y"
  ),
  "\n",
  sep = ""
)
cat(
  "Duration seconds: ",
  round(
    as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    2
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")

q("no")
