#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "Usage: Rscript diagnose_metaxcan_alignment.R ",
      "<target> <brain_region> [base_input_dir]"
    )
  )
}

TARGET <- args[1]
brain_region <- args[2]

if (length(args) >= 3) {
  base_input_dir <- args[3]
} else {
  base_input_dir <- Sys.getenv("FDRREG_RESULTS_DIR", "")
}

input_dir <- file.path(base_input_dir, TARGET, "06.metaxcan")

extract_trait_from_filename <- function(filename, region) {
  trait_name <- basename(filename)
  trait_name <- sub("^gtex_v8_", "", trait_name)

  suffix_pattern <- paste0(
    "\\.overlap\\.4magma_in_",
    region,
    "\\.csv$"
  )

  trait_name <- sub(suffix_pattern, "", trait_name)
  return(trait_name)
}

pattern <- paste0("_in_", brain_region, "\\.csv$")
filenames <- list.files(
  path = input_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(filenames) == 0) {
  stop(
    "No files found in: ",
    input_dir,
    " for region: ",
    brain_region
  )
}

trait_names <- vapply(
  filenames,
  extract_trait_from_filename,
  character(1),
  region = brain_region
)

target_idx <- which(trait_names == TARGET)

if (length(target_idx) != 1) {
  cat("Detected trait names:\n")
  print(data.frame(
    trait = trait_names,
    file = filenames,
    stringsAsFactors = FALSE
  ))

  stop(
    "Expected exactly one target file, but found: ",
    length(target_idx)
  )
}

cat("============================================================\n")
cat("MetaXcan alignment diagnostic\n")
cat("============================================================\n")
cat("Target:", TARGET, "\n")
cat("Brain region:", brain_region, "\n")
cat("Input directory:", input_dir, "\n")
cat("Number of files:", length(filenames), "\n\n")

data_list <- vector("list", length(filenames))
file_summary <- vector("list", length(filenames))

for (i in seq_along(filenames)) {
  dt <- fread(filenames[i], data.table = FALSE)

  required_columns <- c("gene", "zscore", "pvalue")
  missing_columns <- setdiff(required_columns, names(dt))

  if (length(missing_columns) > 0) {
    file_summary[[i]] <- data.frame(
      trait = trait_names[i],
      file = filenames[i],
      raw_rows = nrow(dt),
      valid_rows = NA_integer_,
      unique_genes = NA_integer_,
      duplicated_gene_rows = NA_integer_,
      invalid_gene_rows = NA_integer_,
      invalid_zscore_rows = NA_integer_,
      invalid_pvalue_rows = NA_integer_,
      missing_columns = paste(missing_columns, collapse = ";"),
      stringsAsFactors = FALSE
    )

    data_list[[i]] <- NULL
    next
  }

  gene <- as.character(dt$gene)
  zscore <- suppressWarnings(as.numeric(dt$zscore))
  pvalue <- suppressWarnings(as.numeric(dt$pvalue))

  invalid_gene <- is.na(gene) | !nzchar(gene)
  invalid_zscore <- !is.finite(zscore)
  invalid_pvalue <- !is.finite(pvalue)

  valid_rows <- (
    !invalid_gene &
    !invalid_zscore &
    !invalid_pvalue
  )

  clean_dt <- dt[valid_rows, , drop = FALSE]
  clean_dt$gene <- gene[valid_rows]
  clean_dt$zscore <- zscore[valid_rows]
  clean_dt$pvalue <- pvalue[valid_rows]

  file_summary[[i]] <- data.frame(
    trait = trait_names[i],
    file = filenames[i],
    raw_rows = nrow(dt),
    valid_rows = nrow(clean_dt),
    unique_genes = length(unique(clean_dt$gene)),
    duplicated_gene_rows = sum(duplicated(clean_dt$gene)),
    invalid_gene_rows = sum(invalid_gene),
    invalid_zscore_rows = sum(invalid_zscore),
    invalid_pvalue_rows = sum(invalid_pvalue),
    missing_columns = "",
    stringsAsFactors = FALSE
  )

  data_list[[i]] <- clean_dt
}

file_summary <- rbindlist(
  file_summary,
  fill = TRUE
)

cat("File-level summary:\n")
print(file_summary)

summary_output <- file.path(
  getwd(),
  paste0(
    TARGET,
    ".",
    brain_region,
    ".file_summary.csv"
  )
)

fwrite(file_summary, summary_output)

if (any(vapply(data_list, is.null, logical(1)))) {
  stop(
    "At least one file is missing required columns. See: ",
    summary_output
  )
}

target_data <- data_list[[target_idx]]
trait_indices <- setdiff(seq_along(data_list), target_idx)

pairwise_summary <- lapply(trait_indices, function(i) {
  target_genes <- target_data$gene
  trait_genes <- data_list[[i]]$gene
  common_genes <- intersect(target_genes, trait_genes)

  data.frame(
    trait = trait_names[i],
    target_genes = length(target_genes),
    trait_genes = length(trait_genes),
    common_genes = length(common_genes),
    missing_from_trait = length(setdiff(target_genes, trait_genes)),
    extra_in_trait = length(setdiff(trait_genes, target_genes)),
    exactly_identical = identical(target_genes, trait_genes),
    stringsAsFactors = FALSE
  )
})

pairwise_summary <- rbindlist(pairwise_summary)

cat("\nPairwise comparison against target:\n")
print(pairwise_summary)

pairwise_output <- file.path(
  getwd(),
  paste0(
    TARGET,
    ".",
    brain_region,
    ".pairwise_summary.csv"
  )
)

fwrite(pairwise_summary, pairwise_output)

all_gene_lists <- lapply(
  data_list,
  function(dt) unique(dt$gene)
)

global_common_genes <- Reduce(
  intersect,
  all_gene_lists
)

cat("\n============================================================\n")
cat("Global intersection results\n")
cat("============================================================\n")
cat("Target valid genes:", nrow(target_data), "\n")
cat(
  "Common genes across all",
  length(data_list),
  "files:",
  length(global_common_genes),
  "\n"
)

cat("\n============================================================\n")
cat("Simulation of the current alignment code\n")
cat("============================================================\n")

simulated_target <- target_data
simulated_traits <- data_list[trait_indices]
simulated_trait_names <- trait_names[trait_indices]

target_genes <- simulated_target$gene
first_mismatch <- NA_integer_

for (i in seq_along(simulated_traits)) {
  if (!identical(target_genes, simulated_traits[[i]]$gene)) {
    first_mismatch <- i

    common_genes <- intersect(
      target_genes,
      simulated_traits[[i]]$gene
    )

    simulated_target <- simulated_target[
      simulated_target$gene %in% common_genes,
      ,
      drop = FALSE
    ]

    for (j in seq_along(simulated_traits)) {
      simulated_traits[[j]] <- simulated_traits[[j]][
        simulated_traits[[j]]$gene %in% common_genes,
        ,
        drop = FALSE
      ]
    }

    break
  }
}

if (is.na(first_mismatch)) {
  cat("No mismatch was detected by the current code.\n")
} else {
  cat(
    "First mismatched trait:",
    simulated_trait_names[first_mismatch],
    "\n"
  )
}

simulated_lengths <- vapply(
  simulated_traits,
  nrow,
  integer(1)
)

simulated_summary <- data.frame(
  trait = simulated_trait_names,
  rows_after_current_alignment = simulated_lengths,
  stringsAsFactors = FALSE
)

cat("\nRows after simulating the current alignment logic:\n")
print(simulated_summary)

cat(
  "\nNumber of distinct trait lengths after current alignment:",
  length(unique(simulated_lengths)),
  "\n"
)

simulated_features <- sapply(
  simulated_traits,
  function(x) abs(x$zscore)
)

cat("\nResult returned by sapply:\n")
cat("Class:", paste(class(simulated_features), collapse = ", "), "\n")
cat("Type:", typeof(simulated_features), "\n")
cat(
  "Dimensions:",
  if (is.null(dim(simulated_features))) {
    "NULL"
  } else {
    paste(dim(simulated_features), collapse = " x ")
  },
  "\n"
)

cat("\nDiagnostic outputs:\n")
cat("File summary:", summary_output, "\n")
cat("Pairwise summary:", pairwise_output, "\n")

if (is.null(dim(simulated_features))) {
  cat(
    "\nConclusion: sapply returned a non-matrix object. ",
    "This directly explains the colnames error.\n",
    sep = ""
  )
}

if (length(unique(simulated_lengths)) > 1) {
  cat(
    "Conclusion: trait files still have different row counts ",
    "after the current alignment logic.\n",
    sep = ""
  )
}

if (length(global_common_genes) < 2) {
  cat(
    "Conclusion: fewer than two genes are shared across all files. ",
    "The upstream MetaXcan files require investigation.\n",
    sep = ""
  )
}
