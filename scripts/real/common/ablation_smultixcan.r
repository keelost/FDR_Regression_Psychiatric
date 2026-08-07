#!/usr/bin/env Rscript
#*********************************************************#
# Ablation sensitivity analysis - S-MultiXcan FDRreg
#   SA1: 13 curated bio annotations only (no disorders)
#   SA2: all disorder traits in folder + 13 curated bio
#   Theoretical null only. No empirical / contribution / lasso.
#*********************************************************#
# R-4.0.2 or later

rm(list = ls())

suppressPackageStartupMessages({
  library(FDRreg)
  library(data.table)
  library(stringr)
  library(HelpersMG)
  library(optparse)
})

BIO_KEEP_COLS <- c("scz.drug", "scz.literature", "scz.animal",
                   "depre.drug", "depre.literature",
                   "bpd.drug", "bpd.literature",
                   "adhd.drug", "adhd.literature", "adhd.animal",
                   "asd.drug", "asd.literature", "asd.animal")

# Non-annotation key column in the ensembl bio file (excluded from "all bio")
BIO_KEY_COLS <- c("ENSEMBL_GENE_ID")

option_list <- list(
  make_option(c("--base_dir"), type = "character",
              default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
              help = "Base directory holding per-target folders. [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "10.ablation", "smultixcan"),
              help = "Ablation output directory. [default: %default]"),
  make_option(c("--bio_file"), type = "character",
              default = Sys.getenv("FDRREG_BIO_ENSEMBL", ""),
              help = "Ensembl-keyed biological annotation file. [default: %default]"),
  make_option(c("--targets"), type = "character", default = NULL,
              help = "Optional comma-separated target list. Default: auto-detect. [default: NULL]"),
  make_option(c("--seed"), type = "integer", default = 100,
              help = "Random seed. [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Safe p-value to z-score conversion (mirrors main pipeline)
p_to_z_safe <- function(p, cap = 8.2) {
  p_clamped <- pmax(pmin(p, 1 - .Machine$double.eps), .Machine$double.eps)
  z <- qnorm(1 - p_clamped / 2)
  pmin(pmax(z, -cap), cap)
}

run_fdrreg_the <- function(z, covars, label) {
  covars <- as.matrix(covars)
  keep <- apply(covars, 2, function(x) length(unique(x[!is.na(x)])) > 1)
  if (any(!keep)) {
    cat(sprintf("    [%s] dropping %d zero-variance covariate(s)\n", label, sum(!keep)))
  }
  covars <- covars[, keep, drop = FALSE]
  if (ncol(covars) == 0) {
    cat(sprintf("    [%s] no usable covariates, skipping\n", label))
    return(NULL)
  }
  tryCatch(
    FDRreg(z, covars, nulltype = "theoretical", method = "pr"),
    error = function(e) {
      cat(sprintf("    [%s] FDRreg failed: %s\n", label, conditionMessage(e)))
      NULL
    }
  )
}

build_bio_matrix <- function(bio_dt, keys, cols) {
  keep <- intersect(cols, names(bio_dt))
  idx <- match(as.character(keys), as.character(bio_dt$ENSEMBL_GENE_ID))
  mat <- as.matrix(bio_dt[idx, ..keep])
  mat[is.na(mat)] <- 0
  storage.mode(mat) <- "double"
  mat
}

detect_targets <- function(base_dir) {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  dirs[vapply(dirs, function(d)
    dir.exists(file.path(base_dir, d, "08.smultixcan_output")), logical(1))]
}

if (!is.null(opt$targets)) {
  targets <- trimws(strsplit(opt$targets, ",")[[1]])
  targets <- targets[targets != ""]
} else {
  targets <- detect_targets(opt$base_dir)
}
cat(sprintf("S-MultiXcan ablation: %d target(s): %s\n",
            length(targets), paste(targets, collapse = ", ")))

bio_dt <- fread(opt$bio_file, showProgress = FALSE)

for (target in targets) {
  cat(sprintf("\n========== S-MultiXcan target: %s ==========\n", target))
  input_dir <- file.path(opt$base_dir, target, "08.smultixcan_output")

  all_files <- list.files(input_dir, pattern = "\\.allbrain\\.txt$", full.names = TRUE)
  if (length(all_files) == 0) {
    cat("  No .allbrain.txt files found, skipping.\n")
    next
  }
  target_basename <- paste0(target, ".allbrain.txt")
  target_idx <- which(basename(all_files) == target_basename)
  if (length(target_idx) == 0) {
    cat(sprintf("  Target file not found: %s\n", target_basename))
    next
  }
  var_idx <- setdiff(seq_along(all_files), target_idx)
  var_names <- str_remove(basename(all_files[var_idx]), "\\.allbrain\\.txt$")

  read_one <- function(x) {
    dt <- fread(x, showProgress = FALSE)
    cols <- intersect(c("gene", "gene_name", "pvalue"), names(dt))
    dt <- dt[, ..cols]
    dt <- dt[complete.cases(pvalue)]
    dt <- dt[order(gene)]
    dt[, zscore := p_to_z_safe(pvalue)]
    dt
  }

  target_data <- read_one(all_files[target_idx])
  var_data <- if (length(var_idx) > 0) lapply(all_files[var_idx], read_one) else list()

  # Align on common genes
  gene_sets <- c(list(target_data$gene), lapply(var_data, `[[`, "gene"))
  common_genes <- Reduce(intersect, gene_sets)
  if (length(common_genes) == 0) {
    cat("  No common genes across files, skipping.\n")
    next
  }
  common_genes <- sort(common_genes)

  target_data <- target_data[match(common_genes, gene)]

  # Target z: only p-value available, assign random sign (seed per target)
  set.seed(opt$seed)
  target_z <- target_data$zscore * sign(rnorm(nrow(target_data)))

  # Disorder covariates (absolute z)
  if (length(var_data) > 0) {
    z_disorder <- sapply(var_data, function(d) {
      d <- d[match(common_genes, gene)]
      abs(d$zscore)
    })
    z_disorder <- as.matrix(z_disorder)
    colnames(z_disorder) <- var_names
  } else {
    z_disorder <- NULL
  }

  # Bio matrices (ensembl key, version stripped)
  gene_ids_clean <- str_split_fixed(common_genes, "\\.", 2)[, 1]
  bio_all_cols <- setdiff(names(bio_dt), BIO_KEY_COLS)   # SA1: all bio annotations
  bio_mat_all  <- build_bio_matrix(bio_dt, gene_ids_clean, bio_all_cols)
  bio_mat_sub  <- build_bio_matrix(bio_dt, gene_ids_clean, BIO_KEEP_COLS)  # SA2 subset

  # ---- SA1: all bio annotations only ---- #
  cat("  SA1: all bio annotations only...\n")
  fdr_sa1 <- run_fdrreg_the(target_z, bio_mat_all, "SA1")

  # ---- SA2: disorders + curated bio subset ---- #
  if (!is.null(z_disorder)) {
    cat("  SA2: selected disorders + curated bio subset...\n")
    covars_sa2 <- cbind(abs(z_disorder), bio_mat_sub)
    fdr_sa2 <- run_fdrreg_the(target_z, covars_sa2, "SA2")
  } else {
    fdr_sa2 <- NULL
  }

  # ---- Save ---- #
  res <- data.frame(
    gene_id = gene_ids_clean,
    gene_name = target_data$gene_name,
    pvalue = target_data$pvalue,
    qval = p.adjust(target_data$pvalue, method = "fdr"),
    SA1_FDR_the = if (!is.null(fdr_sa1)) fdr_sa1$FDR else NA_real_,
    SA2_FDR_the = if (!is.null(fdr_sa2)) fdr_sa2$FDR else NA_real_,
    stringsAsFactors = FALSE
  )

  target_out_dir <- file.path(opt$out_dir, target)
  dir.create(target_out_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(res, file.path(target_out_dir, paste0(target, ".ablation.fdrreg.txt")), sep = " ")
  cat(sprintf("  Saved: %s\n",
              file.path(target_out_dir, paste0(target, ".ablation.fdrreg.txt"))))
}

cat("\n========== S-MultiXcan ablation complete ==========\n")
quit(save = "no", status = 0, runLast = FALSE)
