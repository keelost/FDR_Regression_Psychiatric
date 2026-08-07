#!/usr/bin/env Rscript
#*********************************************************#
# Ablation sensitivity analysis - MAGMA FDRreg
#   SA1: 13 curated bio annotations only (no disorders)
#   SA2: all disorder traits in folder + 13 curated bio
#   Theoretical null only. No empirical / contribution / lasso.
#*********************************************************#
# R-4.0.2 or later

rm(list = ls())

suppressPackageStartupMessages({
  library(FDRreg)
  library(data.table)
  library(HelpersMG)
  library(optparse)
})

# ---- Curated bio annotation columns for SA2 (13-column subset) ---- #
BIO_KEEP_COLS <- c("scz.drug", "scz.literature", "scz.animal",
                   "depre.drug", "depre.literature",
                   "bpd.drug", "bpd.literature",
                   "adhd.drug", "adhd.literature", "adhd.animal",
                   "asd.drug", "asd.literature", "asd.animal")

# Non-annotation key columns in the entrez bio file (excluded from "all bio")
BIO_KEY_COLS <- c("ID", "GeneName")

# ---- Argument parsing ---- #
option_list <- list(
  make_option(c("--base_dir"), type = "character",
              default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
              help = "Base directory holding per-target folders. [default: %default]"),
  make_option(c("--out_dir"), type = "character",
              default = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "10.ablation", "magma"),
              help = "Ablation output directory. [default: %default]"),
  make_option(c("--bio_file"), type = "character",
              default = Sys.getenv("FDRREG_BIO_ENTREZ", ""),
              help = "Entrez-keyed biological annotation file. [default: %default]"),
  make_option(c("--targets"), type = "character", default = NULL,
              help = "Optional comma-separated target list. Default: auto-detect. [default: NULL]"),
  make_option(c("--seed"), type = "integer", default = 100,
              help = "Random seed. [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Helper: run theoretical FDRreg with zero-variance guard ---- #
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

# ---- Helper: build bio matrix aligned to entrez gene keys ---- #
#   cols = character vector of annotation columns to keep
build_bio_matrix <- function(bio_dt, keys, cols) {
  keep <- intersect(cols, names(bio_dt))
  idx <- match(as.character(keys), as.character(bio_dt$ID))
  mat <- as.matrix(bio_dt[idx, ..keep])
  mat[is.na(mat)] <- 0
  storage.mode(mat) <- "double"
  mat
}

# ---- Auto-detect targets ---- #
detect_targets <- function(base_dir) {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  dirs[vapply(dirs, function(d)
    dir.exists(file.path(base_dir, d, "04.magma_output")), logical(1))]
}

if (!is.null(opt$targets)) {
  targets <- trimws(strsplit(opt$targets, ",")[[1]])
  targets <- targets[targets != ""]
} else {
  targets <- detect_targets(opt$base_dir)
}
cat(sprintf("MAGMA ablation: %d target(s): %s\n",
            length(targets), paste(targets, collapse = ", ")))

# ---- Load bio library once ---- #
bio_dt <- fread(opt$bio_file, showProgress = FALSE)

# ---- Process each target ---- #
for (target in targets) {
  cat(sprintf("\n========== MAGMA target: %s ==========\n", target))

  magma_dir <- file.path(opt$base_dir, target, "04.magma_output")
  target_file <- file.path(magma_dir, paste0(target, ".genes.out"))
  if (!file.exists(target_file)) {
    cat(sprintf("  Target file not found, skipping: %s\n", target_file))
    next
  }

  # All .genes.out files; variables = everything except target
  all_files <- list.files(magma_dir, pattern = "\\.genes\\.out$", full.names = TRUE)
  trait_names_all <- sub("\\.genes\\.out$", "", basename(all_files))
  var_files <- all_files[trait_names_all != target]
  var_names <- trait_names_all[trait_names_all != target]

  if (length(var_files) == 0) {
    cat("  No disorder trait files found; SA2 will be skipped.\n")
  }

  # Read target and variables
  target_gene <- fread(target_file, showProgress = FALSE)
  var_data <- lapply(var_files, fread, showProgress = FALSE)

  # Align on common genes across target + all variables
  gene_sets <- c(list(target_gene$GENE), lapply(var_data, `[[`, "GENE"))
  common_genes <- Reduce(intersect, gene_sets)
  if (length(common_genes) == 0) {
    cat("  No common genes across files, skipping.\n")
    next
  }
  common_genes <- sort(common_genes)

  target_gene <- target_gene[match(common_genes, target_gene$GENE), ]

  # Target z: MAGMA has no direction, assign random sign (seed per target)
  set.seed(opt$seed)
  target_z <- abs(qnorm(target_gene$P / 2)) * sign(rnorm(nrow(target_gene)))

  # Disorder covariate matrix (absolute z), aligned to common genes
  if (length(var_files) > 0) {
    z_disorder <- sapply(var_data, function(d) {
      d <- d[match(common_genes, d$GENE), ]
      abs(qnorm(d$P / 2))
    })
    z_disorder <- as.matrix(z_disorder)
    colnames(z_disorder) <- var_names
  } else {
    z_disorder <- NULL
  }

  # Bio matrices aligned to common genes (entrez key = GENE)
  bio_all_cols <- setdiff(names(bio_dt), BIO_KEY_COLS)   # SA1: all bio annotations
  bio_mat_all  <- build_bio_matrix(bio_dt, common_genes, bio_all_cols)
  bio_mat_sub  <- build_bio_matrix(bio_dt, common_genes, BIO_KEEP_COLS)  # SA2: 13-col subset

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

  # ---- Assemble and save ---- #
  res <- data.frame(
    gene_id = common_genes,
    P = target_gene$P,
    qval = p.adjust(target_gene$P, method = "fdr"),
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

cat("\n========== MAGMA ablation complete ==========\n")
quit(save = "no", status = 0, runLast = FALSE)
