#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function(args) {
  if (any(args %in% c("-h", "--help"))) {
    cat("Usage: prepare_dist_input.R --input FILE --output FILE [--clear-info-threshold 0.6]\n")
    quit(status = 0)
  }
  opt <- list(clear_info_threshold = 0.6)
  i <- 1L
  while (i <= length(args)) {
    if (i == length(args) || !startsWith(args[[i]], "--")) stop("Invalid command-line arguments")
    key <- gsub("-", "_", substring(args[[i]], 3), fixed = TRUE)
    opt[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opt$clear_info_threshold <- as.numeric(opt$clear_info_threshold)
  opt
}
opt <- parse_cli(commandArgs(trailingOnly = TRUE))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("--input and --output are required")
}
if (!file.exists(opt$input)) stop("Input does not exist: ", opt$input)

gwas <- fread(opt$input, fill = TRUE)
setnames(gwas, names(gwas), tolower(names(gwas)))

aliases <- list(
  snpid = c("snpid", "snp", "rsid", "markername"),
  chr = c("chr", "chrom", "chromosome"),
  bpos = c("bpos", "bp", "pos", "position"),
  a1 = c("a1", "effect_allele", "ea"),
  a2 = c("a2", "other_allele", "nea")
)
for (canonical in names(aliases)) {
  found <- aliases[[canonical]][aliases[[canonical]] %in% names(gwas)]
  if (!canonical %in% names(gwas) && length(found)) {
    setnames(gwas, found[[1]], canonical)
  }
}

required <- c("snpid", "chr", "bpos", "a1", "a2")
missing <- setdiff(required, names(gwas))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

if (!"z" %in% names(gwas)) {
  if (!all(c("beta", "se") %in% names(gwas))) {
    stop("Input must contain z, or both beta and se")
  }
  gwas[, z := beta / se]
}

n_input <- nrow(gwas)
if ("info" %in% names(gwas)) {
  gwas <- gwas[is.na(info) | info >= opt$clear_info_threshold]
}

gwas[, `:=`(
  snpid = as.character(snpid),
  chr = sub("^chr", "", as.character(chr), ignore.case = TRUE),
  bpos = as.integer(bpos),
  a1 = toupper(as.character(a1)),
  a2 = toupper(as.character(a2)),
  z = as.numeric(z)
)]

# DIST expects biallelic SNPs. The cleaned inputs should already satisfy this,
# but enforcing it here makes the upstream contract explicit and auditable.
gwas <- gwas[
  !is.na(snpid) & snpid != "" &
    chr %in% as.character(1:22) & !is.na(bpos) & is.finite(z) &
    a1 %chin% c("A", "C", "G", "T") & a2 %chin% c("A", "C", "G", "T") &
    a1 != a2
]
gwas <- unique(gwas, by = "snpid")
gwas[, chr_order__ := as.integer(chr)]
setorder(gwas, chr_order__, bpos)
gwas[, chr_order__ := NULL]

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(gwas[, .(snpid, chr, bpos, a1, a2, z)], opt$output, sep = " ", na = "NA")
cat(sprintf("Prepared %s: %d of %d rows retained\n", opt$output, nrow(gwas), n_input))
