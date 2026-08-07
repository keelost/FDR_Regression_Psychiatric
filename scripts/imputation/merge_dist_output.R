#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function(args) {
  if (any(args %in% c("-h", "--help"))) {
    cat(paste0(
      "Usage: merge_dist_output.R --trait NAME --input-dir DIR ",
      "--variant-map FILE --output FILE [--imputed-info-threshold 0.5]\n"
    ))
    quit(status = 0)
  }
  opt <- list(imputed_info_threshold = 0.5)
  i <- 1L
  while (i <= length(args)) {
    if (i == length(args) || !startsWith(args[[i]], "--")) stop("Invalid command-line arguments")
    key <- gsub("-", "_", substring(args[[i]], 3), fixed = TRUE)
    opt[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opt$imputed_info_threshold <- as.numeric(opt$imputed_info_threshold)
  opt
}
opt <- parse_cli(commandArgs(trailingOnly = TRUE))

required_args <- c("trait", "input_dir", "variant_map", "output")
missing_args <- required_args[vapply(required_args, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_args)) stop("Missing required options: ", paste0("--", gsub("_", "-", missing_args), collapse = ", "))
if (!file.exists(opt$variant_map)) stop("Variant map does not exist: ", opt$variant_map)

files <- file.path(opt$input_dir, sprintf("%s.chr%d.imputation.txt", opt$trait, 1:22))
missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop("Missing DIST chromosome output(s): ", paste(basename(missing_files), collapse = ", "))
}

read_chr <- function(path) {
  x <- fread(path, fill = TRUE)
  if (!"snpid" %in% names(x)) setnames(x, names(x)[1], "snpid")
  setnames(x, names(x), tolower(names(x)))
  if ("bpos" %in% names(x) && !"bp" %in% names(x)) setnames(x, "bpos", "bp")
  needed <- c("snpid", "chr", "bp", "a1", "a2", "z", "info", "pval")
  absent <- setdiff(needed, names(x))
  if (length(absent)) stop(basename(path), " is missing: ", paste(absent, collapse = ", "))
  x
}

imputed <- rbindlist(lapply(files, read_chr), fill = TRUE, use.names = TRUE)
imputed <- imputed[!is.na(info) & info >= opt$imputed_info_threshold]
imputed[, `:=`(
  snpid = as.character(snpid),
  chr = sub("^chr", "", as.character(chr), ignore.case = TRUE),
  bp = as.integer(bp)
)]

variant_map <- fread(opt$variant_map, fill = TRUE)
setnames(variant_map, names(variant_map), tolower(names(variant_map)))
pick <- function(candidates) {
  hit <- candidates[candidates %in% names(variant_map)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}
chr_col <- pick(c("chr", "chrom", "chromosome"))
pos_col <- pick(c("pos", "bp", "bpos", "position"))
rsid_col <- pick(c("rsid", "snpid", "snp"))
if (anyNA(c(chr_col, pos_col, rsid_col))) {
  stop("Variant map must contain chromosome, position, and rsID columns")
}
variant_map <- unique(variant_map[, .(
  chr = sub("^chr", "", as.character(get(chr_col)), ignore.case = TRUE),
  bp = as.integer(get(pos_col)),
  mapped_rsid = as.character(get(rsid_col))
)], by = c("chr", "bp"))

has_rsid <- grepl("^rs[0-9]+$", imputed$snpid, ignore.case = TRUE)
with_id <- imputed[has_rsid]
without_id <- imputed[!has_rsid]
if (nrow(without_id)) {
  without_id <- merge(without_id, variant_map, by = c("chr", "bp"), all = FALSE)
  without_id[, snpid := mapped_rsid]
  without_id[, mapped_rsid := NULL]
}

result <- rbindlist(list(with_id, without_id), fill = TRUE, use.names = TRUE)
result <- result[!is.na(snpid) & snpid != ""]
result <- unique(result, by = "snpid")
result[, chr_order__ := as.integer(chr)]
setorder(result, chr_order__, bp)
result[, chr_order__ := NULL]

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(result, opt$output, sep = " ", na = "NA")
cat(sprintf("Merged %s: %d variants retained at INFO >= %.3g\n", opt$output, nrow(result), opt$imputed_info_threshold))
