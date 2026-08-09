# Project: FDRreg Analysis Pipeline - Post-hoc SNP Audit
# File: audit_snp_overlap.R
# Description: For every target already processed by the FDRreg pipeline:
#   (1) Report the SNP count taken directly from the generated result files.
#   (2) Recompute the overlapping-SNP count when the input sources are
#       switched to the "targets" (target SNPs) and "library" (trait SNPs,
#       stored as .txt.gz) directories.
#   (3) Recover the trait list for each target directly from the finished
#       FDRreg outputs (model-info RDS, with MAGMA filenames as fallback).
#
#   IMPORTANT (v1.1): The recomputed overlap now mirrors the updated
#   FDRreg_pipeline.R, which performs allele harmonisation BEFORE overlap
#   identification. Both the target and every trait are passed through the
#   same harmonisation filter (allele-orientation classification and, when
#   enabled, dropping of strand-ambiguous palindromic SNPs). SNPs whose
#   alleles are incompatible with the target - or palindromic when
#   DROP_PALINDROMIC is TRUE - are removed before the intersection is taken.
#   This keeps the audit's "new overlap" count comparable to what the
#   pipeline actually feeds into FDRreg.
#
# Note: This script only inspects/derives; it does NOT re-run FDRreg.
# Version: 1.1

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Transparent reader for optionally-gzipped files. Some data.table versions do
# not decompress .gz natively (they require the R.utils package); piping the
# file through "gzip -dc" removes that dependency and works everywhere.
fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(fread(cmd = paste("gzip -dc", shQuote(path)), ...))
  }
  fread(path, ...)
}

#-------------------------------------------------#
#---- Harmonisation configuration (mirror FDR) ---#
#-------------------------------------------------#
# These MUST match FDRreg_pipeline.R so the recomputed overlap is comparable.
# In the targets/library sources the effect/other alleles are stored as a1/a2.
TARGET_EA <- "a1"; TARGET_OA <- "a2"
TRAIT_EA  <- "a1"; TRAIT_OA  <- "a2"

# Drop strand-ambiguous palindromic SNPs (A/T, C/G) when TRUE. Keep in sync
# with the DROP_PALINDROMIC value used when the pipeline was run.
DROP_PALINDROMIC <- TRUE

#-------------------------------------------------#
#---- Harmonisation helpers (mirror FDR logic) ---#
#-------------------------------------------------#

# Reverse-complement of allele strings (vectorized).
complement_allele <- function(a) chartr("ACGT", "TGCA", a)

# Strand-ambiguous palindromic SNPs: A/T or C/G pairs.
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
  (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# Classify allele orientation between target (e1/e2) and trait (o1/o2).
# Returns a numeric sign vector:
#    1  -> alleles aligned (direct or strand-flipped)
#   -1  -> alleles swapped (swap or strand-flipped + swap)
#   NA  -> incompatible, or palindromic when DROP_PALINDROMIC is TRUE (dropped)
# Only the "is it usable" information (non-NA) is needed for the audit count.
classify_sign <- function(e1, e2, o1, o2) {
  co1 <- complement_allele(o1)
  co2 <- complement_allele(o2)

  aligned <- (o1 == e1 & o2 == e2) | (co1 == e1 & co2 == e2)
  swapped <- (o1 == e2 & o2 == e1) | (co1 == e2 & co2 == e1)

  s <- rep(NA_real_, length(e1))
  s[aligned] <- 1
  s[swapped] <- -1

  if (DROP_PALINDROMIC) {
    s[is_palindromic(o1, o2)] <- NA_real_
  }
  s
}

#-------------------------------#
#---- Fixed Path Parameters ----#
#-------------------------------#

# Directory that holds one sub-directory per already-processed target.
results_root <- "/path/to/SO_Lab/18.fdrreg_rebuild"

# New input sources for requirement (2).
targets_dir <- "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets"   # target SNPs
library_dir <- "/exeh_3/rstao/dr.so/002.meta/002.clear.data/library"   # trait SNPs (.txt.gz)

# Where to write the audit summary (dedicated tmp/summary directory).
audit_dir <- "/path/to/SO_Lab/18.fdrreg_rebuild/01.extra.analysis/07.dist"
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)
audit_out <- file.path(audit_dir, "snp_overlap_audit.csv")

# Number of cores for trait-file reads.
n_cores <- max(1, min(detectCores() - 1, 4))

#-------------------------------#
#---- Small Helper Functions ---#
#-------------------------------#

# Locate a data file for a given base name inside a directory, probing a set
# of plausible suffixes (gzipped variants included). Returns NA if not found.
# Actual layout:
#   library/<name>.clear.txt.gz   e.g. library/alco2018.clear.txt.gz
#   targets/<name>.clear.txt      e.g. targets/ptsd2017.clear.txt
find_data_file <- function(dir_path, base_name) {
  if (!dir.exists(dir_path)) return(NA_character_)
  suffixes <- c(".clear.txt.gz", ".clear.txt",   # primary layout
                ".txt.gz", ".txt",
                ".impute.map.txt.gz", ".impute.map.txt",
                ".csv.gz", ".csv", ".gz")
  for (sfx in suffixes) {
    cand <- file.path(dir_path, paste0(base_name, sfx))
    if (file.exists(cand)) return(cand)
  }
  NA_character_
}

# Generic case-insensitive column resolver against a header vector.
resolve_col <- function(hdr, candidates) {
  hit <- hdr[tolower(hdr) %in% tolower(candidates)]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

# Detect the SNP-ID column from a header vector. Falls back to the first column.
detect_snpid_from_hdr <- function(hdr) {
  if (length(hdr) == 0) return(NA_character_)
  candidates <- c("snpid", "snp", "rsid", "rs_id", "rsids",
                  "markername", "marker", "id", "variant_id", "variant")
  hit <- resolve_col(hdr, candidates)
  if (!is.na(hit)) return(hit)
  hdr[1]  # fallback: assume the first column holds the identifier
}

# Read snpid + effect/other alleles from a file and normalise column names to
# snpid / a1 / a2. Alleles are upper-cased so that classification is robust to
# case differences between the targets and library directories. Returns NULL
# when the file is missing/unreadable or the required columns are absent.
read_snp_alleles <- function(path, ea = TRAIT_EA, oa = TRAIT_OA) {
  if (is.na(path) || !file.exists(path)) return(NULL)

  hdr <- tryCatch(names(fread_any(path, nrows = 0L)), error = function(e) character(0))
  if (length(hdr) == 0) return(NULL)

  snp_col <- detect_snpid_from_hdr(hdr)
  ea_col  <- resolve_col(hdr, ea)
  oa_col  <- resolve_col(hdr, oa)
  if (is.na(snp_col) || is.na(ea_col) || is.na(oa_col)) return(NULL)

  dt <- tryCatch(fread_any(path, select = c(snp_col, ea_col, oa_col)),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  setnames(dt, c(snp_col, ea_col, oa_col), c("snpid", "a1", "a2"))
  dt[, snpid := as.character(snpid)]
  dt[, a1 := toupper(as.character(a1))]
  dt[, a2 := toupper(as.character(a2))]

  # Deduplicate on snpid to keep the subsequent merge well-defined.
  unique(dt, by = "snpid")
}

# Read target alleles and expose them as snpid / e1 / e2 (reference alleles).
read_target_alleles <- function(path) {
  dt <- read_snp_alleles(path, TARGET_EA, TARGET_OA)
  if (is.null(dt)) return(NULL)
  setnames(dt, c("a1", "a2"), c("e1", "e2"))
  dt
}

# Count rows in a result file cheaply (row count = SNP count for per-SNP files).
count_rows <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_integer_)
  dt <- tryCatch(fread(path, select = 1L), error = function(e) NULL)
  if (is.null(dt)) return(NA_integer_)
  nrow(dt)
}

#-------------------------------------------------#
#---- Requirement 3: recover traits per target ---#
#-------------------------------------------------#

# Try the model-info RDS first; fall back to MAGMA-input filenames.
get_traits_for_target <- function(target_dir, target_name) {
  fdr_dir   <- file.path(target_dir, "02.fdrreg_results")
  magma_dir <- file.path(target_dir, "01.magma_input")

  # (a) Preferred source: fdrreg_model_info_<mode>.rds
  rds_files <- list.files(fdr_dir, pattern = "^fdrreg_model_info_.*\\.rds$",
                          full.names = TRUE)
  if (length(rds_files) > 0) {
    info <- tryCatch(readRDS(rds_files[1]), error = function(e) NULL)
    if (!is.null(info)) {
      traits <- unique(c(info$traits_with_overlap, info$traits_without_overlap))
      traits <- traits[nchar(traits) > 0]
      if (length(traits) > 0) {
        return(list(traits = traits, source = "rds"))
      }
    }
  }

  # (b) Fallback: derive from MAGMA-input filenames minus the target itself.
  magma_files <- list.files(magma_dir, pattern = "\\.overlap\\.4magma\\.txt$")
  if (length(magma_files) > 0) {
    names_only <- sub("\\.overlap\\.4magma\\.txt$", "", magma_files)
    traits <- setdiff(names_only, target_name)
    if (length(traits) > 0) {
      return(list(traits = traits, source = "magma_filenames"))
    }
  }

  list(traits = character(0), source = "none")
}

#-------------------------------------------------#
#---- Requirement 1: SNP count from results ------#
#-------------------------------------------------#

get_result_snp_count <- function(target_dir, target_name) {
  fdr_dir     <- file.path(target_dir, "02.fdrreg_results")
  overlap_dir <- file.path(target_dir, "00.overlap_data")

  # (a) Preferred: per-SNP FDR result file (one row per analyzed SNP).
  per_snp <- list.files(fdr_dir, pattern = "^fdr_values_per_snp_.*\\.csv$",
                        full.names = TRUE)
  if (length(per_snp) > 0) {
    n <- count_rows(per_snp[1])
    if (!is.na(n)) return(list(count = n, source = basename(per_snp[1])))
  }

  # (b) Fallback: overlapping-SNP list.
  snp_list <- file.path(overlap_dir, "overlapping_snps_list.csv")
  if (file.exists(snp_list)) {
    n <- count_rows(snp_list)
    if (!is.na(n)) return(list(count = n, source = "overlapping_snps_list.csv"))
  }

  # (c) Fallback: target overlap CSV.
  tgt_overlap <- file.path(overlap_dir, paste0(target_name, ".overlap.csv"))
  if (file.exists(tgt_overlap)) {
    n <- count_rows(tgt_overlap)
    if (!is.na(n)) return(list(count = n, source = basename(tgt_overlap)))
  }

  list(count = NA_integer_, source = "none")
}

#-------------------------------------------------#
#---- Requirement 2: recompute overlap (new src) -#
#---- with allele harmonisation (mirror FDR) -----#
#-------------------------------------------------#

recompute_overlap <- function(target_name, traits) {
  missing <- character(0)

  # Target SNPs + alleles from the "targets" directory (reference alleles).
  tgt_file       <- find_data_file(targets_dir, target_name)
  target_alleles <- read_target_alleles(tgt_file)
  if (is.null(target_alleles) || nrow(target_alleles) == 0) {
    return(list(count = NA_integer_,
                missing = paste0("target:", target_name)))
  }
  setkey(target_alleles, snpid)

  # Trait files from the "library" directory (.txt.gz).
  trait_files <- vapply(traits, function(tr) find_data_file(library_dir, tr),
                        character(1))

  # Per-trait harmonisation pass (Pass A equivalent): keep only SNPs whose
  # alleles are compatible with the target after orientation classification.
  #   NULL         -> file missing / unreadable (flagged, excluded)
  #   character(0) -> file present but no shared or no usable SNP
  #   character(>0)-> usable SNP ids
  kept_list <- mclapply(seq_along(traits), function(k) {
    f <- trait_files[k]
    if (is.na(f) || !file.exists(f)) return(NULL)

    ta <- read_snp_alleles(f, TRAIT_EA, TRAIT_OA)
    if (is.null(ta) || nrow(ta) == 0) return(NULL)

    m <- merge(target_alleles, ta, by = "snpid")
    if (nrow(m) == 0L) return(character(0))

    s <- classify_sign(m$e1, m$e2, m[["a1"]], m[["a2"]])
    m$snpid[!is.na(s)]
  }, mc.cores = min(length(traits), n_cores))
  names(kept_list) <- traits

  # Flag traits whose source file could not be read.
  is_missing <- vapply(kept_list, is.null, logical(1))
  if (any(is_missing)) {
    missing <- c(missing, paste0("library:", traits[is_missing]))
  }

  # Intersect over all readable traits. A readable trait with zero usable SNPs
  # correctly forces the overlap to zero (it is a genuine harmonisation result,
  # not a missing input).
  present <- kept_list[!is_missing]
  if (length(present) == 0) {
    return(list(count = NA_integer_, missing = missing))
  }

  overlap_ids <- Reduce(intersect, c(list(target_alleles$snpid), present))
  list(count = length(overlap_ids), missing = missing)
}

#-------------------------------#
#---- Main Audit Loop ----------#
#-------------------------------#

cat("============================================================\n")
cat("FDRreg post-hoc SNP overlap audit (harmonisation-aware)\n")
cat("  Results root     : ", results_root, "\n")
cat("  Targets dir      : ", targets_dir, "\n")
cat("  Library dir      : ", library_dir, "\n")
cat("  Drop palindromic : ", DROP_PALINDROMIC, "\n")
cat("============================================================\n\n")

# Target directories = immediate sub-directories of the results root that
# contain the expected FDRreg output structure.
candidate_dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
target_dirs <- candidate_dirs[
  dir.exists(file.path(candidate_dirs, "02.fdrreg_results")) |
  dir.exists(file.path(candidate_dirs, "00.overlap_data"))
]

if (length(target_dirs) == 0) {
  stop("No processed target directories found under: ", results_root)
}

audit_rows <- vector("list", length(target_dirs))

for (i in seq_along(target_dirs)) {
  target_dir  <- target_dirs[i]
  target_name <- basename(target_dir)

  cat(sprintf("[%d/%d] Target: %s\n", i, length(target_dirs), target_name))

  # Requirement 3: traits from finished results.
  trait_info <- get_traits_for_target(target_dir, target_name)
  traits     <- trait_info$traits

  # Requirement 1: SNP count from generated files.
  res_count <- get_result_snp_count(target_dir, target_name)

  # Requirement 2: recomputed overlap using new input sources + harmonisation.
  if (length(traits) == 0) {
    new_overlap <- list(count = NA_integer_, missing = "no_traits_found")
  } else {
    new_overlap <- recompute_overlap(target_name, traits)
  }

  cat(sprintf("    traits (%s, n=%d): %s\n",
              trait_info$source, length(traits),
              ifelse(length(traits) == 0, "(none)", paste(traits, collapse = ", "))))
  cat(sprintf("    SNPs from results : %s  [%s]\n",
              ifelse(is.na(res_count$count), "NA", res_count$count), res_count$source))
  cat(sprintf("    SNPs new overlap  : %s  (post-harmonisation)\n",
              ifelse(is.na(new_overlap$count), "NA", new_overlap$count)))
  if (length(new_overlap$missing) > 0) {
    cat(sprintf("    missing inputs    : %s\n",
                paste(new_overlap$missing, collapse = ", ")))
  }
  cat("\n")

  audit_rows[[i]] <- data.table(
    target                 = target_name,
    n_traits               = length(traits),
    traits                 = paste(traits, collapse = ";"),
    traits_source          = trait_info$source,
    snp_count_from_results = res_count$count,
    snp_count_source       = res_count$source,
    snp_count_new_overlap  = new_overlap$count,
    drop_palindromic       = DROP_PALINDROMIC,
    missing_inputs         = paste(new_overlap$missing, collapse = ";")
  )
}

audit_dt <- rbindlist(audit_rows, use.names = TRUE, fill = TRUE)

fwrite(audit_dt, audit_out)

cat("============================================================\n")
cat("Audit complete.\n")
cat("  Targets audited : ", nrow(audit_dt), "\n")
cat("  Summary written : ", audit_out, "\n")
cat("============================================================\n")

print(audit_dt)

q("no")
