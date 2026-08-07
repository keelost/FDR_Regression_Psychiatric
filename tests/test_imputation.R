suppressPackageStartupMessages(library(data.table))

root <- normalizePath(".", winslash = "/")
tmp <- normalizePath(tempfile("fdrreg-imputation-"), winslash = "/", mustWork = FALSE)
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

clear <- data.table(
  snpid = c("rs1", "rs1", "rs2", "rs_indel", "rs_low_info"),
  chr = c(1, 1, 2, 2, 2),
  bpos = c(100, 100, 200, 250, 300),
  a1 = c("a", "a", "C", "AT", "G"),
  a2 = c("g", "g", "T", "A", "A"),
  beta = c(0.2, 0.2, -0.4, 0.1, 0.1),
  se = c(0.1, 0.1, 0.2, 0.1, 0.1),
  info = c(0.9, 0.9, 0.8, 0.9, 0.5)
)
clear_file <- file.path(tmp, "mini.clear.txt")
dist_input <- file.path(tmp, "mini.4impute.txt")
fwrite(clear, clear_file, sep = " ")

status <- system2(
  "Rscript",
  c(file.path(root, "scripts/imputation/prepare_dist_input.R"),
    "--input", clear_file, "--output", dist_input),
  stdout = TRUE, stderr = TRUE
)
stopifnot(is.null(attr(status, "status")) || attr(status, "status") == 0)
prepared <- fread(dist_input)
stopifnot(nrow(prepared) == 2L)
stopifnot(identical(prepared$snpid, c("rs1", "rs2")))
stopifnot(all.equal(prepared$z, c(2, -2)))
stopifnot(identical(prepared$a1, c("A", "C")))

for (chr in 1:22) {
  rows <- if (chr == 1) {
    data.table(marker = c("rs1", "1:150:A:G"), chr = 1, bp = c(100, 150),
      a1 = "A", a2 = "G", af1 = 0.2, z = c(2, 1),
      info = c(0.9, 0.8), pval = c(0.04, 0.2), type = "imputed")
  } else if (chr == 2) {
    data.table(marker = c("rs2", "rs_low"), chr = 2, bp = c(200, 300),
      a1 = "C", a2 = "T", af1 = 0.3, z = c(-2, 0.5),
      info = c(0.5, 0.49), pval = c(0.04, 0.6), type = "imputed")
  } else {
    data.table(marker = character(), chr = integer(), bp = integer(),
      a1 = character(), a2 = character(), af1 = numeric(), z = numeric(),
      info = numeric(), pval = numeric(), type = character())
  }
  fwrite(rows, file.path(tmp, sprintf("mini.chr%d.imputation.txt", chr)), sep = " ")
}

variant_map <- data.table(chr = 1, pos = 150, rsid = "rs150")
variant_map_file <- file.path(tmp, "variant_map.txt")
fwrite(variant_map, variant_map_file, sep = " ")
final_file <- file.path(tmp, "mini.impute.map.txt")

status <- system2(
  "Rscript",
  c(file.path(root, "scripts/imputation/merge_dist_output.R"),
    "--trait", "mini", "--input-dir", tmp, "--variant-map", variant_map_file,
    "--output", final_file),
  stdout = TRUE, stderr = TRUE
)
stopifnot(is.null(attr(status, "status")) || attr(status, "status") == 0)
final <- fread(final_file)
stopifnot(identical(final$snpid, c("rs1", "rs150", "rs2")))
stopifnot(!"1:150:A:G" %in% final$snpid)
stopifnot(!"rs_low" %in% final$snpid)

cat("Imputation fixture tests passed.\n")
