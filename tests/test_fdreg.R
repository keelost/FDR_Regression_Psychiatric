#!/usr/bin/env Rscript

all_args <- commandArgs(FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
file_arg <- sub("^--file=", "", file_arg[1])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = FALSE)
if (!dir.exists(root)) root <- normalizePath(".")

targets <- read.delim(file.path(root, "config", "targets.tsv"), check.names = FALSE)
stopifnot(nrow(targets) == 16L)
stopifnot(length(unique(targets$target)) == 16L)
stopifnot(all(nzchar(targets$target)))

regions <- readLines(file.path(root, "config", "regions.txt"), warn = FALSE)
regions <- regions[nzchar(trimws(regions))]
stopifnot(length(regions) == 13L)
stopifnot(length(unique(regions)) == 13L)

cat("R configuration checks passed\n")
