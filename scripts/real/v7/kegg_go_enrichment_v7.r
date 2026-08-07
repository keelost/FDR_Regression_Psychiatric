#-----------------------------#
#---- KEGG/GO  ----#
#-----------------------------#

# SoftWare: R-4.5.3
# Platform: Server 203
# Date: 2026/07/08
# Author: Jinghong QIU
# Description: KEGG/GO enrichment on FDRreg gene lists, using each file's own
#              gene set as the enrichment background. Results are written per
#              target (dataset) as tab-separated files.
# Guide to del: install.packages('name') BiocManager::install('name')
# export NUM_THREADS=80
# Usage:
#   Rscript kegg_go_enrichment.r --p 0.01 --bio Y
#     --p   significance threshold (sig gene = FDR column < p)
#     --bio Y -> use bio.FDR.the ; N -> use FDR.the

#---- load packages and set dir ----#
rm(list = ls())
chooseCRANmirror(ind = 22) # for HK mirror setting, getCRANmirrors() can check the number
options(stringsAsFactors = FALSE)
start_time <- Sys.time() # time tracking
set.seed(100)
# load packages
packages <- c("data.table", "gprofiler2", "dplyr", "optparse")
invisible(lapply(packages, library, character.only = TRUE))

#---- Pipeline ----#
# ============================================================
# 0. Parse command-line arguments
# ============================================================
option_list <- list(
  make_option(c("-p", "--p"), type = "double", default = 0.05,
              help = "Significance threshold; a gene is significant when the chosen FDR column < p [default %default]"),
  make_option(c("-b", "--bio"), type = "character", default = "N",
              help = "Y -> use bio.FDR.the column; N -> use FDR.the column [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

SIG_FDR_THRESHOLD <- opt$p
bio_flag <- toupper(opt$bio)
if (!bio_flag %in% c("Y", "N")) {
  stop("--bio must be 'Y' or 'N'")
}
SIG_FDR_COLUMN <- if (bio_flag == "Y") "bio.FDR.the" else "FDR.the"
significance_label <- paste0(SIG_FDR_COLUMN, " < ", SIG_FDR_THRESHOLD)
message("Significance rule: ", significance_label)

# ============================================================
# 1. Configuration: input/output paths
# ============================================================
# Root directory that holds one sub-directory per phenotype/target.
input_root_dir <- Sys.getenv("FDRREG_RESULTS_DIR", "")
if (!nzchar(input_root_dir)) stop("Set FDRREG_RESULTS_DIR before running v7 KEGG/GO enrichment.")

# Sub-directory markers that distinguish the two analysis methods.
magma_subdir_marker      <- "05.magma_fdrreg"
smultixcan_subdir_marker <- "09.smultixcan_fdrreg_v7"

# Common file suffix for both methods.
result_file_pattern <- "\\.gene\\.bio\\.fdrreg\\.txt$"

# Gene ID column names differ between methods.
magma_gene_column      <- "ID"                # MAGMA reports Entrez gene IDs
smultixcan_gene_column <- "ENSEMBL_GENE_ID"   # SMultiXCan reports Ensembl gene IDs

# Output directory (a per-target sub-folder is created under this path).
output_dir <- file.path(input_root_dir, "01.extra.analysis", "04.kegg_go_v7")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Enrichment source sets.
pathway_sources <- c("REAC", "KEGG", "WP")
go_sources      <- c("GO:BP", "GO:MF", "GO:CC")

# Tag used in output file names to record the significance rule.
sig_tag <- paste0(gsub("\\.", "_", SIG_FDR_COLUMN), "_lt_", SIG_FDR_THRESHOLD)

# ============================================================
# 2. Discover result files and split by method
# ============================================================
all_result_files <- list.files(
  path       = input_root_dir,
  pattern    = result_file_pattern,
  full.names = TRUE,
  recursive  = TRUE
)

magma_files      <- grep(paste0("/", magma_subdir_marker, "/"),      all_result_files, value = TRUE)
smultixcan_files <- grep(paste0("/", smultixcan_subdir_marker, "/"), all_result_files, value = TRUE)

message("Found ", length(magma_files), " MAGMA files and ",
        length(smultixcan_files), " SMultiXCan files.")

# ============================================================
# 3. Derive dataset (target) name from file path
#    Path layout: <input_root_dir>/<dataset>/<method_subdir>/<dataset>.gene.bio.fdrreg.txt
#    so the dataset name is the grandparent directory.
# ============================================================
get_dataset_name <- function(file_path) {
  basename(dirname(dirname(file_path)))
}

# ============================================================
# 4. Detect gene ID namespace (Ensembl vs Entrez)
# ============================================================
detect_gene_id_type <- function(gene_ids) {
  gene_ids <- unique(gene_ids)
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]

  if (length(gene_ids) == 0) return(NA_character_)

  if (all(grepl("^ENSG", gene_ids))) {
    return("ensembl")
  } else if (all(grepl("^[0-9]+$", gene_ids))) {
    return("entrez")
  } else {
    return("mixed")
  }
}

# ============================================================
# 5. Core enrichment runner
#    Runs gprofiler2::gost with the file's own gene set as the
#    custom background (domain_scope = "custom_annotated").
# ============================================================
run_gost_enrichment <- function(query_genes, background_genes, enrichment_sources,
                                 dataset, method_label, filter_label) {
  query_genes <- unique(as.character(query_genes))
  query_genes <- query_genes[!is.na(query_genes) & query_genes != ""]

  background_genes <- unique(as.character(background_genes))
  background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]

  if (length(query_genes) < 2) {
    warning(paste("Skip", dataset, method_label, filter_label, ": fewer than 2 query genes"))
    return(NULL)
  }

  query_gene_id_type <- detect_gene_id_type(query_genes)
  if (is.na(query_gene_id_type) || query_gene_id_type == "mixed") {
    warning(paste("Skip", dataset, method_label, filter_label,
                  ": cannot resolve a single gene ID namespace"))
    return(NULL)
  }

  gost_args <- list(
    query             = query_genes,
    organism          = "hsapiens",
    sources           = enrichment_sources,
    correction_method = "fdr",
    significant       = FALSE,
    evcodes           = TRUE,
    custom_bg         = background_genes,   # use the file's genes as background
    domain_scope      = "custom_annotated"  # restrict universe to annotated custom genes
  )
  if (query_gene_id_type == "entrez") {
    gost_args$numeric_ns <- "ENTREZGENE_ACC"
  }

  gost_result <- tryCatch(
    do.call(gost, gost_args),
    error = function(e) {
      warning(paste("gost failed for", dataset, method_label, filter_label, ":", e$message))
      return(NULL)
    }
  )

  if (is.null(gost_result) || is.null(gost_result$result) || nrow(gost_result$result) == 0) {
    return(NULL)
  }

  list(
    result_table       = as.data.frame(gost_result$result),
    query_gene_id_type = query_gene_id_type
  )
}

# ============================================================
# 6. Formatters for pathway and GO results
# ============================================================
format_pathway_result <- function(gost_output, dataset, method_label, filter_label) {
  if (is.null(gost_output)) return(NULL)

  result_table <- gost_output$result_table
  result_table$q_value <- p.adjust(result_table$p_value, method = "fdr")

  result_table %>%
    transmute(
      dataset                       = dataset,
      source_method                 = method_label,
      filter_type                   = filter_label,
      gene_id_type                  = gost_output$query_gene_id_type,
      p_value                       = p_value,
      q_value                       = q_value,
      pathway                       = term_name,
      source                        = source,
      external_id                   = term_id,
      members_input_overlap         = intersection,
      members_input_overlap_geneids = intersection,
      size                          = term_size,
      effective_size                = effective_domain_size
    ) %>%
    arrange(p_value)
}

format_go_result <- function(gost_output, dataset, method_label, filter_label) {
  if (is.null(gost_output)) return(NULL)

  result_table <- gost_output$result_table
  result_table$q_value <- p.adjust(result_table$p_value, method = "fdr")

  result_table %>%
    transmute(
      dataset                       = dataset,
      source_method                 = method_label,
      filter_type                   = filter_label,
      gene_id_type                  = gost_output$query_gene_id_type,
      `p-value`                     = p_value,
      `q-value`                     = q_value,
      term_goid                     = term_id,
      term_category                 = source,
      term_level                    = NA_integer_,
      term_name                     = term_name,
      members_input_overlap         = intersection,
      members_input_overlap_geneids = intersection,
      size                          = term_size,
      effective_size                = effective_domain_size
    ) %>%
    arrange(`p-value`)
}

# ============================================================
# 7. Process a single result file
#    background_genes = all genes in the file
#    query_genes      = genes passing the significance rule
# ============================================================
process_file <- function(file_path, method_label, gene_column) {
  gene_stat_table <- fread(file_path)
  dataset <- get_dataset_name(file_path)

  if (!(SIG_FDR_COLUMN %in% names(gene_stat_table)) || !(gene_column %in% names(gene_stat_table))) {
    warning(paste("Skip", method_label, "file (missing", SIG_FDR_COLUMN, "or",
                  gene_column, "column):", file_path))
    return(NULL)
  }

  background_genes <- as.character(gene_stat_table[[gene_column]])

  # Coerce the FDR column to numeric so that any "NA" strings or non-numeric
  # entries become proper NA values before the threshold comparison.
  fdr_values <- suppressWarnings(as.numeric(gene_stat_table[[SIG_FDR_COLUMN]]))
  significant_rows <- fdr_values < SIG_FDR_THRESHOLD
  significant_rows[is.na(significant_rows)] <- FALSE
  query_genes <- as.character(gene_stat_table[[gene_column]][significant_rows])

  pathway_output <- run_gost_enrichment(
    query_genes, background_genes, pathway_sources,
    dataset, method_label, significance_label
  )
  go_output <- run_gost_enrichment(
    query_genes, background_genes, go_sources,
    dataset, method_label, significance_label
  )

  list(
    pathway = format_pathway_result(pathway_output, dataset, method_label, significance_label),
    go      = format_go_result(go_output, dataset, method_label, significance_label)
  )
}

# ============================================================
# 8. Sanitise gene-list fields and write tab-separated output
#    Gene lists contain commas; converting to ';' and writing TSV
#    avoids any clash with field separators.
# ============================================================
sanitise_gene_list_fields <- function(result_dt) {
  if ("members_input_overlap" %in% names(result_dt)) {
    result_dt[, members_input_overlap := gsub(",", ";", members_input_overlap)]
  }
  if ("members_input_overlap_geneids" %in% names(result_dt)) {
    result_dt[, members_input_overlap_geneids := gsub(",", ";", members_input_overlap_geneids)]
  }
  result_dt
}

write_target_result <- function(result_list, dataset, analysis_label) {
  combined <- if (length(result_list) > 0) rbindlist(result_list, fill = TRUE) else data.table()
  if (nrow(combined) == 0) {
    message("No ", analysis_label, " results for target ", dataset, "; nothing written.")
    return(invisible(NULL))
  }
  combined <- sanitise_gene_list_fields(combined)

  target_output_dir <- file.path(output_dir, dataset)
  if (!dir.exists(target_output_dir)) dir.create(target_output_dir, recursive = TRUE)

  output_file <- file.path(
    target_output_dir,
    paste0(dataset, "_", analysis_label, "_enrichment_", sig_tag, ".tsv")
  )
  fwrite(combined, output_file, sep = "\t", quote = TRUE)
  message("Wrote ", output_file)
  invisible(output_file)
}

# ============================================================
# 9. Batch process, grouped per target (dataset)
# ============================================================
# Map each target to its MAGMA and SMultiXCan files.
file_registry <- rbindlist(list(
  if (length(magma_files) > 0)
    data.table(file_path = magma_files, method_label = "MAGMA",
               gene_column = magma_gene_column),
  if (length(smultixcan_files) > 0)
    data.table(file_path = smultixcan_files, method_label = "SMultiXCan",
               gene_column = smultixcan_gene_column)
), fill = TRUE)

if (nrow(file_registry) == 0) {
  stop("No result files found under ", input_root_dir)
}

file_registry[, dataset := vapply(file_path, get_dataset_name, character(1))]

for (current_dataset in sort(unique(file_registry$dataset))) {
  message("==== Target: ", current_dataset, " ====")
  target_files <- file_registry[dataset == current_dataset]

  pathway_result_list <- list()
  go_result_list      <- list()

  for (i in seq_len(nrow(target_files))) {
    file_path    <- target_files$file_path[i]
    method_label <- target_files$method_label[i]
    gene_column  <- target_files$gene_column[i]

    message("Processing ", method_label, ": ", basename(file_path))
    file_result <- process_file(file_path, method_label, gene_column)
    if (is.null(file_result)) next
    if (!is.null(file_result$pathway)) {
      pathway_result_list[[length(pathway_result_list) + 1]] <- file_result$pathway
    }
    if (!is.null(file_result$go)) {
      go_result_list[[length(go_result_list) + 1]] <- file_result$go
    }
  }

  write_target_result(pathway_result_list, current_dataset, "pathway")
  write_target_result(go_result_list,      current_dataset, "GO")
}

cat("Done!\n")

#--- quit and save result ----#
# Time track
end_time <- Sys.time()
cat('Time Start:', format(start_time, "%a %b %d %X %Y"), '\n')
cat('Time End:', format(end_time, "%a %b %d %X %Y"), '\n')
cat('Time consuming:', end_time - start_time, '\n')

q('no')
