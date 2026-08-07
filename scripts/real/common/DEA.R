#*********************************************#
#-----Drug Enrichment Analysis Pipeline-------#
#-----Unified: enrichment + ATC + summaries---#
#-----R-4.0.2 / fdrreg_rebuild version--------#
#*********************************************#

## ============================================================
## 1. CONFIG
## ============================================================
rm(list = ls())
suppressMessages({
  library(data.table); library(stringr); library(dplyr)
  library(Hotelling); library(mppa); library(ICSNP); library(atc); library(readr)
})

cfg <- list(
  ## input root: one sub-dir per target
  base_dir  = Sys.getenv("FDRREG_RESULTS_DIR", ""),

  ## output root + sub-dirs
  out_dir    = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "03.drug_enrichment"),

  ## shared resources (global.files)
  gdir       = Sys.getenv("FDRREG_GLOBAL_DIR", ""),
  lib_entrez = Sys.getenv("FDRREG_BIO_ENTREZ", ""),
  lib_ensbl  = Sys.getenv("FDRREG_BIO_ENSEMBL", ""),
  drug_mat   = Sys.getenv("FDRREG_DRUG_MATRIX", ""),
  atc_lists  = Sys.getenv("FDRREG_ATC_LISTS", ""),
  fn_atc     = Sys.getenv("FDRREG_ATC_FUNCTION", ""),

  ## per-target relative paths
  magma_sub  = "05.magma_fdrreg",
  smx_sub    = "09.smultixcan_fdrreg",

  ## loop dimensions
  targets = c("adhd2016","adhd2019","asd2015.pgc","asd2019","bd2012","bd2018",
              "mdd2013","mdd2019","mddco","sa.bpd.pgc","sa.ipsych","sa.scz.pgc",
              "scz2012","scz2014","scz.eas2019","wal.scz2018"),
  levels  = c("magma","smultixcan"),

  ## type label -> source column used as the enrichment p-value
  type_map = c(qval = "qval", fdr = "FDR.the", biofdr = "bio.FDR.the"),

  ## drug categories of interest for summary subsets
  drug7 = data.frame(
    level3_codes = c("N03A","N06A","N05A","N04B","N06B","N05B","N05C"),
    NAME = c("Anti-epileptics","Antidepressants","Antipsychotics",
             "Dopaminergic agents","Psychostimulants","Anxiolytics","Hypnotics"),
    stringsAsFactors = FALSE),

  ## stage switches
  run_stage_a = TRUE,   # enrichment + ATC category
  run_stage_b = TRUE    # summary tables
)
cfg$drug5 <- cfg$drug7[1:5, ]   # first five categories

## output sub-dirs
dir_enrich <- file.path(cfg$out_dir, "01.enrich")
dir_cat    <- file.path(cfg$out_dir, "02.category")
dir_summ   <- file.path(cfg$out_dir, "03.summary")
for (d in c(dir_enrich, dir_cat, dir_summ))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## ============================================================
## 2. FUNCTIONS
## ============================================================

## Normalised path builders: {target}_{level}_{type}.<suffix> ---------------
enrich_path <- function(target, level, type, suffix)
  file.path(dir_enrich, paste0(target, "_", level, "_", type, ".enrich.", suffix))

cat_path <- function(target, level, type)
  file.path(dir_cat, paste0(target, "_", level, "_", type, ".category.csv"))

## Load gene-based result for one (target, level); returns df with $Gene -----
load_gene <- function(target, level, lib_entrez) {
  if (level == "magma") {
    f <- file.path(cfg$base_dir, target, cfg$magma_sub,
                   paste0(target, ".gene.bio.fdrreg.txt"))
    g <- fread(f)                                   # space/auto delimited
    if (!"GeneName" %in% colnames(g))
      g <- left_join(g, lib_entrez[, c("ID", "GeneName")], by = "ID")
    g$Gene <- g$GeneName
  } else {
    f <- file.path(cfg$base_dir, target, cfg$smx_sub,
                   paste0(target, ".gene.bio.fdrreg.txt"))
    g <- fread(f)
    g$Gene <- g$ENSEMBL_GENE_ID_name                # symbol already present
  }
  as.data.frame(g)
}

## Pick the source column as "Pvalue" (only per-type difference) -------------
resolve_pvalue <- function(dt, type_col) {
  if (!type_col %in% colnames(dt))
    stop(sprintf("column '%s' not found in gene file", type_col))
  dt$Pvalue <- as.numeric(dt[[type_col]])
  dt
}

## Drug enrichment test (explicit args; robust -Inf / duplicate handling) ----
drug_enrichment_test <- function(trait_gene, mat.drug) {
  pv <- trait_gene$Pvalue
  pv[pv == 1] <- 0.9999                             # guard qnorm(1) = Inf
  mat.disease <- data.frame(Gene = trait_gene$Gene, zval_trait = qnorm(pv))

  common <- intersect(mat.drug$Gene, mat.disease$Gene)
  mat.disease2 <- arrange(mat.disease[mat.disease$Gene %in% common, ], Gene)
  mat.disease2$zval_trait[mat.disease2$zval_trait == -Inf] <- -5e-08
  mat.drug2 <- arrange(mat.drug[mat.drug$Gene %in% common, ], Gene)
  if (length(mat.disease2$Gene) != length(mat.drug2$Gene))
    mat.disease2 <- mat.disease2[!duplicated(mat.disease2$Gene), ]

  no.drugs <- ncol(mat.drug2) - 1
  pvals <- numeric(no.drugs)
  for (i in seq_len(no.drugs)) {
    genes.for.drug <- mat.drug2[, i + 1]
    if (sum(genes.for.drug == 1) < 5) {             # need >= 5 associated genes
      pvals[i] <- 999
    } else {
      pvals[i] <- t.test(mat.disease2$zval_trait ~ genes.for.drug,
                         alternative = "greater")$p.value
    }
  }
  data.frame(drug = colnames(mat.drug2)[-1],
             t.test.p.one.sided_trait = pvals, stringsAsFactors = FALSE)
}

## BH-adjust, skipping the 999 placeholders ---------------------------------
add_fdr <- function(x, col_in, col_out) {
  x[[col_out]] <- 999
  keep <- which(x[[col_in]] != 999)
  if (length(keep) > 0)
    x[[col_out]][keep] <- p.adjust(x[[col_in]][keep], method = "fdr")
  x
}

## ============================================================
## 3. STAGE A : enrichment + ATC category
## ============================================================
if (cfg$run_stage_a) {
  source(cfg$fn_atc)                 # ATC_enrichment()
  lib_entrez <- fread(cfg$lib_entrez)
  load(cfg$drug_mat)                 # -> mat.drug
  load(cfg$atc_lists)                # -> ATC_drug_lists (+ level3_codes)

  for (target in cfg$targets) {
    for (level in cfg$levels) {
      gene_raw <- tryCatch(load_gene(target, level, lib_entrez), error = function(e) {
        cat(sprintf("[SKIP] %s | %s : %s\n", target, level, conditionMessage(e))); NULL })
      if (is.null(gene_raw)) next

      for (tl in names(cfg$type_map)) {
        tcol <- cfg$type_map[[tl]]
        cat(sprintf("[STAGE A] %s | %s | %s (%s)\n", target, level, tl, tcol))

        trait_gene <- resolve_pvalue(gene_raw, tcol)

        ## ---- per-drug enrichment (with BH) ----
        DEA <- drug_enrichment_test(trait_gene, mat.drug)
        DEA <- DEA[order(DEA$t.test.p.one.sided_trait), ]
        DEA <- add_fdr(DEA, "t.test.p.one.sided_trait", "p.adj")
        save(DEA, file = enrich_path(target, level, tl, "rdata"))
        fwrite(DEA, enrich_path(target, level, tl, "csv"), sep = ",")

        ## ---- ATC category enrichment (BH on BOTH one/two sample) ----
        atc_res <- ATC_enrichment(enrich_path(target, level, tl, "csv"),
                                  "t.test.p.one.sided_trait", "drug")
        atc_res <- add_fdr(atc_res, "pval.oneSamp.t", "p.adj.one")
        atc_res <- add_fdr(atc_res, "pval.twoSamp.t", "p.adj.two")
        atc_res <- arrange(atc_res, pval.twoSamp.t)
        fwrite(atc_res, cat_path(target, level, tl), sep = ",")
      }
    }
  }
}

## ============================================================
## 4. STAGE B : summary tables (tidy long + focused subsets)
## ============================================================
if (cfg$run_stage_b) {

  ## ---- 4.1 collect all category results into one long table ----
  long <- vector("list", 0L)
  for (target in cfg$targets)
    for (level in cfg$levels)
      for (tl in names(cfg$type_map)) {
        fp <- cat_path(target, level, tl)
        if (!file.exists(fp)) next
        d <- fread(fp, data.table = FALSE)
        d$target <- target; d$level <- level; d$type <- tl
        long[[length(long) + 1L]] <- d
      }
  long <- rbindlist(long, use.names = TRUE, fill = TRUE)
  fwrite(long, file.path(dir_summ, "category_all_long.csv"), sep = ",")

  ## ---- 4.2 focused subsets: 5-drug and 7-drug tables ----
  make_subset <- function(codes_df, fname) {
    sub <- long[long$level3_codes %in% codes_df$level3_codes, ]
    sub <- left_join(sub, codes_df, by = "level3_codes")   # attach friendly NAME
    fwrite(sub, file.path(dir_summ, fname), sep = ",")
  }
  make_subset(cfg$drug5, "drug5_detail.csv")
  make_subset(cfg$drug7, "drug7_detail.csv")
}

q("no")
