#*********************************************#
#-----DEA Pipeline - SINGLE TARGET SMOKE TEST-#
#-----R-4.0.2 / fdrreg_rebuild version--------#
#*********************************************#

## ============================================================
## 1. CONFIG (test)
## ============================================================
rm(list = ls())
suppressMessages({
  library(data.table); library(stringr); library(dplyr)
  library(Hotelling); library(mppa); library(ICSNP); library(atc); library(readr)
})

cfg <- list(
  test_target = "scz2014",              # <-- single target under test

  base_dir  = "/path/to/SO_Lab/18.fdrreg_rebuild/",
  out_dir   = "/path/to/SO_Lab/18.fdrreg_rebuild/01.extra.analysis/03.drug_enrichment/TEST/",

  lib_entrez = "/path/to/global.files/magma-library-uniq-entrez.csv",
  drug_mat   = "/path/to/global.files/mat.drug_DSigDB.Rdata",
  atc_lists  = "/path/to/global.files/ATC_drug_lists_ALL.Rdata",
  fn_atc     = "/path/to/global.files/005.drug_ATC.category_asFunc.R",

  magma_sub  = "05.magma_fdrreg",
  smx_sub    = "09.smultixcan_fdrreg",

  levels   = c("magma","smultixcan"),
  type_map = c(qval = "qval", fdr = "FDR.the", biofdr = "bio.FDR.the"),

  drug7 = data.frame(
    level3_codes = c("N03A","N06A","N05A","N04B","N06B","N05B","N05C"),
    NAME = c("Anti-epileptics","Antidepressants","Antipsychotics",
             "Dopaminergic agents","Psychostimulants","Anxiolytics","Hypnotics"),
    stringsAsFactors = FALSE)
)
cfg$drug5 <- cfg$drug7[1:5, ]

dir_enrich <- file.path(cfg$out_dir, "01.enrich")
dir_cat    <- file.path(cfg$out_dir, "02.category")
dir_summ   <- file.path(cfg$out_dir, "03.summary")
for (d in c(dir_enrich, dir_cat, dir_summ))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## small helper for check-point logging
chk <- function(pass, msg) {
  cat(sprintf("  [%s] %s\n", ifelse(pass, "OK  ", "FAIL"), msg))
  if (!pass) stop(msg)
}

## ============================================================
## 2. FUNCTIONS (identical to the full pipeline)
## ============================================================
enrich_path <- function(target, level, type, suffix)
  file.path(dir_enrich, paste0(target, "_", level, "_", type, ".enrich.", suffix))
cat_path <- function(target, level, type)
  file.path(dir_cat, paste0(target, "_", level, "_", type, ".category.csv"))

load_gene <- function(target, level, lib_entrez) {
  if (level == "magma") {
    f <- file.path(cfg$base_dir, target, cfg$magma_sub,
                   paste0(target, ".gene.bio.fdrreg.txt"))
    if (!file.exists(f)) stop(paste("missing file:", f))
    g <- fread(f)
    if (!"GeneName" %in% colnames(g))
      g <- left_join(g, lib_entrez[, c("ID", "GeneName")], by = "ID")
    g$Gene <- g$GeneName
  } else {
    f <- file.path(cfg$base_dir, target, cfg$smx_sub,
                   paste0(target, ".gene.bio.fdrreg.txt"))
    if (!file.exists(f)) stop(paste("missing file:", f))
    g <- fread(f)
    g$Gene <- g$ENSEMBL_GENE_ID_name
  }
  as.data.frame(g)
}

resolve_pvalue <- function(dt, type_col) {
  if (!type_col %in% colnames(dt))
    stop(sprintf("column '%s' not found", type_col))
  dt$Pvalue <- as.numeric(dt[[type_col]])
  dt
}

drug_enrichment_test <- function(trait_gene, mat.drug) {
  pv <- trait_gene$Pvalue
  pv[pv == 1] <- 0.9999
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
    if (sum(genes.for.drug == 1) < 5) {
      pvals[i] <- 999
    } else {
      pvals[i] <- t.test(mat.disease2$zval_trait ~ genes.for.drug,
                         alternative = "greater")$p.value
    }
  }
  list(res = data.frame(drug = colnames(mat.drug2)[-1],
                        t.test.p.one.sided_trait = pvals, stringsAsFactors = FALSE),
       n_common = length(common))
}

add_fdr <- function(x, col_in, col_out) {
  x[[col_out]] <- 999
  keep <- which(x[[col_in]] != 999)
  if (length(keep) > 0)
    x[[col_out]][keep] <- p.adjust(x[[col_in]][keep], method = "fdr")
  x
}

## ============================================================
## 3. LOAD SHARED RESOURCES (with checks)
## ============================================================
cat("== Loading shared resources ==\n")
chk(file.exists(cfg$fn_atc),    "ATC function script exists")
source(cfg$fn_atc)
chk(exists("ATC_enrichment"),   "ATC_enrichment() loaded")

chk(file.exists(cfg$lib_entrez),"entrez library exists")
lib_entrez <- fread(cfg$lib_entrez)
chk(all(c("ID","GeneName") %in% colnames(lib_entrez)), "entrez library has ID/GeneName")

chk(file.exists(cfg$drug_mat),  "drug matrix Rdata exists")
load(cfg$drug_mat)
chk(exists("mat.drug"),         "mat.drug loaded")
cat(sprintf("       mat.drug dim = %d x %d\n", nrow(mat.drug), ncol(mat.drug)))

chk(file.exists(cfg$atc_lists), "ATC lists Rdata exists")
load(cfg$atc_lists)
chk(exists("ATC_drug_lists"),   "ATC_drug_lists loaded")
cat(sprintf("       level3_codes present: %s\n", exists("level3_codes")))

## ============================================================
## 4. STAGE A on the single target
## ============================================================
target <- cfg$test_target
cat(sprintf("\n== STAGE A on target: %s ==\n", target))

for (level in cfg$levels) {
  cat(sprintf("\n-- level: %s --\n", level))
  gene_raw <- load_gene(target, level, lib_entrez)
  chk(nrow(gene_raw) > 0, sprintf("gene file loaded (%d rows)", nrow(gene_raw)))
  chk("Gene" %in% colnames(gene_raw), "Gene column present")
  chk(sum(!is.na(gene_raw$Gene)) > 0, "Gene column non-empty")

  ## quick match check against drug matrix
  n_match <- length(intersect(mat.drug$Gene, gene_raw$Gene))
  cat(sprintf("       genes matched to mat.drug = %d\n", n_match))
  chk(n_match >= 100, "at least 100 genes match mat.drug")

  for (tl in names(cfg$type_map)) {
    tcol <- cfg$type_map[[tl]]
    cat(sprintf("   type %s (col %s):\n", tl, tcol))
    chk(tcol %in% colnames(gene_raw), sprintf("source column '%s' exists", tcol))

    trait_gene <- resolve_pvalue(gene_raw, tcol)
    rng <- range(trait_gene$Pvalue, na.rm = TRUE)
    cat(sprintf("       Pvalue range = [%.4g, %.4g], NA = %d\n",
                rng[1], rng[2], sum(is.na(trait_gene$Pvalue))))

    ## per-drug enrichment
    de <- drug_enrichment_test(trait_gene, mat.drug)
    DEA <- de$res
    DEA <- DEA[order(DEA$t.test.p.one.sided_trait), ]
    DEA <- add_fdr(DEA, "t.test.p.one.sided_trait", "p.adj")
    n_tested <- sum(DEA$t.test.p.one.sided_trait != 999)
    cat(sprintf("       drugs tested (>=5 genes) = %d / %d\n", n_tested, nrow(DEA)))
    chk(n_tested > 0, "at least one drug tested")
    save(DEA, file = enrich_path(target, level, tl, "rdata"))
    fwrite(DEA, enrich_path(target, level, tl, "csv"), sep = ",")

    ## ATC category enrichment
    atc_res <- ATC_enrichment(enrich_path(target, level, tl, "csv"),
                              "t.test.p.one.sided_trait", "drug")
    chk(nrow(atc_res) > 0, "ATC_enrichment returned rows")
    atc_res <- add_fdr(atc_res, "pval.oneSamp.t", "p.adj.one")
    atc_res <- add_fdr(atc_res, "pval.twoSamp.t", "p.adj.two")
    atc_res <- arrange(atc_res, pval.twoSamp.t)
    fwrite(atc_res, cat_path(target, level, tl), sep = ",")
    cat(sprintf("       category rows = %d ; wrote %s\n",
                nrow(atc_res), basename(cat_path(target, level, tl))))
  }
}

## ============================================================
## 5. STAGE B on the single target
## ============================================================
cat("\n== STAGE B (summary) ==\n")
long <- list()
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
chk(nrow(long) > 0, sprintf("long table built (%d rows)", nrow(long)))

make_subset <- function(codes_df, fname) {
  sub <- long[long$level3_codes %in% codes_df$level3_codes, ]
  sub <- left_join(sub, codes_df, by = "level3_codes")
  fwrite(sub, file.path(dir_summ, fname), sep = ",")
  cat(sprintf("       %s : %d rows\n", fname, nrow(sub)))
}
make_subset(cfg$drug5, "drug5_detail.csv")
make_subset(cfg$drug7, "drug7_detail.csv")

cat("\n== SMOKE TEST FINISHED OK ==\n")
cat(sprintf("Outputs under: %s\n", cfg$out_dir))
q("no")
