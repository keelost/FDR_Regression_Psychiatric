#*********************************************#
#-----Drug Enrichment Analysis Pipeline-------#
#-----Unified: enrichment + ATC + summaries---#
#-----R-4.0.2 / fdrreg_rebuild v7-------------#
#*********************************************#

## ============================================================
## 1. CONFIG
## ============================================================

rm(list = ls())

suppressMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(Hotelling)
  library(mppa)
  library(ICSNP)
  library(atc)
  library(readr)
})

cfg <- list(

  ## Input root: one subdirectory per target.
  base_dir = Sys.getenv("FDRREG_RESULTS_DIR", ""),

  ## The "_v7" suffix distinguishes these results from old versions.
  out_dir = file.path(Sys.getenv("FDRREG_RESULTS_DIR", ""), "01.extra.analysis", "03.drug_enrichment_v7"),

  ## Shared resources.
  gdir = Sys.getenv("FDRREG_GLOBAL_DIR", ""),

  lib_entrez = paste0(
    Sys.getenv("FDRREG_BIO_ENTREZ", "")
  ),

  lib_ensbl = paste0(
    Sys.getenv("FDRREG_BIO_ENSEMBL", "")
  ),

  drug_mat = paste0(
    Sys.getenv("FDRREG_DRUG_MATRIX", "")
  ),

  atc_lists = paste0(
    Sys.getenv("FDRREG_ATC_LISTS", "")
  ),

  fn_atc = paste0(
    Sys.getenv("FDRREG_ATC_FUNCTION", "")
  ),

  ## Per-target relative paths.
  magma_sub = "05.magma_fdrreg",
  smx_sub = "09.smultixcan_fdrreg_v7",

  ## MetaXcan path retained for future tissue-level analyses.
  metaxcan_sub = "07.metaxcan_fdrreg_v7/01.fdrreg_results",

  ## Analysis dimensions.
  targets = c(
    "adhd2016",
    "adhd2019",
    "asd2015.pgc",
    "asd2019",
    "bd2012",
    "bd2018",
    "mdd2013",
    "mdd2019",
    "mddco",
    "sa.bpd.pgc",
    "sa.ipsych",
    "sa.scz.pgc",
    "scz2012",
    "scz2014",
    "scz.eas2019",
    "wal.scz2018"
  ),

  levels = c(
    "magma",
    "smultixcan"
  ),

  ## Type label and preferred source column.
  type_map = c(
    qval = "qval",
    fdr = "FDR.the",
    biofdr = "bio.FDR.the"
  ),

  ## Drug categories used in focused summaries.
  drug7 = data.frame(
    level3_codes = c(
      "N03A",
      "N06A",
      "N05A",
      "N04B",
      "N06B",
      "N05B",
      "N05C"
    ),
    NAME = c(
      "Anti-epileptics",
      "Antidepressants",
      "Antipsychotics",
      "Dopaminergic agents",
      "Psychostimulants",
      "Anxiolytics",
      "Hypnotics"
    ),
    stringsAsFactors = FALSE
  ),

  ## Stage switches.
  run_stage_a = TRUE,
  run_stage_b = TRUE
)

cfg$drug5 <- cfg$drug7[1:5, , drop = FALSE]

## Output subdirectories.
dir_enrich <- file.path(cfg$out_dir, "01.enrich")
dir_cat <- file.path(cfg$out_dir, "02.category")
dir_summ <- file.path(cfg$out_dir, "03.summary")

for (d in c(dir_enrich, dir_cat, dir_summ)) {
  if (!dir.exists(d)) {
    dir.create(
      d,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

## ============================================================
## 2. PATH BUILDERS
## ============================================================

enrich_path <- function(target, level, type, suffix) {
  file.path(
    dir_enrich,
    paste0(
      target,
      "_",
      level,
      "_",
      type,
      ".enrich.",
      suffix
    )
  )
}

cat_path <- function(target, level, type) {
  file.path(
    dir_cat,
    paste0(
      target,
      "_",
      level,
      "_",
      type,
      ".category.csv"
    )
  )
}

## ============================================================
## 3. INPUT FUNCTIONS
## ============================================================

load_gene <- function(target, level, lib_entrez) {

  if (level == "magma") {

    f <- file.path(
      cfg$base_dir,
      target,
      cfg$magma_sub,
      paste0(
        target,
        ".gene.bio.fdrreg.txt"
      )
    )

    if (!file.exists(f)) {
      stop(sprintf(
        "MAGMA result file does not exist: %s",
        f
      ))
    }

    ## fill = TRUE is required because the v7 MAGMA header contains
    ## two lasso columns that are absent from the data records.
    g <- fread(
      f,
      header = TRUE,
      fill = TRUE,
      data.table = FALSE,
      showProgress = FALSE
    )

    cat(sprintf(
      paste0(
        "[INFO] Loaded MAGMA file: %s | ",
        "rows=%d | columns=%d\n"
      ),
      f,
      nrow(g),
      ncol(g)
    ))

    if (!"ID" %in% colnames(g)) {
      stop(sprintf(
        "MAGMA file does not contain the required ID column: %s",
        f
      ))
    }

    ## Report completely empty columns.
    empty_columns <- colnames(g)[
      vapply(
        g,
        function(x) all(is.na(x)),
        logical(1L)
      )
    ]

    if (length(empty_columns) > 0L) {
      cat(sprintf(
        "[INFO] Completely empty MAGMA columns: %s\n",
        paste(empty_columns, collapse = ", ")
      ))
    }

    ## Add gene symbols if they are not already present.
    if (!"GeneName" %in% colnames(g)) {

      annotation <- as.data.frame(lib_entrez)[
        ,
        c("ID", "GeneName"),
        drop = FALSE
      ]

      g$ID <- as.character(g$ID)
      annotation$ID <- as.character(annotation$ID)
      annotation$GeneName <- trimws(
        as.character(annotation$GeneName)
      )

      valid_annotation <- !is.na(annotation$ID) &
        annotation$ID != "" &
        !is.na(annotation$GeneName) &
        annotation$GeneName != ""

      annotation <- annotation[
        valid_annotation,
        ,
        drop = FALSE
      ]

      duplicated_annotation_ids <- sum(
        duplicated(annotation$ID)
      )

      if (duplicated_annotation_ids > 0L) {
        cat(sprintf(
          paste0(
            "[INFO] Removing %d duplicated annotation IDs ",
            "before MAGMA gene mapping.\n"
          ),
          duplicated_annotation_ids
        ))

        annotation <- annotation[
          !duplicated(annotation$ID),
          ,
          drop = FALSE
        ]
      }

      ## Use match instead of left_join to guarantee unchanged row count.
      annotation_index <- match(
        g$ID,
        annotation$ID
      )

      g$GeneName <- annotation$GeneName[
        annotation_index
      ]
    }

    g$Gene <- trimws(
      as.character(g$GeneName)
    )

    missing_gene <- is.na(g$Gene) |
      g$Gene == "" |
      g$Gene == "NA"

    duplicated_gene_n <- sum(
      duplicated(g$Gene[!missing_gene])
    )

    cat(sprintf(
      paste0(
        "[INFO] MAGMA gene annotation: total=%d, ",
        "missing_gene=%d, duplicated_gene=%d\n"
      ),
      nrow(g),
      sum(missing_gene),
      duplicated_gene_n
    ))

  } else if (level == "smultixcan") {

    f <- file.path(
      cfg$base_dir,
      target,
      cfg$smx_sub,
      paste0(
        target,
        ".gene.bio.fdrreg.txt"
      )
    )

    if (!file.exists(f)) {
      stop(sprintf(
        "SMultiXcan v7 result file does not exist: %s",
        f
      ))
    }

    g <- fread(
      f,
      header = TRUE,
      fill = TRUE,
      data.table = FALSE,
      showProgress = FALSE
    )

    cat(sprintf(
      paste0(
        "[INFO] Loaded SMultiXcan file: %s | ",
        "rows=%d | columns=%d\n"
      ),
      f,
      nrow(g),
      ncol(g)
    ))

    if ("ENSEMBL_GENE_ID_name" %in% colnames(g)) {

      g$Gene <- trimws(
        as.character(g$ENSEMBL_GENE_ID_name)
      )

    } else if ("gene_name" %in% colnames(g)) {

      g$Gene <- trimws(
        as.character(g$gene_name)
      )

    } else {
      stop(sprintf(
        paste0(
          "SMultiXcan file does not contain ",
          "'ENSEMBL_GENE_ID_name' or 'gene_name': %s"
        ),
        f
      ))
    }

    missing_gene <- is.na(g$Gene) |
      g$Gene == "" |
      g$Gene == "NA"

    cat(sprintf(
      paste0(
        "[INFO] SMultiXcan gene annotation: total=%d, ",
        "missing_gene=%d, duplicated_gene=%d\n"
      ),
      nrow(g),
      sum(missing_gene),
      sum(duplicated(g$Gene[!missing_gene]))
    ))

  } else {
    stop(sprintf(
      "Unsupported analysis level: %s",
      level
    ))
  }

  g <- as.data.frame(g)

  ## Remove records without usable gene symbols.
  valid_gene <- !is.na(g$Gene) &
    g$Gene != "" &
    g$Gene != "NA"

  g <- g[
    valid_gene,
    ,
    drop = FALSE
  ]

  rownames(g) <- NULL

  g
}

## ============================================================
## 4. P-VALUE/FDR RESOLUTION
## ============================================================

resolve_pvalue <- function(dt, type_col) {

  dt <- as.data.frame(dt)

  ## Column aliases across MAGMA, SMultiXcan and MetaXcan.
  column_aliases <- list(

    qval = c(
      "qval"
    ),

    FDR.the = c(
      "FDR.the",
      "FDR_theoretical"
    ),

    bio.FDR.the = c(
      "bio.FDR.the",
      "bio_FDR_theoretical"
    ),

    FDR.emp = c(
      "FDR.emp",
      "FDR_empirical"
    ),

    bio.FDR.emp = c(
      "bio.FDR.emp",
      "bio_FDR_empirical"
    ),

    FDR.the.lasso = c(
      "FDR.the.lasso",
      "bio.FDR.the.lasso",
      "lasso_FDR_theoretical"
    ),

    FDR.emp.lasso = c(
      "FDR.emp.lasso",
      "bio.FDR.emp.lasso",
      "lasso_FDR_empirical"
    )
  )

  if (type_col %in% names(column_aliases)) {
    candidate_columns <- column_aliases[[type_col]]
  } else {
    candidate_columns <- type_col
  }

  available_columns <- candidate_columns[
    candidate_columns %in% colnames(dt)
  ]

  if (length(available_columns) == 0L) {
    stop(sprintf(
      paste0(
        "None of the candidate columns [%s] ",
        "was found in the gene file."
      ),
      paste(candidate_columns, collapse = ", ")
    ))
  }

  selected_column <- available_columns[1L]

  raw_value <- dt[[selected_column]]

  dt$Pvalue <- suppressWarnings(
    as.numeric(raw_value)
  )

  conversion_failed <- !is.na(raw_value) &
    trimws(as.character(raw_value)) != "" &
    is.na(dt$Pvalue)

  if (any(conversion_failed)) {
    warning(sprintf(
      paste0(
        "%d non-empty values in column '%s' ",
        "could not be converted to numeric."
      ),
      sum(conversion_failed),
      selected_column
    ))
  }

  invalid_value <- is.na(dt$Pvalue) |
    !is.finite(dt$Pvalue) |
    dt$Pvalue < 0 |
    dt$Pvalue > 1

  cat(sprintf(
    paste0(
      "[INFO] Column '%s': total=%d, ",
      "invalid=%d, valid=%d\n"
    ),
    selected_column,
    nrow(dt),
    sum(invalid_value),
    sum(!invalid_value)
  ))

  if (all(invalid_value)) {
    stop(sprintf(
      "No valid values remain in column '%s'.",
      selected_column
    ))
  }

  dt <- dt[
    !invalid_value,
    ,
    drop = FALSE
  ]

  valid_range <- range(
    dt$Pvalue,
    na.rm = TRUE
  )

  cat(sprintf(
    "[INFO] Column '%s' range: [%s, %s]\n",
    selected_column,
    format(valid_range[1L], scientific = TRUE),
    format(valid_range[2L], scientific = TRUE)
  ))

  attr(
    dt,
    "selected_pvalue_column"
  ) <- selected_column

  dt
}

## ============================================================
## 5. DRUG ENRICHMENT
## ============================================================

drug_enrichment_test <- function(trait_gene, mat.drug) {

  trait_gene <- as.data.frame(trait_gene)

  required_trait_columns <- c(
    "Gene",
    "Pvalue"
  )

  missing_trait_columns <- setdiff(
    required_trait_columns,
    colnames(trait_gene)
  )

  if (length(missing_trait_columns) > 0L) {
    stop(sprintf(
      "trait_gene is missing required columns: %s",
      paste(missing_trait_columns, collapse = ", ")
    ))
  }

  if (!"Gene" %in% colnames(mat.drug)) {
    stop(
      "mat.drug does not contain the required Gene column."
    )
  }

  trait_gene$Gene <- trimws(
    as.character(trait_gene$Gene)
  )

  trait_gene$Pvalue <- suppressWarnings(
    as.numeric(trait_gene$Pvalue)
  )

  drug_gene <- trimws(
    as.character(mat.drug$Gene)
  )

  valid_trait <- !is.na(trait_gene$Gene) &
    trait_gene$Gene != "" &
    !is.na(trait_gene$Pvalue) &
    is.finite(trait_gene$Pvalue) &
    trait_gene$Pvalue >= 0 &
    trait_gene$Pvalue <= 1

  trait_gene <- trait_gene[
    valid_trait,
    c("Gene", "Pvalue"),
    drop = FALSE
  ]

  if (nrow(trait_gene) == 0L) {
    stop(
      "No valid disease genes remain after filtering."
    )
  }

  valid_drug_gene <- !is.na(drug_gene) &
    drug_gene != ""

  if (!all(valid_drug_gene)) {
    mat.drug <- mat.drug[
      valid_drug_gene,
      ,
      drop = FALSE
    ]

    drug_gene <- drug_gene[
      valid_drug_gene
    ]
  }

  if (length(drug_gene) == 0L) {
    stop(
      "No valid genes remain in the drug matrix."
    )
  }

  ## Collapse duplicated disease gene symbols.
  ## The minimum value retains the strongest evidence for each symbol.
  duplicated_trait_n <- sum(
    duplicated(trait_gene$Gene)
  )

  if (duplicated_trait_n > 0L) {
    cat(sprintf(
      paste0(
        "[INFO] Collapsing %d duplicated disease-gene rows ",
        "using the minimum Pvalue.\n"
      ),
      duplicated_trait_n
    ))

    trait_gene <- aggregate(
      Pvalue ~ Gene,
      data = trait_gene,
      FUN = function(x) {
        min(x, na.rm = TRUE)
      }
    )
  }

  ## The tested drug matrix contains no duplicated genes.
  ## Stop with a clear diagnostic if another matrix is used later.
  duplicated_drug_n <- sum(
    duplicated(drug_gene)
  )

  if (duplicated_drug_n > 0L) {
    stop(sprintf(
      paste0(
        "The drug matrix contains %d duplicated gene rows. ",
        "The current input matrix is expected to contain unique genes."
      ),
      duplicated_drug_n
    ))
  }

  common_genes <- intersect(
    trait_gene$Gene,
    drug_gene
  )

  common_genes <- sort(
    unique(common_genes)
  )

  cat(sprintf(
    paste0(
      "[INFO] Gene overlap: disease=%d, ",
      "drug_matrix=%d, common=%d\n"
    ),
    nrow(trait_gene),
    length(drug_gene),
    length(common_genes)
  ))

  if (length(common_genes) < 10L) {
    stop(sprintf(
      paste0(
        "Only %d common genes were found between ",
        "disease and drug data."
      ),
      length(common_genes)
    ))
  }

  ## Use match to guarantee identical ordering.
  disease_index <- match(
    common_genes,
    trait_gene$Gene
  )

  drug_index <- match(
    common_genes,
    drug_gene
  )

  if (anyNA(disease_index) || anyNA(drug_index)) {
    stop(
      "Internal gene matching error after defining common genes."
    )
  }

  aligned_disease_gene <- trait_gene$Gene[
    disease_index
  ]

  aligned_drug_gene <- drug_gene[
    drug_index
  ]

  if (!all(
    aligned_disease_gene == aligned_drug_gene
  )) {
    mismatch <- which(
      aligned_disease_gene != aligned_drug_gene
    )

    stop(sprintf(
      paste0(
        "Gene alignment failed after match; ",
        "the first mismatch is at position %d."
      ),
      mismatch[1L]
    ))
  }

  ## Protect qnorm from values exactly equal to zero or one.
  pv <- trait_gene$Pvalue[
    disease_index
  ]

  pv <- pmax(
    pmin(pv, 1 - 1e-12),
    1e-300
  )

  zval_trait <- qnorm(pv)

  if (any(!is.finite(zval_trait))) {
    stop(
      "Non-finite z-scores remain after numerical protection."
    )
  }

  drug_columns <- setdiff(
    colnames(mat.drug),
    "Gene"
  )

  no_drugs <- length(drug_columns)

  if (no_drugs == 0L) {
    stop(
      "No drug-membership columns were found in mat.drug."
    )
  }

  cat(sprintf(
    "[INFO] Starting enrichment tests for %d drugs.\n",
    no_drugs
  ))

  pvals <- rep(
    999,
    no_drugs
  )

  for (i in seq_len(no_drugs)) {

    raw_membership <- mat.drug[[drug_columns[i]]][drug_index]

    if (is.factor(raw_membership)) {
      raw_membership <- as.character(
        raw_membership
      )
    }

    membership_numeric <- suppressWarnings(
      as.numeric(raw_membership)
    )

    genes_for_drug <- as.integer(
      !is.na(membership_numeric) &
        membership_numeric > 0
    )

    n_associated <- sum(
      genes_for_drug == 1L
    )

    n_unassociated <- sum(
      genes_for_drug == 0L
    )

    ## At least five associated genes are required.
    if (n_associated < 5L ||
        n_unassociated < 2L) {
      pvals[i] <- 999
      next
    }

    associated_z <- zval_trait[
      genes_for_drug == 1L
    ]

    unassociated_z <- zval_trait[
      genes_for_drug == 0L
    ]

    ## This is equivalent to the original formula:
    ## t.test(zval_trait ~ drug_membership,
    ##        alternative = "greater")
    ##
    ## Because smaller FDR values produce smaller qnorm values,
    ## enrichment corresponds to higher values in the non-member group.
    pvals[i] <- tryCatch(
      t.test(
        unassociated_z,
        associated_z,
        alternative = "greater"
      )$p.value,
      error = function(e) {
        warning(sprintf(
          "t-test failed for drug '%s': %s",
          drug_columns[i],
          conditionMessage(e)
        ))
        999
      }
    )

    if (i %% 1000L == 0L) {
      cat(sprintf(
        "[INFO] Completed %d of %d drug tests.\n",
        i,
        no_drugs
      ))
    }
  }

  cat(sprintf(
    "[INFO] Completed all %d drug tests.\n",
    no_drugs
  ))

  data.frame(
    drug = drug_columns,
    t.test.p.one.sided_trait = pvals,
    stringsAsFactors = FALSE
  )
}

## ============================================================
## 6. MULTIPLE-TESTING CORRECTION
## ============================================================

add_fdr <- function(x, col_in, col_out) {

  if (!col_in %in% colnames(x)) {
    stop(sprintf(
      "Input column '%s' was not found.",
      col_in
    ))
  }

  x[[col_out]] <- 999

  keep <- which(
    !is.na(x[[col_in]]) &
      is.finite(x[[col_in]]) &
      x[[col_in]] != 999
  )

  if (length(keep) > 0L) {
    x[[col_out]][keep] <- p.adjust(
      x[[col_in]][keep],
      method = "fdr"
    )
  }

  x
}

## ============================================================
## 7. LOAD AND VALIDATE SHARED RESOURCES
## ============================================================

if (!file.exists(cfg$lib_entrez)) {
  stop(sprintf(
    "MAGMA annotation file does not exist: %s",
    cfg$lib_entrez
  ))
}

if (!file.exists(cfg$drug_mat)) {
  stop(sprintf(
    "Drug matrix file does not exist: %s",
    cfg$drug_mat
  ))
}

if (cfg$run_stage_a &&
    !file.exists(cfg$fn_atc)) {
  stop(sprintf(
    "ATC function file does not exist: %s",
    cfg$fn_atc
  ))
}

if (cfg$run_stage_a &&
    !file.exists(cfg$atc_lists)) {
  stop(sprintf(
    "ATC list file does not exist: %s",
    cfg$atc_lists
  ))
}

lib_entrez <- fread(
  cfg$lib_entrez,
  data.table = FALSE,
  showProgress = FALSE
)

if (!all(
  c("ID", "GeneName") %in%
    colnames(lib_entrez)
)) {
  stop(
    "The MAGMA annotation library must contain ID and GeneName."
  )
}

loaded_drug_objects <- load(
  cfg$drug_mat
)

if (!"mat.drug" %in% loaded_drug_objects ||
    !exists("mat.drug")) {
  stop(
    "The drug RData file did not create an object named mat.drug."
  )
}

if (!"Gene" %in% colnames(mat.drug)) {
  stop(
    "The loaded mat.drug object does not contain a Gene column."
  )
}

cat(sprintf(
  paste0(
    "[INFO] Drug matrix loaded: rows=%d, ",
    "columns=%d, drugs=%d\n"
  ),
  nrow(mat.drug),
  ncol(mat.drug),
  ncol(mat.drug) - 1L
))

cat(sprintf(
  "[INFO] Duplicated drug-matrix genes: %d\n",
  sum(
    duplicated(
      trimws(
        as.character(mat.drug$Gene)
      )
    )
  )
))

## ============================================================
## 8. STAGE A: ENRICHMENT AND ATC CATEGORY
## ============================================================

if (cfg$run_stage_a) {

  source(cfg$fn_atc)

  loaded_atc_objects <- load(
    cfg$atc_lists
  )

  cat(sprintf(
    "[INFO] Loaded ATC objects: %s\n",
    paste(
      loaded_atc_objects,
      collapse = ", "
    )
  ))

  for (target in cfg$targets) {

    for (level in cfg$levels) {

      gene_raw <- tryCatch(
        load_gene(
          target,
          level,
          lib_entrez
        ),
        error = function(e) {
          cat(sprintf(
            "[SKIP] %s | %s : %s\n",
            target,
            level,
            conditionMessage(e)
          ))
          NULL
        }
      )

      if (is.null(gene_raw)) {
        next
      }

      for (tl in names(cfg$type_map)) {

        tcol <- cfg$type_map[[tl]]

        cat(sprintf(
          "\n[STAGE A] %s | %s | %s (%s)\n",
          target,
          level,
          tl,
          tcol
        ))

        trait_gene <- tryCatch(
          resolve_pvalue(
            gene_raw,
            tcol
          ),
          error = function(e) {
            cat(sprintf(
              "[ERROR] %s | %s | %s | resolve_pvalue: %s\n",
              target,
              level,
              tl,
              conditionMessage(e)
            ))
            NULL
          }
        )

        if (is.null(trait_gene)) {
          next
        }

        selected_pvalue_column <- attr(
          trait_gene,
          "selected_pvalue_column"
        )

        cat(sprintf(
          "[INFO] Using p-value/FDR column: %s\n",
          selected_pvalue_column
        ))

        ## Per-drug enrichment.
        DEA <- tryCatch(
          drug_enrichment_test(
            trait_gene,
            mat.drug
          ),
          error = function(e) {
            cat(sprintf(
              "[ERROR] %s | %s | %s | enrichment: %s\n",
              target,
              level,
              tl,
              conditionMessage(e)
            ))
            NULL
          }
        )

        if (is.null(DEA)) {
          next
        }

        DEA <- DEA[
          order(
            DEA$t.test.p.one.sided_trait
          ),
          ,
          drop = FALSE
        ]

        DEA <- add_fdr(
          DEA,
          "t.test.p.one.sided_trait",
          "p.adj"
        )

        rdata_output <- enrich_path(
          target,
          level,
          tl,
          "rdata"
        )

        csv_output <- enrich_path(
          target,
          level,
          tl,
          "csv"
        )

        save(
          DEA,
          file = rdata_output
        )

        fwrite(
          DEA,
          csv_output,
          sep = ","
        )

        cat(sprintf(
          "[INFO] DEA output written: %s\n",
          csv_output
        ))

        ## ATC category enrichment.
        atc_res <- tryCatch(
          ATC_enrichment(
            csv_output,
            "t.test.p.one.sided_trait",
            "drug"
          ),
          error = function(e) {
            cat(sprintf(
              "[ERROR] %s | %s | %s | ATC enrichment: %s\n",
              target,
              level,
              tl,
              conditionMessage(e)
            ))
            NULL
          }
        )

        if (is.null(atc_res)) {
          next
        }

        atc_res <- add_fdr(
          atc_res,
          "pval.oneSamp.t",
          "p.adj.one"
        )

        atc_res <- add_fdr(
          atc_res,
          "pval.twoSamp.t",
          "p.adj.two"
        )

        atc_res <- atc_res[
          order(atc_res$pval.twoSamp.t),
          ,
          drop = FALSE
        ]

        category_output <- cat_path(
          target,
          level,
          tl
        )

        fwrite(
          atc_res,
          category_output,
          sep = ","
        )

        cat(sprintf(
          "[INFO] ATC category output written: %s\n",
          category_output
        ))
      }
    }
  }
}

## ============================================================
## 9. STAGE B: SUMMARY TABLES
## ============================================================

if (cfg$run_stage_b) {

  long_list <- vector(
    "list",
    0L
  )

  for (target in cfg$targets) {

    for (level in cfg$levels) {

      for (tl in names(cfg$type_map)) {

        fp <- cat_path(
          target,
          level,
          tl
        )

        if (!file.exists(fp)) {
          next
        }

        d <- fread(
          fp,
          data.table = FALSE,
          showProgress = FALSE
        )

        d$target <- target
        d$level <- level
        d$type <- tl

        long_list[[length(long_list) + 1L]] <- d
      }
    }
  }

  if (length(long_list) == 0L) {

    cat(
      "[WARNING] No category result files were found for Stage B.\n"
    )

  } else {

    long <- rbindlist(
      long_list,
      use.names = TRUE,
      fill = TRUE
    )

    all_long_output <- file.path(
      dir_summ,
      "category_all_long.csv"
    )

    fwrite(
      long,
      all_long_output,
      sep = ","
    )

    cat(sprintf(
      "[INFO] Combined category output written: %s\n",
      all_long_output
    ))

    make_subset <- function(codes_df, fname) {

      if (!"level3_codes" %in% colnames(long)) {
        warning(
          "The combined category table does not contain level3_codes."
        )
        return(invisible(NULL))
      }

      sub <- long[
        long$level3_codes %in%
          codes_df$level3_codes,
        ,
        drop = FALSE
      ]

      sub <- left_join(
        as.data.frame(sub),
        codes_df,
        by = "level3_codes"
      )

      output_file <- file.path(
        dir_summ,
        fname
      )

      fwrite(
        sub,
        output_file,
        sep = ","
      )

      cat(sprintf(
        "[INFO] Focused summary output written: %s\n",
        output_file
      ))

      invisible(NULL)
    }

    make_subset(
      cfg$drug5,
      "drug5_detail.csv"
    )

    make_subset(
      cfg$drug7,
      "drug7_detail.csv"
    )
  }
}

cat("\n[INFO] DEA v7 pipeline completed.\n")

q("no")
