#!/usr/bin/env Rscript
#******************************#
#-----FDR_REG: MAGMA Analysis---#
#******************************#
# R-4.0.2 or later
# Optimized version with argument parsing

#----Load packages and set directories----#
rm(list = ls())

# Record start time for execution time calculation
start_time <- Sys.time()

# Load required packages
suppressPackageStartupMessages({
  library(FDRreg)
  library(data.table)
  library(powerplus)
  library(stringr)
  library(dplyr)
  library(glmnet)
  library(HelpersMG)
  library(optparse)
})

# Set global seed once
set.seed(100)

# ---- Argument Parsing ---- #
option_list <- list(
  make_option(c("-t", "--trait"), type = "character", default = NULL,
              help = "Trait name to analyze (required). [REQUIRED]"),
  make_option(c("-v", "--variables"), type = "character", default = NULL,
              help = "Path to text file with variable traits (one per line) or comma-separated list. [REQUIRED]"),
  make_option(c("-b", "--bio"), type = "character", 
              default = Sys.getenv("FDRREG_BIO_ENTREZ", ""),
              help = "Path to biological annotation file. [default: %default]"),
  make_option(c("-l", "--lasso"), type = "logical", default = TRUE,
              help = "Perform Lasso selection on biological annotations. [default: %default]"),
  make_option(c("-o", "--output"), type = "character", 
              default = Sys.getenv("FDRREG_RESULTS_DIR", ""),
              help = "Base output directory. [default: %default]"),
  make_option(c("-s", "--seed"), type = "integer", default = 100,
              help = "Random seed for reproducibility. [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$trait)) {
  stop("--trait is required. See: Rscript script.R --help")
}

if (is.null(opt$variables)) {
  stop("--variables is required. See: Rscript script.R --help")
}

# Update seed if provided
if (!is.null(opt$seed)) {
  set.seed(opt$seed)
}

# ---- Configuration ---- #
trait_name <- opt$trait
bio_annotation_path <- opt$bio
perform_lasso <- opt$lasso
base_output_path <- opt$output

# Parse variables argument (file or comma-separated string)
if (file.exists(opt$variables)) {
  # Read from file
  variable_traits <- fread(opt$variables, header = FALSE)$V1
  cat(sprintf("Loaded %d variable traits from file: %s\n", 
              length(variable_traits), opt$variables))
} else {
  # Parse as comma-separated string
  variable_traits <- strsplit(opt$variables, ",")[[1]]
  variable_traits <- trimws(variable_traits[variable_traits != ""])
  cat(sprintf("Parsed %d variable traits from string\n", length(variable_traits)))
}

# ---- Directory Setup ---- #
trait_dir <- file.path(base_output_path, trait_name)
magma_dir <- file.path(trait_dir, "04.magma_output")
output_dir <- file.path(trait_dir, "05.magma_fdrreg")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n========================================\n"))
cat(sprintf("Processing trait: %s\n", trait_name))
cat(sprintf("Variables: %s\n", paste(variable_traits, collapse = ", ")))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("========================================\n"))

# ---- Load MAGMA Data ---- #
# Load target trait's MAGMA output
target_file <- file.path(magma_dir, paste0(trait_name, ".genes.out"))
if (!file.exists(target_file)) {
  stop(sprintf("Target MAGMA file not found: %s", target_file))
}

target_gene <- fread(target_file)
target_gene_order <- target_gene[order(target_gene$GENE), ]

# Load variable traits' MAGMA outputs
variable_files <- file.path(magma_dir, paste0(variable_traits, ".genes.out"))
missing_files <- variable_files[!file.exists(variable_files)]

if (length(missing_files) > 0) {
  warning(sprintf("Missing %d variable files:\n%s", 
                  length(missing_files), paste(missing_files, collapse = "\n")))
  # Remove missing traits from list
  existing_indices <- file.exists(variable_files)
  variable_traits <- variable_traits[existing_indices]
  variable_files <- variable_files[existing_indices]
}

if (length(variable_traits) == 0) {
  stop("No valid variable traits found")
}

# Load variable data
files <- lapply(variable_files, fread)
names(files) <- variable_traits

# Prepare combined z-scores
fileout_order <- files[[1]][order(files[[1]]$GENE), ]
z_combine <- abs(qnorm(fileout_order$P / 2))

for (i in 2:length(files)) {
  geneout <- files[[i]]
  geneout_order <- geneout[order(geneout$GENE), ]
  z_combine <- cbind(z_combine, abs(qnorm(geneout_order$P / 2)))
}
colnames(z_combine) <- variable_traits

# ---- FDRreg Analysis ---- #
# Extract target z-scores with random sign
target_z <- abs(qnorm(target_gene_order$P / 2)) * sign(rnorm(length(target_gene_order$P)))

# Generate gene list for annotation
gene_list <- target_gene_order$GENE
fwrite(as.data.frame(gene_list), 
       file.path(output_dir, paste0(trait_name, ".gene.list.txt")))

# Perform FDR regression (both null types)
fdr_target_theore <- FDRreg(target_z, z_combine, nulltype = 'theoretical', method = 'pr')
fdr_target_empiri <- FDRreg(target_z, z_combine, nulltype = 'empirical', method = 'pr')

# ---- Assessment of FDRreg ---- #
features_se <- SEfromHessian(fdr_target_theore$model$hessian)
features_coef <- fdr_target_theore$model$coef
features_z_score <- features_coef[c(2:length(features_coef))] / features_se[c(2:length(features_se))]
features_pvalue <- 2 * pnorm(abs(features_z_score), lower.tail = FALSE)

assessment <- data.frame(
  feature = colnames(z_combine),
  p1 = features_pvalue,
  beta1 = features_coef[c(2:length(features_coef))],
  se1 = features_se[c(2:length(features_se))],
  stringsAsFactors = FALSE
)

# ---- Biological Annotation Integration ---- #
# Load biological annotation data
bio_annotation <- fread(bio_annotation_path)

# Merge target gene data with biological annotation
colnames(target_gene_order) <- gsub("GENE", "ID", colnames(target_gene_order))
target_bio <- left_join(target_gene_order[, 1:3], bio_annotation, by = 'ID')
target_bio[is.na(target_bio)] <- 0

# Prepare bio features
bio_cols <- grep("^(dise|bio|expre|pathway|tfbs|ASDDenovo|SCZDenovo|scz|depre|bpd|adhd|asd)", 
                 names(target_bio), value = TRUE)

if (length(bio_cols) > 0) {
  bio_features <- target_bio[, bio_cols, with = FALSE]
  z_combine_bio <- as.matrix(cbind(abs(z_combine), abs(bio_features)))
} else {
  z_combine_bio <- as.matrix(abs(z_combine))
  bio_features <- NULL
}

# FDRreg with biological annotations
fdr_target_theore_bio <- FDRreg(target_z, z_combine_bio, nulltype = 'theoretical', method = 'pr')
fdr_target_empiri_bio <- FDRreg(target_z, z_combine_bio, nulltype = 'empirical', method = 'pr')

# ---- Lasso Selection (if requested) ---- #
lasso_performed <- FALSE
z_combine_lasso_bio <- NULL

# Initialize Lasso model objects to prevent save errors
fdr_target_theore_bio_lasso <- NULL
fdr_target_empiri_bio_lasso <- NULL

if (perform_lasso && !is.null(bio_features) && ncol(bio_features) > 0) {
  tryCatch({
    cat("Performing Lasso selection on biological features...\n")
    
    # Set parallel = FALSE to avoid backend warning
    lasso_cv <- cv.glmnet(as.matrix(bio_features), abs(target_z), 
                         family = 'gaussian', nlambda = 50, alpha = 1, 
                         standardize = TRUE, parallel = FALSE)
    coef_list_min <- as.matrix(coef(lasso_cv, s = "lambda.min"))
    coef_list_min <- unlist(coef_list_min)[-1]
    
    # Filter features with non-zero coefficients
    non_zero_features <- which(coef_list_min != 0)
    if (length(non_zero_features) > 0) {
      z_combine_lasso_bio <- as.matrix(cbind(abs(z_combine), 
                                             abs(bio_features[, non_zero_features, with = FALSE])))
      lasso_performed <- TRUE
      cat(sprintf("Lasso selected %d biological features\n", length(non_zero_features)))
      
      # Save selected features
      selected_features <- data.frame(
        feature = names(bio_features)[non_zero_features],
        coefficient = coef_list_min[non_zero_features]
      )
      fwrite(selected_features, 
             file.path(output_dir, paste0(trait_name, ".lasso_selected_features.csv")))
      
      # Run FDR regression with Lasso-selected features
      fdr_target_theore_bio_lasso <- FDRreg(target_z, z_combine_lasso_bio, 
                                            nulltype = 'theoretical', method = 'pr')
      fdr_target_empiri_bio_lasso <- FDRreg(target_z, z_combine_lasso_bio, 
                                            nulltype = 'empirical', method = 'pr')
    } else {
      cat("Lasso selected no biological features\n")
    }
  }, error = function(e) {
    cat(sprintf("Lasso failed: %s\n", conditionMessage(e)))
  })
}

# ---- Assessment for Biological Annotations ---- #
if (fdr_target_theore_bio$FDR[1] != "NA") {
  features_se_bio <- SEfromHessian(fdr_target_theore_bio$model$hessian)
  features_coef_bio <- fdr_target_theore_bio$model$coef
  features_z_score_bio <- features_coef_bio[c(2:length(features_coef_bio))] / 
                         features_se_bio[c(2:length(features_se_bio))]
  features_pvalue_bio <- 2 * pnorm(abs(features_z_score_bio), lower.tail = FALSE)
  
  assessment2 <- data.frame(
    feature = colnames(z_combine_bio),
    p2 = features_pvalue_bio,
    beta2 = features_coef_bio[c(2:length(features_coef_bio))],
    se2 = features_se_bio[c(2:length(features_se_bio))],
    stringsAsFactors = FALSE
  )
  
  assessment_sum <- merge(assessment, assessment2, by = 'feature', all = TRUE)
} else {
  assessment_sum <- assessment
  assessment_sum$p2 <- NA
  assessment_sum$beta2 <- NA
  assessment_sum$se2 <- NA
}

# ---- Complete FDR Results ---- #
target_gene_order$FDR.the <- fdr_target_theore$FDR
target_gene_order$FDR.emp <- fdr_target_empiri$FDR
target_gene_order$qval <- p.adjust(target_gene_order$P, method = "fdr", n = length(target_gene_order$P))

if (fdr_target_theore_bio$FDR[1] != "NA") {
  target_gene_order$bio.FDR.the <- fdr_target_theore_bio$FDR
  target_gene_order$bio.FDR.emp <- fdr_target_empiri_bio$FDR
} else {
  target_gene_order$bio.FDR.the <- NA
  target_gene_order$bio.FDR.emp <- NA
}

# ---- Summary Statistics ---- #
thresholds <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 
                0.001, 5e-04, 5e-06, 5e-08)

# Helper function to count significant genes
count_sig <- function(fdr_vec, th_vec) {
  vapply(th_vec, function(th) sum(fdr_vec < th, na.rm = TRUE), integer(1))
}

qval <- count_sig(target_gene_order$qval, thresholds)
fdr_the <- count_sig(target_gene_order$FDR.the, thresholds)
fdr_emp <- count_sig(target_gene_order$FDR.emp, thresholds)

bio_fdr_the <- bio_fdr_emp <- bio_fdr_the_lasso <- bio_fdr_emp_lasso <- 
  rep(NA, length(thresholds))

if (fdr_target_theore_bio$FDR[1] != "NA") {
  bio_fdr_the <- count_sig(target_gene_order$bio.FDR.the, thresholds)
  bio_fdr_emp <- count_sig(target_gene_order$bio.FDR.emp, thresholds)
}

# ---- Lasso FDR Analysis (if applicable) ---- #
if (lasso_performed && !is.null(z_combine_lasso_bio)) {
  # Assessment for lasso model
  features_se_lasso <- SEfromHessian(fdr_target_theore_bio_lasso$model$hessian)
  features_coef_lasso <- fdr_target_theore_bio_lasso$model$coef
  features_z_score_lasso <- features_coef_lasso[c(2:length(features_coef_lasso))] / 
                           features_se_lasso[c(2:length(features_se_lasso))]
  features_pvalue_lasso <- 2 * pnorm(abs(features_z_score_lasso), lower.tail = FALSE)
  
  assessment3 <- data.frame(
    feature = colnames(z_combine_lasso_bio),
    p3 = features_pvalue_lasso,
    beta3 = features_coef_lasso[c(2:length(features_coef_lasso))],
    se3 = features_se_lasso[c(2:length(features_se_lasso))],
    stringsAsFactors = FALSE
  )
  
  assessment_sum <- merge(assessment_sum, assessment3, by = 'feature', all = TRUE)
  
  target_gene_order$bio.FDR.the.lasso <- fdr_target_theore_bio_lasso$FDR
  target_gene_order$bio.FDR.emp.lasso <- fdr_target_empiri_bio_lasso$FDR
  
  bio_fdr_the_lasso <- count_sig(target_gene_order$bio.FDR.the.lasso, thresholds)
  bio_fdr_emp_lasso <- count_sig(target_gene_order$bio.FDR.emp.lasso, thresholds)
} else {
  assessment_sum$p3 <- NA
  assessment_sum$beta3 <- NA
  assessment_sum$se3 <- NA
  
  target_gene_order$bio.FDR.the.lasso <- NA
  target_gene_order$bio.FDR.emp.lasso <- NA
}

# ---- Save Results ---- #
# Save gene-level results (with Lasso columns if available)
fwrite(as.data.frame(target_gene_order), 
       file.path(output_dir, paste0(trait_name, ".gene.bio.fdrreg.txt")), sep = ' ')

# Save assessment/contribution table
fwrite(as.data.frame(assessment_sum), 
       file.path(output_dir, paste0(trait_name, ".contribution.csv")), sep = ',')

# Save summary statistics
sum_data <- data.frame(
  threshold = thresholds,
  qval = qval,
  fdr_the = fdr_the,
  fdr_emp = fdr_emp,
  bio_fdr_the = bio_fdr_the,
  bio_fdr_emp = bio_fdr_emp,
  bio_fdr_the_lasso = bio_fdr_the_lasso,
  bio_fdr_emp_lasso = bio_fdr_emp_lasso
)
fwrite(as.data.frame(sum_data), 
       file.path(output_dir, paste0(trait_name, ".summary.csv")), sep = ',')

# Save FDRreg models - Conditionally include Lasso models if they exist
objects_to_save <- c("fdr_target_theore", "fdr_target_empiri", 
                     "fdr_target_theore_bio", "fdr_target_empiri_bio")
if (!is.null(fdr_target_theore_bio_lasso) && !is.null(fdr_target_empiri_bio_lasso)) {
  objects_to_save <- c(objects_to_save, 
                       "fdr_target_theore_bio_lasso", "fdr_target_empiri_bio_lasso")
}
save(list = objects_to_save, 
     file = file.path(output_dir, paste0(trait_name, ".fdrreg_models.RData")))

# ---- Summary Output ---- #
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\n---------- Summary ----------\n")
cat("All gene number:", length(gene_list), "\n")
cat("Lasso performed:", perform_lasso, "\n")
cat("Lasso applied:", lasso_performed, "\n\n")

cat("Significant genes at various thresholds:\n")
print(sum_data)

cat(sprintf("\nResults saved to: %s\n", output_dir))
cat(sprintf("Total execution time: %.2f seconds\n", execution_time))

# Exit gracefully
q("no")
