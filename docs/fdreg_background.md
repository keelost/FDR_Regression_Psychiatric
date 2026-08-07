# FDRreg background

FDRreg is an independently developed statistical method and R package. The
authors of this repository did not create or publish that package. This project
uses FDRreg as a third-party statistical dependency and contributes only the
psychiatric GWAS application workflow and its study-specific preprocessing,
integration, validation, and reporting code. Users should cite the upstream
FDRreg method/package separately from this repository.

FDRreg models the probability that each test is non-null as a function of
external covariates. In this project, signed GWAS z-scores are the responses;
genetically correlated traits and biological annotations are the covariates.

The SNP branch harmonises alleles first. For traits with sample overlap, LDSC
genetic-covariance intercepts are assembled into a covariance matrix and the
corresponding z-scores are decorrelated before FDRreg. Traits without sample
overlap retain their harmonised z-scores.

Gene-level branches use MAGMA, MetaXcan, or SMultiXcan gene statistics and can
fit theoretical-null, empirical-null, LASSO, and biological-annotation models
according to the validated script options. The source scripts save model
metadata, FDR values, q-values, and contribution summaries for auditing.

The validated R files under `scripts/real` call the upstream FDRreg package and
add the study-specific data preparation and analysis logic. The ordinary Bash
runner supplies paths and controls stage ordering.
