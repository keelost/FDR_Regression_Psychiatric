# Pipeline overview

This repository starts from cleaned GWAS summary statistics or from prepared
target directories. Initial download and cleaning remain external because the
original GWAS cannot be redistributed. DIST imputation from cleaned inputs is
implemented; later harmonisation and overlap construction occur inside the SNP
FDRreg code.

## Execution flow

The ordinary Bash entry point is `scripts/run_pipeline.sh`. It reads
`config/project.env`, `config/targets.tsv`, `config/regions.txt`, and
`config/versions.tsv`.

1. `imputation`: convert cleaned target/library GWAS files to DIST inputs, run
   chromosomes 1-22, retain imputed variants at INFO >= 0.5, map non-rs IDs,
   remove duplicates, and create `*.impute.map.txt`.
2. `ldsc`: calculate heritability and pairwise genetic correlation from
   already munged `*.sumstats.gz` inputs, then build the source-labelled LDSC
   CSV used for decorrelation.
3. `snp`: run SNP-level FDRreg and construct each target's overlap data.
4. `magma`: run gene-level FDRreg from prepared MAGMA inputs and outputs.
5. `metaxcan`: run tissue-level FDRreg for v7, v8, or both.
6. `smultixcan`: run multi-tissue FDRreg for v7, v8, or both.
7. `sensitivity`: run rg and absolute-z sensitivity analyses, keeping v7 and
   v8 results separate.
8. `enrichment`: run versioned drug and KEGG/GO enrichment analyses.
9. `validation`: run PPV and temporal validation.
10. `summary`: generate result summaries and risk-locus tables.

The `all` stage runs stages 3 through 10. Imputation and LDSC are intentionally
separate because they are expensive upstream preparations that are commonly
run once and reused.

## Imputation contract

`clear.data` means a downloaded GWAS that has already been normalized to the
project column names, filtered to INFO >= 0.6 when INFO is available, stripped
of duplicate SNPs and InDels, and assigned validated SNP identifiers. The
preparation script checks and enforces these conditions where possible and
computes `z = beta / se` when no z-score column exists.

DIST 1.0.0 is run against the ancestry-matched 1000 Genomes reference listed in
`config/imputation_refs.tsv`. Per-chromosome results are merged at INFO >= 0.5.
Variants without an rsID are mapped by chromosome and position using
`FDRREG_DIST_VARIANT_MAP`, then the final data are deduplicated by rsID.

## Version isolation

MetaXcan and SMultiXcan versions never share output directories:

| Analysis | v8 | v7 |
|---|---|---|
| MetaXcan input | `06.metaxcan` | `06.metaxcan_v7` |
| MetaXcan FDRreg | `07.metaxcan_fdrreg` | `07.metaxcan_fdrreg_v7` |
| SMultiXcan input | `08.smultixcan_output` | `08.smultixcan_output_v7` |
| SMultiXcan FDRreg | `09.smultixcan_fdrreg` | `09.smultixcan_fdrreg_v7` |

## LDSC decorrelation

For traits with sample overlap, LDSC covariance intercepts form a correlation
matrix. The SNP scripts apply its inverse square root to signed z-scores before
FDRreg. Traits listed in `traits_no` retain their harmonised signed z-scores.
The exact expected files and columns are documented in `input_inventory.md`.

## Analysis branches

- SNP-level FDRreg identifies associated variants from target and covariate
  z-scores.
- MAGMA FDRreg aggregates variant evidence at gene level and can include
  biological annotations.
- MetaXcan FDRreg operates independently across the 13 configured brain
  regions.
- SMultiXcan FDRreg combines evidence across brain regions.
- Sensitivity, ablation, enrichment, validation, diagnostics, and summary
  scripts operate on these versioned outputs.

Historical scripts with inconsistent paths or superseded logic are retained
under `scripts/real/legacy` and `scripts/ldsc/legacy`; the runner does not call
them.
