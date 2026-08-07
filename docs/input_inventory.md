# Server Input and Upload Inventory

This project does not contain GWAS data or licensed reference resources. Copy
`config/project.env.example` to `config/project.env`, then point each variable
to the corresponding server file or directory.

## 1. GWAS access and cleaned inputs

The original GWAS summary statistics are not distributed with this repository
because their providers impose data-use and redistribution restrictions. Apply
for or download each dataset from the source listed in the study's
supplementary tables.

Before DIST, create one `<trait>.clear.txt` (optionally gzip-compressed) for
each required target and auxiliary trait. A `clear.data` file is the downloaded
GWAS after:

- retaining INFO >= 0.6 where an INFO field is supplied;
- removing duplicate SNP identifiers and InDels;
- assigning or validating the correct SNP identifier;
- normalizing columns to `snpid`, `chr`, `bpos`, `a1`, `a2`, plus either `z`
  or both `beta` and `se`.

Place these files in `FDRREG_CLEAR_TARGET_DIR` and
`FDRREG_CLEAR_LIBRARY_DIR`. They must be available on the server but must not
be committed to Git.

## 2. DIST resources and generated SNP inputs

To regenerate `*.impute.map.txt`, provide:

| Configuration | Required resource |
|---|---|
| `FDRREG_DIST_BIN` | DIST 1.0.0 executable |
| `FDRREG_DIST_REF_EUR` | EUR files named `chr1.1kg.eur.gz` through `chr22.1kg.eur.gz` |
| `FDRREG_DIST_REF_EAS` | EAS files named `chr1.1kg.eas.gz` through `chr22.1kg.eas.gz` |
| `FDRREG_DIST_REF_PATTERN_EUR/EAS` | Override reference filenames with `{chr}` and `{panel}` placeholders |
| `FDRREG_DIST_VARIANT_MAP` | Variant table with chromosome, position, and rsID columns |
| `FDRREG_IMPUTE_WORK_DIR` | Writable location for DIST inputs and chromosome outputs |
| `config/imputation_refs.tsv` | Trait-to-reference-panel mapping; verify it against each GWAS ancestry |

The runner applies INFO >= 0.6 to cleaned inputs when an INFO column is
present and INFO >= 0.5 to DIST outputs. Results are mapped to rsIDs,
deduplicated, and written as `<trait>.impute.map.txt` under
`FDRREG_INPUT_DIR`.

Generated files need not be uploaded when DIST will be rerun. If the server
will start directly at SNP FDRreg, upload the final `*.impute.map.txt` files.

## 3. Required for SNP FDRreg

Upload or make available:

- One `<trait>.impute.map.txt` file for every target and predictor listed in
  `config/targets.tsv`.
- Required columns: `snpid`, `a1`, `a2`, `z`, and `pval`. A signed `beta`
  column is retained when present.
- The prepared files are placed together under `FDRREG_INPUT_DIR`.
- `genetic_correlation_results.with_source.csv`, configured through
  `FDRREG_LDSC_RESULTS`.

The LDSC CSV must contain:

```text
p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval
```

`p1` and `p2` must use the source prefixes produced by the included LDSC
parser, for example `T:scz2012` and `L:adhd2019`. SNP decorrelation uses
`gcov_intercept`; the other fields are retained for auditing.

## 4. Files needed to build the LDSC CSV

Set the following paths if the LDSC stage will be rerun:

| Configuration | Required resource |
|---|---|
| `FDRREG_LDSC_TARGET_DIR` | Target `*.clear.txt.sumstats.gz` files |
| `FDRREG_LDSC_LIBRARY_DIR` | Predictor/library `*.sumstats.gz` files |
| `FDRREG_LDSC_BIN` | LDSC `ldsc.py` |
| `FDRREG_LDSC_REF` | Per-chromosome reference LD scores |
| `FDRREG_LDSC_WEIGHTS` | Per-chromosome regression weights, usually the same EUR LD directory |
| `FDRREG_LDSC_OUT_DIR` | Writable directory for logs and CSV results |

If raw summary statistics must be munged again, also provide
`munge_sumstats.py` and the HapMap3 `w_hm3.snplist`. The current ordinary LDSC
entry point expects already munged files; the original missing-target munging
scripts are retained under `scripts/ldsc/legacy` for reference.

Generated LDSC files do not need to be uploaded when they can be rebuilt:

- `heritability_results.csv`
- `genetic_correlation_results.csv`
- `genetic_correlation_results.with_source.csv`
- `logs/*.log`

Only `genetic_correlation_results.with_source.csv` is consumed directly by the
SNP FDRreg stage.

## 5. MAGMA and biological annotations

For MAGMA execution, provide the MAGMA executable, a population-matched PLINK
reference panel, gene annotation file, and a trait/sample-size table. The
gene-level FDRreg stage expects, per target:

```text
<base>/<target>/04.magma_output/<trait>.genes.out
```

The exact frozen matrices used in the study are supplied in
`data/annotations`:

- `magma-library-all.csv`: Entrez, Ensembl, gene name, and all annotations
- `magma-library-uniq-entrez.csv`: MAGMA-ready Entrez-keyed table; set as
  `FDRREG_BIO_ENTREZ`
- `magma-library-uniq-ensembl.csv`: TWAS-ready Ensembl-keyed table; set as
  `FDRREG_BIO_ENSEMBL`

These are manually frozen research annotations derived from DAVID 6.8,
Open Targets (February 2020 extraction), and denovo-db 1.6.1. Because DAVID and
Open Targets change over time, current queries will not reproduce the reported
matrix exactly.

## 6. MetaXcan and SMultiXcan resources

Keep the two model releases separate. Download MetaXcan code, GTEx prediction
models, covariance files, and annotations from the official project:
<https://github.com/hakyimlab/MetaXcan>.

### GTEx v8

- Mashr prediction model `.db` files
- Model covariance files required by S-PrediXcan
- `gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz`
- MetaXcan `SPrediXcan.py` and `SMulTiXcan.py`
- An rsID-to-v8-varID mapping generated from the model databases

Outputs use `06.metaxcan`, `07.metaxcan_fdrreg`,
`08.smultixcan_output`, and `09.smultixcan_fdrreg`.

### GTEx v7

- GTEx v7 brain model `.db` files and covariance files
- `snp_covariance_v7.txt.gz`
- Compatible v7 MetaXcan scripts/environment

Outputs use the same directory names with `_v7` suffixes.

## 7. Downstream enrichment resources

Drug enrichment additionally requires:

- `mat.drug_DSigDB.Rdata`
- `ATC_drug_lists_ALL.Rdata`
- `005.drug_ATC.category_asFunc.R`
- Entrez and Ensembl biological annotation CSV files

KEGG/GO enrichment calls `gprofiler2`; the server therefore needs network
access to g:Profiler when those results are regenerated.

## 8. What belongs in Git

Commit the three frozen annotation CSVs because they reproduce the reported
models and are small enough for Git. Do not commit raw GWAS files, DIST or LDSC
reference panels, MetaXcan models, MAGMA reference panels, logs, or result
directories.

Readers seeking a compact package interface can use
<https://github.com/keelost/fdrregGenomics>. This repository intentionally
contains the detailed study workflow and intermediate-stage code.
