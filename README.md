# FDRreg Psychiatric GWAS Analysis

This repository contains the validated R, Python, and Bash analysis code used
for FDR regression across 16 psychiatric GWAS datasets. It uses an ordinary
Bash command-line entry point with independently runnable R and Python scripts.

## Scope and attribution

This repository does **not** introduce, implement, or claim authorship of the
FDRreg statistical method or the third-party `FDRreg` R package. FDRreg was
developed and published by its original authors. This repository contributes
the study-specific psychiatric genomics workflow: GWAS preparation, covariate
selection, decorrelation, orchestration across SNP/MAGMA/MetaXcan/SMultiXcan
analyses, sensitivity analyses, enrichment, validation, and reporting.

Users must install and cite the upstream `FDRreg` package and method according
to its documentation, in addition to citing this psychiatric workflow and its
associated study. See [NOTICE](NOTICE) for the ownership boundary.

The pipeline accepts cleaned GWAS summary statistics and can reproduce the
DIST summary-statistic imputation that generates `*.impute.map.txt`. Because
the source GWAS are access controlled, users must obtain them from the original
providers and perform the documented initial cleaning before running this
project.

## Complete analysis workflow

The scientific workflow begins with the original GWAS summary statistics. The
source GWAS cannot be distributed with this repository because of data-use
restrictions; users must obtain each dataset from its original provider. The
main dependency chain is:

```text
raw data
  -> clear data
  -> LDSC selection
  -> imputation
  -> harmonisation
  -> decorrelation
       |-> SNP FDRreg -> risk loci
       |-> MAGMA conversion -> MAGMA FDRreg
       |-> MetaXcan conversion -> MetaXcan FDRreg
       `-> MetaXcan conversion -> SMultiXcan conversion -> SMultiXcan FDRreg
```

The four FDRreg branches use the independently developed third-party `FDRreg`
R package. This repository supplies the psychiatric-genomics data preparation,
covariate construction, orchestration, sensitivity analyses, and reporting
around that package.

### 1. Raw data and cleaning

`raw data` means GWAS summary statistics obtained from the cited consortium or
repository. Cleaning produces one normalized `*.clear.txt` file per trait. The
cleaning procedure standardizes column names and alleles, filters variants at
INFO >= 0.6 when INFO is available, removes duplicate SNPs and InDels, checks
or restores the correct SNP identifier, and calculates a signed z-score from
`beta / se` when required. Cleaning is documented here but remains an upstream
user responsibility because the licensed source data are not included.

### 2. LDSC selection

LDSC is run on munged target and candidate-covariate summary statistics. Its
heritability and pairwise genetic-correlation results are used to select the
trait covariates listed for each of the 16 targets in `config/targets.tsv`.
For datasets with sample overlap, LDSC cross-trait intercepts also supply the
covariance matrix used later for decorrelation. The reconstructed,
source-labelled LDSC CSV is therefore both a covariate-selection record and a
required decorrelation input.

### 3. Imputation with DIST

The imputation workflow converts each cleaned target and covariate GWAS to
DIST input, runs chromosomes 1-22 against the ancestry-matched 1000 Genomes
reference, keeps imputed variants at INFO >= 0.5, maps non-rs variant IDs by
chromosome and position, removes duplicates, and writes the final
`*.impute.map.txt` files. Reference ancestry is configured in
`config/imputation_refs.tsv`.

### 4. Harmonisation and decorrelation

Target and covariate files are intersected by SNP ID. Effect alleles are
aligned so that all signed z-scores have the same orientation (`harmo`). For
traits with sample overlap, the LDSC intercept covariance matrix is converted
to an inverse square-root transform and applied to the harmonised z-scores
(`decorrelation`). Traits without sample overlap retain their harmonised
scores. These aligned data feed every downstream analysis branch.

### 5. SNP FDRreg and risk loci

SNP FDRreg uses the decorrelated target z-score as the response and the
selected trait z-scores as covariates. The workflow saves theoretical-null,
empirical-null, LASSO, and biological-covariate results where configured. The
significant SNP results are then clumped with PLINK to derive independent risk
loci and risk-locus summary tables.

### 6. MAGMA branch

The decorrelated SNP data are converted to MAGMA input, mapped to genes, and
analysed with MAGMA. Gene-level statistics are joined to the frozen,
manually curated biological annotations in `data/annotations`, then passed to
MAGMA FDRreg. This produces the gene-level FDRreg results used by ablation,
validation, and enrichment analyses.

### 7. MetaXcan branch

The decorrelated GWAS statistics are converted to the format expected by
MetaXcan and analysed separately for the 13 configured GTEx brain regions.
Each region's gene-level output is then passed to MetaXcan FDRreg. GTEx v7 and
v8 models, covariance files, annotations, inputs, and outputs remain strictly
separate throughout this branch.

### 8. SMultiXcan branch

SMultiXcan conversion combines the decorrelated GWAS data with the regional
MetaXcan conversion outputs. S-MultiXcan then integrates evidence across the
13 brain regions, after which the multi-tissue gene statistics are analysed by
SMultiXcan FDRreg. This branch is also run independently for GTEx v7 and v8.

## Extra analyses

The supplementary analyses consume the following main-workflow products:

| Input | Extra analysis | Purpose |
|---|---|---|
| Imputation | DIST audit | Check overlap, mapping, and variant retention in the DIST-imputed files |
| SNP FDRreg | Decorrelation diagnosis (`deco_diagnosis`) | Compare pre- and post-decorrelation correlations and validate the transform |
| SNP, MAGMA, and SMultiXcan FDRreg | Genetic-correlation sensitivity (`rg`) | Refit or compare results under alternative genetic-correlation covariates |
| SNP and MetaXcan FDRreg | Absolute-z sensitivity (`zabs`) | Test whether conclusions are robust to using absolute rather than signed z-scores |
| SNP, MAGMA, MetaXcan, and SMultiXcan FDRreg | Sensitivity/PPV (`sen_ppv`) | Perform temporal validation and compare sensitivity and positive predictive value |
| MAGMA, MetaXcan, and SMultiXcan FDRreg | Ablation | Quantify the contribution of covariate groups and model components |
| MAGMA and SMultiXcan FDRreg | Drug enrichment | Test significant genes against drug and ATC resources |
| MAGMA and SMultiXcan FDRreg | KEGG/GO enrichment | Test pathway, GO, KEGG, and Reactome enrichment using analysis-specific backgrounds |

Summary and diagnostic utilities aggregate these versioned outputs without
mixing GTEx releases. The v7 branch is retained because it more closely
reproduces the earlier study results; v8 provides results using the newer
models.

## Code organization

Validated source scripts are organized under:

```text
scripts/real/common     release-independent and GTEx v8 analysis code
scripts/real/v7         GTEx v7-specific analysis code
scripts/real/utilities  conversion, audit, and summary utilities
scripts/real/legacy     superseded or path-inconsistent historical versions
scripts/ldsc            LDSC execution and result reconstruction
scripts/imputation      DIST input preparation, execution, and output mapping
```

When duplicate source files existed, the version in the original test
directory was selected. GTEx v7 and v8 scripts are intentionally both kept.

The frozen biological-annotation matrices used in the study are provided in
`data/annotations`. They contain the manually curated DAVID-derived scores and
the February 2020 Open Targets and denovo-db-derived features described in the
manuscript and supplementary methods.

## Configuration

```bash
cp config/project.env.example config/project.env
vim config/project.env
```

Important configuration files:

- `config/project.env`: server paths and executable locations
- `config/targets.tsv`: the 16 target datasets and sample-overlap groups
- `config/regions.txt`: 13 GTEx brain regions
- `config/versions.tsv`: v7/v8 directory and filename conventions

Run the input checker before starting:

```bash
bash scripts/check_inputs.sh config/project.env
```

See [docs/input_inventory.md](docs/input_inventory.md) for the exact files that
must be uploaded or made available on the server, including the LDSC CSV.

## Running

```bash
# Show available options
bash scripts/run_pipeline.sh --help

# List configured targets
bash scripts/run_pipeline.sh --stage list --config config/project.env

# Build LDSC results and the source-labelled correlation CSV
bash scripts/run_pipeline.sh --stage ldsc --config config/project.env

# Build *.impute.map.txt for one target and its configured covariates
bash scripts/run_pipeline.sh --stage imputation --target scz2012 --jobs 8 \
  --config config/project.env

# SNP FDRreg for one target
bash scripts/run_pipeline.sh --stage snp --target scz2012 \
  --config config/project.env

# MetaXcan FDRreg for both model releases
bash scripts/run_pipeline.sh --stage metaxcan --target scz2012 \
  --version both --jobs 4 --config config/project.env

# SMultiXcan FDRreg using GTEx v7 only
bash scripts/run_pipeline.sh --stage smultixcan --target scz2012 \
  --version v7 --config config/project.env
```

Each validated R/Python script can still be executed directly. Use its
`--help` output for analysis-specific options.

## Expected directory layout

```text
<base>/
  <target>/
    00.overlap_data/
    01.magma_input/
    02.fdrreg_results/
    03.contribution/
    04.magma_output/
    05.magma_fdrreg/
    06.metaxcan/                 # GTEx v8
    06.metaxcan_v7/              # GTEx v7
    07.metaxcan_fdrreg/
    07.metaxcan_fdrreg_v7/
    08.smultixcan_output/
    08.smultixcan_output_v7/
    09.smultixcan_fdrreg/
    09.smultixcan_fdrreg_v7/
  01.extra.analysis/
```

Versioned result directories must not be merged. The v7 branch reproduces the
earlier analysis more closely; v8 uses the current GTEx mashr models.

MetaXcan code, prediction models, covariance files, and annotations are not
vendored here. Download release-compatible resources from the official
[MetaXcan repository](https://github.com/hakyimlab/MetaXcan). Readers who only
need a compact interface to the method can use
[fdrregGenomics](https://github.com/keelost/fdrregGenomics); this repository is
the detailed, study-specific reproducibility workflow.

## Installation

Create the base environment:

```bash
conda env create -f environment.yml
conda activate fdrreg-psychiatric
```

Some research packages and external tools are not consistently available from
Conda and may need installation from their original projects. See
[docs/installation.md](docs/installation.md).

## Testing

Local verification is intentionally limited to syntax, configuration, help
commands, and small fixtures. Full analyses require the server data and
reference resources.

```bash
bash tests/run_checks.sh
```

## Citation

See [CITATION.cff](CITATION.cff) for the project citation and the associated
FDRreg psychiatric GWAS study.
