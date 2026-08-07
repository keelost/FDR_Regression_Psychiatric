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

## Analysis branches

- SNP-level FDRreg with allele harmonisation and LDSC intercept decorrelation
- DIST summary-statistic imputation from cleaned GWAS inputs
- MAGMA gene-level FDRreg with biological annotations
- MetaXcan and SMultiXcan using both GTEx v7 and GTEx v8 models
- Decorrelation diagnostics, rg and z-score sensitivity analyses
- Ablation analysis, temporal sensitivity/PPV validation, and risk-locus clumping
- Drug, ATC, KEGG, GO, and Reactome enrichment
- Result auditing, summaries, and gene/trait annotation utilities

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
