# Usage

## 1. Configure the server

```bash
cp config/project.env.example config/project.env
vim config/project.env
bash scripts/check_inputs.sh config/project.env
```

The complete upload and reference-resource checklist is in
`docs/input_inventory.md`.

## 2. LDSC

The LDSC runner consumes already munged target and library summary statistics:

```bash
bash scripts/run_pipeline.sh --stage ldsc --config config/project.env
```

It writes heritability results, pairwise genetic correlations, LDSC logs, and
`genetic_correlation_results.with_source.csv`. Set `FDRREG_LDSC_RESULTS` to
that final CSV before SNP FDRreg.

## 3. DIST imputation

Configure the cleaned target/library directories, DIST executable,
ancestry-specific reference directories, mapping table, and writable work
directory in `config/project.env`. Then run one target plus its configured
covariates:

```bash
bash scripts/run_pipeline.sh --stage imputation --target scz2012 --jobs 8 --config config/project.env
```

The output is one `<trait>.impute.map.txt` under `FDRREG_INPUT_DIR`. Use
`--target all` to prepare the union of every target and covariate in
`config/targets.tsv`.

## 4. Main analyses

```bash
bash scripts/run_pipeline.sh --stage snp --target scz2012 --config config/project.env
bash scripts/run_pipeline.sh --stage magma --target scz2012 --config config/project.env
bash scripts/run_pipeline.sh --stage metaxcan --target scz2012 --version both --jobs 4 --config config/project.env
bash scripts/run_pipeline.sh --stage smultixcan --target scz2012 --version both --jobs 4 --config config/project.env
```

Use `--target all` to expand all 16 configured targets. A comma-separated list
is also accepted.

## 5. Downstream analyses

Validated downstream scripts are under `scripts/real/common` and
`scripts/real/v7`. They retain their original command-line interfaces and
output conventions. Run `Rscript <script> --help` where supported.

The top-level runner currently orchestrates the LDSC and four main FDRreg
branches. Downstream sensitivity, validation, enrichment, and reporting remain
explicit standalone steps because they have different resource and completion
requirements.

## 6. GTEx versions

`--version v8` uses the directories without suffix. `--version v7` uses the
`_v7` directories. `--version both` executes them independently and never
copies results between releases.
