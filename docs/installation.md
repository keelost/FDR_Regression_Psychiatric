# Installation

## Base environment

```bash
conda env create -f environment.yml
conda activate fdrreg-psychiatric
```

The analysis was developed for Linux servers. Bash 4+, R 4.x, Python 3.9+
and GNU core utilities are expected.

## R packages

Required across the included scripts:

```r
install.packages(c(
  "data.table", "dplyr", "stringr", "readr", "optparse", "glmnet",
  "ggplot2", "scales", "doParallel", "Hotelling", "ICSNP", "mppa",
  "gprofiler2"
))
```

Install `FDRreg`, `powerplus`, `HelpersMG`, and `atc` from the same sources or
package snapshots used on the analysis server. These packages are central to
the validated statistical behavior and should not be silently substituted.
`FDRreg` is a third-party package developed by its original authors; it is not
authored or redistributed by this project and must be cited separately.

## External software

- LDSC, including `ldsc.py` and optionally `munge_sumstats.py`
- MAGMA
- PLINK 1.9/2
- MetaXcan/S-PrediXcan and SMultiXcan from
  <https://github.com/hakyimlab/MetaXcan>
- GTEx v7 and/or v8 prediction models, covariance resources, and annotations
  from the official MetaXcan downloads
- DIST 1.0.0 and ancestry-matched 1000 Genomes reference files when upstream
  imputation is being reproduced

Record every executable and reference path in `config/project.env`.

## Verification

```bash
bash scripts/check_inputs.sh config/project.env
bash tests/run_checks.sh
```
