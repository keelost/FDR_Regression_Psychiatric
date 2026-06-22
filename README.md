<!-- README.md -->
# FDRreg Psychiatric GWAS Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.9-blue.svg)](https://www.python.org)

## Overview

This repository implements a reproducible bioinformatics pipeline that applies
an **FDR regression (FDRreg)** empirical-Bayes framework to Genome-Wide
Association Study (GWAS) summary statistics for six psychiatric traits
(schizophrenia, bipolar disorder, major depressive disorder, autism spectrum
disorder, ADHD, and anxiety). Using 42 genetically correlated traits as
covariates and biological annotation features from pathway databases, the
pipeline identifies significantly associated SNPs, genes, and enriched
pathways across SNP-level, gene-level (MAGMA), tissue-specific transcriptomic
(S-PrediXcan across 13 GTEx brain regions), and multi-tissue transcriptomic
(SMultiXcan) analyses. Internal validation compares discoveries against
larger-scale GWAS, and simulation studies benchmark FDRreg performance.

---

## Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         INPUT DATA                                       │
│       GWAS Targets (6 psychiatric traits)                                │
│       GWAS Covariates (42 genetically correlated traits)                 │
└────────────────────────────────┬─────────────────────────────────────────┘
                                 │
                        ┌────────▼─────────┐
                        │  STAGE 0:        │
                        │  Quality Control │
                        │  Imputation (DIST│
                        │  / 1000G)        │
                        │  LDSC Genetic    │
                        │  Correlation     │
                        │  Harmonization   │
                        │  SNP Overlap     │
                        │  Variable Select.│
                        │  Decorrelation   │
                        └────────┬─────────┘
                                 │
         ┌───────────────┬───────┼───────┬───────────────┐
         │               │       │       │               │
  ┌──────▼──────┐ ┌──────▼─────┐ │ ┌─────▼──────┐ ┌──────▼───────┐
  │ BRANCH 1    │ │ BRANCH 2   │ │ │ BRANCH 3   │ │ BRANCH 4     │
  │ SNP-level   │ │ Gene-level │ │ │ Transcript.│ │ Multi-Tissue │
  │ FDRreg      │ │ MAGMA +    │ │ │ MetaXcan + │ │ SMultiXcan + │
  │ (with LDSC  │ │ FDRreg     │ │ │ FDRreg     │ │ FDRreg       │
  │ decorr.)    │ │ + bio ann. │ │ │ (13 brain  │ │ (combined    │
  │             │ │            │ │ │  regions)  │ │  brain)      │
  └──────┬──────┘ └──────┬─────┘ │ └─────┬──────┘ └──────┬───────┘
         │               │       │       │               │
         └───────────────┴───────┴───────┴───────────────┘
                                 │
                        ┌────────▼─────────┐
                        │  DOWNSTREAM:     │
                        │  Drug Enrichment │
                        │  KEGG / GO /     │
                        │  Reactome Enrich.│
                        │  Internal Valid. │
                        │  Simulation      │
                        │  In vivo valid.  │
                        └──────────────────┘
```

---

## Quick Start

### 1. Install

```bash
git clone https://github.com/username/FDRreg-psychiatric-GWAS.git
cd FDRreg-psychiatric-GWAS
conda env create -f environment.yml
conda activate fdreg-gwas
```

See [docs/installation.md](docs/installation.md) for full setup instructions.

### 2. Configure

```bash
# Edit configuration files to point to your data
vim config/config.yaml     # Paths, parameters, compute resources
vim config/traits.yaml     # Target and covariate trait definitions
vim config/fdreg_params.yaml  # FDRreg hyperparameters
```

### 3. Run

```bash
# Dry run (validate the DAG)
snakemake -s workflow/Snakefile --cores 1 --use-conda -n

# Full run
snakemake -s workflow/Snakefile --cores all --use-conda

# Run with custom config
snakemake -s workflow/Snakefile --cores 32 --use-conda \
  --config targets='["SCZ","BIP"]' fdr_threshold=0.01
```

---

## External Tools Required

| Tool | Version | Description |
|------|---------|-------------|
| [LDSC](https://github.com/bulik/ldsc) | ≥ 1.0 | LD Score Regression for genetic correlation |
| [MAGMA](https://ctg.cncr.nl/software/magma) | ≥ 1.08 | Gene-based GWAS analysis |
| [MetaXcan / S-PrediXcan](https://github.com/hakyimlab/MetaXcan) | Latest | Transcriptome-wide association (tissue-specific) |
| [SMultiXcan](https://github.com/hakyimlab/MetaXcan) | Latest | Multi-tissue TWAS |
| [FDRreg](https://github.com/jgscott/FDRreg) | Latest | FDR regression R package |
| [PLINK2](https://www.cog-genomics.org/plink/2.0/) | ≥ 2.0 | Genetic data manipulation |
| [DIST](https://github.com/omerwe/dist) | Latest | GWAS summary statistics imputation |
| [Snakemake](https://snakemake.readthedocs.io) | ≥ 7.0 | Workflow management |

---

## Input Data Requirements

| Data | Format | Required Fields | Description |
|------|--------|-----------------|-------------|
| Target GWAS | Space/tab-delimited | `snpid, chr, bpos, a1, a2, beta, se, pval, n` | Summary stats for 6 psychiatric traits |
| Covariate GWAS | Space/tab-delimited | `snpid, chr, bpos, a1, a2, beta (or z), se, pval, n` | Summary stats for 42 correlated traits |
| 1000 Genomes | PLINK binary | `.bed, .bim, .fam` | European reference panel (for imputation & MAGMA) |
| LDSC reference | LD score files | Per-chromosome `.l2.ldscore.gz` | EUR LD scores + `w_hm3.snplist` |
| MAGMA annotation | `.genes.annot` | Gene–SNP mapping | 10 kb window gene annotation |
| GTEx prediction models | `.db` + covariance | Tissue-specific eQTL weights | v7 imputed European models (13 brain regions) |
| DAVID annotations | CSV | Entrez ID + binary feature columns | Biological annotations for gene-level models |
| Drug database | CSV | Drug ID, ATC code, mechanism | 17,840 drugs → 267 clusters |

---

## Output Files

| Branch | Output File Pattern | Description |
|--------|---------------------|-------------|
| 1 — SNP | `results/snp_fdreg/{trait}.snp.fdrreg.txt` | Per-SNP FDRreg results with decorrelated z, FDR, q-values |
| 2 — Gene | `results/magma_fdreg/{trait}.gene.fdrreg.txt` | Per-gene FDRreg with GWAS + biological annotation covariates |
| 3 — Region | `results/metaxcan_fdreg/{trait}/{region}.fdrreg.txt` | Per-gene FDRreg per brain region |
| 4 — Combined | `results/smultixcan_fdreg/{trait}.smultixcan.fdrreg.txt` | Combined multi-tissue FDRreg |
| Enrichment | `results/enrichment/{trait}.drug_enrichment.csv` | Drug cluster enrichment (267 clusters) |
| Enrichment | `results/enrichment/{trait}.kegg_go_enrichment.csv` | KEGG, Reactome, GO BP enrichment |
| Validation | `results/validation/{trait_pair}.validation.txt` | Small vs. large GWAS overlap statistics |
| Simulation | `results/simulation/sim_results.csv` | FDRreg benchmarking under simulated signals |

---

## Citation

If you use this pipeline, please cite:

> Rao, S. T., Qiu, J. H., Zhi, Y. Q., Lin, Y. P., Zhang, R. Y., Chen, X. T.,
>  ... & So, H. C. (2024). Discovering additional genetic loci associated
> with six psychiatric disorders/traits via FDR regression model leveraging
> external genetic and biological data. medRxiv, 2024-01.](https://doi.org/10.1101/2024.01.29.24301912)

```bibtex
@article {Rao2024.01.29.24301912,
	author = {Rao, Shi-tao and Qiu, Jing-hong and Zhi, Yi-qiang and Lin, Yu-ping and Zhang, Ruo-yu and Chen, Xiao-tong and Xu, Dan and So, Hon-Cheong},
	title = {Discovering additional genetic loci associated with six psychiatric disorders/traits via FDR regression model leveraging external genetic and biological data},
	elocation-id = {2024.01.29.24301912},
	year = {2024},
	doi = {10.1101/2024.01.29.24301912},
	publisher = {Cold Spring Harbor Laboratory Press},
	abstract = {Background Common psychiatric disorders have substantial heritability influenced by multiple genes. While a number of susceptibility variants have been identified, many associated variants remain undiscovered. This study aimed to identify additional genetic loci associated with common psychiatric disorders/traits by leveraging correlated traits and biological annotations.Methods We proposed application of the false discovery rate (FDR) regression model to uncover additional genetic loci for six psychiatric disorders/traits. To enhance the likelihood of discovering additional significant genetic loci and genes, we utilized a set of 42 correlated traits and 21 biological annotations as covariates. Internal validation analysis and drug cluster enrichment analysis were conducted to validate the biological significance of the additional genetic loci/genes uncovered. We also experimentally validated two additional genes revealed for autism spectrum disorder (ASD).Results The FDR regression (FDRreg) analysis strategy revealed hundreds of additional significant genes (FDR\&lt;0.01) in gene-level analyses, surpassing the number of significant genes found in the original studies. Specifically, in 11/16 trait analyses, FDRreg identified more significant genes based on gene-based analysis with MAGMA, and in 12/16 analyses, FDRreg identified more significant genes based on imputed expression in the brain. In SNP-level results, the majority of analyses (13/16) identified an equal or higher number of genomic risk loci (FDR\&lt;0.01). We found that FDRreg is able to reveal genes that are later known to be significant in subsequent larger-scale GWAS. Drug cluster enrichment analysis demonstrated a stronger enrichment in psychiatry-related drug clusters. In utero electroporation (IUE) experiments provided evidence to support two additional genes identified for ASD in critical embryonic brain development processes.Conclusions By integrating genetically correlated traits and biological annotations, the FDRreg strategy enables the identification of a greater number of additional significant genes and risk loci. Moreover, the new associated genes exhibited meaningful biological and clinical implications. This study presents a valuable approach for uncovering the genetic basis of psychiatric disorders and gaining insights into their underlying biology.Competing Interest StatementThe authors have declared no competing interest.Funding StatementThis study was supported by Fujian Provincial Natural Science Foundation Youth Innovation Project (grant no. 2021J05050), Fujian Province Joint Innovation Project (grant no. 2021Y9030), Research start-up funds for high-level talents from Fujian Medical University (grant no. XRCZX2021009).Author DeclarationsI confirm all relevant ethical guidelines have been followed, and any necessary IRB and/or ethics committee approvals have been obtained.YesI confirm that all necessary patient/participant consent has been obtained and the appropriate institutional forms have been archived, and that any patient/participant/sample identifiers included were not known to anyone (e.g., hospital staff, patients or participants themselves) outside the research group so cannot be used to identify individuals.YesI understand that all clinical trials and any other prospective interventional studies must be registered with an ICMJE-approved registry, such as ClinicalTrials.gov. I confirm that any such study reported in the manuscript has been registered and the trial registration ID is provided (note: if posting a prospective study registered retrospectively, please provide a statement in the trial ID field explaining why the study was not registered in advance).YesI have followed all appropriate research reporting guidelines, such as any relevant EQUATOR Network research reporting checklist(s) and other pertinent material, if applicable.YesAll data produced in the present study are available upon reasonable request to the authors.},
	URL = {https://www.medrxiv.org/content/early/2024/01/30/2024.01.29.24301912},
	eprint = {https://www.medrxiv.org/content/early/2024/01/30/2024.01.29.24301912.full.pdf},
	journal = {medRxiv}
}
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
