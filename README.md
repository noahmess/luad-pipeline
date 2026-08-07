# Nuclear envelope remodelling in lung adenocarcinoma

Analysis pipeline supporting the manuscript *Nuclear envelope remodelling in lung adenocarcinoma: prognostic value of LMNB2, TMPO, NDC1 and LBR, and identification of NUP88 as a Ki-67-uncoupled transcript associated with reduced myeloid infiltration*.

This repository contains the R scripts and output tables needed to reproduce every figure, table and statistic reported in the paper, including the revision analyses (cross-signature benchmark, leave-one-gene-out sensitivity, random-signature permutation test, and the GSE72094 replication of the NUP88-immune correlation).

**Scope of this repository.** Only code (`R/`) and output CSVs (`outputs/tables/`) are version-controlled. Raw and cached data, figures and manuscript documents are excluded (see *Distribution notes*) and are either regenerable from the code and the public accessions listed below, or are not part of the reproducibility chain.

---

## Overview

The pipeline takes TCGA-LUAD raw data and three independent cohorts (GSE31210, GSE50081, GSE72094) and produces:

- a coordinated expression signature of twenty nuclear envelope genes,
- a four-gene prognostic score (Score A: LMNB2, TMPO, NDC1, LBR),
- a comparative benchmark of Score A against four published signatures,
- a leave-one-gene-out sensitivity analysis,
- a random-signature permutation benchmark (Score A vs. 1000 random gene panels, two null pools),
- a replication of the NUP88-immune correlation (xCell) in an independent cohort,
- all manuscript figures and tables.

Cohort sizes used throughout the manuscript:

| Cohort | n | Role |
|---|---|---|
| TCGA-LUAD expression cohort | 540 tumours, 59 normals | DEG, co-expression, immune, ROC |
| TCGA-LUAD harmonised survival cohort | **394** (129 events) | Cox, Score A, KM, permutation benchmark discovery |
| TCGA held-out (30% partition, seed 42) | 118 | Internal validation |
| GSE31210 | 204 (Japanese, 30 events) | External validation |
| GSE50081 | 181 (Canadian, 75 events) | External validation |
| GSE72094 | 442 (398 with complete survival) | External replication, NUP88-immune correlation |

---

## Repository structure

```
luad-pipeline/
├── R/                                   Pipeline scripts (numbered execution order)
│   ├── 00_setup.R                       Seed, paths, gene sets, palettes, helpers
│   ├── 01_download_tcga.R               TCGAbiolinks download
│   ├── 02_differential_expression.R     DESeq2
│   ├── 03_functional_enrichment.R       clusterProfiler (GO, KEGG)
│   ├── 04_cohort_harmonization.R        Builds n = 540 and n = 394 cohorts
│   ├── 05_survival_cox.R                Univariate and multivariate Cox, Score A
│   ├── 05b_table1_clinical.R            Manuscript Table 1
│   ├── 05c_incremental_analyses.R       Incremental value, stage I, C-index
│   ├── 06_external_validation.R         Score A in GSE31210 and GSE50081
│   ├── 07_immune_analysis.R             Immune marker correlations (TCGA)
│   ├── 08_figures.R                     Figures 1-9 and S1
│   ├── 09_figure10_nup88.R              Figure 10 (NUP88 integrated profile)
│   ├── 10_signatures.R                  Definitions of the five benchmarked signatures
│   ├── 11_benchmark.R                   Cross-signature benchmark on four cohorts
│   ├── 12_sensitivity_LOO.R             Leave-one-gene-out + nested LRT for NDC1
│   ├── 13_random_signature_benchmark.R  Score A vs. 1000 random gene panels
│   ├── 14_methylation_7genes.R          Figure 14, promoter methylation (7 genes)
│   ├── 15_gse72094_xcell.R              GSE72094 QC + xCell deconvolution
│   └── 16_nup88_immune_correlation_gse72094.R   NUP88-immune replication in GSE72094
│
├── outputs/
│   ├── tables/                          Manuscript and supplementary tables (CSV) - tracked
│   ├── figures/                         Figures (PNG + TIFF, 300 dpi) - NOT tracked, regenerate locally
│   └── sessionInfo.txt                  R and package versions this pipeline was run against
│
├── run_all.R                            Master script, steps 00-09 (see note below)
├── run_revision.R                       Master script, steps 10-12 (see note below)
├── .gitignore
├── LICENSE                              MIT (code only, see below)
└── README.md                            This file
```

`data/`, `GDCdata/`, `archive/`, `*.RData`, `*.rds`, `*.docx`, `*.xlsx` and `outputs/figures/*.tiff|png` exist in the working copy but are excluded from version control (see *Distribution notes*).

**Not currently wired into a master script.** `run_all.R` covers steps 00-09 only and does not call `05b`/`05c`; `run_revision.R` covers steps 10-12 only. Scripts `05b`, `05c`, `13`, `14`, `15`, `16` must each be run with `source("R/<script>.R")` individually, in numeric order, after their prerequisites exist. This is accurate as of this archive - it is not an oversight to be silently fixed, and extending the master scripts is a reasonable follow-up if you want single-command reproduction of the full paper.

---

## Requirements

- **R** 4.3 or later (this archive was run on 4.6.0, macOS aarch64 - see `outputs/sessionInfo.txt` for the exact versions and full dependency tree)
- **Bioconductor 3.18+**

Required CRAN packages:

```
survival, glmnet, timeROC, ggplot2, ggrepel, ggpubr, survminer, pROC,
RColorBrewer, pheatmap, cowplot, patchwork, dplyr, tidyr, readxl
```

Required Bioconductor packages:

```
TCGAbiolinks, DESeq2, clusterProfiler, org.Hs.eg.db, GEOquery,
Biobase, SummarizedExperiment, GSVA (>= 1.50, for the ssgseaParam interface)
```

Required GitHub package:

```
xCell (dviraran/xCell) - install with remotes::install_github("dviraran/xCell")
```

Install in one step:

```r
install.packages(c("survival", "glmnet", "timeROC", "ggplot2", "ggrepel", "ggpubr",
                   "survminer", "pROC", "RColorBrewer", "pheatmap",
                   "cowplot", "patchwork", "dplyr", "tidyr", "readxl", "remotes"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "clusterProfiler",
                       "org.Hs.eg.db", "GEOquery", "Biobase",
                       "SummarizedExperiment", "GSVA"))

remotes::install_github("dviraran/xCell", upgrade = "never")
```

---

## Reproducing the main pipeline

The full pipeline from raw TCGA download to figures runs in roughly 30 to 60 minutes depending on the network and on whether the TCGA data has already been cached.

```r
setwd("~/Desktop/luad-pipeline")
source("run_all.R")
```

`run_all.R` sources scripts 00 to 09 in order, in one R session, so objects persist between steps. To also reproduce Table 1 and the incremental analyses, run `05b` and `05c` right after (they only need objects already in memory from steps 00-05):

```r
source("R/05b_table1_clinical.R")
source("R/05c_incremental_analyses.R")
```

**Output**

- `outputs/figures/` populates with Figures 1 to 10 and Figure S1 (not tracked in this repository - regenerate locally).
- `outputs/tables/` populates with Table 1, Cox univariate and multivariate tables, Score A external validation and immune correlation tables.

---

## Reproducing the revision analyses

```r
setwd("~/Desktop/luad-pipeline")
source("run_revision.R")
```

Requires `data/TCGA_LUAD_cohort.RData` and `data/TCGA_LUAD_data.RData` (produced by the main pipeline). Runtime is about one minute once GEO ExpressionSets are cached in `data/GEO_cache/`.

Steps performed:

1. **`R/10_signatures.R`** defines the five benchmarked signatures:
   - **Score A**: equal-weight Z-score average of LMNB2, TMPO, NDC1, LBR (coefficient-free).
   - **LASSO-Cox**: refit on a 70/30 partition of TCGA (seed = 42), candidate pool of nine genes.
   - **Li 2023 TME6**: PLK1, LDHA, FURIN, FSCN1, RAB27B, MS4A1 (Li et al., 2023, *Sci Rep* 13:13854).
   - **Zhou 2022**: COL1A1, GPX3, PLEK2 (Zhou et al., 2022, *Transl Lung Cancer Res* 11:1827).
   - **CCP-like**: 31-gene cell-cycle progression score (Cuzick 2011, *Lancet Oncol* 12:245; Wistuba 2013).

2. **`R/11_benchmark.R`** (the canonical benchmark script - see *Note on 11_benchmark* below) runs each signature across four cohorts and reports HR per 1 SD, Cox p-value, Harrell C-index with 95% CI, and time-dependent AUC at 12, 36 and 60 months (timeROC, Blanche 2013). Also writes `BENCHMARK_ENV.RData`, a cache of the four cohorts and helper functions reused by `R/12`, `R/13`.

3. **`R/12_sensitivity_LOO.R`** runs four Score A variants with one gene removed at a time, plus a nested likelihood-ratio test on TCGA full comparing `Surv ~ LMNB2 + TMPO + LBR` against `Surv ~ LMNB2 + TMPO + LBR + NDC1`.

**Output**

- `outputs/tables/Benchmark_main_3cohorts.csv` (held-out, GSE31210, GSE50081).
- `outputs/tables/Benchmark_suppl_TCGAfull.csv` (TCGA n = 394).
- `outputs/tables/Supplementary_Table_S3_LOO.csv` (leave-one-gene-out).
- `BENCHMARK_ENV.RData` (cached cohorts and helper functions; not redistributed, regenerable, ~210 MB).

**Note on `11_benchmark`.** Two near-duplicate scripts existed during development, `11_benchmark.R` and `11_benchmark-2.R`. Only `11_benchmark.R` is complete and current - `run_revision.R` sources it, and it matches the `BENCHMARK_ENV.RData` on disk. `11_benchmark-2.R` was missing one assignment (`gene_symbols <- as.character(rowData(data)$gene_name)`) and would fail if run; it has been moved to `archive/old_scripts/11_benchmark-2_broken_duplicate.R` (not tracked in this repository) rather than silently deleted.

---

## Reproducing the additional revision analyses (not in `run_revision.R`)

These three analyses were added after `run_revision.R` was written and are run independently.

### Random-signature permutation benchmark

```r
source("R/13_random_signature_benchmark.R")
```

Requires `BENCHMARK_ENV.RData` (from `R/11_benchmark.R`) and `data/TCGA_LUAD_data.RData` (for GENCODE gene-biotype annotation). Compares Score A against 1000 random four-gene panels on four cohorts (TCGA discovery + three independent validation cohorts), against two null pools (any protein-coding gene; protein-coding genes correlated with MKI67, |rho| > 0.40). Seed: `SEED` (42, from `R/00_setup.R`).

Output: `outputs/tables/Random_signature_benchmark_detail.csv`, `outputs/tables/Random_signature_benchmark_summary.csv`.

Result: Score A beats random panels only in the TCGA discovery cohort (empirical p = 0.004 on C-index); this does not replicate as statistically significant in the three independent cohorts (empirical p 0.11-0.54), though direction is consistent. See manuscript Discussion for the power-based interpretation.

### Methylation (Figure 14)

```r
source("R/14_methylation_7genes.R")
```

Requires the TCGA 450k methylation level-3 beta values cached under `GDCdata/` (not tracked - several GB; rebuild with a GDC methylation query via `TCGAbiolinks`, the same tool used in `R/01_download_tcga.R`, but no dedicated download script for this specific query is included in this archive). No random seed required (deterministic Wilcoxon/t-tests).

Output: `outputs/tables/Methylation_7genes_stats.csv`, `outputs/figures/Figure14_methylation_7genes.tiff`.

### NUP88-immune correlation, replicated in GSE72094

```r
source("R/15_gse72094_xcell.R")                        # QC + xCell deconvolution
source("R/16_nup88_immune_correlation_gse72094.R")      # correlation + replication verdict
```

`R/15` downloads GSE72094 (cached in `data/GEO_cache/`), runs three pre-analysis checks (NUP88 probe presence, confirmed LUAD histology for all 442 samples, survival completeness), then runs xCell. `R/16` correlates NUP88 (mean of its 3 probes - concordant, pairwise Spearman rho 0.76-0.94) against four xCell scores (`aDC`, `Macrophages`, `Macrophages M1`, `B-cells`) and compares against the existing TIMER3.0 reference values in `S2 TABLE TIMER3 NUP88 LMBN2.xlsx` (not re-derived - TIMER3.0 was not re-run on TCGA, to avoid a second, divergent TCGA number; see the script header for the full rationale). No random seed required (deterministic correlation tests).

The mapping from xCell's raw score names to TIMER3.0's display labels (`aDC` -> "Myeloid dendritic cell activated") was confirmed against the `immunedeconv` source code (`inst/extdata/cell_type_mapping.xlsx`, github.com/omnideconv/immunedeconv), not assumed.

**Methodological asymmetry to note when citing this result:** GSE72094 correlations are unadjusted Spearman correlations; the TIMER3.0 reference values are purity-adjusted partial Spearman correlations (no tumour-purity estimate is available for GSE72094).

Replication criterion (fixed before computing): replicated if >= 3/4 populations show the same (negative) direction and >= 2/4 are significant after Benjamini-Hochberg correction. Result: 4/4 negative, 4/4 significant.

Output: `outputs/tables/NUP88_immune_correlation_GSE72094_vs_TIMER3.csv`.

---

## Data sources

All data analysed in this study are publicly available.

| Source | Accession / access | Used by | Cached at |
|---|---|---|---|
| TCGA-LUAD | NCI Genomic Data Commons (`TCGAbiolinks`) | `R/01`-`R/09`, `R/13` | `data/TCGA_LUAD_data.RData`, `data/TCGA_LUAD_cohort.RData` |
| GSE31210 | NCBI GEO, Okayama et al. (2012) | `R/11`, `R/12` | `data/GEO_cache/GSE31210_eset.rds` |
| GSE50081 | NCBI GEO, Der et al. (2014) | `R/11`, `R/12` | `data/GEO_cache/GSE50081_eset.rds` |
| GSE72094 | NCBI GEO (442-patient lung adenocarcinoma cohort; consult the GEO record for the associated publication) | `R/15`, `R/16` | `data/GEO_cache/GSE72094_eset.rds` |
| TIMER3.0 | Cui et al. (2025), web portal, purity-adjusted partial Spearman correlations | `R/07`, `R/16` (reference only) | `S2 TABLE TIMER3 NUP88 LMBN2.xlsx` (not tracked - see *Distribution notes*) |
| Human Protein Atlas, cBioPortal (PanCancer Atlas LUAD), UALCAN, KM-Plotter | web portals, no programmatic download | manuscript cross-checks | - |

The TCGA methylation level-3 beta values used by `R/14` (`GDCdata/`) are excluded from this repository because of size; re-querying via `TCGAbiolinks` will rebuild what is needed.

---

## Reproducibility notes (seeds and parameters)

Every script that calls a random-number-consuming function (`sample()`, `geom_jitter`, `cv.glmnet`'s internal CV fold assignment) sets `set.seed()` beforehand, using the shared `SEED <- 42L` constant from `R/00_setup.R`:

| Script | Stochastic step | Seed |
|---|---|---|
| `R/08_figures.R` | `geom_jitter` (point placement) | `set.seed(SEED)`, re-seeded immediately before the plot |
| `R/11_benchmark.R` | 70/30 TCGA train/held-out partition; `cv.glmnet` fold assignment (LASSO-Cox) | `set.seed(42)` at script top |
| `R/13_random_signature_benchmark.R` | Random gene-panel sampling (1000 permutations x 2 null pools x 4 cohorts) | `set.seed(SEED)` at script top |

Scripts with no stochastic step (`00`-`07` except noted, `05b`, `05c`, `09`, `10`, `12`, `14`, `15`, `16`) do not set a seed - they don't need one; this was checked by grepping for `sample(`, `rnorm(`, `runif(`, `boot(` and jitter/repel calls across every script in `R/`, not assumed.

Key non-seed parameters, declared as named constants near the top of their script rather than hard-coded inline:

| Parameter | Value | Script |
|---|---|---|
| `N_PERM` | 1000 | `R/13_random_signature_benchmark.R` |
| `MKI67_RHO_CUTOFF` | 0.40 | `R/13_random_signature_benchmark.R` |
| `CONCORDANCE_RHO_CUTOFF` (NUP88 probe agreement) | 0.5 | `R/16_nup88_immune_correlation_gse72094.R` |
| LASSO-Cox: `alpha`, `nfolds`, `lambda` | 1, 10, `lambda.min` | `R/11_benchmark.R` |
| DEG significance threshold | `padj < 0.05` | `R/02_differential_expression.R` |

Caveat: reproducing exact values requires running each script in full, from the top, in a fresh state - the seed fixes the RNG stream from that point forward, so sourcing only part of a script, or re-running a code block out of order in an existing session, will not reproduce the published numbers.

---

## Distribution notes

Only `R/*.R` and `outputs/tables/*.csv` are version-controlled (see `.gitignore`). Excluded, and why:

- **`data/`, `GDCdata/`** - raw and cached data (TCGA RangedSummarizedExperiment, DESeq2 results, GEO ExpressionSets, 450k methylation betas). Several hundred MB to multiple GB; fully regenerable from the public accessions above via `R/01`, `R/15`, or the manual GDC methylation query used for `R/14`.
- **`BENCHMARK_ENV.RData` (~210 MB), `SESSION_BACKUP.RData` (~410 MB), `.RData`** - working snapshots/caches, regenerable, not required for reproduction.
- **`*.rds` at repository root (`df.rds`, `results.rds`), `gene_table-*.csv`** - orphaned artefacts not referenced by any current script; excluded rather than deleted.
- **`outputs/figures/*.tiff`, `*.png`** - binary, regenerable by `R/08`, `R/09`, `R/14`.
- **`*.docx`, `*.xlsx`** - manuscript drafts and the TIMER3.0 reference table (`S2 TABLE TIMER3 NUP88 LMBN2.xlsx`). Not code, not a pipeline-generated CSV, per this archive's scope. If you want the TIMER3.0 reference table bundled with the Zenodo deposit for readers who don't have the manuscript, add it back explicitly (`git add -f "S2 TABLE TIMER3 NUP88 LMBN2.xlsx"`) - it is small (24 KB).
- **`archive/`** - superseded drafts and one broken duplicate script, kept locally for the author's own traceability only.

---

## Manuscript outputs at a glance

| Output | Produced by | Used in |
|---|---|---|
| Figure 1 (volcano) | `08_figures.R` | Section 3.1 |
| Figure 2 (heatmap, 20 genes) | `08_figures.R` | Section 3.1 |
| Figure 3 (GO BP, CC) | `08_figures.R` | Section 3.1 |
| Figure 4 (Kaplan-Meier, four genes, n=394) | `08_figures.R` | Section 3.3 |
| Figure 5 (Kaplan-Meier, Score A, n=394) | `08_figures.R` | Section 3.3 |
| Figure 6 (Kaplan-Meier, Score A, stage I) | `08_figures.R` | Section 3.3 |
| Figure 7 (external validation) | `08_figures.R` | Section 3.4 |
| Figure 8 (ROC, tumour vs normal) | `08_figures.R` | Section 3.4 |
| Figure 9 (co-expression network) | `08_figures.R` | Section 3.5 |
| Figure 10 (NUP88 integrated profile) | `09_figure10_nup88.R` | Sections 3.5-3.6 |
| Figure 14 (promoter methylation, 7 genes) | `14_methylation_7genes.R` | Section 3.7 |
| Figure S1 (Ki-67 correlation) | `08_figures.R` | Section 3.5 |
| Table 1 (clinical characteristics, n=394) | `05b_table1_clinical.R` | Section 3.2 |
| Table 2 (univariate Cox) | `05_survival_cox.R` | Section 3.2 |
| Table 3 (multivariate Cox) | `05_survival_cox.R` | Section 3.2 |
| Table 4 (Score A summary) | `05_survival_cox.R` | Section 3.3 |
| Table 6 (cross-signature benchmark) | `11_benchmark.R` | Section 3.9 |
| Supplementary Table S2 (TIMER3.0, NUP88 and LMNB2) | external (TIMER3.0) | Section 3.6 |
| Supplementary Table S3 (leave-one-gene-out, NDC1 LRT) | `12_sensitivity_LOO.R` | Section 3.9 |
| Random-signature permutation benchmark | `13_random_signature_benchmark.R` | Discussion (revision) |
| NUP88-immune correlation, GSE72094 replication | `15_gse72094_xcell.R`, `16_nup88_immune_correlation_gse72094.R` | Results (revision) |

---

## Citation

If you use this pipeline or its derivatives, please cite the manuscript (citation block will be added once the DOI is assigned).

A Zenodo archive of this repository will be deposited at the time of acceptance, with a permanent DOI.

---

## License

The code in this repository is released under the **MIT License** (see `LICENSE`). The data are governed by their original providers' terms (TCGA, GEO, Human Protein Atlas, UALCAN, KM-Plotter, TIMER3.0, cBioPortal).

---

## Contact

For questions about the code or the analyses, please open an issue on the GitHub repository.
