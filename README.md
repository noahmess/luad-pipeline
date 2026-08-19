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
│   ├── 08_figures.R                     Manuscript Figures 1, 2, 3, 4, 5, 8, 9, 10, 13, S1 (script's internal figure numbers differ - see "Manuscript outputs at a glance")
│   ├── 09_figure10_nup88.R              Manuscript Figure 11, NUP88 integrated profile (script's internal name says "Figure10" - see "Manuscript outputs at a glance")
│   ├── 10_signatures.R                  Definitions of the five benchmarked signatures
│   ├── 11_benchmark.R                   Cross-signature benchmark on four cohorts
│   ├── 12_sensitivity_LOO.R             Leave-one-gene-out + nested LRT for NDC1
│   ├── 13_random_signature_benchmark.R  Score A vs. 1000 random gene panels
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

**Known issue: `R/01_download_tcga.R` can fail during extraction.** `GDCdownload()` fetches the STAR-Counts query (601 files, ~2.5 GB - larger than the "several hundred MB" the script's own header suggests) in chunked tar archives and extracts each with R's `untar()`. On at least one tested configuration (macOS, system `tar` = `bsdtar` 3.5.3/libarchive, TCGAbiolinks 2.40.0) this reproducibly failed on every attempt, always partway through chunk extraction, with:

```
tar: Error exit delayed from previous errors.
Erreur dans if (ret == 1) break : l'argument est de longueur nulle
```

Reading TCGAbiolinks' source (`GDCdownload.aux`, `GDCdownload.by.chunk`) shows `untar()`'s exit status isn't propagated cleanly when the underlying `tar` reports this kind of non-fatal, delayed error (a known `bsdtar` behaviour on macOS for archives containing entries it partially can't restore) - `GDCdownload.by.chunk`'s retry loop ends up comparing `NULL` instead of a real return code, and the whole download aborts instead of retrying or failing with the library's own (fairly informative) intended error message. Two things this rules out, and one that could plausibly help but is untested here:

- `GDCdownload(query, method = "client")` uses the external `gdc-client` binary instead of R's `untar()`, which would sidestep this specific bug - but `gdc-client` is a separate program that has to be installed and on `PATH`; it is not part of this project's R dependencies and wasn't available where this was tested.
- Leftover partial downloads from a failed attempt land as numbered `.tar.gz` files and per-file UUID directories **at the repository root** (`GDCdownload`'s default working directory when a chunk fails before cleanup), not inside `GDCdata/` - delete those manually if a download attempt is interrupted.
- Passing a smaller `files.per.chunk` to `GDCdownload()` would produce smaller tar archives per extraction and might avoid whatever triggers the delayed `tar` error - plausible, but not verified: testing it requires re-attempting the full ~2.5 GB download, which wasn't safe to do again in the disk conditions this was investigated under.

**Workaround that is verified to work**, for anything downstream that only needs `data/TCGA_LUAD_data.RData` for its gene-symbol/biotype annotation (`gene_symbols`, or `rowData(data)$gene_name` / `$gene_type`) rather than the raw count matrix itself - `R/11_benchmark.R` and `R/13_random_signature_benchmark.R` are exactly this case, since the actual counts they use come from `norm_tumor` in `data/TCGA_LUAD_cohort.RData`, cached separately by `R/04`. If `data/TCGA_LUAD_DEG_results.RData` already exists (produced by `R/02_differential_expression.R`), it carries the same annotation as a `gene_info` object (`DFrame`, 60660 rows, `gene_name` and `gene_type` columns, GENCODE-derived) that can stand in for `data/TCGA_LUAD_data.RData` without re-running `R/01`:

```r
load("data/TCGA_LUAD_DEG_results.RData")   # gene_info, gene_symbols
data <- SummarizedExperiment::SummarizedExperiment(rowData = gene_info)
save(data, file = "data/TCGA_LUAD_data.RData")
```

This was cross-checked against the real thing: `gene_info$gene_name` is identical to the `gene_symbols` vector `R/01` itself produces (60660/60660 match), and re-running `R/11` and `R/13` with this substitute reproduced the existing `Random_signature_benchmark_detail.csv` / `_summary.csv` byte-for-byte (identical MD5). It does not help scripts that need the actual count matrix from `data` (`R/01`-`R/03`, `R/08`) - those still need a real download.

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

### Methylation (Table S5, Section 3.7) - not a pipeline step

**Note:** the manuscript's Figure 14 is the cBioPortal OncoPrint of genomic alterations (`R/08_figures.R`/manuscript Section 3.6-3.7, PanCancer Atlas, 566 samples), not methylation. Promoter methylation has no dedicated figure in the manuscript; it is reported as text and in **Supplementary Table S5**.

**The promoter-methylation values and p-values reported in Table S5 come entirely from the UALCAN web portal** (Chandrashekar et al., 2022; ualcan.path.uab.edu; TCGA-LUAD, 473 primary tumour vs. 32 normal), accessed on **[UALCAN access date - to confirm]**. UALCAN is the sole declared source - see *Data sources* below. The beta values themselves were read manually off UALCAN's boxplots (UALCAN does not expose them as downloadable numbers), with an estimated reading uncertainty of about 0.005 on the beta scale; the beta differences cited in the main text come from that manual reading. Because of this uncertainty, **Table S5 reports only direction (hyper-/hypomethylation) and p-value, not the beta values themselves.**

A script that recomputed promoter methylation locally from TCGA 450k data as a consistency check against these UALCAN values (`14_methylation_7genes.R`) has been moved to `archive/old_scripts/` (gitignored, not part of this repository) - it never reproduced UALCAN's probe selection/direction for all seven genes and was not retained as a source of any published number. It is not part of the reproduction instructions below.

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
| TCGA-LUAD | NCI Genomic Data Commons via `TCGAbiolinks` (STAR-Counts RNA-seq + clinical), **downloaded 30 May 2026** | `R/01`-`R/09`, `R/13` | `data/TCGA_LUAD_data.RData`, `data/TCGA_LUAD_cohort.RData` |
| GSE31210 | NCBI GEO, Okayama et al. (2012) | `R/11`, `R/12` | `data/GEO_cache/GSE31210_eset.rds` |
| GSE50081 | NCBI GEO, Der et al. (2014) | `R/11`, `R/12` | `data/GEO_cache/GSE50081_eset.rds` |
| GSE72094 | NCBI GEO (442-patient lung adenocarcinoma cohort; consult the GEO record for the associated publication) | `R/15`, `R/16` | `data/GEO_cache/GSE72094_eset.rds` |
| TIMER3.0 | Cui et al. (2025), web portal, purity-adjusted partial Spearman correlations | `R/07`, `R/16` (reference only) | `S2 TABLE TIMER3 NUP88 LMBN2.xlsx` (not tracked - see *Distribution notes*) |
| **UALCAN** | Chandrashekar et al. (2022), web portal (ualcan.path.uab.edu), TCGA-LUAD, 473 tumour / 32 normal. **Sole declared source for the promoter-methylation values and p-values in Supplementary Table S5** (manually read off the portal's boxplots and interface, not programmatically downloaded - see the *Methylation* note above). Accessed **[UALCAN access date - to confirm]**. | manuscript Section 3.7, Table S5 | - (values hard-coded in the archived consistency-check script, see below) |
| **KM-Plotter** | Gyorffy B. (2024), *The Innovation* 5(3):100625, web portal (kmplot.com), OS, JetSet best probe, adenocarcinoma histology, univariate Cox. **Sole declared source for the individual-gene external validation in Table 5.** Accessed **8 August 2026**. | manuscript Section 3.4, Table 5 | `data/kmplotter_permanent_links.txt` (settings, values and permanent `pa_id` links for LMNB2, TMPO, NDC1, LBR) |
| Human Protein Atlas, cBioPortal (PanCancer Atlas LUAD) | web portals, no programmatic download | manuscript cross-checks | - |

A local recomputation of promoter methylation from TCGA 450k beta values (`GDCdata/`, excluded from this repository because of size) was run as a consistency check against the UALCAN values above; it is not part of this repository's pipeline (see the *Methylation* note above - the script now lives in `archive/old_scripts/`, gitignored). `data/met_450k_LUAD.rds`, an untracked interactive `GDCprepare()` cache kept to avoid recomputation, was removed from the working copy - it was never written by a script and fed into no published number.

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

- **`data/`, `GDCdata/`** - raw and cached data (TCGA RangedSummarizedExperiment, DESeq2 results, GEO ExpressionSets, 450k methylation betas). Several hundred MB to multiple GB; fully regenerable from the public accessions above via `R/01`, `R/15`, or the manual GDC methylation query used by the archived methylation consistency-check script.
- **`BENCHMARK_ENV.RData` (~210 MB), `SESSION_BACKUP.RData` (~410 MB), `.RData`** - working snapshots/caches, regenerable, not required for reproduction.
- **`*.rds` at repository root (`df.rds`, `results.rds`), `gene_table-*.csv`** - orphaned artefacts not referenced by any current script; excluded rather than deleted.
- **`outputs/figures/*.tiff`, `*.png`** - binary, regenerable by `R/08`, `R/09`.
- **`*.docx`, `*.xlsx`** - manuscript drafts and the TIMER3.0 reference table (`S2 TABLE TIMER3 NUP88 LMBN2.xlsx`). Not code, not a pipeline-generated CSV, per this archive's scope. If you want the TIMER3.0 reference table bundled with the Zenodo deposit for readers who don't have the manuscript, add it back explicitly (`git add -f "S2 TABLE TIMER3 NUP88 LMBN2.xlsx"`) - it is small (24 KB).
- **`archive/`** - superseded drafts and old/superseded scripts (including the methylation consistency-check script, not retained as a source of published numbers), kept locally for the author's own traceability only.

---

## Manuscript outputs at a glance

Rebuilt from the final manuscript PDF's actual figure/table numbering (`1 final.pdf`), not from script filenames - several scripts use an older internal numbering left over from before a revision reordered and added figures (e.g. `08_figures.R` internally calls the co-expression network `"Figure4_..."` and the ROC curves `"Figure6_..."`; both are correct in content, just numbered differently inside the script than in the manuscript). Where a script's own output filename disagrees with the manuscript number, the script's internal name is given in parentheses so you can locate the code. **Bold** rows have no producing script in `R/` - either genuinely external, or genuinely missing.

### Figures

| Manuscript figure | Produced by | Used in |
|---|---|---|
| Figure 1 (volcano) | `08_figures.R` (`Figure1_Volcano_Plot`) | Section 3.1 |
| Figure 2 (heatmap, 20 genes) | `08_figures.R` (`Figure2_Heatmap_20genes`) | Section 3.1 |
| Figure 3 (GO BP, CC) | `08_figures.R` (`Figure3A_GO_BP`, `Figure3B_GO_CC`) | Section 3.1 |
| Figure 4 (Kaplan-Meier, four genes, n=394) | `08_figures.R` (internal name `Figure5_KaplanMeier_n394`) | Section 3.3 |
| Figure 5 (Kaplan-Meier, Score A, n=394) | `08_figures.R` (internal name `Figure_KM_ScoreA_n394`) | Section 3.3 |
| Figure 6 (Kaplan-Meier, Score A, stage I) | `05c_incremental_analyses.R` (`ScoreA_StageI_KM`) | Section 3.3 |
| **Figure 7 (external validation, GSE31210 + GSE50081)** | **not found** - `06_external_validation.R` computes the underlying statistics (`External_validation_GEO.csv`) but contains no plotting code; no script writes a Figure 7 image anywhere in `R/` | Section 3.4 |
| Figure 8 (ROC, tumour vs normal) | `08_figures.R` (internal name `Figure6_ROC_Curves`) | Section 3.4 |
| Figure 9 (co-expression network) | `08_figures.R` (internal name `Figure4_Correlation_Matrix`) | Section 3.5 |
| Figure 10 (correlation with tumour stage) | `08_figures.R` (internal name `Figure7_Expression_by_Stage_n394`) | Section 3.5-3.6 |
| Figure 11 (NUP88 integrated profile) | `09_figure10_nup88.R` (internal name `Figure10_NUP88_GoldenTriangle`) | Sections 3.5-3.6 |
| **Figure 12 (NUP88/LMNB2 immune infiltration, TIMER3.0)** | **not found** - purity-adjusted partial Spearman values come from the external TIMER3.0 reference (`S2 TABLE TIMER3 NUP88 LMBN2.xlsx`, see *Data sources*); no script renders this figure | Section 3.6 |
| Figure 13 (PD-L1/CD274 correlation) | `08_figures.R` (internal name `Figure9_PDL1_Scatterplots`) | Section 3.6 |
| **Figure 14 (OncoPrint, genomic alterations)** | **external** - cBioPortal web portal, PanCancer Atlas, 566 samples; not regenerated by any script | Section 3.6-3.7 |
| Figure S1 (Ki-67 correlation) | `08_figures.R` (`FigureS1_Ki67_Correlation`) | Section 3.5 |

**Orphan pipeline output, not in the manuscript:** `08_figures.R` also renders `"Figure8_Immune_Correlation"` (a 9-gene x 18-marker local Spearman correlation heatmap, fed by `07_immune_analysis.R`'s `cor_immune`). Its own header comment still says "feeds Figure 8" - that was true before the revision that replaced it with the TIMER3.0-based, purity-adjusted Figure 12 above (2 genes only, partial correlation). This local heatmap does not correspond to any figure in the final manuscript; it is still generated by `08_figures.R` and written to `outputs/figures/` on every run, but nothing in the manuscript cites it.

### Tables

| Manuscript table | Produced by | Used in |
|---|---|---|
| Table 1 (clinical characteristics, n=394) | `05b_table1_clinical.R` (`Table1_Clinical.csv`) | Section 3.2 |
| Table 2 (univariate Cox) | `05_survival_cox.R` (`Cox_univarie_n394.csv`) | Section 3.2 |
| Table 3 (multivariate Cox) | `05_survival_cox.R` (`Cox_multivarie_n394.csv`) | Section 3.2 |
| Table 4 (Score A summary, 5-yr OS) | `05_survival_cox.R` (internal name `Comparison_Individual_vs_Signature.csv`) | Section 3.3 |
| Table 5 (external validation, KM-Plotter) | external (KM-Plotter) - `data/kmplotter_permanent_links.txt` (settings, values, permanent `pa_id` links) | Section 3.4 |
| Table 6 (cross-signature benchmark) | `11_benchmark.R` (internal name `Benchmark_main_3cohorts.csv`) | Section 3.9 |
| Supplementary Table S1 (benchmark, full TCGA cohort) | `11_benchmark.R` (internal name `Benchmark_suppl_TCGAfull.csv`) | Section 3.9 |
| Supplementary Table S2 (TIMER3.0, NUP88 and LMNB2) | external (TIMER3.0) | Section 3.6 |
| Supplementary Table S3 (leave-one-gene-out) | `12_sensitivity_LOO.R` (`Supplementary_Table_S3_LOO.csv`) | Section 3.9 |
| **Supplementary Table S4 (4- vs 5-gene score, LMNB1 sensitivity)** | **not found** - no script in `R/` computes the 4-gene-vs-5-gene (+LMNB1) HR/C-index/AUC comparison reported in this table | Discussion |
| Supplementary Table S5 (promoter methylation, direction + p-value only) | external (UALCAN, manually transcribed - see *Methylation* note above) | Section 3.7 |

**Pipeline outputs that support the text but aren't a numbered table:** `05_survival_cox.R` also writes `Schoenfeld_results.csv` (proportional-hazards check, cited in Methods prose, no table number) and `05c_incremental_analyses.R` writes `Incremental_analyses_summary.csv` (LR test / C-index gain from adding Score A to the clinical model, supports Section 3.3 narrative, no table number).

| Other output | Produced by | Used in |
|---|---|---|
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
