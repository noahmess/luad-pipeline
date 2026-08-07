# =============================================================================
# run_all.R
# Run the whole pipeline end to end, from the project root:
#     Rscript run_all.R
# or, inside an R session started at the project root:
#     source("run_all.R")
#
# Steps 01-04 cache their .RData outputs in data/ and are skipped on re-runs;
# delete the relevant data/*.RData file to force recomputation.
# =============================================================================

PROJECT_ROOT <- getwd()
source("R/00_setup.R")

steps <- c(
  "R/01_download_tcga.R",          # GDC download (slow, cached)
  "R/02_differential_expression.R",# DESeq2, 20 NE genes
  "R/03_functional_enrichment.R",  # GO / KEGG
  "R/04_cohort_harmonization.R",   # n = 540 and n = 394 cohorts
  "R/05_survival_cox.R",           # Cox uni/multi, Score A, Schoenfeld
  "R/06_external_validation.R",    # GEO validation (reconstructed)
  "R/07_immune_analysis.R",        # immune / Ki-67 / PD-L1
  "R/08_figures.R",                # Figures 1-9 + S1
  "R/09_figure10_nup88.R"          # Figure 10 composite
)

for (s in steps) {
  message("\n========== ", s, " ==========")
  t0 <- Sys.time()
  source(s)
  message(sprintf("---------- done in %.1f s", as.numeric(Sys.time() - t0, "secs")))
}
message("\nPipeline complete. Figures in outputs/figures/, tables in outputs/tables/.")
