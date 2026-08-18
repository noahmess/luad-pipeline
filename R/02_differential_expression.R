# =============================================================================
# 02_differential_expression.R
# Differential expression between tumour and adjacent-normal tissue with DESeq2.
#
# CANONICAL SPECIFICATION (transcriptome-wide):
#   - Input          : TCGA-LUAD STAR-Counts, assay "unstranded" (raw integer
#                      counts), 60660 genes x 601 samples.
#   - Samples        : 540 Primary Tumour + 59 Solid Tissue Normal = 599.
#                      The 2 Recurrent Tumour samples are excluded.
#   - Pre-filtering  : none applied by hand; DESeq2's default independent
#                      filtering is used when computing padj.
#   - Design         : ~ condition   (factor levels Normal, Tumor).
#   - Contrast       : Tumor vs Normal.
#   - Size factors   : median-of-ratios, estimated on the FULL transcriptome
#                      (all 60660 genes), NOT on the 20-gene panel.
#   - log2FC         : maximum-likelihood estimate (MLE), NO lfcShrink.
#   - padj           : Benjamini-Hochberg across all tested genes.
#   - DESeq2 1.52.0, R 4.6.0 (see outputs/sessionInfo.txt).
#   The 20 nuclear-envelope panel genes are EXTRACTED from the full-transcriptome
#   result AFTER fitting, so their fold-changes and padj reflect whole-genome
#   normalisation. This is the standard DESeq2 specification.
#
# Outputs:
#   data/TCGA_LUAD_DEG_results.RData      (results_deg, results_df)
#   outputs/tables/Differential_expression_transcriptomewide.csv  (20-gene table)
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("DESeq2", "SummarizedExperiment"))

if (file.exists(RDATA$deg)) {
  message("02: cache found (", basename(RDATA$deg), ") - skipping DESeq2.")
} else {
  load(RDATA$tcga)  # data, gene_symbols, tumor_idx, normal_idx

  # --- Build the DESeq2 dataset (tumour vs normal) ---------------------------
  counts <- SummarizedExperiment::assay(data, "unstranded")
  cols   <- c(tumor_idx, normal_idx)
  coldata <- data.frame(
    row.names = colnames(counts)[cols],
    condition = factor(c(rep("Tumor",  length(tumor_idx)),
                         rep("Normal", length(normal_idx))),
                       levels = c("Normal", "Tumor"))
  )
  dds <- DESeqDataSetFromMatrix(countData = counts[, cols],
                                colData   = coldata,
                                design    = ~ condition)

  # Map rownames to gene symbols and record the row indices of the 20 NE-panel
  # genes. NOTE: dds still holds the FULL transcriptome here - no panel
  # subsetting is done before fitting. DESeq() below estimates size factors and
  # dispersions on all 60660 genes; the 20 panel rows are EXTRACTED only
  # afterwards (results_deg[ne_rows, ]), so their log2FC and padj reflect
  # whole-transcriptome normalisation.
  rownames(dds) <- make.unique(gene_symbols)
  ne_rows <- match(NE_GENES_20, gene_symbols)        # one Ensembl row per symbol
  stopifnot(!anyNA(ne_rows))

  dds <- DESeq(dds)
  results_deg <- results(dds, contrast = c("condition", "Tumor", "Normal"))

  # --- Assemble the 20-gene results table ------------------------------------
  results_df <- as.data.frame(results_deg[ne_rows, ])
  results_df$gene <- NE_GENES_20
  # Order by significance (most significant first) for readable tables/labels.
  results_df <- results_df[order(results_df$padj), ]

  save(results_deg, results_df, file = RDATA$deg)

  # CSV export (canonical transcriptome-wide differential-expression table)
  csv_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat",
                "pvalue", "padj")
  write.csv(results_df[, csv_cols],
            file.path(PATHS$tables,
                      "Differential_expression_transcriptomewide.csv"),
            row.names = FALSE)

  message("02: ", sum(results_df$padj < 0.05, na.rm = TRUE),
          "/20 genes differentially expressed (padj < 0.05, transcriptome-wide). ",
          "Saved ", basename(RDATA$deg),
          " and outputs/tables/Differential_expression_transcriptomewide.csv")
}
