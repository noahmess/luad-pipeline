# =============================================================================
# 02_differential_expression.R
# Differential expression between tumour and adjacent-normal tissue with DESeq2,
# restricted to the 20 nuclear-envelope genes selected a priori.
#
# Cohort: all tumour + normal samples (the "n = 540" expression side; DESeq2
# itself uses every tumour and every normal column - no survival filtering).
#
# Output: data/TCGA_LUAD_DEG_results.RData containing
#   results_deg : full DESeqResults object
#   results_df  : data.frame of the 20 NE genes (log2FoldChange, padj, gene, ...)
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

  # Collapse to gene symbols: keep the 20 NE genes (first matching row each).
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
  message("02: ", sum(results_df$padj < 0.05, na.rm = TRUE),
          "/20 genes differentially expressed (padj < 0.05). Saved ",
          basename(RDATA$deg))
}
