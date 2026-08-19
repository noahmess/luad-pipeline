# =============================================================================
# 18_dedup_sensitivity_DE.R  (SENSITIVITY ANALYSIS - does not touch main results)
# Deduplicated differential-expression sensitivity check for the pseudo-
# replication flagged by R/01: 12 tumour patients contribute multiple samples.
#
# Method (mirrors the canonical transcriptome-wide spec of R/02, only the
# tumour sample set changes):
#   - Keep ONE tumour sample per patient: for patients with several samples,
#     keep the FIRST in barcode order -> 517 unique-patient tumours.
#   - Same 59 Solid Tissue Normal samples.
#   - Full-transcriptome DESeq2, size factors on all genes, Wald test,
#     contrast Tumor vs Normal, MLE log2FC, BH padj. Extract the 20 NE genes.
#   - Compare against the main 16-DEG set (padj < 0.05 in R/02).
#
# Output (separate sensitivity file; main results untouched):
#   outputs/tables/Differential_expression_dedup517_sensitivity.csv
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("DESeq2", "SummarizedExperiment"))
load(RDATA$tcga)  # data, gene_symbols, tumor_idx, normal_idx

counts <- SummarizedExperiment::assay(data, "unstranded")
bc     <- colnames(counts)

# --- Deduplicate tumours: one sample per patient, first by barcode order ------
tum_bc  <- bc[tumor_idx]
ord     <- order(tum_bc)                       # barcode order
tum_ord <- tumor_idx[ord]
pat     <- substr(bc[tum_ord], 1, 12)
keep    <- tum_ord[!duplicated(pat)]           # first sample per patient
message("18: tumours ", length(tumor_idx), " -> ", length(keep),
        " after deduplication (", length(tumor_idx) - length(keep),
        " duplicate samples dropped).")

cols    <- c(keep, normal_idx)
coldata <- data.frame(
  row.names = bc[cols],
  condition = factor(c(rep("Tumor",  length(keep)),
                       rep("Normal", length(normal_idx))),
                     levels = c("Normal", "Tumor")))

dds <- DESeqDataSetFromMatrix(countData = counts[, cols],
                              colData = coldata, design = ~ condition)
rownames(dds) <- make.unique(gene_symbols)
ne_rows <- match(NE_GENES_20, gene_symbols)
stopifnot(!anyNA(ne_rows))
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("condition", "Tumor", "Normal"))[ne_rows, ])
res$gene <- NE_GENES_20

# --- Compare with the main analysis (committed transcriptome-wide table) -------
main <- read.csv("outputs/tables/Differential_expression_transcriptomewide.csv")
main16 <- main$gene[!is.na(main$padj) & main$padj < 0.05]     # the 16 DEGs

cmp <- merge(
  data.frame(gene = main$gene,
             log2FC_main = round(main$log2FoldChange, 3),
             padj_main = main$padj),
  data.frame(gene = res$gene,
             log2FC_dedup = round(res$log2FoldChange, 3),
             padj_dedup = res$padj),
  by = "gene")
cmp$in_main_16DEG      <- cmp$gene %in% main16
cmp$sig_dedup          <- !is.na(cmp$padj_dedup) & cmp$padj_dedup < 0.05
cmp$direction_preserved <- sign(cmp$log2FC_main) == sign(cmp$log2FC_dedup)
cmp <- cmp[order(cmp$padj_dedup), ]

write.csv(cmp, "outputs/tables/Differential_expression_dedup517_sensitivity.csv",
          row.names = FALSE)

# --- Verdict on the 16 main DEGs ----------------------------------------------
d16 <- cmp[cmp$in_main_16DEG, ]
cat("\n=== 16 main DEGs under deduplication (517 tumours) ===\n")
print(d16[, c("gene", "log2FC_main", "log2FC_dedup", "padj_dedup",
              "sig_dedup", "direction_preserved")], row.names = FALSE)
cat(sprintf("\n16 DEGs still significant (padj<0.05): %d/16 | direction preserved: %d/16\n",
            sum(d16$sig_dedup), sum(d16$direction_preserved)))
message("18: wrote Differential_expression_dedup517_sensitivity.csv")
