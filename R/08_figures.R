# =============================================================================
# 08_figures.R
# Render Figures 1-9 and Figure S1 in their final, corrected versions.
# Figure 10 (NUP88 composite) is built separately in 09_figure10_nup88.R.
#
# Each figure is written to outputs/figures/ as TIFF (300 dpi, LZW) + PNG.
#
# Two fixes carried by this script (see README.md "Corrections"):
#   * Figure 1: the LMNB1 label was clipped to "LMN" because the right-most
#     point's text ran off the panel. Fixed with ggrepel + an expanded x-axis.
#   * Figures 5 and 7 are the n = 394, four-gene versions (not the old n = 540 /
#     9-gene drafts).
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("ggplot2", "ggrepel", "patchwork", "survival", "survminer",
        "pheatmap", "corrplot", "SummarizedExperiment",
        "clusterProfiler", "enrichplot"))
set.seed(SEED)

load(RDATA$tcga)    # data, gene_symbols, tumor_idx, normal_idx
load(RDATA$deg)     # results_df
load(RDATA$cohort)  # norm_tumor, expr9, expr_all, expr4_v, time_v, status_v, stage_v
load(RDATA$immune)  # cor_immune, sig_matrix, ki67_cors, immune_label_map

# Save a base-graphics / pheatmap / corrplot figure to TIFF + PNG.
save_base <- function(draw, name, width, height) {
  for (dev in c("tiff", "png")) {
    path <- file.path(PATHS$figures, paste0(name, ".", dev))
    if (dev == "tiff") tiff(path, width, height, "in", res = 300, compression = "lzw")
    else               png(path, width, height, "in", res = 300)
    draw(); dev.off()
  }
  message("  saved: ", name, ".tiff + .png")
}

# ============================================================================
# FIGURE 1 - Volcano plot (20 NE genes, tumour vs normal)
# ============================================================================
message("Figure 1: volcano plot")
volc <- results_df
volc$neg_log10_p <- -log10(volc$padj)
volc$Status <- factor(
  ifelse(volc$padj >= 0.05, "NS",
         ifelse(volc$log2FoldChange > 0, "Upregulated", "Downregulated")),
  levels = c("Upregulated", "Downregulated", "NS"))

set.seed(SEED)  # ggrepel is stochastic
p1 <- ggplot(volc, aes(log2FoldChange, neg_log10_p, color = Status)) +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = gene), size = 3, show.legend = FALSE,
                           max.overlaps = Inf, min.segment.length = 0,
                           box.padding = 0.4, seed = SEED) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = PAL$deg,
                     labels = c("Upregulated", "Downregulated", "Not significant")) +
  # Expand the x-axis so the right-most label (LMNB1) is no longer clipped.
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(x = expression(log[2]~"(Fold Change)"),
       y = expression(-log[10]~"(adjusted p-value)"), color = "Expression") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"))
save_fig(p1, "Figure1_Volcano_Plot", 8, 6)

# ============================================================================
# FIGURE 2 - Heatmap of the 20 NE genes (50 tumours + all normals)
# ============================================================================
message("Figure 2: heatmap")
idx20 <- match(NE_GENES_20, gene_symbols)
counts <- SummarizedExperiment::assay(data, "unstranded")
heat <- log2(counts[idx20, c(tumor_idx[1:50], normal_idx)] + 1)
rownames(heat) <- NE_GENES_20
annot <- data.frame(Group = c(rep("Tumor", 50), rep("Normal", length(normal_idx))))
rownames(annot) <- colnames(heat)
save_base(function()
  pheatmap(heat, scale = "row", annotation_col = annot,
           annotation_colors = list(Group = c(Tumor = "#E74C3C", Normal = "#2E86C1")),
           show_colnames = FALSE, color = PAL$diverging,
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize = 10, fontsize_row = 9),
  "Figure2_Heatmap_20genes", 10, 6)

# ============================================================================
# FIGURE 3 - GO enrichment (A: Biological Process, B: Cellular Component)
# ============================================================================
message("Figure 3: GO enrichment")
load(RDATA$enrichment)  # go_bp, go_cc
if (nrow(as.data.frame(go_bp)) > 0) {
  p3a <- barplot(go_bp, showCategory = 10) + ggtitle("GO Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  save_fig(p3a, "Figure3A_GO_BP", 9, 6)
}
if (nrow(as.data.frame(go_cc)) > 0) {
  p3b <- barplot(go_cc, showCategory = 10) + ggtitle("GO Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  save_fig(p3b, "Figure3B_GO_CC", 9, 6)
}

# ============================================================================
# FIGURE 4 - Spearman co-expression matrix of the 9-gene panel (n = 540)
# ============================================================================
message("Figure 4: co-expression matrix")
cor_mat <- cor(expr9, method = "spearman")
save_base(function()
  corrplot(cor_mat, method = "color", type = "upper", addCoef.col = "black",
           number.cex = 0.8, tl.col = "black", tl.srt = 45, tl.cex = 1,
           cl.cex = 0.9, col = PAL$diverging, mar = c(0, 0, 1, 0)),
  "Figure4_Correlation_Matrix", 7, 7)

# ============================================================================
# FIGURE 5 - Kaplan-Meier of the four signature genes (n = 394)
# p-values are dichotomised log-rank (median split); continuous-Cox HRs are in
# Table 2. TMPO and LBR are NS by log-rank but significant in continuous Cox.
# ============================================================================
message("Figure 5: Kaplan-Meier, 4 genes (n = 394)")
km_plots <- list()
for (g in SIGNATURE_4) {
  grp <- factor(ifelse(expr4_v[, g] > median(expr4_v[, g]), "High", "Low"),
                levels = c("Low", "High"))
  df  <- data.frame(time = time_v, status = status_v, grp = grp)
  fit <- survfit(Surv(time, status) ~ grp, data = df)
  km_plots[[g]] <- ggsurvplot(
    fit, data = df, pval = TRUE, pval.size = 4,
    risk.table = TRUE, risk.table.height = 0.28, title = g,
    legend.labs = c("Low expression", "High expression"), legend.title = "",
    xlab = "Time (months)", ylab = "Overall survival",
    palette = unname(PAL$surv), ggtheme = theme_classic(base_size = 12),
    font.main = c(13, "bold"), surv.median.line = "hv")
}
fig5 <- arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2, print = FALSE)
save_fig(fig5, "Figure5_KaplanMeier_n394", 11, 11)

# Bonus: Score A Kaplan-Meier (n = 394), referenced by the manuscript.
scoreA <- rowMeans(scale(expr4_v))
grpA <- factor(ifelse(scoreA > median(scoreA), "High", "Low"), c("Low", "High"))
dfA  <- data.frame(time = time_v, status = status_v, grp = grpA)
pA   <- ggsurvplot(
  survfit(Surv(time, status) ~ grp, data = dfA), data = dfA, pval = TRUE,
  pval.size = 4.5, risk.table = TRUE, risk.table.height = 0.25,
  title = "Score A : Mean Z-score (LMNB2 + TMPO + NDC1 + LBR)",
  legend.labs = c("Low Risk", "High Risk"), legend.title = "",
  xlab = "Time (months)", ylab = "Overall survival probability",
  palette = unname(PAL$surv), ggtheme = theme_classic(base_size = 12),
  surv.median.line = "hv", tables.theme = theme_cleantable())
save_fig(print(pA), "Figure_KM_ScoreA_n394", 7, 8)

# ============================================================================
# FIGURE 6 - ROC curves, tumour vs normal (9 panel genes)
# ============================================================================
message("Figure 6: ROC curves")
.need("DESeq2")
dds_roc <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[, c(tumor_idx, normal_idx)],
  colData = data.frame(cond = c(rep("T", length(tumor_idx)),
                                rep("N", length(normal_idx)))), design = ~ 1)
norm_full <- DESeq2::counts(DESeq2::estimateSizeFactors(dds_roc), normalized = TRUE)
expr_roc  <- t(log2(norm_full[match(PANEL_9, gene_symbols), ] + 1))
colnames(expr_roc) <- PANEL_9
labels_roc <- c(rep(1, length(tumor_idx)), rep(0, length(normal_idx)))
roc_cols <- c("#E74C3C", "#2E86C1", "#27AE60", "#F39C12", "#8E44AD",
              "#1ABC9C", "#E67E22", "#3498DB", "#95A5A6")
auc_vals <- numeric(length(PANEL_9))
draw_roc <- function() {
  par(mar = c(5, 5, 2, 2))
  plot(0, 0, type = "n", xlim = 0:1, ylim = 0:1,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       cex.lab = 1.2, cex.axis = 1.1)
  abline(0, 1, lty = 2, col = "grey60")
  for (i in seq_along(PANEL_9)) {
    pred <- expr_roc[, PANEL_9[i]]
    th  <- sort(unique(pred), decreasing = TRUE)
    tpr <- sapply(th, function(t) sum(pred >= t & labels_roc == 1) / sum(labels_roc == 1))
    fpr <- sapply(th, function(t) sum(pred >= t & labels_roc == 0) / sum(labels_roc == 0))
    auc_vals[i] <<- round(abs(sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)), 3)
    lines(fpr, tpr, col = roc_cols[i], lwd = 2)
  }
  legend("bottomright", legend = paste0(PANEL_9, " (AUC = ", auc_vals, ")"),
         col = roc_cols, lwd = 2, cex = 0.85, bty = "n")
}
save_base(draw_roc, "Figure6_ROC_Curves", 8, 7)
# Persist the AUC values shown on the figure (computed in draw_roc) so the
# 9 numbers are traceable without regenerating the image.
write.csv(data.frame(gene = PANEL_9, AUC = auc_vals),
          file.path(PATHS$tables, "Figure_ROC_AUC_9genes.csv"), row.names = FALSE)

# ============================================================================
# FIGURE 7 - Expression by tumour stage, four signature genes (n = 394)
# Kruskal-Wallis p annotated per panel.
# ============================================================================
message("Figure 7: expression by stage, 4 genes (n = 394)")
set.seed(SEED)  # geom_jitter
stage_plots <- list()
for (g in SIGNATURE_4) {
  df_g <- data.frame(expr = expr4_v[, g], stage = stage_v)
  df_g <- df_g[!is.na(df_g$stage), ]
  kw   <- kruskal.test(expr ~ stage, data = df_g)
  stage_plots[[g]] <- ggplot(df_g, aes(stage, expr, fill = stage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(aes(color = stage), width = 0.15, size = 0.8, alpha = 0.5) +
    annotate("text", x = 1, y = max(df_g$expr) * 1.01, hjust = 0, size = 3.5,
             label = paste0("Kruskal-Wallis, p = ",
                            formatC(kw$p.value, format = "f", digits = 3))) +
    scale_fill_manual(values = PAL$stage) +
    scale_color_manual(values = PAL$stage) +
    labs(title = g, x = "Tumor Stage", y = "Expression (log2)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "none")
}
fig7 <- (stage_plots[["LMNB2"]] | stage_plots[["LBR"]]) /
        (stage_plots[["TMPO"]]  | stage_plots[["NDC1"]])
save_fig(fig7, "Figure7_Expression_by_Stage_n394", 10, 9)

# ============================================================================
# FIGURE 8 - Immune-marker correlation heatmap (9 genes x 18 markers, n = 540)
# ============================================================================
message("Figure 8: immune correlation heatmap")
cor_disp <- cor_immune
colnames(cor_disp) <- immune_label_map[colnames(cor_disp)]
save_base(function()
  corrplot(cor_disp, method = "color", type = "full", addCoef.col = "black",
           number.cex = 0.6, tl.col = "black", tl.srt = 45, tl.cex = 0.85,
           cl.cex = 0.8, col = PAL$diverging, mar = c(0, 0, 1, 0)),
  "Figure8_Immune_Correlation", 10, 6)

# ============================================================================
# FIGURE 9 - PD-L1 (CD274) scatterplots, four signature genes (n = 540)
# ============================================================================
message("Figure 9: PD-L1 scatterplots")
pdl1_plots <- list()
for (g in SIGNATURE_4) {
  df <- data.frame(x = expr_all[, g], y = expr_all[, "CD274"])
  rr <- suppressWarnings(cor.test(df$x, df$y, method = "spearman"))
  pdl1_plots[[g]] <- ggplot(df, aes(x, y)) +
    geom_point(size = 1, alpha = 0.3, color = "#2E86C1") +
    geom_smooth(method = "lm", color = "#E74C3C", se = FALSE, linewidth = 1) +
    annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.5, size = 3.5,
             label = paste0("rho = ", round(rr$estimate, 3),
                            "\np = ", formatC(rr$p.value, format = "e", digits = 1))) +
    labs(x = paste0(g, " expression (log2)"),
         y = "PD-L1 / CD274 expression (log2)", title = paste(g, "vs PD-L1")) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}
fig9 <- (pdl1_plots[["LMNB2"]] | pdl1_plots[["TMPO"]]) /
        (pdl1_plots[["NDC1"]]  | pdl1_plots[["LBR"]])
save_fig(fig9, "Figure9_PDL1_Scatterplots", 10, 8)

# ============================================================================
# FIGURE S1 - Ki-67 (MKI67) correlation barplot (9 genes, n = 540)
# NUP88 (rho ~ 0.005) is the proliferation-uncoupled outlier.
# ============================================================================
message("Figure S1: Ki-67 correlation")
ki <- ki67_cors[order(-ki67_cors$rho), ]
ki$gene <- factor(ki$gene, levels = ki$gene)
pS1 <- ggplot(ki, aes(gene, rho, fill = rho)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = rho), vjust = -0.5, size = 3.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0.4) +
  labs(x = "", y = "Spearman correlation with Ki-67 (MKI67)") +
  ylim(0, 1) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none")
save_fig(pS1, "FigureS1_Ki67_Correlation", 7, 5)

message("08: Figures 1-9 + S1 done.")
