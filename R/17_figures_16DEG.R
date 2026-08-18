# =============================================================================
# 17_figures_16DEG.R  (Path A, transcriptome-wide canonical DE)
# Regenerates Figure 1 (volcano) and Figure 3 (GO enrichment) for the 16-DEG
# set, reading ONLY the committed output tables on this branch:
#   outputs/tables/Differential_expression_transcriptomewide.csv   (Figure 1)
#   outputs/tables/Enrichment_16DEG_GO_KEGG.csv                    (Figure 3)
# Visual style is kept identical to R/08_figures.R (same PAL$deg colours,
# ggrepel labels, theme_classic; clusterProfiler-style GO barplots).
#
# Outputs (outputs/figures/, PNG + TIFF, 300 dpi):
#   Figure1_volcano_16DEG.{png,tiff}
#   Figure3_enrichment_16DEG.{png,tiff}
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork)
})
SEED <- 42
FIG_DIR <- "outputs/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Palette identical to R/00_setup.R PAL$deg
PAL_DEG <- c(Upregulated = "#E74C3C", Downregulated = "#2E86C1", NS = "grey60")

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".tiff")), plot,
         width = width, height = height, dpi = 300, compression = "lzw")
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300)
  message("  saved: ", name, ".tiff + .png")
}

# ---------------------------------------------------------------------------
# FIGURE 1 - Volcano plot (20 NE genes, tumour vs normal; 16 significant)
# ---------------------------------------------------------------------------
volc <- read.csv("outputs/tables/Differential_expression_transcriptomewide.csv")
volc$neg_log10_p <- -log10(volc$padj)
volc$Status <- factor(
  ifelse(volc$padj >= 0.05, "NS",
         ifelse(volc$log2FoldChange > 0, "Upregulated", "Downregulated")),
  levels = c("Upregulated", "Downregulated", "NS"))

set.seed(SEED)
p1 <- ggplot(volc, aes(log2FoldChange, neg_log10_p, color = Status)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = gene), size = 3, show.legend = FALSE,
                  max.overlaps = Inf, min.segment.length = 0,
                  box.padding = 0.4, seed = SEED) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = PAL_DEG,
                     labels = c("Upregulated", "Downregulated", "Not significant")) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(x = expression(log[2]~"(Fold Change)"),
       y = expression(-log[10]~"(adjusted p-value)"), color = "Expression") +
  theme_classic(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"))
save_fig(p1, "Figure1_volcano_16DEG", 8, 6)

# ---------------------------------------------------------------------------
# FIGURE 3 - GO enrichment (Panel A: Biological Process, B: Cellular Component)
# clusterProfiler-style barplot: x = gene count, fill = adjusted p-value.
# ---------------------------------------------------------------------------
enr <- read.csv("outputs/tables/Enrichment_16DEG_GO_KEGG.csv")

go_bar <- function(df, title, n = 10) {
  d <- df[order(df$padj), ]
  d <- head(d, n)
  d$term <- factor(d$term, levels = rev(d$term[order(d$count)]))
  ggplot(d, aes(x = count, y = term, fill = padj)) +
    geom_col() +
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust",
                          guide = guide_colorbar(reverse = TRUE)) +
    labs(x = "Count", y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

p3a <- go_bar(enr[enr$category == "GO_BP", ], "GO Biological Process")
p3b <- go_bar(enr[enr$category == "GO_CC", ], "GO Cellular Component")
p3  <- (p3a / p3b) + plot_annotation(tag_levels = "A")
save_fig(p3, "Figure3_enrichment_16DEG", 9, 9)

message("17: figures written to ", FIG_DIR)
