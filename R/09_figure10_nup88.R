# =============================================================================
# 09_figure10_nup88.R
# Figure 10 - the NUP88 "golden triangle" composite, regenerated ENTIRELY from
# R (no image editing). Four panels:
#   (a) Spearman rho with Ki-67       (NUP88 ~ 0.005 vs 0.49-0.82 for the rest)
#   (b) Spearman rho with tumour stage(NUP88 ~ 0.028 vs 0.10-0.19 for the rest)
#   (c) Spearman rho with immune markers (NUP88 negative vs mixed for others)
#   (d) Summary comparison table
#
# CORRECTION carried here: panel (d), row "Cox univariate", previously read
#   "Yes (HR = 1.37)"  ->  now "No (HR = 1.159, p = 0.090)".
# NUP88 is NOT significant in univariate Cox on the harmonised n = 394 cohort.
#
# The per-gene Spearman values below are the published values; they are
# reproducible from 07_immune_analysis.R (Ki-67 / immune) and from the stage
# correlation in the survival data. They are listed explicitly so the figure
# renders identically and independently of cached objects.
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("ggplot2", "ggpubr"))
set.seed(SEED)

genes <- c("LMNB1", "LMNB2", "TMPO", "NDC1", "NUP107", "NUP62", "LBR", "BANF1", "NUP88")
ki67_rho  <- c(0.822, 0.747, 0.642, 0.610, 0.605, 0.528, 0.492, 0.160, 0.005)
stage_rho <- c(0.152, 0.193, 0.111, 0.108, 0.101, 0.112, 0.021, 0.150, 0.028)
grp <- ifelse(genes == "NUP88", "NUP88", "Other genes")

# A reusable horizontal bar panel (Ki-67 / stage).
bar_panel <- function(rho, ylab, title, ymax) {
  df <- data.frame(Gene = factor(genes, levels = genes[order(rho)]),
                   rho = rho, group = grp)
  lev <- levels(df$Gene)
  ggplot(df, aes(Gene, rho, fill = group)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = round(rho, 3)), hjust = -0.2, size = 2.8) +
    scale_fill_manual(values = c(NUP88 = PAL$nup88[["NUP88"]],
                                 "Other genes" = PAL$nup88[["Other"]])) +
    coord_flip() + ylim(min(0, min(rho)) - 0.05, ymax) +
    labs(x = "", y = ylab, title = title) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14, hjust = 0),
          axis.text.y = element_text(
            face  = ifelse(lev == "NUP88", "bold", "plain"),
            color = ifelse(lev == "NUP88", PAL$nup88[["NUP88"]], "grey30")))
}

pA <- bar_panel(ki67_rho,  expression("Spearman " * rho * " with Ki-67"), "a", 1.0)
pB <- bar_panel(stage_rho, expression("Spearman " * rho * " with stage"), "b", 0.25)

# Panel C: immune markers, NUP88 vs mean of the other eight genes.
immune_markers <- c("CD4", "FOXP3", "CTLA4", "PD-L1", "LAG3",
                    "CD8A", "NK", "CD68", "IFN-gamma")
nup88_immune <- c(-0.297, -0.319, -0.269, -0.254, -0.222, -0.080, -0.189, -0.182, -0.106)
mean_other   <- c(-0.128,  0.057,  0.012,  0.081,  0.168,  0.100, -0.091, -0.053,  0.179)
df_imm <- data.frame(
  Marker = factor(rep(immune_markers, 2), levels = immune_markers),
  rho = c(nup88_immune, mean_other),
  Group = factor(rep(c("NUP88", "Mean (8 other genes)"), each = 9),
                 levels = c("NUP88", "Mean (8 other genes)")))
pC <- ggplot(df_imm, aes(Marker, rho, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(NUP88 = PAL$nup88[["NUP88"]],
                               "Mean (8 other genes)" = PAL$nup88[["Other"]])) +
  ylim(-0.4, 0.3) +
  labs(x = "", y = expression("Spearman " * rho), title = "c", fill = "") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        legend.position = c(0.8, 0.9), legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# Panel D: summary table - CORRECTED Cox row.
comparison <- data.frame(
  Feature = c("Upregulated in LUAD", "Cox univariate (p < 0.05)",
              "Ki-67 correlation", "Stage correlation",
              "Module co-expression", "Immune markers (17/18 negative)"),
  NUP88 = c("Yes",
            "No (HR = 1.159, p = 0.090)",   # corrected (was "Yes (HR = 1.37)")
            "0.005 (NS)", "0.028 (NS)", "Uncorrelated", "Yes"),
  Other_8_genes = c("Yes", "8/8 Yes", "0.49 to 0.82", "0.10 to 0.19",
                    "Tightly coupled", "No (mixed)"))
pD <- ggtexttable(comparison, rows = NULL,
                  cols = c("Feature", "NUP88", "Other 8 genes"),
                  theme = ttheme(
                    colnames.style = colnames_style(fill = "grey90", color = "grey20",
                                                    face = "bold", size = 9),
                    tbody.style = tbody_style(fill = c("white", "grey97"),
                                              color = "grey20", size = 8,
                                              hjust = 0, x = 0.05))) +
  labs(title = "d") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

# Assemble the 2 x 2 composite and save.
top    <- ggarrange(pA, pB, ncol = 2, widths = c(1, 1))
bottom <- ggarrange(pC, pD, ncol = 2, widths = c(1.2, 0.8))
fig10  <- ggarrange(top, bottom, nrow = 2, heights = c(1, 1.1))

for (dev in c("tiff", "png")) {
  path <- file.path(PATHS$figures, paste0("Figure10_NUP88_GoldenTriangle.", dev))
  if (dev == "tiff") tiff(path, 12, 10, "in", res = 300, compression = "lzw")
  else               png(path, 12, 10, "in", res = 300)
  print(fig10); dev.off()
}
message("09: Figure 10 written with corrected panel (d): ",
        "NUP88 Cox univariate = No (HR = 1.159, p = 0.090).")
