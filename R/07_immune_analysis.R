# =============================================================================
# 07_immune_analysis.R
# Microenvironment analysis on the n = 540 tumour cohort:
#   * Spearman correlation of the 9 panel genes with 18 immune / proliferation
#     markers (matrix + significance stars) -> feeds Figure 8.
#   * Correlation with Ki-67 (MKI67) -> feeds Figure S1; NUP88 is the outlier.
#   * Correlation with PD-L1 (CD274) for the four signature genes -> Figure 9.
#   * Nuclear-envelope score vs cytotoxic and checkpoint scores.
#
# Output: data/TCGA_LUAD_immune.RData containing
#   cor_immune, p_immune : 9 x 18 Spearman rho and p
#   sig_matrix           : 9 x 18 significance stars ("*", "**", "***")
#   ki67_cors            : per-gene rho with MKI67
#   immune_label_map     : marker -> human-readable label
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
load(RDATA$cohort)  # expr_all, found_markers, expr9
set.seed(SEED)

if (file.exists(RDATA$immune)) {
  message("07: cache found (", basename(RDATA$immune), ") - skipping.")
} else {
  markers <- found_markers

  # ---- Spearman rho and p between each panel gene and each marker -----------
  cor_immune <- p_immune <- matrix(NA_real_, length(PANEL_9), length(markers),
                                   dimnames = list(PANEL_9, markers))
  for (g in PANEL_9) for (m in markers) {
    tt <- suppressWarnings(cor.test(expr_all[, g], expr_all[, m],
                                    method = "spearman"))
    cor_immune[g, m] <- tt$estimate
    p_immune[g, m]   <- tt$p.value
  }

  # Significance stars matrix (same shape).
  sig_matrix <- matrix("", nrow(p_immune), ncol(p_immune),
                       dimnames = dimnames(p_immune))
  sig_matrix[p_immune < 0.05]  <- "*"
  sig_matrix[p_immune < 0.01]  <- "**"
  sig_matrix[p_immune < 0.001] <- "***"

  # Human-readable marker labels for the figure.
  immune_label_map <- c(
    CD8A = "CD8+ T cells", CD4 = "CD4+ T cells", FOXP3 = "Tregs",
    CD68 = "Macrophages", CD163 = "M2 Macrophages", NOS2 = "M1 Macrophages",
    CD19 = "B cells", NCAM1 = "NK cells", CD274 = "PD-L1", PDCD1 = "PD-1",
    CTLA4 = "CTLA-4", HAVCR2 = "TIM-3", LAG3 = "LAG-3", IFNG = "IFN-gamma",
    GZMA = "Granzyme A", PRF1 = "Perforin", MKI67 = "Ki-67", PECAM1 = "CD31")

  # ---- Ki-67 (proliferation) correlation, all 9 genes -----------------------
  ki67_cors <- NULL
  if ("MKI67" %in% colnames(expr_all)) {
    for (g in PANEL_9) {
      tt <- suppressWarnings(cor.test(expr_all[, g], expr_all[, "MKI67"],
                                      method = "spearman"))
      ki67_cors <- rbind(ki67_cors,
                         data.frame(gene = g, rho = round(tt$estimate, 3),
                                    pvalue = tt$p.value))
    }
    ki67_cors <- ki67_cors[order(-ki67_cors$rho), ]
  }

  # ---- Nuclear-envelope score vs cytotoxic / checkpoint scores --------------
  score_env  <- rowMeans(scale(expr_all[, PANEL_9]))
  cyto       <- intersect(c("CD8A", "GZMA", "PRF1", "IFNG"), colnames(expr_all))
  checkpoint <- intersect(c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"),
                          colnames(expr_all))
  rho_cyto  <- cor(score_env, rowMeans(scale(expr_all[, cyto])),       method = "spearman")
  rho_check <- cor(score_env, rowMeans(scale(expr_all[, checkpoint])), method = "spearman")
  message(sprintf("07: NE score vs cytotoxic rho = %.3f | vs checkpoint rho = %.3f",
                  rho_cyto, rho_check))

  save(cor_immune, p_immune, sig_matrix, ki67_cors, immune_label_map,
       file = RDATA$immune)
  write.csv(ki67_cors, file.path(PATHS$tables, "Ki67_correlation.csv"),
            row.names = FALSE)
  message("07: immune objects saved -> ", basename(RDATA$immune))
}
