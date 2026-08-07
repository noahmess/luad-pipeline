# =============================================================================
# 03_functional_enrichment.R
# Gene Ontology (BP / CC / MF) and KEGG over-representation analysis of the
# dysregulated nuclear-envelope genes, with clusterProfiler.
#
# Output: data/TCGA_LUAD_enrichment.RData containing
#   go_bp, go_cc, go_mf : enrichResult objects
#   kegg                : enrichResult object (may be empty)
#   entrez              : SYMBOL -> ENTREZID mapping used
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("clusterProfiler", "org.Hs.eg.db"))

if (file.exists(RDATA$enrichment)) {
  message("03: cache found (", basename(RDATA$enrichment), ") - skipping.")
} else {
  # Enrichment is run on the 9-gene up-regulated panel (the coherent
  # over-expressed module that drives the biological story).
  entrez <- bitr(PANEL_9, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  message("03: ", nrow(entrez), "/", length(PANEL_9), " genes mapped to Entrez.")

  go_args <- list(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  go_bp <- do.call(enrichGO, c(go_args, ont = "BP"))
  go_cc <- do.call(enrichGO, c(go_args, ont = "CC"))
  go_mf <- do.call(enrichGO, c(go_args, ont = "MF"))
  kegg  <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                      pvalueCutoff = 0.05)

  save(go_bp, go_cc, go_mf, kegg, entrez, file = RDATA$enrichment)
  message("03: enrichment saved -> ", basename(RDATA$enrichment))
}
