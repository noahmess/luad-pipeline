# =============================================================================
# 03_functional_enrichment.R
# Gene Ontology (BP / CC / MF) and KEGG over-representation analysis of the
# differentially expressed nuclear-envelope genes, with clusterProfiler.
#
# GENE SET (canonical, Voie A): the differentially expressed genes from the
# transcriptome-wide DESeq2 analysis (R/02), i.e. results_df rows with
# padj < 0.05. Under the transcriptome-wide specification this is the 16-gene
# set (LMNB1, NUP155, NDC1, LMNB2, SUN2, NUP62, NUP107, TMPO, NUP88, SYNE1,
# LBR, TOR1AIP1, BANF1, SUN1, ANKLE2, LMNA). Enrichment is run on the full DE
# set, not on a hand-picked sub-panel.
#
# STRING network (nodes / edges / clustering coefficient / PPI-enrichment p) is
# obtained separately from the STRING v12 REST API on the same DE set:
#   https://string-db.org/api/tsv/ppi_enrichment?identifiers=<genes>&species=9606
# (default medium confidence 0.4). Its summary is written to
#   outputs/tables/STRING_network_16DEG.csv
#
# Outputs:
#   data/TCGA_LUAD_enrichment.RData   (go_bp, go_cc, go_mf, kegg, entrez)
#   outputs/tables/Enrichment_16DEG_GO_KEGG.csv
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("clusterProfiler", "org.Hs.eg.db"))

if (file.exists(RDATA$enrichment)) {
  message("03: cache found (", basename(RDATA$enrichment), ") - skipping.")
} else {
  load(RDATA$deg)  # results_df (transcriptome-wide DE, from R/02)
  deg <- results_df$gene[!is.na(results_df$padj) & results_df$padj < 0.05]
  message("03: enrichment on ", length(deg), " differentially expressed genes.")

  entrez <- bitr(deg, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  message("03: ", nrow(entrez), "/", length(deg), " genes mapped to Entrez.")

  go_args <- list(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  go_bp <- do.call(enrichGO, c(go_args, ont = "BP"))
  go_cc <- do.call(enrichGO, c(go_args, ont = "CC"))
  go_mf <- do.call(enrichGO, c(go_args, ont = "MF"))
  kegg  <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                      pvalueCutoff = 0.05)

  save(go_bp, go_cc, go_mf, kegg, entrez, file = RDATA$enrichment)

  # CSV export of all significant terms
  mk <- function(x, ont) {
    if (is.null(x)) return(NULL)
    d <- as.data.frame(x); if (nrow(d) == 0) return(NULL)
    data.frame(category = ont, ID = d$ID, term = d$Description,
               pvalue = signif(d$pvalue, 3), padj = signif(d$p.adjust, 3),
               count = d$Count, geneID = d$geneID)
  }
  all_terms <- do.call(rbind, list(mk(go_cc, "GO_CC"), mk(go_bp, "GO_BP"),
                                   mk(go_mf, "GO_MF"), mk(kegg, "KEGG")))
  write.csv(all_terms,
            file.path(PATHS$tables, "Enrichment_16DEG_GO_KEGG.csv"),
            row.names = FALSE)
  message("03: enrichment saved -> ", basename(RDATA$enrichment),
          " and outputs/tables/Enrichment_16DEG_GO_KEGG.csv")
}
